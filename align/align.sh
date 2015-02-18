#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
# Align input interleaved fastq files to a reference genome
# Input either:
#   1. $interleaved=file.fastq.gz - interleaved paired end Illumina reads
#   2. $paired1, $paired2 - first and second reads of Illumina paired end reads
#
# Output:
# chr1.bam, chr1.bam.bai .. chrN.bam, chrN.bam.bai
#

set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)


if [ "$interleaved" != "" ] # input must be defined
then
    input=$(./swe fetch $interleaved)
    #-p for interleaved sequences
    bwa_opts="-p"
else
    #paired reads are specified
    paired1=$(./swe fetch $paired1)
    paired2=$(./swe fetch $paired2)
    input="$paired1 $paired2"
fi

[ "$sample_id" != "" ] # sample_id must be defined

[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)  # get local copy of reference files


cpu_cores=32
GROUP_ID="@RG\tID:1\tPL:ILLUMINA\tPU:pu\tLB:group1\tSM:$sample_id"



#samtools sort is limited to 62M per thread, so that it won't take more than 2GB total
bwa mem -M $bwa_opts -t $cpu_cores -R "$GROUP_ID"   $gatk_data/hg19/ucsc.hg19.fasta $input \
    | samtools view -@ $cpu_cores -F 0x100 -1 -bt   $gatk_data/hg19/ucsc.hg19.fasta.fai - \
    | samtools sort -@ $cpu_cores -m $[2000/$cpu_cores]M -l 0 - raw

samtools index raw.bam


#split alignments by chromosomes. Emit separate bam files for each contig:  chr1.bam chr2.bam ...
#will require some modofication for 1000g reference build without "chr" prefixes
chr_list=$(samtools idxstats raw.bam| cut -f 1 |grep chr)
for chr in $chr_list
do
	samtools view -@ $cpu_cores -F 4 -b raw.bam  $chr > $chr.bam
	samtools index $chr.bam 
	# run emits in parallel
	{ ./swe emit file $chr.bam     || touch emit.failed & }
	{ ./swe emit file $chr.bam.bai || touch emit.failed & }
done

wait
[ ! -e emit.failed ]


exit 0

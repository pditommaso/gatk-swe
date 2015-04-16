params.bqsr_chr = 'chr22'
params.splits = 10

/* 
 # Copyright (C) 2015 Center for Genomic Regulation 
 # Copyright (C) 2013-2015 Clusterk Inc
 # Align input interleaved fastq files to a reference genome
 # Input either:
 #   1. $interleaved=file.fastq.gz - interleaved paired end Illumina reads
 #   2. $paired1, $paired2 - first and second reads of Illumina paired end reads
 #
 # Output:
 # chr1.bam, chr1.bam.bai .. chrN.bam, chrN.bam.bai
 #
 */


process align {
  input: 
  file genome_fasta 
  file genome_fai 
  file paired1 
  file paired2 
  
  output: 
  file 'chr*.bam' into chr_splits

  """
    #samtools sort is limited to 62M per thread, so that it won't take more than 2GB total
    bwa mem -M -t ${tas.cpus} -R "$GROUP_ID"  $genome_fasta $paired1 $paired2 \
        | samtools view -@ ${tasks.cpus} -F 0x100 -1 -bt $genome_fai - \
        | samtools sort -@ ${tasks.cpus} -l 0 - raw 

    samtools index raw.bam


    #split alignments by chromosomes. Emit separate bam files for each contig:  chr1.bam chr2.bam ...
    #will require some modification for 1000g reference build without "chr" prefixes
    chr_list=\$(samtools idxstats raw.bam| cut -f 1 |grep chr)
    for chr in \$chr_list
    do
        samtools view -@ ${task.cpus} -F 4 -b raw.bam  \$chr > \$chr.bam
        #samtools index \$chr/\$chr.bam 
    done

  """
}


chr_groupped = chr_splits
                  .groupBy { it.name }
                  .flatMap { it.collect { item -> tuple(item.key, item.value) } }

process combine {

  input: 
  set chr, file(bam_list) chr_groupped 

  output: 
  set chr, file('*.bam'), file('*.bai') into chr_combined
  set chr, file('*.bam'), file('*.bai') into chr_combined2 
  
  """
  samtools merge -@ ${task.cpus} ${chr}.bam $bam_list
  samtools index ${chr}.bam
  """

}

chr_bqsr = chr_combined2.filter { chr, bam, bai -> chr == params.bqsr_chr  }

process bqsr {

    input:
    file genome
    file dbsnp
    set chr, file(bam), file(bai) from chr_bqsr 
    
    output:
    file 'bqsr.grp'

  """

    tabix -h $dbsnp $chr:1-300000000 > interval.dbsnp.vcf
    
    gatk.sh \
        -T BaseRecalibrator \
        -I $bam \
        -R $genome \
        -knownSites interval.dbsnp.vcf \
        --covariate ContextCovariate  \
        --covariate ReadGroupCovariate  \
        --covariate QualityScoreCovariate \
        --covariate CycleCovariate  \
        -o bqsr.grp  \
        --disable_indel_quals \
        -nct ${tasks.cpus} 

  """

}

process chr_split {
    input: 
    file genome
    set chr, file(input), file('*.bai') from chr_combined 
    
    output: 
     set chr, file('split_*.bam'), file('split_*.bai'), file('*.interval') into chr_split 

    """ 
    split_chr/equally_spaced_intervals.pl --bam $input --blocks ${params.splits} --chr $chr > es_intervals.txt

    cat es_intervals.txt | parallel -j ${task.cpus} "split_chr/advanced_splitter.pl --reference $genome --bam $input --region {} > {}.breakpoints"


    for interval in `cat es_intervals.txt`
    do
        cat $interval.breakpoints >> breakpoints.lst 
    done

    breakpoints2intervals.pl --bam $input --breakpoints breakpoints.lst --splits=$splits > interval.lst 

    for i in $(seq 1 $splits)
    do
        interval=$(head -n $i interval.lst|tail -n 1)
        echo $interval > $i.interval

        samtools view -@ $cpu_cores -b $input $interval > split_$i.bam 
        samtools index split_$i.bam 

    done
    """

}


process variant_call {
    input: 
    file genome
    file bqsr
    file interval
    
    output: 
    file 'raw.vcf'

"""

tabix -h $gatk_data/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz $interval > VCF1.vcf
tabix -h $gatk_data/1000G_phase1.indels.hg19.vcf.gz $interval                  > VCF2.vcf
tabix -h $gatk_data/dbsnp_137.hg19.vcf.gz $interval         > interval.dbsnp_137.hg19.vcf


java -Xmx2g -jar gatk/MarkDuplicates.jar INPUT=$input  OUTPUT=deduplicated.bam REMOVE_DUPLICATES=true METRICS_FILE=duplication.metrics
samtools index deduplicated.bam

gatk.sh \
    -I deduplicated.bam \
    -R $genome \
    -T RealignerTargetCreator \
    -o realigned.intervals \
    --known VCF1.vcf \
    --known VCF2.vcf \
    -L $interval 

gatk.sh \
    -I deduplicated.bam \
    -R $genome \
    -T IndelRealigner \
    -targetIntervals realigned.intervals \
    -o realigned.bam \
    -known VCF1.vcf \
    -known VCF2.vcf \
    --consensusDeterminationModel KNOWNS_ONLY \
    -L $interval  \
    -LOD 0.4  \
    --downsample_to_coverage 250 \
    -compress 0

gatk.sh \
    -R $genome \
    -I realigned.bam  \
    -T PrintReads  \
    -o recalibrated.bam  \
    --disable_indel_quals  \
    -BQSR $bqsr \
    -compress 0


#full featured GATK, run HaplotypeCaller
gatk.sh \
    -R $genome \
    -T HaplotypeCaller \
    -I recalibrated.bam \
    --dbsnp interval.dbsnp_137.hg19.vcf \
    -stand_call_conf 30.0 \
    -stand_emit_conf 30.0 \
    -o raw.hc.vcf \
    -L $interval \
    -nct ${task.cpus} \
    -rf BadCigar

gatk.sh \
    -R $genome \
    -T VariantAnnotator \
    -I recalibrated.bam \
    --variant raw.hc.vcf \
    -A MappingQualityZero \
    --dbsnp interval.dbsnp_137.hg19.vcf \
    -o raw.vcf \
    -L $interval \
    -nt ${task.cpus} \
    -rf BadCigar
"""
}

process combine_vcf {


file 'raw.vcf'

"""
#
# Combine multiple VCFs for the same sample into a single VCF
# Input: $input space separated list of input files (1234:chr2.bam 1235:chr2.bam 1236:chr2.bam) 
# Output: raw.vcf 

for vcf in $input
do
cat $(./swe fetch $vcf) >> raw.vcf.tmp
done


#obtain header from the first VCF in the list
set -- $input
first_vcf=$1

cat $(./swe fetch $first_vcf) |grep -P '^\#' >header.vcf


cat header.vcf > raw.vcf
cat raw.vcf.tmp | grep -vP "^\#" | vcf-sort -c >> raw.vcf
"""
}


process vqsr {

}
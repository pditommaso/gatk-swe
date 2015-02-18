#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
#
# Split chromsome into N subregions by finding safe split locations
# Input:
#        $input  - input bam.file (chr4.bam)
#        $splits - number of splits
#        $chr    - chromosome, we expect that input bam file contains only one chromosome
#

set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)

[ "$input"  != "" ] && input=$(./swe fetch $input)
[ "$splits" != "" ]
[ "$chr" != "" ]
[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)

cpu_cores=32

samtools index $input #fix it 
split_chr/equally_spaced_intervals.pl --bam $input --blocks $splits --chr $chr > es_intervals.txt

cat es_intervals.txt | parallel -j $cpu_cores "split_chr/advanced_splitter.pl --reference $gatk_data/hg19/ucsc.hg19.fasta --bam $input --region {} > {}.breakpoints"


[ ! -e breakpoints.lst ] || rm breakpoints.lst

for interval in `cat es_intervals.txt`
do
    cat $interval.breakpoints >> breakpoints.lst 
done

split_chr/breakpoints2intervals.pl --bam $input --breakpoints breakpoints.lst --splits=$splits > interval.lst 

for i in $(seq 1 $splits)
do
    interval=$(head -n $i interval.lst|tail -n 1)
    echo $interval > $i.interval

    { samtools view -@ $cpu_cores -b $input $interval > $i.bam && samtools index $i.bam && ./swe emit file $i.interval $i.bam $i.bam.bai && rm $i.bam || touch emit.failed & }

    #run 5 processes in parallel at most
    [ "$(($i % 5))" != "0" ]  || wait

    [ ! -e emit.failed ]

done
wait


[ ! -e emit.failed ]

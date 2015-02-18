#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
#
# Run base quality recalibration
# Input:  $input  - input bam file
# Output: bqsr.grp
#
set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)


[ "$input" != "" ]
[ "$chr" != "" ]
[ "$gatk_jar" != "" ]
[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)

bai=$(./swe fetch $input.bai)
input=$(./swe fetch $input)


gatk_jar=$(./swe fetch $gatk_jar)

cpu_cores=32


samtools index $input #fix me

tabix -h $gatk_data/dbsnp_137.hg19.vcf.gz $chr:1-300000000 > interval.dbsnp_137.hg19.vcf
java -Xmx6g -jar $gatk_jar \
    -T BaseRecalibrator \
    -I $input \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
    -knownSites interval.dbsnp_137.hg19.vcf \
    --covariate ContextCovariate  \
    --covariate ReadGroupCovariate  \
    --covariate QualityScoreCovariate \
    --covariate CycleCovariate  \
    -o bqsr.grp  \
    --disable_indel_quals \
    -nct $cpu_cores 

./swe emit file bqsr.grp
exit 0


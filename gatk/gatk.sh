#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
#
# Run Variant calling for a chromosomal region (Deduplications, Indel realignment, Variant Calling, and Variant Annotation to add mq0)
# Input:
#      $input    - input bam file (e.g. chr3:1-10000.bam)
#      $bqsr     - base quality recalibration file (bqsr.grp)
#      $interval - chromosomal interval (chr3:1-10000)
#
# Output: raw.vcf

set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)

[ "$input"     != "" ] && input=$(./swe fetch $input)
[ "$gatk_jar"  != "" ] && gatk_jar=$(./swe fetch $gatk_jar)
[ "$bqsr"      != "" ] && bqsr=$(./swe fetch $bqsr)
[ "$interval"  != "" ] && interval_file=$(./swe fetch $interval)
[ "$GATK_REFERENCE" != "" ] && gatk_data=$(./swe dc $GATK_REFERENCE)

interval=$(cat $interval_file)
cpu_cores=32

tabix -h $gatk_data/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz $interval > VCF1.vcf
tabix -h $gatk_data/1000G_phase1.indels.hg19.vcf.gz $interval                  > VCF2.vcf
tabix -h $gatk_data/dbsnp_137.hg19.vcf.gz $interval         > interval.dbsnp_137.hg19.vcf


java -Xmx2g -jar gatk/MarkDuplicates.jar INPUT=$input  OUTPUT=deduplicated.bam REMOVE_DUPLICATES=true METRICS_FILE=duplication.metrics
samtools index deduplicated.bam

java -Xmx2g -jar $gatk_jar \
    -I deduplicated.bam \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
    -T RealignerTargetCreator \
    -o realigned.intervals \
    --known VCF1.vcf \
    --known VCF2.vcf \
    -L $interval 

java -Xmx2g -jar $gatk_jar \
    -I deduplicated.bam \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
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

java -Xmx2g -jar $gatk_jar \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
    -I realigned.bam  \
    -T PrintReads  \
    -o recalibrated.bam  \
    --disable_indel_quals  \
    -BQSR $bqsr \
    -compress 0

if [[ $gatk_jar =~ GenomeAnalysisTKLite ]] ;
then
    java -Xmx6g -jar $gatk_jar \
        -R $gatk_data/hg19/ucsc.hg19.fasta \
        -T UnifiedGenotyper \
        -I recalibrated.bam \
        --dbsnp interval.dbsnp_137.hg19.vcf \
        -o raw.vcf \
        -glm BOTH \
        -L $interval \
        -stand_call_conf 30.0 \
        -stand_emit_conf 30.0 \
        -nct $cpu_cores \
        -rf BadCigar
else 
    #full featured GATK, run HaplotypeCaller
    java -Xmx6g -jar $gatk_jar \
        -R $gatk_data/hg19/ucsc.hg19.fasta \
        -T HaplotypeCaller \
        -I recalibrated.bam \
        --dbsnp interval.dbsnp_137.hg19.vcf \
        -stand_call_conf 30.0 \
        -stand_emit_conf 30.0 \
        -o raw.hc.vcf \
        -L $interval \
        -nct $cpu_cores \
        -rf BadCigar

    java -Xmx6g -jar $gatk_jar \
        -R $gatk_data/hg19/ucsc.hg19.fasta \
        -T VariantAnnotator \
        -I recalibrated.bam \
        --variant raw.hc.vcf \
        -A MappingQualityZero \
        --dbsnp interval.dbsnp_137.hg19.vcf \
        -o raw.vcf \
        -L $interval \
        -nt $cpu_cores \
        -rf BadCigar

fi

./swe emit file raw.vcf
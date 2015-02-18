#!/bin/bash
set -e
set -x
set -o pipefail

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


java -Xmx2g -jar ./MarkDuplicates.jar INPUT=$input  OUTPUT=deduplicated.bam REMOVE_DUPLICATES=true METRICS_FILE=duplication.metrics
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

freebayes -f $gatk_data/hg19/ucsc.hg19.fasta recalibrated.bam > fb.vcf

java -Xmx6g -jar $gatk_jar \
    -R $gatk_data/hg19/ucsc.hg19.fasta \
    -T VariantAnnotator \
    -I recalibrated.bam \
    --variant fb.vcf \
    -A QualByDepth \
    -A MappingQualityRankSumTest \
    -A FisherStrand \
    -A ReadPosRankSumTest \
    -A MappingQualityZero \
    --dbsnp interval.dbsnp_137.hg19.vcf \
    -o raw.vcf \
    -L $interval \
    -nt $cpu_cores \
    -rf BadCigar
	    					    	    							    
./swe emit file raw.vcf
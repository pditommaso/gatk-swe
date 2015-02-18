set -e -x

# GATK reference file, should have trailing slash
export GATK_REFERENCE=s3://gapp-east/gatk-reference/

# Sample exome data
export SAMPLE_DATA=s3://gapp-east/sample

#temporary storage bucket, we recommend setting up lifecycle policy that would erase files after 7 days after their creation time 
export SWE_S3_STORAGE=s3://east-temp-0
# Genome FAI 
export GENOME_FAI=./bin/ucsc.hg19.fasta.fai


VariantCaller=UG
VariantCaller=HC
if [ "$GATK" == "UG" ]
then
    #Gatk Jar file, if ends with Lite.jar we use Unified genotyper, and HaplotypeCaller otherwise
    export GATK_JAR=s3://gapp-east/GenomeAnalysisTKLite.jar
else
    #Requires commercial license from Appistry unless you're in academia
    export GATK_JAR=s3://gapp-east/GenomeAnalysisTK.jar
fi

#CHROMOSOMES=chr22
    
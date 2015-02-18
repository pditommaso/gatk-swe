#!/bin/bash
set -e -x 
#Download reference files and save them to S3
# best way to launch it ksub -s -u general_settings.sh -u download_reference.sh -de 20000 -dm 20000 --wrap="bash download_reference.sh"
#

. general_settings.sh

#es3 testhg19/ucsc.hg19.fasta && es3 test ${GATK_REFERENCE}1000G_phase1.indels.hg19.vcf.gz && touch re

if [ ! -e ./gatk-reference.tar.gz ]
then

    wget -O ./gatk-reference.tar.gz http://gapp-west.s3.amazonaws.com/gatk-reference.tar.gz
    tar xvf ./gatk-reference.tar.gz
    echo Reference files downloaded successfully
fi

es3 sync ./gatk-reference/ ${GATK_REFERENCE}




if [ ! -e ./gcat_set_025_2.fastq.gz ]
then
    wget	-O ./gcat_set_025_1.fastq.gz http://gapp-west.s3.amazonaws.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_1.fastq.gz
    wget	-O ./gcat_set_025_2.fastq.gz http://gapp-west.s3.amazonaws.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025_2.fastq.gz
    wget	-O ./gcat_set_025.bed.gz     http://gapp-west.s3.amazonaws.com/sample_exome/025_Bioplanet_GCAT_30x/gcat_set_025.bed.gz
fi

es3 sync ./gcat_set_025.bed.gz ./gcat_set_025_1.fastq.gz ./gcat_set_025_2.fastq.gz $SAMPLE_DATA/

echo completed successfully

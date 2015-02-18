#!/bin/bash
set -e -x 
#Download reference files and save them to S3
# best way to launch it ksub -s -u general_settings.sh -u download_reference.sh -de 20000 -dm 20000 --wrap="bash download_reference.sh"
#

. general_settings.sh




#for url in ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U0a ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_006_AH81VLADXX/Project_RM8398/Sample_U0a ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131223_D00360_007_BH88WKADXX/Project_RM8398/Sample_U0a ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131223_D00360_008_AH88U0ADXX/Project_RM8398/Sample_U0a ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140313_D00360_0014_AH8GGVADXX/Project_RM8398/Sample_U0a ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140313_D00360_0015_BH9258ADXX/Project_RM8398/Sample_U0a ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140407_D00360_0016_AH948VADXX/Project_RM8398/Sample_U0a ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140407_D00360_0017_BH947YADXX/Project_RM8398/Sample_U0a 
#do

#ksub -e 1000 -de 20000 -dm 20000 --wrap="wget -r $url && es3 sync ./ftp-trace.ncbi.nih.gov/ $SAMPLE_DATA/"
#done



for url in ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U0b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_006_AH81VLADXX/Project_RM8398/Sample_U0b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131223_D00360_007_BH88WKADXX/Project_RM8398/Sample_U0b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/131223_D00360_008_AH88U0ADXX/Project_RM8398/Sample_U0b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140313_D00360_0014_AH8GGVADXX/Project_RM8398/Sample_U0b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140313_D00360_0015_BH9258ADXX/Project_RM8398/Sample_U0b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140407_D00360_0016_AH948VADXX/Project_RM8398/Sample_U0b ftp://ftp-trace.ncbi.nih.gov/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140407_D00360_0017_BH947YADXX/Project_RM8398/Sample_U0b 
do

ksub -e 1000 -de 20000 -dm 20000 --wrap="wget -r $url && es3 sync ./ftp-trace.ncbi.nih.gov/ $SAMPLE_DATA/"
done




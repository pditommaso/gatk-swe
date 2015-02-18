#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
set -e -x
set -o pipefail


. ./general_settings.sh

. $1


export PATH=$PATH:./bin
export SHELL=/bin/bash
export SWE_ENGINE

if [ "$SWE_ENGINE" != "local" ]
then
    #Using cirrus for scheduling and S3 as a backend
    [ "$SWE_ENGINE" == "clusterk" ]
    [ "$GATK_REFERENCE" != "" ] && es3 test $GATK_REFERENCE # check that $GATK_REFERENCE is defined and exists



    #check reference file consistency    
    es3 test ${GATK_REFERENCE}hg19/ucsc.hg19.fasta
    es3 test ${GATK_REFERENCE}1000G_phase1.indels.hg19.vcf.gz

    [ "$SWE_S3_STORAGE" != "" ] && es3 touch $SWE_S3_STORAGE/test # SWE s3 storage is deefined and we can write to it
    [ -e $GENOME_FAI ] # reference genome file 

    [ "$GATK_JAR" != "" ] #gatk jar must be defined
    es3 test $GATK_JAR


else
    # A simple hack that allows local execution, or execution in absense of Cirrus and es3

    [ -e $GATK_JAR ]
    [ -e ${GATK_REFERENCE}hg19/ucsc.hg19.fasta ]
    [ -e ${GATK_REFERENCE}1000G_phase1.indels.hg19.vcf.gz ]
    [ -e $GENOME_FAI ] # reference genome file 

fi


#create bin.tar.gz which will contain most of the executables needed for analysis
[ -e bin.tar.gz ] || tar czvf bin.tar.gz ./bin ./align ./bqsr ./combine ./combine_vcf ./gatk ./split_chr ./split_fastq ./vqsr ./publish


#chromosomes should start with BQSR_CHROMOSOME
[ "$CHROMOSOMES" != "" ] || export CHROMOSOMES="chr22 chr1 chr21 chr20 chr2 chr19 chr18 chr3 chr17 chr4 chr16 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chrX chrY chrM"

BQSR_CHR=$(echo $CHROMOSOMES|cut -f 1 -d ' ' )
echo Will be using $BQSR_CHR for base quality recalibration. Consider putting smaller chromosome at the beggining of the list e.g. chr22


#create new job for restarts and sample management
job_id=$(kjob add -n $NAME-`date "+%Y-%m-%d-%H:%M:%S"`)
export K_JOB_ID=$job_id

#set this job to batch mode, do not actually start any tasks until the entire job is commited
kjob batch $job_id


#oldest jobs get higher priority
priority=$[10000000-$job_id]

export CLUSTERK_QUEUE=default
# -q $CLUSTERK_QUEUE          - override default task queue
# -j $job_id                  - set job ID
# -p $priority                - set priority for all tasks
# -u ./swe                    - attach swe to tasks and provide it in task working directory
# -u bin.tar.gz               - attach bin.tar.gz and extract it to task working directory
# -t ANALYSIS=$ANALYSIS       - set ANALYSIS tag to $ANALYSIS (WGS or exome)
# -c auto:ANALYSIS,STAGE      - set default number of cores to be determined automatically based on previous history
# -e auto:ANALYSIS,STAGE      - set default allocated RAM to be determined automatically based on previous history
# -de auto:ANALYSIS,STAGE,CHR - set default allocated disk space to be determined automatically based on previous history
# -th auto:ANALYSIS,STAGE,CHR - set expected task time to be determined automatically based on previous history
# -m 15000                    - set maximum RAM usage to 15GB
# -dm 30000                   - set maximum scracth space usage to 30GB
export common_params=" -q $CLUSTERK_QUEUE -j $job_id -p $priority -u ./swe -u bin.tar.gz -t ANALYSIS=$ANALYSIS -c auto:ANALYSIS,STAGE -e auto:ANALYSIS,STAGE -m 15000 -dm 30000 -de auto:ANALYSIS,STAGE,CHR -th auto:ANALYSIS,STAGE,CHR -t CHR=NA  "


NAME_PREFIX="$NAME:$ANALYSIS";

#####################################################################################################
# Process input files:
# 1. Split each input file into chunks of 1GB
# 2. Launch alignment job per each chunk
# 3. Each alignment will return one sorted file per chromosome, which will later be combined.
#
# supported input formats:   paired[s3_path_for_read1,s3_path_for_read2]
#                            local[local_path_for_read1,local_path_for_read2]
#
#

#approximate size of input splits in bytes. 500MB 
input_split_size=500000000

#list of alignment job IDs
align_job_ids=""
for input in $INPUT_FASTQ
do
	#check input string format
	#compute number of alignment splits. one split per GB of input data.
	if [[ $input =~ paired\[(.*),(.*)\] ]] ;
	then
		 file1=${BASH_REMATCH[1]}
		 file2=${BASH_REMATCH[2]}
		 file1_size=$(es3 ls $file1 | head -n 1| cut -f 2)

	elif [[ $input =~ local\[(.*),(.*)\] ]] ;
	then
		file1=${BASH_REMATCH[1]}
		file2=${BASH_REMATCH[2]}
		file1_size=$(stat -c%s ${BASH_REMATCH[1]})
	
	else
		echo Not implemented
		false
	fi

	splits=$[$file1_size/$input_split_size+1]

	echo Processing $input, will split it in $splits chunks

	# input.sh:  accepts input fastq, splits in N chunks, saves them as N.fastq.gz

	if [ "$splits" == "1" ]
	then
	#if only one split required - start alignment directly
			align_job_id=$(ksub $common_params  \
							 -v sample_id=SAMPLE \
							 -v paired1=$file1 \
							 -v paired2=$file2 \
							 -t NAME=$NAME_PREFIX:align -t STAGE=align\
							 --wrap="bash align/align.sh")
			align_job_ids="$align_job_ids $align_job_id"

	else
	
		split_job_id=$(ksub $common_params \
						 -v splits=$splits \
						 -v input1=$file1 \
						 -v input2=$file2 \
						 -t NAME=$NAME_PREFIX:interleave -t STAGE=interleave\
						 -c 8 -dm 50000 \
						 --wrap="bash split_fastq/split.sh")
	
	    #align each split, reference them via $split_job_id

		for split in $(seq 1 $splits)
		do
			# align.sh: accepts input split file, produces alignment, split by chromosome
			align_job_id=$(ksub $common_params  \
							 -v sample_id=SAMPLE \
							 -v interleaved=$split_job_id:$split.fastq.gz \
							 -d $split_job_id \
							 -t NAME=$NAME_PREFIX:align -t STAGE=align\
							 --wrap="bash align/align.sh")
			align_job_ids="$align_job_ids $align_job_id"
		done
	fi

done


########### GATK  stage
# 1. combine chromosome files from each alignment job into single BAM file per chromosome
# 2. Use one of the chromosomes for base quality recalibration (chr22)
# 3. For each chromosome, run split_chr to find safe spots where chromosome can be split
# 4. For each resulting interval run GATK HaplotypeCaller and the rest of Best Practices GATK pipeline
#
#

if [ "$ANALYSIS" == "exome" ]
then
    chr_split_size=30000000
else

    #for GATK LITE use larger split size, since Unified Genotyper is much faster
    if [[ $GATK_JAR =~ GenomeAnalysisTKLite ]] 
    then
	chr_split_size=15000000
    else
	chr_split_size=5000000
    fi
    

fi


for chr in $CHROMOSOMES
do
	#get chromosome size and compute optimal number of splits
	echo Processing chromosome $chr
	chr_size=$(grep "^$chr	" $GENOME_FAI |cut -f 2)
	gatk_splits=$[$chr_size/$chr_split_size+1]
	[ "$gatk_splits" != "0" ]

	#create comma separated list of alignment jobs
	align_job_list=$(echo $align_job_ids |tr " " ",")

	#create list of input alignment files for current chromosome
	input_array=""
	for align_job in $align_job_ids
	do
		input_array="$input_array $align_job:$chr.bam"
	done

	#submit comine jobs, and pass list of alignent files as input, and all alignment jobs as prerequsite
	#combine.sh: accepts a list of aligned bam files, produced combined file for a given chromsome
	# output: $combine_job_id:$chr.bam


	combine_job_id=$(ksub $common_params \
							-v chr=$chr \
							-v input="$input_array" \
							-d $align_job_list \
							-c 8 \
							-t NAME=$NAME_PREFIX:combine:$chr -t STAGE=combine -t CHR=$chr \
							--wrap="bash combine/combine.sh" )

	#for one of the chromosomes submit BQSR job to produce Base Quality Recalibration files
	if [ "$chr" == "$BQSR_CHR" ]
	then
		#bqsr.sh: runs Base Quality recalibration on chr22
		# output is $bqsr_job_id:bqsr.grp
		bqsr_job_id=$(ksub $common_params \
						     -v chr=$chr \
						     -v gatk_jar=$GATK_JAR \
						     -v input=$combine_job_id:$chr.bam \
						     -d $combine_job_id \
						     -t NAME=$NAME_PREFIX:bqsr:$chr -t STAGE=bqsr \
						     --wrap="bash bqsr/bqsr.sh")
	fi
	
	# for each combined chromosomes, run split_chr.sh to obtain list of 
	#chr_split.sh: finds genomic locations where it is safe to split a chromosomes
	#               returns list of bam files:  split.$chr.$split_id.bam
	
	chr_split_id=$(ksub  $common_params \
						  -v input=$combine_job_id:$chr.bam \
						  -v splits=$gatk_splits \
						  -v chr=$chr \
						  -d $combine_job_id \
						  -t NAME=$NAME_PREFIX:chr_split:$chr -t STAGE=chr_split  -t CHR=$chr \
						  -c 8 \
						  --wrap="bash split_chr/split_chr.sh")

	#for each chromosome split, start a GATK analysis job
	#assign lower priority to gatk tasks
	gatk_priority=$[$priority-10000]

	gatk_job_ids=""
	for split_id in $(seq 1 $gatk_splits)
	do

		#gatk.sh runs gatk on a sub-interval and applies BQSR
		# output: $gatk_job_id:raw.vcf
		#submit GATKs with lower priority to let chr_splits finish faster
		gatk_job_id=$(ksub $common_params -p $gatk_priority\
							-v input=$chr_split_id:$split_id.bam \
							-v interval=$chr_split_id:$split_id.interval \
							-v bqsr=$bqsr_job_id:bqsr.grp \
							-v gatk_jar=$GATK_JAR \
							-d $chr_split_id,$bqsr_job_id \
							-t NAME=$NAME_PREFIX:gatk:$chr:$split_id -t STAGE=gatk \
							--wrap="bash gatk/gatk.sh")
		gatk_job_ids="$gatk_job_ids $gatk_job_id"
		
	done
	# collect all GATK job ids, and create comma separated list for dependendencies (-d)
	gatk_job_list=$(echo $gatk_job_ids |tr " " ",")
	
	input_array=""
	for gatk_job in $gatk_job_ids
	do
	    input_array="$input_array $gatk_job:raw.vcf"
	done


	combine_vcf_job_id=$(ksub  $common_params \
			    -v input="$input_array" \
			    -d $gatk_job_list \
			    -t NAME=$NAME_PREFIX:combine_vcf:$chr -t STAGE=combine_vcf \
			    --wrap="bash combine_vcf/combine_vcf.sh" )
	
	combine_vcf_job_ids="$combine_vcf_job_ids $combine_vcf_job_id"


done

# collect all GATK job ids, and create comma separated list for dependendencies (-d)
combine_vcf_job_list=$(echo $combine_vcf_job_ids |tr " " ",")

input_array=""
for combine_vcf_job in $combine_vcf_job_ids
do
    input_array="$input_array $combine_vcf_job:raw.vcf"
done

#submit a job that combines all sub-region vcf files, into one sorted VCF
#combine_vcf.sh concatenate sub-region VCFs into final VCF

combine_vcf_job_id=$(ksub  $common_params \
			    -v input="$input_array" \
			    -d $combine_vcf_job_list \
			    -t NAME=$NAME_PREFIX:combine_vcf -t STAGE=combine_vcf \
			    --wrap="bash combine_vcf/combine_vcf.sh" )

#submit a job that run 
#run variant quality recalibration
# output: $vqsr_job_id:recalibrated.filtered.vcf.gz


vqsr_job_id=$(ksub $common_params \
		    -v gatk_jar=$GATK_JAR \
		    -v input=$combine_vcf_job_id:raw.vcf \
		    -v ANALYSIS=$ANALYSIS \
		    -d $combine_vcf_job_id \
		    -t NAME=$NAME_PREFIX:vqsr -t STAGE=vqsr \
		    --wrap="bash vqsr/vqsr.sh")


# save results to S3 and make them publicly accessible over HTTP
publish_job_id=$(ksub $common_params \
		    -v input=$vqsr_job_id:recalibrated.filtered.vcf.gz \
		    -v path="$SAMPLE_DATA/$NAME" \
		    -d $vqsr_job_id \
		    -t NAME=$NAME_PREFIX:publish -t STAGE=publish \
		    --wrap="bash publish/publish.sh")


echo Commiting all tasks to job $job_id
#submit all tasks to the server
kjob commit $job_id

#wait for all jobs to finish

#kwait $vqsr_job_id
# output files can be fetched from S3 using command below. Alternatively you can add another stage that will publish them to another S3 bucket
#./swe fetch $vqsr_job_id:recalibrated.filtered.vcf.gz

exit 0



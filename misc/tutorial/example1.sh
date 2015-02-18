#!/bin/bash

set -e # abort on any errors
set -x # display commands being run

# This sample script will submit two tasks - one that compresses file and another one that decompresses file.

#S3 buckes for storing temporary files emitted by tasks
export SWE_S3_STORAGE=s3://gapp-east-temp
export SWE_ENGINE=clusterk


#create a sample input file
head -c 1024 /dev/urandom > input.txt

#upload input file to S3, return S3 path 
input_s3_name=$(./swe store input.txt)

echo input.txt is stored as $input_s3_name


#submit compression task, in receives "input" variable and emits output.gz
#Task will be submitted to "default" queue
# -u swe                   - upload swe script, it will be available in task working directory
# -v input=$input_s3_name  - set environmental variable $input 
# --wrap                   - execute script using bash -c
#ksub returs numerical task ID on success
compress_task_id=$(ksub \
			-u swe \
			-v     input=$input_s3_name \
			--wrap="local_name=\$(./swe fetch \$input) && \
				cat \$local_name | gzip -c > output.gz &&
				./swe emit file output.gz " 
		  )

echo Compression task $compress_task_id submitted 		
#check that previous task completed successfully
#if task failed we can view the log file via the the web interface or klog $compress_task_id
#kwait $compress_task_id && echo $compress_task_id finished successfully

#we can now submit decompression tasks that will depend on compress_task_id and use it's output as input
# -d $compress_task_id   - wait for $compress_task_id to finish
decompress_task_id=$(ksub \
			    -u swe \
			    -d $compress_task_id \
			    -v input=$compress_task_id:output.gz \
			    --wrap="local_name=\$(./swe fetch \$input) && \
				    zcat \$local_name > output.txt &&
				    ./swe emit file output.txt " 
		  )
echo deompression task $decompress_task_id submitted 		

#wait for decompression task to finish
kwait $decompress_task_id 

echo $decompress_task_id finished successfully

#download results locally 
output_local_name=$(./swe fetch $decompress_task_id:output.txt)

echo decompression task $decompress_task_id output is saved as $output_local_name

#check that original file matches
diff input.txt $output_local_name && echo files match, test is successfull


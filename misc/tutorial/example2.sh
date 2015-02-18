#!/bin/bash

set -e # abort on any errors
set -x # display commands being run

# This sample script will submit two tasks - one that compresses file and another one that decompresses file.

#S3 buckes for storing temporary files emitted by tasks
export SWE_S3_STORAGE=s3://gapp-east-temp
export SWE_ENGINE=clusterk


#input file
cat > lowercase.txt << END
file
with
lines
in
lower
case
that
will
be
converted
to
upper
case
END



#upload input file to S3, return S3 path 
lowercase_s3_path=$(./swe store lowercase.txt)

echo lowercase is stored as $lowercase_s3_path



#submit split task 
#Task will be submitted to "default" queue
# -u swe                   - upload swe script, it will be available in task working directory
# -u split.sh		   - attach split.sh
# -v input=$input_s3_name  - set environmental variable $input 
# --wrap                   - execute script using bash -c
#ksub returs numerical task ID on success
split_task_id=$(ksub \
			-u swe \
			-u split.sh \
			-v lines=2 \
			-v input=$lowercase_s3_path \
			--wrap="bash split.sh split"
		  )


echo Split task $split_task_id submitted 		
#check that previous task completed successfully
#if task failed we can view the log file via the the web interface or klog $split_task_id
kwait $split_task_id 

echo $split_task_id finished successfully

file_list=$(./swe fetch $split_task_id:file.list)
echo split names are: $(cat $file_list)

for split_name in `cat $file_list`
do
    echo submitting a task for $split_name

    upcase_task_id=$(ksub \
			    -u swe \
			    -u split.sh \
			    -v input=$split_task_id:$split_name \
			    --wrap="bash split.sh upcase"
			)    
    upcase_ids="$upcase_ids $upcase_task_id"
    upcase_outputs="$upcase_outputs $upcase_task_id:output.txt"
done


#replace spaces with commas for lowercase_ids to use in dependency list
upcase_ids=$(echo $upcase_ids |tr " " ",")

echo upcase task IDS: $upcase_ids
echo upcase outputs IDS: $upcase_outputs


#submit combine task that will put it all together

combine_task_id=$(ksub \
			    -u swe \
			    -u split.sh \
			    -d $upcase_ids \
			    -v input="$upcase_outputs" \
			    --wrap="bash split.sh combine"
		  )

echo Combine task $combine_task_id submitted 		

kwait $combine_task_id

echo combine task_id $combine_task_id finished successfully

#download results locally 
echo output:
cat $(./swe fetch $combine_task_id:output.txt)


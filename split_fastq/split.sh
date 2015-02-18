#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
#
# Split input fastq files into N chunks, output interleaved fastq
#
# Input:
#       $input1 - first  read file
#       $input2 - second read file
#       $splits - number of splits to produce
#
# Output: 1.fastq.gz, 2.fastq.gz, ... , N.fastq.gz
#      
set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)


[ "$splits" != "" ]
[ "$input1" != "" ] && input1=$(./swe fetch $input1)
[ "$input2" != "" ] && input2=$(./swe fetch $input2)


./split_fastq/split_input_fastq.pl --input paired[$input1,$input2] --splits $splits
# split_input_fastq.pl will output 1.fastq.gz 2.fastq.gz .... N.fastq.gz  - interleaved FASTQ files

for i in $(seq 1 $splits)
do
    ./swe emit file $i.fastq.gz
done

exit 0
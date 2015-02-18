#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
#
# Combine multiple bams for a given chromosome into a single bam
#
# Input : $input - space separated list of input files
# Output: $chr.bam, $chr.bam.bai
#

set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)


[ "$input" != "" ]
[ "$chr" != "" ]


cpu_cores=32


# compute number of input files
words=( $input )

if [ "${#words[@]}" == "1" ]
then
    #single input file, nothing to merge
    cp $(./swe fetch $input) $chr.bam
    cp $(./swe fetch $input.bai ) $chr.bam.bai
else
    # multiple input files

    local_bams=""
    for bam in $input
    do
	local_bam=$(./swe fetch $bam)
	local_bams="$local_bams $local_bam"
    done

    samtools merge -@ $cpu_cores $chr.bam $local_bams
    samtools index $chr.bam
fi

    ./swe emit file $chr.bam
    ./swe emit file $chr.bam.bai

exit 0

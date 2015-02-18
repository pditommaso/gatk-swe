#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
#
# Combine multiple VCFs for the same sample into a single VCF
# Input: $input space separated list of input files (1234:chr2.bam 1235:chr2.bam 1236:chr2.bam) 
# Output: raw.vcf 
#
set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)


[ "$input" != "" ]

[ ! -e raw.vcf.tmp ] || rm raw.vcf.tmp

#run_in_parallel=1

# cat all input VCFs into raw.vcf.tmp
if [ "$run_in_parallel" != "" ]
then
    # save input files into input.XXXX.vcf in parallel
    parallel -j 10 -i bash -c "mv \`./swe fetch {}\` input.\$\$.vcf " -- $input

    cat input.*.vcf > raw.vcf.tmp
else
    for vcf in $input
    do
	cat $(./swe fetch $vcf) >> raw.vcf.tmp
    done
fi


#obtain header from the first VCF in the list
set -- $input
first_vcf=$1

cat $(./swe fetch $first_vcf) |grep -P '^\#' >header.vcf


cat header.vcf > raw.vcf
cat raw.vcf.tmp | grep -vP "^\#" | vcf-sort -c >> raw.vcf

./swe emit file raw.vcf

exit 0


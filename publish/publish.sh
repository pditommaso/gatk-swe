#!/bin/bash
# Copyright (C) 2013-2015 Clusterk Inc
#
# Publish output file to S3. This wont work unless Cirrus CLI (es3) tools are installed.
# Input:
#        $input  - input file
#        $path   - S3 path where this file should be published

set -e
set -x
set -o pipefail
# Abort execution if any dependencies failed
[ "$CLUSTERK_FAILED_DEPS" == "" ] || ( echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting && exit 1)

input=$(./swe fetch $input)

es3 sync $input $path/
es3 publish $path

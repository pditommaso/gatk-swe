#!/bin/bash
#this script will install necessary ubuntu packages for gatk-pipeline

set -e -x
apt-get update

apt="apt-get install -qq "

for dep in "basename" "realpath" "dirname" "pigz" "zcat" "aria2" "wget" "curl"
do
which $dep >/dev/null 2>/dev/null || $apt $dep
done

which java || $apt openjdk-7-jre

#needed for parallel
$apt moreutils
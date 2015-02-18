#!/bin/bash

set -e -x















exit
#!/usr/bin/perl
#
# Applies BED filter to a VCF file
$filter=$ARGV[0];
die "Filter file is not specified" unless -e $filter;
die "Incorrect filter format" unless $filter =~/\.(bed|BED)/;


#hash based coarse filter
open(FL,"$filter")||die;
while (<FL>)
{
    @a=split("\t");
    $a[1]=int($a[1]/10000);
    $a[2]=int($a[2]/10000);

    for ($i=$a[1]-1;$i<=$a[2]+1;$i++)
    {
	$OK{"$a[0]-$i"}=1;
    }
}
$tmp=int(rand()*1000000)."bed";
$tmp="out.bed";

open(OUT,">$tmp");

while (<STDIN>)
{
    if (/^\#/)
    {$header.=$_;next;}
    @a=split("\t");

    $red=int($a[1]/10000);

    next unless $OK{"$a[0]-$red"};
    print OUT "$a[0]\t$a[1]\t$a[1]\t$_";
}

close OUT;

print $header;
open(FL,"bedtools intersect -u -a $tmp -b $filter|") ||die;
while (<FL>)
{
    @a=split("\t");
    for ($i=3;$i<scalar(@a);$i++)
    {
	print "\t" unless $i==3;
	print $a[$i];
    }

}

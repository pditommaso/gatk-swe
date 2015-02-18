#!/usr/bin/perl

# expects samtools view to be piped into it
#open(FL,"samtools mpileup -l test.bed -f ~/gatk-data/hg19/ucsc.hg19.fasta chr10.bam|")||die;
use strict;

use Getopt::Long;

my $bam; my $ref; my $region;
GetOptions ('bam=s' => \$bam, "reference=s"=>\$ref, "region=s"=>\$region);

die unless $region;

my $gap=800;
my $max_splits=2;
my $old_pos;
my $old_chr;
my $safe;
my $last_safe;
my $midpoint;
my @breakpoints;


die unless -e $bam; die unless -e $ref;
open(PILEUP, "samtools view $bam -1 -b -\@ 5 $region | samtools mpileup -f $ref -|")||die;
my $tm=time;
while (<PILEUP>)
{
    my ($chr,$pos,$ref,$depth,$calls)=split("\t");


    if ($tm) {print STDERR "samtools seek took ".(time-$tm)." seconds for region $region\n";$tm=0;}
    $safe=1;
    
    if ($old_pos and ($pos-$old_pos>$gap)) # large coverage gap
    {
	
	$midpoint=int( ($pos + $old_pos) /2 );
	push(@breakpoints, $midpoint);
	print STDERR "Coverage gap at $midpoint\n";
	$safe=0;
    }

    $old_pos=$pos;
    if (scalar(@breakpoints)>=$max_splits)
    {	
	foreach (@breakpoints) {print "$_\n";}
	exit;
    }

    unless ($chr eq $old_chr)
    {
	die "Currently only one chromosome is supported" if $old_chr;
	$safe=0;
	$old_chr=$chr;
    }

    
    # reference must be uppercase, meaning that it's not a repeat region
    $safe=0 unless $ref=~/[ACGT]/;
    
    # call string must not have  more than 1 nonreference bases, if depth is high enough

    if ($depth<15)
    {
	$safe=0 if $calls =~ /[ACGTNacgtn]/;
	# call string must not have any insertions or deletions
	$safe=0 if $calls=~/(\+|\-)/;

    }
    else
    {
	#for higher depth up to 2 mismatches are allowed
	$calls =~ s/[^ACGTNacgtn\+\-]//g;
	$safe=0 if length($calls)>2;
    }

    
    if ($safe)
    {
	if ($last_safe)
	{
	    if ($pos-$last_safe > $gap)
	    {
		$midpoint=int( ($pos + $last_safe) /2 );
		push(@breakpoints, $midpoint);
		$last_safe=$pos;
		print STDERR "we can break at $midpoint $last_safe\n";
	    }
	}
	else
	{
	    $last_safe=$pos;;
	}

    }
    else
    {
	$last_safe=0;

    }

}


foreach (@breakpoints) {print "$_\n";}
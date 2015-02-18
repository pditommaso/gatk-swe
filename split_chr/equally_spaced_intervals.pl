#!/usr/bin/perl
#generate simple 1mb.intervals
use strict;

use Getopt::Long;

my $bam; my $ref;
my $block_length;
my $blocks;
my $ref_chr;
GetOptions ('bam=s' => \$bam, "blocks=i"=>\$blocks,"chr=s"=>\$ref_chr);



die unless -e $bam;
die unless -e "$bam.bai";

my $region;
my $old_chr;
my $max_coord;
open(FL,"samtools idxstats $bam|")||die;
while (<FL>)
{
    my($chr,$len,$reads)=split(/\s+/);
    next unless $chr eq $ref_chr;
    die unless $reads;

    $block_length=int($len/$blocks);
#    die "Multiple chromosomes per bam file!\n" if ($old_chr and ($old_chr ne $chr));
    $old_chr=$chr;

    
    for (my $i=1;$i<=$len;$i+=$block_length)
    {
	$max_coord=$i+$block_length;
	$max_coord=$len if $max_coord>$len;

	$region++;
	print "$chr:$i-$max_coord\n";
	
    }

}


#!/usr/bin/perl

use Getopt::Long;

my $brkp;
GetOptions ('bam=s' => \$bam, "breakpoints=s"=>\$brkp,"splits=i"=>\$splits);

die unless -e "$bam.bai";

die unless $splits;

my $chr; my $len; my $reads;
my $chr_length=0; my $chromosome;
open(FL,"samtools idxstats $bam|")||die;
while (<FL>)
{
    ($chr,$len,$reads)=split(/\s+/);
    next unless $reads;
    

    die "Multiple chromosomes per bam file!\n" if ($old_chr and ($old_chr ne $chr));
    $old_chr=$chr;

    $chr_length=$len;
    $chromosome=$chr;
}
die "can't establish chr length" unless $chr_length;

print STDERR "chromosome length $chr_length\n";

open(FL,$brkp) || die;
@BREAKPOINTS=<FL>;
chomp @BREAKPOINTS;

if ($chr_length>1000000)
{
    warn "No breakpoints found for sufficiently large chromosome!" unless scalar @BREAKPOINTS;
}

$block=int($chr_length/$splits);
my $start=1;
my $cnt=0;
my $pos;
foreach (sort {$a <=> $b} @BREAKPOINTS)
{
    $pos=$_;
    if ($chr_length-$pos<$block)
    {
	print "$chromosome:$start-$chr_length\n";
	$cnt++;
	last;
    }
    if ($pos-$start> $block) 
    {
	print "$chromosome:$start-$pos\n";
	$cnt++;
	$start=$pos+1;
    }
}

unless($chr_length-$pos<$block)
{
	print "$chromosome:$start-$chr_length\n";
	$cnt++;
}
die "$cnt $splits" if $cnt > $splits;
for ($i=0;$i<$splits-$cnt;$i++){print "$chromosome:1-1\n";}

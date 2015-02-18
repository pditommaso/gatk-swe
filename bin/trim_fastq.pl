#!/usr/bin/perl
# Trim fastq reads to certain length
# e,g,   cat input.fastq | trim_fastq.pl 150 > output.fastq
use strict;

my $trim_length=$ARGV[0] || die "please specify nonzero trim length";


while (<STDIN>)
{
    print ; #read name goes unchanged
    
    #tream read
    my $a=<STDIN>;
    chomp($a);
    print substr($a,0,$trim_length)."\n";

    #plus goes unchanged
    my $a=<STDIN>;
    print $a;

    #trim quality string
    my $a=<STDIN>;
    chomp($a);
    print substr($a,0,$trim_length)."\n";
    
}
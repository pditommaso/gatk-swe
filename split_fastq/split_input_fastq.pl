#!/usr/bin/perl

use strict;
use File::Basename;
use Digest::MD5 qw(md5 md5_hex md5_base64);


use Getopt::Long;

# TODO: check that filenames are valid for security purposes

# input is either:
#   interleaved[file_on_s3.fastq];
# or
#  paired[pair_1.fastq,pair_2.fastq]

die "pigz not installed " if system("pigz --help  >/dev/null 2>/dev/null");
die "zcat not installed " if system("zcat --help  >/dev/null 2>/dev/null");

my $i;
my $input; 
my $splits;
GetOptions ('input=s' => \$input, "splits=i"=>\$splits);

die "Incorrect input format $input " unless $input=~/^(paired)\[(\S+)\]$/;
my ($type, $files)=($1,$2);

#compute pseudo-random prefix so that we can we use multiple file names

my $pairs_written=0;
my $bases_written=0;

if ($type eq "paired")
{
    print STDERR "Determined file type to be paired FASTQ ($type, $files)\n";
    my @input_files=split(/\,/,$files);
    die "we expect two input files for paired  fastq: $files" unless scalar (@input_files)==2;

    my $fn_prefix;
    #die "$input_files[0] | $input_files[1]";
    foreach my $input_file (@input_files)
    {
	die "$input_file must end in either .fastq or .fastq.gz" unless $input_file=~/.*([^\/]+)(\_1|\_2|\_R1|\_R2|\_r1|\_r2)\S*(\.fastq|\.fastq\.gz|\.fq|\.fq.gz)$/;
    }


    my $local_fn= $input_files[0];
    die unless -e $local_fn;
    if ($local_fn =~ /\.gz$/ )
    {
	open(L,"zcat $local_fn|") ||die;
    }
    else
    {
	open(L,"$local_fn") ||die;
    }

    my $local_fn=$input_files[1];
    die unless -e $local_fn;
    if ($local_fn =~ /\.gz$/ )
    {
	open(R,"zcat $local_fn|") ||die;
    }
    else
    {
	open(R,"$local_fn") ||die;
    }

    die "number of splits is not defined" unless $splits;

    my @FH;
    for ($i=1;$i<=$splits;$i++)
    {
	local *O;
	open(O,"|pigz > $i.fastq.gz.tmp")||die;
	push (@FH,*O)
    }
    my $read_cnt=0;
    my $fh;
    my $chunk=0;
    while (<L>)
    {
	my $name_L=$_;
	my $seq_L=<L>;
	my $plus_L=<L>;
	my $qual_L=<L>;
 
	my $name_R=<R>;
	my $seq_R=<R>;
	my $plus_R=<R>;
	my $qual_R=<R>;


	$fh=$FH[$read_cnt % $splits];

	print $fh $name_L.$seq_L.$plus_L.$qual_L.$name_R.$seq_R.$plus_R.$qual_R;

	$read_cnt++;
    }

    foreach (@FH) {close $_;}

    for ($i=1;$i<=$splits;$i++)
    {
	die if system("mv $i.fastq.gz.tmp $i.fastq.gz");
    }
}


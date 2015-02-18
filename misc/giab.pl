#!/usr/bin/perl


open(FL,"es3  lsr s3://gapp-east/sample/giab/|");
while (<FL>)
{
    chomp;
    #s3://gapp-east/sample/giab/ftp/technical/NIST_NA12878_HG001_HiSeq_300x/140407_D00360_0017_BH947YADXX/Project_RM8398/Sample_U0a/U0a_CGATGT_L002_R2_004.fastq.gz
    next unless /(\S+_)(R1)(_\d+\.fastq.gz)/;
    my ($pre,$r,$post)=($1,$2,$3);
#    print "$_\n";
    die if $K{$_};
    $K{$_}=1;
    
    push(@input,"paired[$1"."R1"."$3,$1"."R2"."$3]");


}
die unless scalar (@input);

while (<DATA>)
{
    print;
}
print "INPUT_FASTQ=\"".join(" ",@input)."\"\n";


__DATA__
NAME=GIAB_NA12878
ANALYSIS=WGS-giab



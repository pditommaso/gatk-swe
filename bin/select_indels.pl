#!/usr/bin/perl


while (<STDIN>)
{
    if (/^\#/) {print;next;}
    @a=split("\t");
    
    $snp=1;
    foreach ($a[3],split(",",$a[4]))
    {
	$snp=0 unless length($_)==1;
    }

    print unless $snp;
}
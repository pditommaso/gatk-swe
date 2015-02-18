#!/usr/bin/perl

while (<STDIN>)
{
    if (/^\#/) { print; next;}

    @a=split("\t");
    print if $a[5]>30;
}
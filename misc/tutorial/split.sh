#!/bin/bash
# This file defines 3 functions: split, upcase, combine
# Split  takes  "input" and "lines" and produces file.list and some number of splits
# upcase takes  "input" and produces "output.txt"
# combined takes a list of files as "input" and produces output.txt
set -e -x



if [ "$1" == "split" ]
then

    #check input parameters
    [ "$input" != "" ] # input must be defined
    [ "$lines" != "" ] # lines must be defined
    
    #save input locally
    input=$(./swe fetch $input)
    
    #split input in chunks with $lines lines each
    split -l $lines $input output.

    #create list of resulting splits
    ls output.* > file.list
    
    #save file.list and splits
    ./swe emit file file.list output.*

exit 0
fi


if [ "$1" == "upcase" ]
then

    #check input parameters
    [ "$input" != "" ] # input must be defined
    #save input locally
    input=$(./swe fetch $input)
    
    tr '[:lower:]' '[:upper:]'  < $input > output.txt

    ./swe emit file output.txt


exit 0
fi


if [ "$1" == "combine" ]
then

    #check input parameters
    [ "$input" != "" ] # input must be defined


    #loop over input files and save them as output.txt
    for input_file in $input
    do
	cat `./swe fetch $input_file` >> output.txt
    done

    ./swe emit file output.txt

exit 0
fi




echo Incorrect command
exit 1
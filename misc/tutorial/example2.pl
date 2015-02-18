#!/usr/bin/perl


#S3 buckes for storing temporary files emitted by tasks
$ENV{SWE_S3_STORAGE}="s3://gapp-east-temp";
$ENV{SWE_ENGINE}="clusterk";


#create lowercase.txt
open(OUT,">lowercase.txt");
while (<DATA>)
{
    print OUT;
}


#upload input file to S3, return S3 path 
$lowercase_s3_path=run("./swe store lowercase.txt");

print "lowercase is stored as $lowercase_s3_path\n";

$split_task_id=run("ksub \
			-u swe \
			-u split.sh \
			-v lines=2 \
			-v input=$lowercase_s3_path \
			--wrap=\"bash split.sh split\" ");

print "Split task $split_task_id submitted\n";

run("kwait $split_task_id");

print "$split_task_id finished successfully\n";


$file_list=run("./swe fetch $split_task_id:file.list");

open(FL,"$file_list") || die;
while(<FL>)
{
    chomp;
    $split_name=$_;

    print "submitting a task for $split_name\n";

    $upcase_task_id=run("ksub 
			    -u swe 
			    -u split.sh 
			    -v input=$split_task_id:$split_name 
			    --wrap=\"bash split.sh upcase\" ");
    push (@upcase_ids, $upcase_task_id);
}


print "upcase task IDS: ".join(" ",@upcase_ids)."\n";


$combine_task_id=run("ksub 
			    -u swe 
			    -u split.sh
			    -d ".join(",",@upcase_ids)." 
			    -v input=\"".join(" ", map {"$_:output.txt"} @upcase_ids)."\" 
			    --wrap=\"bash split.sh combine\" ");

print "Combine task $combine_task_id submitted\n";


run("kwait $combine_task_id");
print "combine task_id $combine_task_id finished successfully\n";

#download results locally 
$local_result=run("./swe fetch $combine_task_id:output.txt");
system("cat $local_result");

exit;


sub run
{

    my $cmd=shift;
    $cmd=~s/\n/ /g;
    my $out=`$cmd`;
    
    die "$cmd exited with code $? " if $?;
    print STDERR "running $cmd\n";
    chomp $out;

    
    return $out;
}

__DATA__
file
with
lines
in
lower
case
that
will
be
converted
to
upper
case
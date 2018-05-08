#! /usr/bin/perl
use strict;
use warnings;

#open the fastq file
my $FileName = $ARGV[0];
open FASTQ, "$FileName" or die("Cannot open the file!\n");

#define some variables
my $criteria = 0;
my $line;
my $counter = 1;
my $subcounter;
my $subsub;
my @temp;

#two dimensional array to store sequence and the score
my @reads;

#the variable which stores the result 
my @AverageQ;
my @AverageQRanged;
my @Mean;
my @Median;
my @GC;
my @GCRanged;

#the constant to calculate the score
my $bias = 33;

#read the second and forth line of each read
while($line = <FASTQ>){
    chop $line;
    if($counter%4 == 2){
        $reads[($counter-2)/4][0] = $line;
    }
    if($counter%4 == 0){
        $reads[($counter-4)/4][1] = $line;
    }
    $counter++;
}


#calculate the mean and median
for($counter=0;$counter<length($reads[0][1]);$counter++){
    undef @temp;

    $Mean[$counter] = 0;
    $Median[$counter] = 0;
    $temp[0] = 0;

    for($subcounter=0;$subcounter<(scalar @reads);$subcounter++){
        $Mean[$counter] += ord(substr($reads[$subcounter][1],$counter,1));
        for($subsub=0;$subsub<(scalar @temp);$subsub++){
            if(ord(substr($reads[$subcounter][1],$counter,1))>$temp[$subsub]){
                splice(@temp,$subsub,0,ord(substr($reads[$subcounter][1],$counter,1)));
                last;
            }
        }
    }
    $Mean[$counter] = $Mean[$counter]/scalar(@reads) - $bias;

    if(scalar(@reads)%2==1){
        $Median[$counter] = $temp[($#temp-1)/2] - $bias;
    }else{
        $Median[$counter] = ($temp[($#temp-2)/2]+$temp[$#temp/2])/2 - $bias;
    }
}

#calculate average scores and GC content across all reads
for($counter=0;$counter<(scalar @reads);$counter++){
    $AverageQ[$counter] = 0;
    $GC[$counter] = 0;

    for($subcounter=0;$subcounter<length($reads[$counter][1]);$subcounter++){
        $AverageQ[$counter] += ord(substr($reads[$counter][1],$subcounter,1));
        if((substr($reads[$counter][0],$subcounter,1) eq "G")||(substr($reads[$counter][0],$subcounter,1) eq "C")){
            $GC[$counter] ++;
        }
    }

    $AverageQ[$counter] = $AverageQ[$counter]/length($reads[$counter][1]);
    $GC[$counter] = $GC[$counter]/length($reads[$counter][0]);
}

my $SCOREinterval = 5;
my $SCORErange = 40;

#calculate the number of sequence within the specified score range
for($counter=0;$counter<$SCORErange/$SCOREinterval;$counter++){
    $AverageQRanged[$counter] = 0;
    for($subcounter=0;$subcounter<(scalar @AverageQ);$subcounter++){
        if(int(($AverageQ[$subcounter]-$bias)/$SCOREinterval)==$counter){
            $AverageQRanged[$counter] ++;
        }
        if(int(($AverageQ[$subcounter]-$bias)/$SCOREinterval)>$counter&&$counter==$SCORErange/$SCOREinterval-1){
            $AverageQRanged[$counter] ++;
        }
    }
}

my $GCinterval = 10;
my $GCrange = 100;

#calculate the number of sequence within the specified GC content range
for($counter=0;$counter<$GCrange/$GCinterval;$counter++){
    $GCRanged[$counter] = 0;
    for($subcounter=0;$subcounter<(scalar @GC);$subcounter++){
        if(int($GCrange*$GC[$subcounter]/$GCinterval)==$counter){
            $GCRanged[$counter] ++;
        }
    }
}

#print the result
print "There are ".length($reads[0][1])." bases in one read.\nFor each base position, the mean quality is: \n";
for($counter=0;$counter<length($reads[0][1]);$counter++){
    $Mean[$counter] = sprintf("%.2f",$Mean[$counter]); 
    print $Mean[$counter]."\t";
    if(($counter+1)%20==0){
        print "\n";
    }
}

print "\n\nFor each base position, the median quality is: \n";
for($counter=0;$counter<length($reads[0][1]);$counter++){
    $Median[$counter] = sprintf("%.2f",$Median[$counter]); 
    print $Median[$counter]."\t";
    if(($counter+1)%20==0){
        print "\n";
    }
}

print "\n\nQuality score distribution over all sequences: \n";
for($counter=0;$counter<(scalar @AverageQRanged);$counter++){
    print join("",(($counter*$SCOREinterval),"\t~ ",(($counter+1)*$SCOREinterval),"\t:\t",$AverageQRanged[$counter],"\n"));
}

print "\n\nGC distribution over all sequence: \n";
for($counter=0;$counter<(scalar @GCRanged);$counter++){
    print join("",(($counter*$GCinterval)."%\t~ ".(($counter+1)*$GCinterval)."%\t:\t".$GCRanged[$counter]."\n"));
}

print "\n";
close FASTQ;

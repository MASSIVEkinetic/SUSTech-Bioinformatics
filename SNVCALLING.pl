#! /usr/bin/perl
use strict;
use warnings;

#open the fastq file
my $REFName = $ARGV[0];
my $SAMName = $ARGV[1];
my $VCFName = $ARGV[2];
open SAM, "$SAMName" or die("Fail to open file! File name: ".$SAMName." !\n");
open REF, "$REFName" or die("Fail to open file! File name: ".$REFName." !\n");
open VCF, ">$VCFName" or die("Cannot create file ".$VCFName." !\n");

print VCF sprintf("%-10s%-10s%-10s%-10s%-10s\n","#CHROM","POS","ID","REF","ALT");

#define some variables
my $criteria = 0;
my $line;
my $counter;
my $subcounter;
my $subsub;
my @temp;

my %SNV;
my $InfoColNum = 10;
my $CHROM = "Ecoli";

my $REFSeq = "";

while($line=<REF>){
    chomp $line;
    next if substr($line,0,1) eq ">";
    $REFSeq = $REFSeq.$line;
}

#variables related to alignment
my @matrix;
my $x;
my $y;
my $top_bottom = 0;
my $up_down = 0;
my $left_right =0;
my $alignedX = " ";
my $alignedY = " ";
my $gapp = -2;
my $match = 2;
my $mismatch = -1;
my $GOP = 0;
my $LOQ;
my $LOD;
my $threshold = 1;
#E-value is not required.

while($line=<SAM>){
    undef @temp;

    chop $line;
    next if substr($line,0,1) eq "@";

    @temp = split(/\s+/,$line);

    $LOQ = length($temp[9]);
    $LOD = $LOQ/$threshold;
    $matrix[0][0] = 0;
    for($x=0;$x<$LOQ;$x++){
        $matrix[$x+1][0] = ($x+1)*(-2);
    }
    for($y=0;$y<$LOQ;$y++){
        $matrix[0][$y+1] = ($y+1)*(-2);
    }

    my @chainY =  split(//,substr($REFSeq,$temp[3]-1,$LOD));
    my @chainX =  split(//,$temp[9]);

    my $counter;
    my @direction;
    #initialization
    for ($counter = 0; $counter <= scalar(@chainX); $counter++){
        push(@direction,2);
    }

    my $alignedX;
    my $alignedY;
    
    $x = 1;
    $y = 1;

    #assign values to the matrix
    for my $Yaa (@chainY){
        for my $Xaa (@chainX){
            if ($Yaa eq $Xaa){
                $top_bottom = $matrix[$x-1][$y-1] + $match;
            }else{
                $top_bottom = $matrix[$x-1][$y-1] + $mismatch;
            }
            $left_right = $matrix[$x-1][$y] + $gapp + $GOP;
            $left_right = $left_right - $GOP if $direction[$x-1] == 1;
            $up_down = $matrix[$x][$y-1] + $gapp + $GOP;
            $up_down = $up_down - $GOP if $direction[$x] == 0;
            my @result = max($up_down,$left_right,$top_bottom);
            $matrix[$x][$y] = $result[0];
            $direction[$x] = $result[1];
            $x++;
        }
        $y++;
        $x=1;
    }

    #calculate the number of total matches, and find the aligned chains
    $y = scalar(@chainY);
    $x = scalar(@chainX);
    do{
        if(($matrix[$x-1][$y-1]>=$matrix[$x][$y-1] && $matrix[$x-1][$y-1]>=$matrix[$x-1][$y]) || ($chainY[$y-1] eq $chainX[$x-1])){
            if($chainY[$y-1] ne $chainX[$x-1]){
                $SNV{$temp[3]-1+$x} = $chainY[$y-1].$chainX[$x-1];
            } 
            $alignedX = $alignedX.$chainX[$x-1];
            $alignedY = $alignedY.$chainY[$y-1];
            $y--;
            $x--;
        }elsif($matrix[$x-1][$y]>$matrix[$x-1][$y-1] && $matrix[$x-1][$y]>$matrix[$x][$y-1]){
            if($x==scalar(@chainX)){
                $y--;
            }else{
                $SNV{$temp[3]-1+$x} = $chainY[$y-1]."-";
                $alignedY = $alignedY.$chainY[$y-1];
                $alignedX = $alignedX."-";
                $y--;
            }
        }else{
            $SNV{$temp[3]-1+$x} = "-".$chainX[$x-1];
            $alignedY = $alignedY."-";
            $alignedX = $alignedX.$chainX[$x-1];
            $x--;
        }
    }while($x>0 && $y>0);

}

for my $key (keys %SNV){
    print VCF sprintf("%-10s%-10s%-10s%-10s%-10s\n",$CHROM,$key,".",substr($SNV{$key},0,1),substr($SNV{$key},1,1)); 
}

close SAM;
close REF;
close VCF;

sub max{
    my $max = pop(@_);
    my $direction = 0;
    for (my $counter = 0 ; $counter < scalar @_ ; $counter++) {
       if ($_[$counter]>$max){
            $max = $_[$counter];
            $direction = $counter+1;
       }
    }
    return $max,$direction;
}

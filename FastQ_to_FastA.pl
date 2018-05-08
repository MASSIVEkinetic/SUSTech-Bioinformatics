#! /usr/bin/perl
use strict;
use warnings;

my $FileName = $ARGV[0];
open FASTQ, "$FileName" or die("Cannot open the file!\n");
$FileName =~ s/q$/a/;
open FASTA, ">$FileName";

my $criteria = 0;
my $line;
while($line = <FASTQ>){
    chop $line;
    if($line eq "+"){
        $criteria = 1;
    }elsif($criteria == 0){
        print FASTA $line."\n";
    }elsif($criteria == 1){
        $criteria = 0;
    }
}

close FASTA;
close FASTQ;
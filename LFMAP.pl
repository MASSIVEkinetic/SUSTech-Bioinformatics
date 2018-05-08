#! /usr/bin/perl
use strict;
use warnings;
use Storable;

#define some variables
my $counter;
my $subcounter;
my $subsub;
my @temp;
my $temp;
my $IllegalFile = "Illegal file name!\n";
my $SequentialError = "Please run LPMAP.pl index first!\n";
my $usage = "Usage:\n\tLPMAP.pl index file_name.fasta\n\tLPMAP.pl map file1.fasta file2.fastq file3.svf\n";

my %table = (
    A => 6, 6 => "A",
    C => 5, 5 => "C",
    F => 4, 4 => "F",
    G => 3, 3 => "G",
    N => 2, 2 => "N",
    T => 1, 1 => "T",
);

my $AllowedSNV = 2;

if($ARGV[0]){

    if($ARGV[0] eq "index"){

        $ARGV[1] or die $IllegalFile;
        my $SEQref = GetRef($ARGV[1]);

        #Add "F" as the marker at the last of the sequence
        ${$SEQref} = ${$SEQref}."F";
        
        #define @UnsortedString to store unsorted string
        my %UnsortedString;
        #define @SortedString to store sorted string
        my @SortedString;
        #put the input into the array as the first string
        my @suffix;
        #define some other variables
        my $CharCounter;

        #perform the algorithm which repeatly get the first char of the strand and put it at the end.
        for($CharCounter=0;$CharCounter<length(${$SEQref});$CharCounter++){
            $UnsortedString{${$SEQref}} = $CharCounter;
            ${$SEQref} = substr(${$SEQref},1,(length(${$SEQref})-1)).substr(${$SEQref},0,1);
        }

        @SortedString = sort keys %UnsortedString;

        #initialize the first and last column
        my $L = "";
        my $F = "";
        my @LIndex;
        undef @temp;

        for($subcounter=0;$subcounter<6;$subcounter++){
            $temp[$subcounter] = 0;
        }

        #loop every element in the sorted array and obtain their last and first number
        for($subcounter=0;$subcounter<length(${$SEQref});$subcounter++){
            $L = $L.substr($SortedString[$subcounter],-1,1);
            $F = $F.substr($SortedString[$subcounter],0,1);
            $suffix[$subcounter] = $UnsortedString{$SortedString[$subcounter]};
            push(@LIndex,$temp[$table{substr($SortedString[$subcounter],length(${$SEQref})-1,1)}-1]);
            $temp[$table{substr($SortedString[$subcounter],length(${$SEQref})-1,1)}-1] ++;
        }

        push(@LIndex,reverse @temp);

        undef @SortedString;
        undef %UnsortedString;

        store \$L,"LAST.bin";
        store \$F,"FIRST.bin";
        store \@suffix,"SUFFIX.bin";
        store \@LIndex,"LINDEX.bin";

    }elsif($ARGV[0] eq "map"){

        #open the fastq file
        my $REFName = $ARGV[1];
        my $ReadsName = $ARGV[2];
        my $SVFName = $ARGV[3];
        my $F = \${retrieve("FIRST.bin")} or die $SequentialError;
        my $L = \${retrieve("LAST.bin")} or die $SequentialError;
        my $suffix = \@{retrieve("SUFFIX.bin")} or die $SequentialError;
        my $LIndex = \@{retrieve("LINDEX.bin")} or die $SequentialError;
        my $len = length(${$L});
        open REF, "$REFName" or die("Fail to open file! File name: ".$REFName." !\n");
        open READS, "$ReadsName" or die("Fail to open file! File name: ".$ReadsName." !\n");
        open SVF, ">$SVFName" or die("Cannot create file ".$SVFName." !\n");

        print SVF sprintf("%-10s%-10s%-10s%-10s%-10s\n","#CHROM","POS","ID","REF","ALT");

        my %SNV;
        my %SingleReadSNV;
        my @result;
        my $line;
        my $CHROM = "Ecoli";
        my $CurrentLocation;
        my $NextLocation;
        my $mismatch;
        my $exit;

        $temp = 0;

        while($line = <READS>){
            undef @result;
            $line =~ s/\s//g;
            if(substr($line,0,1) eq "@"){
                next;
            }elsif($line eq "+"){
                $temp = 1;
            }elsif($temp == 0){
                for($subcounter=0;$subcounter<$len;$subcounter++){
                    if(substr(${$F},$subcounter,1) eq substr($line,-1,1)){
                        $CurrentLocation = $subcounter;
                        $mismatch = 0;
                        $exit = 0;
                        undef %SingleReadSNV;
                        for($counter=length($line)-2;$counter>=0;$counter--){
                            $NextLocation = ${$LIndex}[$CurrentLocation];                           
                            for($subsub=7-$table{substr(${$L},$CurrentLocation,1)};$subsub>1;$subsub--){
                                $NextLocation += ${$LIndex}[$subsub-8];
                            }
                            
                            if(substr($line,$counter,1) eq substr(${$L},$CurrentLocation,1)){
                                $CurrentLocation = $NextLocation;
                                next;
                            }else{
                                if($mismatch<$AllowedSNV){
                                    
                                    $SingleReadSNV{${$suffix}[$NextLocation]+1} = substr(${$L},$CurrentLocation,1).substr($line,$counter,1);
                                    $mismatch ++;
                                    $CurrentLocation = $NextLocation;
                                    next;
                                }else{
                                    
                                    $exit = 1;
                                    last;
                                }
                            }                            
                        }
                        
                        if($exit==0){
                            for my $key (keys %SingleReadSNV){
                                $SNV{$key} = $SingleReadSNV{$key};
                            }
                        }
                    }
                }
            }elsif($temp == 1){
                $temp = 0;
            }
        }

        for my $key (keys %SNV){
            print SVF sprintf("%-10s%-10s%-10s%-10s%-10s\n",$CHROM,$key,".",substr($SNV{$key},0,1),substr($SNV{$key},1,1)); 
        }

        close REF;
        close SVF;

    }else{

        print $usage;

    }

}else{
    print $usage;
}

# my $IndexInternal = 32;

sub GetRef{
    my $line;
    my $REFSeq = "";
    #open the fasta file
    open REF, "$_[0]" or die("Fail to open file! File name: ".$_[0]." !\n");
    while($line=<REF>){
        $line =~ s/\s//g;
        next if substr($line,0,1) eq ">";
        $line =~ tr/[a-z]/[A-Z]/;
        $REFSeq = $REFSeq.$line;
    }
    return \$REFSeq;
}

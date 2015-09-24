#author         :Weihua Chen
#email          :chenwh550@gmail.com
#last modified  :Dec 15, 2008

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;
use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s","d:s");

if (!$opts{i} or !$opts{o} or !$opts{d}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input gene list file
        -o out file
	-d mgi mouse to go file
----------------------------------------------------------------------\n";
    exit;
    
}
my $ifeblastn=$opts{i};
my $ofgo=$opts{o};
my $database=$opts{d};


open (MGI,$database)||die "Cannot open database : $database!!\n";
my %hGene2GOs = ();
my $id='';
while(<MGI>){
    chomp;
    next if(/^\!/);
    my @arr = split(/\s+/, $_);
    if($arr[3] =~ /GO:/){
        $hGene2GOs{$arr[2]}{$arr[3]} = 1;
    }
}
close MGI;

my %hProcessedGenes = ();
open (OUTGO,">$ofgo")||die "Cannot create output file : $ofgo!!\n" ;
open (IN,$ifeblastn)||die "Cannot open file : $ifeblastn!!\n";
while(<IN>){
    chomp;
    if(/(\S+)/){ ## gene id
        my $gene = $1;
        if(exists $hGene2GOs{$gene}){
            $hProcessedGenes{$gene}=1;
            print OUTGO $gene,"\t",join("\t",sort keys %{$hGene2GOs{$gene}}),"\n";
        }
    }
}
close IN;
close OUTGO;


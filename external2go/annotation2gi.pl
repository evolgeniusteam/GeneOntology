#author         :Weihua Chen
#email          :chenwh550@hotmail.com
#last modified  :

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;

use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input blast.e file
        -o out file
----------------------------------------------------------------------\n";
    exit;
    
}

my %hAnotation2gi=();
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
while (<IN>){
    next if (!/\d+/);
    chomp;
    my ($id,$gi)=/^(\S+).*gi\|(\d+)\|/;
    next if (!$id or !$gi);
    if (exists $hAnotation2gi{$id}){
        $hAnotation2gi{$id}.="\t$gi";
    }else{
        $hAnotation2gi{$id}=$gi;
    }
}
close IN;

while (my ($k,$h)=each %hAnotation2gi) {
	print OUT "$k\t$h\n";
}
close OUT;







#-------------------------------------------------------------------------
#DOCUMENTATION
#-------------------------------------------------------------------------
#this perl script extract gi numbers from blast.e file
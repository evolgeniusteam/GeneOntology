#author         :Weihua Chen
#email          :chenwh550@hotmail.com
#last modified  :June 24, 2008
#

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
        -i input list of gene symbol, one each line
        -o out symbol to GOTerm file
		-d gene_association database from UNIPROT, 'gene_association.goa_human'
----------------------------------------------------------------------\n";
    exit;
    
}
my $ifgenesymbols=$opts{i};
my $ofgo=$opts{o};
my $database=$opts{d};

require '/home/wchen/Dropbox/func.pl';

my %hGeneList = &readfile($ifgenesymbols, 0, 0);

open (UNIPROT,$database)||die "Cannot open database : $database!!\n";
my %storehash;
while(<UNIPROT>)
{
	if (/^UniProt/){
	    my ($db, $acc, $symbol, $goid) = split(/\s+/,$_);
	    $storehash{$symbol}{$goid} = 1 if ($goid =~ /GO/ & exists $hGeneList{$symbol});
	}

}
close UNIPROT;

open (OUTGO,">$ofgo")||die "Cannot create output file : $ofgo!!\n" ;
while(my ($symbol, $refhash) = each %storehash){
	print OUTGO join("\t", $symbol, sort keys %{$refhash}), "\n";
}
close OUTGO;


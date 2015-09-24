#author         :Weihua Chen
#email          :chenwh550@hotmail.com
#last modified  :2006-1-19

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;
use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s","d:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    This is special version for rice only!
----------------------------------------------------------------------
    USAGE: perl $0
        -i input eblast file
        -o out file
----------------------------------------------------------------------\n";
    exit;
    
}
my $ifeblastn=$opts{i};
my $ofgo=$opts{o};
my $database='/database/GO/gene_association.gramene_oryza';
if (defined $opts{d}){
	$database=$opts{d};
}

open (EBLAST,$ifeblastn)||die "Cannot open file : $ifeblastn!!\n";
open (OUTGO,">$ofgo")||die "Cannot create output file : $ofgo!!\n" ;
open (UNIPROT,$database)||die "Cannot open database : $database!!\n";
my %storehash;
while(<UNIPROT>)
{
	if(my ($uniprotID,$goID)=/^GR\t(\S+)\t.*(GO:\S+)\s+/)#$1=prID,$2=GO
	{
		push (@{$storehash{$uniprotID}},$goID);#now every GI has an array containing GOs, all of which stored in a hash:storehash
	}
}
my %hash;
while(<EBLAST>)
{
	chomp;
	next if (!/\d+/);
	next if (!$_);
	my ($q_name,$letter,$queryX,$queryY,$sbjctX,$sbjctY,$length,$score,
	    $e_value,$overlap_total,$identity,$sbject,$anontation)=split(/\t/,$_);
	my ($uniprotID)=$anontation=~/^\((\S+)\)/;
		next if(exists $hash{$q_name});#this ensures every unigene counted only once
		if(exists $storehash{$uniprotID})
		{
			print OUTGO $q_name;
			$hash{$q_name}=1;
			foreach my $goid (@{$storehash{$uniprotID}})
			{
				
				print OUTGO "\t";
				print OUTGO $goid;
			}
			print OUTGO "\n";
		}
}
close EBLAST;
close UNIPROT;
close OUTGO;


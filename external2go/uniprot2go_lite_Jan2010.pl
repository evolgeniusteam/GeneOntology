#author         :Weihua Chen
#email          :chenwh550@hotmail.com
#last modified  :Jan 12, 2010 --
#

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;
use Getopt::Long;

## -- importing the script 'func.pl' depends on host OS, added on June 23, 2010 --
my $os = $^O;
if($os eq 'darwin'){ ## if mac --
    require '/Users/wchen/Dropbox/perl_scripts/func.pl';    
}elsif($os eq 'linux'){ ## -- if linux --
    require '/home/wchen/Dropbox/perl_scripts/func.pl';
}

my %opts=();
GetOptions(\%opts,"i:s","o:s","d:s");

if (!$opts{i} or !$opts{o} or !$opts{d}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
		-i input eblast file
		-o out file
		-d two column acc2go gene association data
			parsed from file 'gene_association.goa_uniprot'
----------------------------------------------------------------------\n";
    exit;
    
}
my $ifeblastn=$opts{i};
my $ofgo=$opts{o};
my $database=$opts{d};

open (EBLAST,$ifeblastn)||die "Cannot open file : $ifeblastn!!\n";
open (OUTGO,">$ofgo")||die "Cannot create output file : $ofgo!!\n" ;
open (UNIPROT,$database)||die "Cannot open database : $database!!\n";
my $hrAcc2GO = &readfile($database, 13, 0, 1);

print STDERR "parse gene_ossociation file done !!\n";
my %hash;
while(<EBLAST>)
{
	chomp;
	next if (!/\d+/);
	next if (!$_);
	my ($q_name,$letter,$queryX,$queryY,$sbjctX,$sbjctY,$length,$score,
	    $e_value,$overlap_total,$identity,$sbject,$anontation)=split(/\t/,$_);
	my ($db, $uniprotID) = split(/\|/, $sbject);
	print $uniprotID, "\n";
	next if(exists $hash{$q_name});#this ensures every unigene counted only once
	if(exists $$hrAcc2GO{$uniprotID}){
		$hash{$q_name}=1;
		print OUTGO join("\t", $q_name, keys %{$$hrAcc2GO{$uniprotID}}), "\n";
	}
}
close EBLAST;
close UNIPROT;
close OUTGO;


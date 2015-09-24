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

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input eblast file
        -o out file
	-d gene_association database, default UNIPROT, optional
      [note]
	last modified : March 2015;
	eblast is in m8 m9 format
----------------------------------------------------------------------\n";
    exit;
    
}
my $ifeblastn=$opts{i};
my $ofgo=$opts{o};
my $database='/database/GO/gene_association.goa_uniprot';
if (defined $opts{d}){
	$database=$opts{d};
}

my %hIDs = ();
open (EBLAST,$ifeblastn)||die "Cannot open file : $ifeblastn!!\n";
while(<EBLAST>)
{
	chomp;
	next if (!/\d+/);
	next if (!$_);
	next if(/^#/); ## if annotation line ...
	my ($q_name,$sbject, $letter,$identity,$overlap_total,$mismatches, $gapopennings, $queryX,$queryY,$sbjctX,$sbjctY,$e_value,$bitscore)=split(/\t/,$_);
	my (undef, $uniprotID, undef)=split(/\|/, $sbject);
	$hIDs{ $uniprotID } ++;
}
close EBLAST;

print STDERR "\teblast done!!\n";

open (EBLAST,$ifeblastn)||die "Cannot open file : $ifeblastn!!\n";
open (OUTGO,">$ofgo")||die "Cannot create output file : $ofgo!!\n" ;
open (UNIPROT,$database)||die "Cannot open database : $database!!\n";
my %storehash;
while(<UNIPROT>)
{
	if (/^UniProt/){
	    my ($db, $id_complete, $id,  $goid) = split(/\s+/,$_);
	    if ($goid =~ /GO/ and exists $hIDs{ $id_complete }){
		$storehash{$id_complete}{$goid} = 1;
		print join("\t", $id_complete, $goid), "\n";
	    }
	}

}

print STDERR "parse gene_ossociation file done !!\n";
my %hash;
while(<EBLAST>)
{
	chomp;
	next if (!/\d+/);
	next if (!$_);
	next if(/^#/); ## if annotation line ...
	my ($q_name,$sbject, $letter,$identity,$overlap_total,$mismatches, $gapopennings, $queryX,$queryY,$sbjctX,$sbjctY,$e_value,$bitscore)=split(/\t/,$_);
	my (undef, $uniprotID, undef)=split(/\|/, $sbject);
		next if(exists $hash{$q_name});#this ensures every unigene counted only once
		if(exists $storehash{$uniprotID})
		{
			print $q_name,"\t",$uniprotID,"\n";
			print OUTGO $q_name;
			$hash{$q_name}=1;
			foreach my $goid (keys %{$storehash{$uniprotID}})
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


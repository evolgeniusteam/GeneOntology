#!/usr/bin/perl
#author:chenweihua
#email:chenwh550@hotmail.com
#for(;;) {Study(); musicit() if (tired);} #:D

#SEE DOCUMENTATION FOR DETAILED INFORMATION
use strict;
my $verion=1.1;
#----------------------------------------------------------------------------------------------------------
#gigo
if (@ARGV!=3){
	print "use:$0 Eblastn output, contig_est file, output file\n";
	exit;
}
my ($ifeblastn,$contig_est,$ofgo)=@ARGV;

open (EBLAST,$ifeblastn)||die;
open (OUTGO,">$ofgo")||die;
open (COMPUGEN,'/data/disk0/zhangb/GO/bin/gene_association.Compugen_GenBank')||die;
my %storehash;
while(<COMPUGEN>)
{
	if(/^CGEN\s+(\S+)\s+GI(\S+)\s+GO:(\S+)\s+/)#$1=prID,$2=GI,$3=GO
	{
		my $array=$2;
		push (@{$storehash{$array}},$3);#now every GI has an array containing GOs, all of which stored in a hash:storehash
	}
}
my %hash;
while(<EBLAST>)
{
	if(/(\S+)[\d\D]*gi\|(\d*)\|/)
	{
		next if(exists $hash{$1});#this ensures every unigene counted only once
		my $term=$1;
		my $gi=$2;
		if(exists $storehash{$gi})
		{
			$hash{$term}=1;
			foreach my $goid (@{$storehash{$gi}})
			{
				print OUTGO $term;
				print OUTGO " ";
				print OUTGO $goid;
				print OUTGO "\n";
			}
		}
	}			
}
close EBLAST;
close COMPUGEN;
close OUTGO;	

#----------------------------------------------------------------------------------------------
#genego.pl
#for GO_UNIGENE

my %idhash=(comp=>'/data/disk0/zhangb/GO/bin/component.id',
	func=>'/data/disk0/zhangb/GO/bin/function.id',
	proc=>'/data/disk0/zhangb/GO/bin/process.id');
	

my %ohash;
my %e_hash;

while (my ($cat, $ifid) = each(%idhash)) {
	%ohash=%e_hash='';
	open(GO,$ofgo)||die;
	my $ofout=$ofgo.'.'.$cat;
	my $oflog=$ofout.'.log';
	open(OUT,">$ofout")||die;
	open(T,">$oflog");
	open(ID,$ifid)||die;
	my @array=<ID>;
	LINE1:while(<GO>)
	{
		if(/(\S+)\s(\d+)/)
		{	
			my $name=$1;
			my $go=$2;	
			my $term;
		
			LINE2:foreach (@array)
			{			
				if($_=~/(\(.*\)[\d\D]*):GO:(.+)/)
				
				{				
					$term=$1;
					my $goid=$2;
					if(!exists $ohash{$term})
					{	
						$ohash{$term}=0;				
					}				
					if($goid=~/($go)/)
					{
						if(exists $e_hash{$name}) 
						{
							foreach (@{$e_hash{$name}})
							{next LINE2 if($_ eq $term);}
						}
						push (@{$e_hash{$name}},$term);
						print T $name."|".$term."\n";
						$ohash{$term}++;
					}
				
				}
				elsif($_=~/GO:($go)/)
				{	
					if(exists $e_hash{$name}) 
					{
						foreach (@{$e_hash{$name}})
						{next LINE2 if($_ eq $term);}
					}
					push (@{$e_hash{$name}},$term);
					print T $name."|".$term."\n";			
					$ohash{$term}++;				
				}
			}
		}		
	}
	my $total=0;
	while( my ($key,$value)=each %ohash )
	{
		next if($value==0);	
		print OUT $key."|";
		print OUT "$value\n";
		$total=$total+$value;
		
	}	
	print "total=$total";
	close GO;
	close OUT;
	close ID;			
}

#----------------------------------------------------------------------------------------------
#contig_est
my %hConest;
open CONEST, $contig_est or die "contig_est file missed\n!!";
while (<CONEST>) {
	if (m/^(\S+)\s+(\d+)\s+.*/) {
		$hConest{$1}=$2;
	}
}
close CONEST;
#----------------------------------------------------------------------------------------------
#genego.pl
#for GO_UNIGENE

%idhash=(comp=>'/data/disk0/zhangb/GO/bin/component.id',
	func=>'/data/disk0/zhangb/GO/bin/function.id',
	proc=>'/data/disk0/zhangb/GO/bin/process.id');


while (my ($cat, $ifid) = each(%idhash)) {
	%ohash=%e_hash='';
	open(GO,$ofgo)||die;
	my $ofout=$ofgo.'.'.$cat.'.EST';
	my $oflog=$ofout.'.log';
	open(OUT,">$ofout")||die;
	open(T,">$oflog");
	open(ID,$ifid)||die;
	my @array=<ID>;
	LINE1:while(<GO>)
	{
		if(/(\S+)\s(\d+)/)
		{	
			my $name=$1;
			my $go=$2;		
			my $size=$hConest{$name};
			my $term;
			LINE2:foreach (@array)
			{			
				if($_=~/(\(.*\)[\d\D]*):GO:(.+)/)
				
				{				
					$term=$1;
					my $goid=$2;
					if(!exists $ohash{$term})
					{	
						$ohash{$term}=0;				
					}				
					if($goid=~/($go)/)
					{
						if(exists $e_hash{$name}) 
						{
							foreach (@{$e_hash{$name}})
							{next LINE2 if($_ eq $term);}
						}
						push (@{$e_hash{$name}},$term);
						print T $name."|".$term."\n";
						$ohash{$term}+=$size;
					}
				
				}
				elsif($_=~/GO:($go)/)
				{	
					if(exists $e_hash{$name}) 
					{
						foreach (@{$e_hash{$name}})
						{next LINE2 if($_ eq $term);}
					}
					push (@{$e_hash{$name}},$term);
					print T $name."|".$term."\n";			
					$ohash{$term}+=$size;				
				}
			}
		}		
	}
	my $total=0;
	while( my ($key,$value)=each %ohash )
	{
		next if($value==0);	
		print OUT $key."|";
		print OUT "$value\n";
		$total=$total+$value;
		
	}	
	print "total=$total";
	close GO;
	close OUT;
	close ID;			
}

#========================================================================================================
#DOCUMENTATION
#========================================================================================================

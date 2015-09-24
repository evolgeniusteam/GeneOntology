if (@ARGV!=2){
	print "use:$0 <Eblast> <out>\n";
	exit;
}
my ($ifeblastn,$ofgo)=@ARGV;

open (EBLAST,$ifeblastn)||die;
open (OUTGO,">$ofgo")||die;
open (COMPUGEN,'/database/GO/gene_association.Compugen_GenBank')||die;
my %storehash;
while(<COMPUGEN>)
{
	if(/^CGEN\s+(\S+)\s+GI(\S+)\s+(GO:\S+)\s+/)#$1=prID,$2=GI,$3=GO
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
			print OUTGO $term;
			$hash{$term}=1;
			foreach my $goid (@{$storehash{$gi}})
			{
				
				print OUTGO "\t";
				print OUTGO $goid;
			}
			print OUTGO "\n";
		}
	}			
}
close EBLAST;
close COMPUGEN;
close OUTGO;	

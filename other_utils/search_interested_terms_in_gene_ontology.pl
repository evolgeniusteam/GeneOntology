#author         :Weihua Chen
#email          :chenwh550@hotmail.com
#last modified  :2005-10-13

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;

use lib qw(/home/chenwh/workspace/lib/perl_lib);
use GO::OntologyProvider::OntologyParser;

#-----------------------------------------------------------------
#geting arguments....

#NOTE: this script searches the three ontology files at specific DIR
#input file should be in native go format
#go-term searching is not case-sensitive...

if (@ARGV!= 3) {
	print "-------------------------------------------------------------------------------------\n";
	print "USAGE: $0 <native_go_file> <outfile> <go_term_interested>";
	print "\n-------------------------------------------------------------------------------------\n\n";
	exit
}


my ($in_go,$out,$term_interested)=@ARGV;

my %hNm_go=();
my %hBioProcess=();
my %hMolFunction=();
my %hCelComponent=();

my $process_ontology='/home/chenwh/pipeline/bin/Gene_Ontology_Tools/process.ontology.txt';
my $function_ontology='/home/chenwh/pipeline/bin/Gene_Ontology_Tools/function.ontology.txt';
my $component_ontology='/home/chenwh/pipeline/bin/Gene_Ontology_Tools/component.ontology.txt';


&read_nm_go_to_hash($in_go,\%hNm_go,);

&retrieve_detailed_infor_by_specific_go_term(\%hBioProcess,\%hNm_go,$process_ontology,$term_interested);
&retrieve_detailed_infor_by_specific_go_term(\%hMolFunction,\%hNm_go,$function_ontology,$term_interested);
&retrieve_detailed_infor_by_specific_go_term(\%hCelComponent,\%hNm_go,$component_ontology,$term_interested);


open OUT, ">$out" or die "Cannot create output file\n";

print OUT "<Cellular Component>\n";
foreach my $cat (sort keys %hCelComponent){
    print OUT "\%$cat\n";
    foreach (@{$hCelComponent{$cat}}){
        print OUT "$_";
    }
}

print OUT "<Molecular function>\n";
foreach my $cat (sort keys %hMolFunction){
    print OUT "\%$cat\n";
    foreach (@{$hMolFunction{$cat}}){
        print OUT "$_";
    }
}

print OUT "<Biological process>\n";
foreach my $cat (sort keys %hBioProcess){
    print OUT "\%$cat\n";
    foreach (@{$hBioProcess{$cat}}){
        print OUT "$_";
    }
}

close OUT;

#------------------------------------------------------------------------
#functions
sub new_array{
	my @arr=@_;
	return \@arr;
}

sub read_nm_go_to_hash{#receive three parameters...1=file_handel, 2=hash_ref, 3=unigene_count;
	my ($in_file,$hStorage)=@_;
	
	print "Reading native nm_go file : $in_file ...";
	
	open NM_GO, $in_file or die "Cannot open file: $in_file!\n";
	while (<NM_GO>){
		chomp;
		my @arr=split(/\t/,$_);
		my $id=shift @arr; #get the first element of the array, the unigene id
		my %hTemp=();
		
		#to elliminate possible redundant go ids for one unigene
		for (@arr){
			if(/\d+/){#if contains validate GO IDs
				$hTemp{$_}=1;
			}
		}
		my $goids=join("\|",keys(%hTemp));
		
		$$hStorage{$id}=$goids;
	}
	
	print "done!!\n\n";
}

sub retrieve_detailed_infor_by_specific_go_term{#receive four parameters, 1=hash_ref_for_catergory_results, 2=hash_ref_of_nm_go, 3=ontology_file_handle, 4=go_term_interested
	my ($hFunctional_catergory,$hReceivedNmGo,$ontology_file,$term_interested)=@_;
	my $go_parser=GO::OntologyProvider::OntologyParser->new(ontologyFile=>$ontology_file);
	#print "now level is $lv\n";
	
	my $errortype='';
	if ($ontology_file=~/process/){
		$errortype='process';
	}elsif($ontology_file=~/component/){
		$errortype='component';
	}elsif($ontology_file=~/function/){
		$errortype='function';
	}
	
	print "Now processing $errortype for current nm_go file ...";
	
	while (my ($unigene_id,$go)=each %{$hReceivedNmGo}){
		my @aGOID=split(/\|/,$go);
		for(@aGOID){
			my $node=$go_parser->nodeFromId($_);
			if ($node){
				my @pathsToRoot=$node->pathsToRoot;
				
				if (@pathsToRoot){
					foreach my $path (@pathsToRoot){
                                            my $msg='';
                                            push @{$path},$node;                                            
 					    for(my $jk=1; $jk<@{$path}; $jk++){
						$msg.="  "x$jk."  ".$path->[$jk]->goid."  ".$path->[$jk]->term."\n";
                                            }
                                            if($msg=~/$term_interested/i){
                                                if (exists $$hFunctional_catergory{$unigene_id}){
                                                    push @{$$hFunctional_catergory{$unigene_id}}, $msg;
                                                }else{
                                                    $$hFunctional_catergory{$unigene_id}=&new_array($msg);
                                                }
                                            }
					}
				}
	
			}
		}

	}

	print "done!!\n\n";
}

#-------------------------------------------------------------------------
#DOCUMENTATION
#-------------------------------------------------------------------------
#
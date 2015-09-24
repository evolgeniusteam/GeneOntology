#author         :Weihua Chen
#email          :chenwh550@hotmail.com
#last modified  :

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;

my $ver=3.0;

print "--------------------------------------------------------------------------------\n";
print "\tversion=$ver\n";
print "--------------------------------------------------------------------------------\n\n";

#require("graph.pm");

use lib qw(/home/chenwh/workspace/lib/perl_lib);
use GO::OntologyProvider::OntologyParser;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"h","i:s@","mark:s@","o:s","p:s","f:s","c:s","level:n","note:s","proc:s","func:s","comp:s","liner","digi");
if (!$opts{i} or !$opts{o}){
    print STDERR <<"    _EOT_";
    Usage: program <options> <specification file> [-i ...[-i ...]]
	   -i       Input nm go file name               [multiple] "nm_go.txt" e.g
           -o       Output file for svg plot            "XXX_go_tree_classify.for plot"
	  [optional]
	   -liner   plot in liner mode
	   -digi    plot in digi instead of percentage mode while using log mode
	   ----------
           -mark    Mark to distinguish multi-column    [multiple] "RICE" e.g
           -h       Help                                
	   ----------
           -c       Component file name                 "component_ontology.txt"   Optional
           -f       Function file name                  "funtion_ontology.txt"     Optional
           -p       Process file name                   "process_ontology.txt"     Optional
	   ----------
	   -level   specify go level that you want retrieve	[2 3]	default=2 Optional
	   ----------
	   -comp    Specify to retrieve the component subcategory	[T/F]	default=True Optional
	   -func
	   -proc
    _EOT_
    exit(1);
}

if (defined $opts{liner} and defined $opts{digi}){
	print STDERR "\n=============================================================\n\tyou cannot use -digi and -liner at the same time!!\n";
	exit(1);
}

my @nm_go_file=@{$opts{i}};
die"\n\tERROR:The number of input file should not bigger than 6!" if (@nm_go_file>6);
my $outfile=$opts{o};

my $process_ontology=defined $opts{p} ? $opts{p}  : '/database/GO/process.ontology';
my $function_ontology=defined $opts{f} ? $opts{f} : '/database/GO/function.ontology';
my $component_ontology=defined $opts{c} ? $opts{c} : '/database/GO/component.ontology';

my $process_search=defined $opts{proc} ? $opts{proc}  : 'T';
my $function_search=defined $opts{func} ? $opts{func} : 'T';
my $component_search=defined $opts{comp} ? $opts{comp} : 'T';


my $level = defined $opts{level} ? $opts{level} : 2;
my $markref=defined $opts{mark} ? $opts{mark} : '';
my @mark=();
@mark=@{$opts{mark}} if ($opts{mark});
my $note=defined $opts{'note'} ? $opts{'note'} : "";

my $isLinerMode=defined $opts{liner} ? 1 : 0;

my $isDigiMode=defined $opts{digi} ? 1 : 0;


#preparing some globel variables...
my $nm_go_file_count=@nm_go_file;
#-----------------------------------------------------------------------
#three hashs to store go information
#see DATA-STRUCTURE for detailed descrptions for their data-structures...
my %hBioProcess=();
my %hMolFunction=();
my %hCelComponent=();

my @maxNumOfUnigenes=();
my %hErrors;
my $maxNumGenes=0;
#-----------------------------------------------------------------------
#set some constant variables that will be used in svg ploting
my @rect_width=("0.6","0.3","0.2","0.15","0.15","0.15");
my @moveper=("0.2","0.2","0.125","0.1","0.1","0.1");
my $note_type="";
#my @color_set=("#9B97C8","#99330C","#0895F7","#F96611");
my @color_set=("#0000FF","#8A2BE2","#A52A2A","#DEB887","#5F9EA0","#7FFF00");

for (my $i=0; $i<@nm_go_file; $i++){
	
	#-----------------------------------------------------------------
	#read the native format of input file, like this
	#unigene	GO:#####	GO:#####
	#one line for one unigene, unigene and its goids seperated by "tabel"
	my %hNm_go=();
	my $total_unigenes;
	
	#-----------------
	&read_nm_go_to_hash($nm_go_file[$i],\%hNm_go,\$total_unigenes);
	
	#please consult the DATA-STRUCTURE section of the documentation for detailed information on the hash %hNm_go
	
	push @maxNumOfUnigenes, $total_unigenes;
	$maxNumGenes=$total_unigenes if ($total_unigenes>$maxNumGenes);
	
	#-----------------------------------------------------------------
	&functional_catergory_of_go(\%hBioProcess,\%hNm_go,\%hErrors,$process_ontology,$total_unigenes,$i,$level) if ($process_search eq 'T');
	&functional_catergory_of_go(\%hMolFunction,\%hNm_go,\%hErrors,$function_ontology,$total_unigenes,$i,$level) if ($function_search eq 'T');
	&functional_catergory_of_go(\%hCelComponent,\%hNm_go,\%hErrors,$component_ontology,$total_unigenes,$i,$level) if $component_search eq 'T';
}

#-----------------------------------------------------------------------
#caculate Yticks for every input file
my %hYTick=();#a hash contains pointers to arrays. like: $hYTick{file1}=[array of YTicks of file1]
for (my $i=0; $i<@nm_go_file; $i++){
	my $max_num=$maxNumOfUnigenes[$i];
	$hYTick{$i}=&new_array(int($max_num/1000),int($max_num/100),int($max_num/10),$max_num);
}


print "Create .list file for svg plot ...";

#-----------------------------------------------------------------------
my %hYData=();#a hash contains pointers to arrays, its data structure is: $hYData{file_num}=@array of frequencies of eash x-element
my @aXTickLabel;#array of labels for x-axis
my @colum_count;#count of each three subcategories 
my $x_end=0;#flag for total elements of x-axis

&get_data_for_GO_plot(\%hCelComponent,\%hMolFunction,\%hBioProcess,\@aXTickLabel,\%hYData,\@colum_count,$isLinerMode,$isDigiMode);

my ($colum_comp, $colum_func, $colum_proc)=@colum_count;
$x_end=$colum_comp+$colum_func+$colum_proc;

my $svg_width=($x_end*45>=1200)?$x_end*45:1500;
my $svg_height=(300+$x_end*3)>400?(300+$x_end*3):400;

if ($x_end eq ''){
	warn"\n!ERROR:\n\tNo one class is generated,please check the input!\n";
	print"\n\tProgram abort!\n\n";
	exit;
}

#-----------------------------------------------------------------------
#creating file for svg plot...
open O, ">$outfile"or die;
my $yStep=33.3;
my $yEnd=100;
my $yStart=0.1;
if ($isLinerMode){
	$yStart=0;
	if ($maxNumGenes<=100){
		$yStep=10;
	}elsif ($maxNumGenes<=500){
		$yStep=50;
	}elsif ($maxNumGenes<=1000){
		$yStep=100;
	}elsif ($maxNumGenes<=5000){
		$yStep=500;
	}elsif ($maxNumGenes<=10000){
		$yStep=1000;
	}elsif ($maxNumGenes<=50000){
		$yStep=5000;
	}else{
		$yStep=10000;
	}
	$yEnd=(int($maxNumGenes/$yStep)+1)*$yStep;	
}elsif ($isDigiMode){
	$yStart=1;
	if ($maxNumGenes<=100){
		$yEnd=100;
		$yStep=50;
	}elsif ($maxNumGenes<=1000){
		$yEnd=1000;
		$yStep=333;
	}elsif ($maxNumGenes<=10000){
		$yEnd=10000;
		$yStep=2500;
	}elsif ($maxNumGenes<=100000){
		$yEnd=100000;
		$yStep=20000;
	}
}
&out($x_end,$rect_width[scalar(@nm_go_file)-1],$moveper[scalar (@nm_go_file)-1],$isLinerMode,$isDigiMode,$yStep,$yEnd,$yStart);
for(@aXTickLabel){
	print O "$_\n";
}
print O ":End\n";
print O "Group:\n";
if ($component_search eq 'T'){
	print O "$colum_comp:Cellular component\n";
}
if ($function_search eq 'T'){
	print O "$colum_func:Molecular function\n";
}
if($process_search eq 'T'){
	print O "$colum_proc:Biological process\n";
}
print O ":End\n\n";

for(my $i=0; $i<@nm_go_file; $i++){
	print O "\nColor:$color_set[$i]\n";
	print O "YMark:r\n";
	print O "Start:0\n";
	print O "End:3\n";
	print O "Step:1\n";
	print O "Scale:\n";
	for(@{$hYTick{$i}}){
		print O "$_\n";
	}
	print O ":End\n";
	print O "Mark:";
	if (@mark){
		print O "$mark[$i]";
	}
	print O "\n";
	for (my $j=0; $j<@{$hYData{$i}};$j++){
		print O "$j:$hYData{$i}[$j]\n";
	}
}
close O;

print "done!!\n\n";

print "Printing out gene_list and error files ...";
#-----------------------------------------------------------------------
#&graphing($list_file,$outfile);
#printing out detailed information for GO functional classifications...
for(my $m=0; $m<@nm_go_file; $m++){
	#the following line ensure the gene.list file output to the current dir
	my @aTemp=split(/\//,$nm_go_file[$m]);
	my $outListFile=$aTemp[-1];
	open OUT, ">$outListFile.gene_list" or warn "Cannot create list file for $outListFile\n";
	
	print OUT "<Cellular component>\n";
	while (my ($k,$h)=each %hCelComponent){
		my $count=$$h{$m}{'count'};
		next if (!$count);
		my $percent= $$h{$m}{'percent'};
		print OUT "  \%$k";
		print OUT "\t$count";
		print OUT "\t$percent\%\n";
		my $aList=$$h{$m}{'list'};
		for(@{$aList}){
			print OUT "    $_\n";
		}
	}

	print OUT "<Molecular Function>\n";
	while (my ($k,$h)=each %hMolFunction){
		my $count=$$h{$m}{'count'};
		next if (!$count);
		my $percent= $$h{$m}{'percent'};
		print OUT "  \%$k\t$count\t$percent\%\n";
		my $aList=$$h{$m}{'list'};
		for(@{$aList}){
			print OUT "    $_\n";
		}
	}

	print OUT "<Biological process>\n";
	while (my ($k,$h)=each %hBioProcess){
		my $count=$$h{$m}{'count'};
		next if (!$count);
		my $percent= $$h{$m}{'percent'};
		print OUT "  \%$k\t$count\t$percent\%\n";
		my $aList=$$h{$m}{'list'};
		for(@{$aList}){
			print OUT "    $_\n";
		}
	}
	close OUT;
	
}

#-----------------------------------------------------------------------
#printing out goids that didn't categorized correctly

open ERR, ">$outfile.errors" or warn "Error messages cannot write into file!!\n";

while (my ($k,$h)=each %hErrors){
	print ERR "Error messages for file $nm_go_file[$k]\n";
	for (@{$h}){
		print ERR $_,"\n";
	}
}

close ERR;

print "done!!\n\n";

print "All jobs done!!\n";

#=========================================================================
#functions...
sub new_array{
	my @arr=@_;
	return \@arr;
}

sub read_nm_go_to_hash{#receive three parameters...1=file_handel, 2=hash_ref, 3=unigene_count;
	my ($in_file,$hStorage,$real_unigenes)=@_;
	
	print "Reading native nm_go file : $in_file ...";
	
	open NM_GO, $in_file or die "Cannot open file: $in_file!\n";
	while (<NM_GO>){
		chomp;
		if (/\d+/){$$real_unigenes++}#count the real unigenes
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

sub functional_catergory_of_go{#receive seven parameters, 1=hash_ref_for_catergory_results, 2=hash_ref_of_nm_go, 3=error_messages, 4=ontology_file_handle, 5=total_num_of_unigenes, 6=file_serial_num, 7=go_levels(go_depth)
	my ($hFunctional_catergory,$hReceivedNmGo,$hErrorMsg,$ontology_file,$total_num,$file_serial,$lv)=@_;
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
	
	print "Now processing $errortype ...";
	
	while (my ($unigene_id,$go)=each %{$hReceivedNmGo}){
		my @aGOID=split(/\|/,$go);
		my %hL2GOforThisUnigene=(); #this hash stores level 2 go infomation teperorily for current unigene
		for(@aGOID){
			my $node=$go_parser->nodeFromId($_);
			if ($node){
				my @pathsToRoot=$node->pathsToRoot;
				
				if (@pathsToRoot){
					foreach my $path (@pathsToRoot){
						#get level 2 go infomation
						push @{$path},$node;
						if (@{$path}>=($lv+1)){#CAUTION: this array does not contain the last node...
							my $tTerm=$path->[$lv]->term;
							#print $path->[2]->goid,"\t",$path->[2]->term,"\n";
							my $termAtLv2=$path->[2]->term;
							if ($termAtLv2!~/obsolete/){
								$hL2GOforThisUnigene{$tTerm}=1;
							}
						}else{
							my $msg=$unigene_id."\n";
							for(my $jk=0; $jk<@{$path}; $jk++){
								$msg.="  "x$jk.$path->[$jk]->goid."  ".$path->[$jk]->term."\n";
							}
							#$msg.="  "x@{$path}.$node->goid."  ".$node->term."\n";
							if ($$hErrorMsg{$file_serial}){
								push @{$$hErrorMsg{$file_serial}}, $msg;
							}else{
								$$hErrorMsg{$file_serial}=&new_array($msg);
							}
						}
		
					}
				}
	
			}
		}

		foreach(keys %hL2GOforThisUnigene){
			next if (/obsolete/);
			if (exists $$hFunctional_catergory{$_}{$file_serial}){
				$$hFunctional_catergory{$_}{$file_serial}{'count'}++;
				push @{$$hFunctional_catergory{$_}{$file_serial}{'list'}},$unigene_id;
			}else{
				$$hFunctional_catergory{$_}{$file_serial}{'count'}=1;
				$$hFunctional_catergory{$_}{$file_serial}{'list'}=&new_array("$unigene_id");
			}
		}
	}
	
	#caculate percentage...
	while (my ($k,$h)=each %{$hFunctional_catergory}){
		if (!$$h{$file_serial}{'count'}){
			$$h{$file_serial}{'count'}=0;
		}
		if (!$$h{$file_serial}{'percent'}){
			$$h{$file_serial}{'percent'}=0;
		}
		$$h{$file_serial}{'percent'}=$$h{$file_serial}{'count'}/$total_num*100;
	}
	
	print "done!!\n\n";
}

sub get_data_for_GO_plot{#accept five parameters, 1,2,3=hash_ref_for_go_catergories, 4=xticklabel, 5=hash_ref_to_receive_Y_data, 6=array_ref_to_receive_num_of_subcatergory;
	my ($go_comp,$go_func,$go_proc,$xticklabel,$y_data,$colum_stat,$linerMode,$digimode)=@_;
	
	#the following codes ensure every sub-category of 
	
	my @comp=();
	foreach my $tHashKey(sort keys %{$go_comp}){
		my $flag=1;#i use this flag to show if all percentages are <0.1
		while (my ($k,$h)=each %{$$go_comp{$tHashKey}}){#$k=file_serial
			if (exists $$h{percent} and $$h{percent}>0.1){
				$flag=0;
			}else {$$h{percent}=0.1;}
		}
		if (!$flag){
			push @comp, $tHashKey;
		}
	}
	push @{$xticklabel}, @comp;
	push @{$colum_stat},scalar(@comp);

	my @func=();
	foreach my $tHashKey(sort keys %{$go_func}){
		my $flag=1;#i use this flag to show if all percentages are <0.1
		while (my ($k,$h)=each %{$$go_func{$tHashKey}}){
			if (exists $$h{percent} and $$h{percent}>0.1){
				$flag=0;
			}else {$$h{percent}=0.1;}
		}
		if (!$flag){
			push @func, $tHashKey;
		}
	}
	push @{$xticklabel}, @func;
	push @{$colum_stat},scalar(@func);

	my @proc=();
	foreach my $tHashKey(sort keys %{$go_proc}){
		my $flag=1;#i use this flag to show if all percentages are <0.1
		while (my ($k,$h)=each %{$$go_proc{$tHashKey}}){
			if (exists $$h{percent} and $$h{percent}>0.1){
				$flag=0;
			}else {$$h{percent}=0.1;}
		}
		if (!$flag){
			push @proc, $tHashKey;
		}
	}	
	push @{$xticklabel}, @proc;
	push @{$colum_stat},scalar(@proc);
	
	my $file_num=scalar(@nm_go_file);
	if (!$linerMode and !$digimode){
		for (my $i=0; $i<$file_num; $i++){
			foreach my $subcat (@comp){
				if (exists $$go_comp{$subcat}{$i}){
					push @{$$y_data{$i}}, $$go_comp{$subcat}{$i}{'percent'};
				}else{
					push @{$$y_data{$i}},0.1;
				}
			}
			foreach my $subcat (@func){
				if (exists $$go_func{$subcat}{$i}){
					push @{$$y_data{$i}}, $$go_func{$subcat}{$i}{'percent'};
				}else{
					push @{$$y_data{$i}},0.1;
				}
			}
			foreach my $subcat (@proc){
				if (exists $$go_proc{$subcat}{$i}){
					push @{$$y_data{$i}}, $$go_proc{$subcat}{$i}{'percent'};
				}else{
					push @{$$y_data{$i}},0.1;
				}
			}
		}
	}else{### liner mode or liner mode
		for (my $i=0; $i<$file_num; $i++){
			foreach my $subcat (@comp){
				if (exists $$go_comp{$subcat}{$i}){
					push @{$$y_data{$i}}, $$go_comp{$subcat}{$i}{'count'};
				}else{
					push @{$$y_data{$i}},0;
				}
			}
			foreach my $subcat (@func){
				if (exists $$go_func{$subcat}{$i}){
					push @{$$y_data{$i}}, $$go_func{$subcat}{$i}{'count'};
				}else{
					push @{$$y_data{$i}},0;
				}
			}
			foreach my $subcat (@proc){
				if (exists $$go_proc{$subcat}{$i}){
					push @{$$y_data{$i}}, $$go_proc{$subcat}{$i}{'count'};
				}else{
					push @{$$y_data{$i}},0;
				}
			}
		}
	}
}

sub out {
	my ($x_end,$width,$moveper,$linerMode,$digimode,$ystep,$yend,$ystart)=@_;
	if ($linerMode){
print O "Type:Simple
Width:$svg_width
Height:$svg_height
BothYAxis:0
MultiRY:0
#MultiY:1
ScaleLen:8
WholeScale:0.9
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
YStep:$ystep
RYStep:$ystep
XScalePos:0.5
XScaleRoate:75
XStart:0
YStart:$ystart
RYStart:$ystart
XEnd:$x_end
YEnd:$yend
RYEnd:$yend
YNeedLog:0
Y:Number of genes
RY:
XUnit:1
Scale:
";
}elsif ($digimode){
print O "Type:Simple
Width:$svg_width
Height:$svg_height
BothYAxis:0
MultiRY:0
#MultiY:1
ScaleLen:8
WholeScale:0.9
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
YStep:$ystep
RYStep:$ystep
XScalePos:0.5
XScaleRoate:75
XStart:0
YStart:$ystart
RYStart:$ystart
XEnd:$x_end
YEnd:$yend
RYEnd:$yend
YNeedLog:10
Y:Number of genes
RY:
XUnit:1
Scale:
";
}else{
print O "Type:Simple
Width:$svg_width
Height:$svg_height
BothYAxis:1
MultiRY:1
#MultiY:1
ScaleLen:8
WholeScale:0.9
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
YStep:$ystep
RYStep:$ystep
XScalePos:0.5
XScaleRoate:75
XStart:0
YStart:$ystart
RYStart:$ystart
XEnd:$x_end
YEnd:$yend
RYEnd:$yend
YNeedLog:10
MarkNoBorder:1
Note:$note_type $note
X:
Y:Percent of genes
RY:Number of genes
XUnit:1
Scale:
";
}
}




#-------------------------------------------------------------------------
#DOCUMENTATION
#-------------------------------------------------------------------------
#i used to use gene_ontology.pl writen by zhk and tools in WEGO at genomics.org.cn to plot my data. unfotunately
#gene_ontology.pl does not provide detailed controls for ploting specific subcatergories and the WEGO
#does not work recently. i decided to write a GO classification script at my wish, surely borrowed much
#from the original gene_ontology.pl and followed its output protocol.

#the original gene_ontology.pl directly outputs a svg file by including a graph.pm and calling its function: graphing
#while in this script only a text file listing elements for ploting as input of the function graphing
#is generated. one of the advantages of doing so is that one can mannually edit the list file, another
#advantage is the some confilts between this script and graph.pm could be avoided, because the 'strict'
#package was used in this script but not in the latter.

#-------------------------------------------------------------------------
#data-structures
#

#------------------------
#%hBioProcess
#$hBioProcess{sub_category_id}=%subHash1;
	#$subHash1{file_serial}=%subHash2;
		#$subHash2{number}=number_of_category_count
		#$subHash2{percent}=percent_of_this_category_to_all_unigenes
		#$subHash2{list}=@, array_of_unigene_list_belonged_to_this_category

#------------------------
#%hNm_go
#$hNm_go{unigene_id}=string
#string='goid|goid|goid|'...

#------------------------
#%hErrors
#$hErrors{file_serial}=%subHash1; 
	#$subHash1{osblete}=%subHash2;	error_type=obsolete|notfound
		#$subHash2{three_GO_categories}=%subHash3;
			#$subHash3{unigeneID}=string. string="goid\tgoid\tgoid". goid seperated by \t.
	#$subHash1{notfound}=string, containing informations from path to root of the GO id.

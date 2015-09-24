#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
my %opts=();
GetOptions(\%opts,"g:s","o:s","i:s@","mark:s@","level:n","title:s","gh:s","ghlevel:n","liner","p:s","om");

#### last modified : Feb 8, 2008
#### Author : Weihua Chen
#### Email : chenwh550@gmail.com

### input file could be :
### geneid  GO:###  GO:###
### or
### geneid  GO:###
### geneid  GO:###

my $ver='3.0';

if (!$opts{g} or !$opts{i} or !$opts{o}){
    print "========================================================================
                    Current version is $ver
=========================================================================
    USAGE: perl $0
        -i input native go file, GENEID GO:######, multiple, no more than 8
        -g input obo_go file
        -o print gene ontology result as text file   
      [optional]
        -p out file for plot
        -om print memmber for each gene ontology classification, require -o option, defalt=FALSE
        -level #th level that you wanna retrieve, default=2
        -mark, multiple
        -title, title of this graph, default='';
        ---------------------
        -gh output go_hierarchy
        -ghlevel go_hierarchy level, default=3
        ---------------------
        -liner plot in liner mode
========================================================================\n";
    exit;
}


my $level=defined $opts{level} ? $opts{level} : 2;
die "EXIT: Level must be >=1\n" if (!$level);

die "EXIT: option -o is required by -om\n" if (exists $opts{om} and !exists $opts{o});

my $ghLevel=defined $opts{ghlevel} ? $opts{ghlevel} : 3;
my $title=defined $opts{title} ? $opts{title} : '';
my $isLinerMode=defined $opts{liner} ? 1 : 0;

#### global variables
my $start_time=time;
my @rect_width=("0.6","0.3","0.2","0.15","0.15","0.15","0.15","0.15");
my @moveper=("0.2","0.2","0.125","0.1","0.1","0.1","0.1","0.1");
#my @color_set=("#99330C","#9B97C8","#A52A2A","#DEB887","#5F9EA0","#7FFF00","#0895F7","#F96611");
my @color_set=("#99330C","#9B97C8","#1E90FF","#8A2BE2","#A52A2A","#DEB887","#5F9EA0","#7FFF00");
###            sky blue, ZiSe     ,brown,    tuhuang,  tianqing, green
my %hGO2Genes=();

my @nm_go_file=@{$opts{i}};
my $fileNumCount=scalar @nm_go_file;
die"\n\tERROR:The number of input file should not bigger than 8!" if (@nm_go_file>8);

my @mark=();
@mark=@{$opts{mark}} if (exists $opts{mark});
my @maxNumOfUnigenes=();
my $maxNumGenes=0;

for (my $i=0; $i<@nm_go_file; $i++){
    my %hUniqGeneIDs = ();
    open IN, $nm_go_file[$i] or die "Cannot open file: $nm_go_file[$i]!\n";
    while(<IN>){
        next if (!/GO:\d+/);
        chomp;
        my ($gene,@aGOIDs)=split(/\s+/,$_);
        $hUniqGeneIDs{$gene} = 1;
        
        #print Dumper (@aGOIDs),"\n";
        foreach my $goid (@aGOIDs){
            $hGO2Genes{$goid}{$i}{$gene}=1; ### add file_serial
        }
    }
    close IN;
    my $geneCount = scalar keys %hUniqGeneIDs;
    $maxNumGenes=$geneCount if ($geneCount>$maxNumGenes);
    push @maxNumOfUnigenes,$geneCount;
}

#print Dumper(@maxNumOfUnigenes);

#### ontology parser
my %hNodesAndRelationship=();
my %hObsoleteTerms=();
&obo_parser($opts{g},\%hNodesAndRelationship,\%hObsoleteTerms);

## searching for putative obsolete go items ..
foreach my $goid (keys %hGO2Genes){
    if (exists $hObsoleteTerms{$goid}){
        my $refhash=$hGO2Genes{$goid};
        foreach my $alt_goid (@{$hObsoleteTerms{$goid}}){
            $hGO2Genes{$alt_goid}=$refhash;
        }
        delete $hGO2Genes{$goid};
    }
}

#### GO search by term
my %hGOResults=();
my %hGOHierarchy=();
while (my ($goid,$refGeneHash)=each %hGO2Genes){
    my $allPath=&paths_to_root($goid);
    #print "now",Dumper(@{$allPath});
    foreach my $local_path(@{$allPath}){
        if (defined $$local_path[$level-1]){
            my $current_node=$$local_path[$level-1]{acc};
            my $name_space=$hNodesAndRelationship{$current_node}{name_space};
            my $name=$hNodesAndRelationship{$current_node}{name};
            #print $name_space,"\t",$name,"\n";
            while (my ($fileSerial,$refHash)=each %{$refGeneHash}){
                #print 'file: ', $fileSerial,"\n";
                foreach my $gene_id (keys %{$refHash}){
                    #$hGOResults{$name_space}{memmber}{$fileSerial}{$gene_id}=1;
                    $hGOResults{$name_space}{$name}{$fileSerial}{$gene_id}=1;
                }
            }
        }

        if (defined $opts{gh} and $fileNumCount==1){
            &gohierarchy(0,$local_path,\%hGOHierarchy,$refGeneHash);
        }
    }
}

if (defined $opts{gh} and $fileNumCount==1){
    open GOH, ">$opts{gh}" or die;
    &printgohierarchy(0,\%hGOHierarchy);
    close GOH;
}

#### print to GRAPH format
#open RECORD,">$opts{rc}" or die "Cannot create file: $opts{rc}!\n";

my %hYData=();#a hash contains pointers to arrays, its data structure is: $hYData{file_num}=@array of frequencies of eash x-element
my @aXTickLabel=();#array of labels for x-axis
my @colum_count=(0,0,0);#count of each three subcategories 

my $yStep=25;
my $yEnd=100;
my $yStart=0.01;
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
}

my $i=0;

if(defined $opts{o}){
    open OUT, ">$opts{o}" or die;
    foreach my $top_level (sort keys %hGOResults){
        my @aSubcat=sort keys %{$hGOResults{$top_level}};
        push @aXTickLabel,@aSubcat;
        my $subcat_count=scalar @aSubcat;
        $colum_count[$i]=$subcat_count;
        $i++;
        print OUT $top_level,"\n";
        foreach my $subcat (@aSubcat){
            print OUT "\t",$subcat,"\n";
            for (my $j=0; $j<@nm_go_file; $j++){
                my $subcat_count=0;
                my @memmber=();
                if (!defined $hGOResults{$top_level}{$subcat}{$j}){
                    if ($isLinerMode){
                        push @{$hYData{$j}},0;
                    }else{
                        push @{$hYData{$j}},$yStart;
                    }
                }else{
                    @memmber=keys %{$hGOResults{$top_level}{$subcat}{$j}};
                    my $count=scalar @memmber;
                    $subcat_count=$count;
                    if ($isLinerMode){
                        push @{$hYData{$j}},$count;
                    }else{
                        if ($count/$maxNumOfUnigenes[$j]*100<=$yStart){
                            push @{$hYData{$j}},$yStart;
                        }else{
                            push @{$hYData{$j}},$count/$maxNumOfUnigenes[$j]*100;
                        }
                    }
                }
                print OUT "\t\t",$subcat_count,"\t",$nm_go_file[$j],"\t";
                print OUT join(",",@memmber) if (exists $opts{om});
                print OUT "\n";
            }
        }
    }
    close OUT;
}
my ($colum_proc,$colum_comp, $colum_func)=@colum_count;
my $xEnd=$colum_comp+$colum_func+$colum_proc;

my $svg_width=($xEnd*45>=1200)?$xEnd*45:1500;
my $svg_height=(300+$xEnd*3)>400?(300+$xEnd*3):400;
if ($xEnd eq ''){
    warn"\n!ERROR:\n\tNo one class is generated,please check the input!\n";
    print"\n\tProgram abort!\n\n";
    exit;
}

#-----------------------------------------------------------------------
#caculate Yticks for every input file
my %hYTick=();#a hash contains pointers to arrays. like: $hYTick{file1}=[array of YTicks of file1]
for (my $i=0; $i<@nm_go_file; $i++){
    my $max_num=$maxNumOfUnigenes[$i];
    @{$hYTick{$i}}=(int($max_num/10000),int($max_num/100),$max_num);
}

&OutputDataForGraph() if (defined $opts{p});

my $used_time=time-$start_time;
print "time used $used_time\n";

sub obo_parser{
    my ($infile,$hash,$obsoleteHash)=@_;

    #### parsing obo_file
    my $backup=$/;
    $/="\n\[";
 
    my $header=1;
    open IN, $infile or die "Cannot open file: $infile!\n";
    while (<IN>){
        if ($header){
            $header=0;
            my ($ver)=/format-version:\s+([0-9.]+)/;
            print "\t===============================================================================\n\t\tParsing OBO file, current version is $ver .\n\t===============================================================================\n";
        }else{
            #### new [item]
            my $acc='';
            my $name_space='';
            my $name='';
            my $is_root=0;
            my $is_obsolete=0;
            
            my @aAlt_id=();
            my @aParents=();
            my @aRecommanded_ids_for_obsolete_id=();
            my @aContents=split(/\n/,$_);
            foreach my $line (@aContents){
                if ($line=~/^id: (\S+)/){
                    $acc=$1;
                    push @aAlt_id,$acc;
                }elsif($line=~/^alt_id: (\S+)/){
                    push @aAlt_id,$1;
                }elsif($line=~/^name: (.+)/){
                    $name=$1;
                }elsif($line=~/^namespace: (.+)/){
                    $name_space=$1;
                }elsif($line=~/^is_a:\s+(GO:\d+)/){
                    push @aParents,$1;
                }elsif($line=~/^relationship: part_of (GO:\d+)/){
                    push @aParents,$1;
                }elsif($line=~/is_obsolete: true/){
                    $is_obsolete=1;
                }elsif($line=~/^comment:/){
                    (@aRecommanded_ids_for_obsolete_id)=$line=~/(GO:\d+)/g;
                }
            }
            #### add new item in hash
            if ($is_obsolete){### if obsolete
                foreach my $local_acc (@aAlt_id){
                    $$obsoleteHash{$local_acc}=\@aRecommanded_ids_for_obsolete_id;
                }
            }else{
                $is_root=1 if !(scalar @aParents);
                my %newHash=(acc=>$acc,parents=>\@aParents,name=>$name,name_space=>$name_space,is_root=>$is_root,is_relationship_build=>0);
                foreach my $local_acc (@aAlt_id){
                    $$hash{$local_acc}=\%newHash;
                }
            }
        }
    }
    close IN;
    
    #### building relationships for all nodes
    foreach my $key (keys %{$hash}){
        if ((!$$hash{$key}{is_root}) and (!$$hash{$key}{is_relationship_build})){ # if not root and relationship not build
            for(my $i=0; $i<@{$$hash{$key}{parents}}; $i++){
                my $parent_go_id=$$hash{$key}{parents}[$i];
                $$hash{$key}{parents}[$i]=$$hash{$parent_go_id};
            }
            $$hash{$key}{is_relationship_build}=1;
        }
    }
}

sub paths_to_root{
    my ($goid)=@_;
    my @aPath=();
    
    #print "goid is $goid\n";
    if (exists $hNodesAndRelationship{$goid}){
        my @aNewPath=($hNodesAndRelationship{$goid});
        push @aPath,\@aNewPath;
        my $has_parents=1;
        while ($has_parents){
            $has_parents=0;
            #### iterate all path
            foreach my $path (@aPath){
                #### check the top node for every path
                if (scalar @{$$path[0]{parents}}){
                    $has_parents++;
                    for (my $i=1; $i<@{$$path[0]{parents}}; $i++){
                        #### duplicate current path
                        my $newPath=&duplicate_path($path);
                        #### add new top node to deplicated path
                        unshift  @{$newPath},$$path[0]{parents}[$i];
                        push @aPath, $newPath;  
                    }
                    #### update current path
                    unshift  @{$path},$$path[0]{parents}[0];
                }else{#### no parents means current node top node
                    #### do nothing
                }
            }
        }
    }
    return \@aPath;
}

sub duplicate_path{
    my ($path)=@_;
    my @aNewPath=();
    foreach my $arrayElement (@{$path}){
        push @aNewPath,$arrayElement;
    }
    return \@aNewPath;
}

sub gohierarchy{
    my ($this_level,$refArrary,$refHash,$refGeneHash)=@_;
    return if ($this_level>=$ghLevel or !defined $$refArrary[$this_level]);
    my $current_node=$$refArrary[$this_level]{acc};
    my $name=$hNodesAndRelationship{$current_node}{name};
    while (my ($fileSerial,$refGeneHash2)=each %{$refGeneHash}){
        foreach my $gene_id (keys %{$refGeneHash2}){
            $$refHash{$name}{memmber}{$gene_id}=1;
        }
    }
    if (!defined $$refHash{$name}{subcat}){
        %{$$refHash{$name}{subcat}}=();
    }
    $this_level++;
    &gohierarchy($this_level,$refArrary,$$refHash{$name}{subcat},$refGeneHash);
}

sub printgohierarchy{
    my ($this_level,$refHash)=@_;
    return if ($this_level>=$ghLevel or (scalar keys %{$refHash})==0);
    foreach my $key (sort keys %{$refHash}){
        my $memmber_count=scalar (keys %{$$refHash{$key}{memmber}});
        print GOH "\t"x$this_level,$key,"\t",$memmber_count,"\t",$memmber_count/$maxNumGenes*100,"\t",join(",",keys %{$$refHash{$key}{memmber}}),"\n";
        &printgohierarchy($this_level+1,$$refHash{$key}{subcat});
    }
}

sub OutputDataForGraph{
    open O,">$opts{p}" or die "Cannot create file: $opts{p}!\n";
print O "Type:Simple
Width:$svg_width
Height:$svg_height
";

print O "BothYAxis:0
MultiRY:0
#MultiY:1
" if ($isLinerMode);

print O "BothYAxis:1
MultiRY:1
#MultiY:1
" if (!$isLinerMode);

print O "ScaleLen:8
WholeScale:0.9
OffsetPer:$rect_width[$fileNumCount-1]
UnitPer:$rect_width[$fileNumCount-1]
MovePer:$moveper[$fileNumCount-1]
XStep:1
YStep:$yStep
RYStep:$yStep
XScalePos:0.5
XScaleRoate:75
XStart:0
YStart:$yStart
RYStart:$yStart
XEnd:$xEnd
YEnd:$yEnd
RYEnd:$yEnd
";

print O "YNeedLog:10
" if (!$isLinerMode);

print O "YNeedLog:0
" if ($isLinerMode);

print O "MarkNoBorder:1
Note:$title
X:
";

print O "Y:Percent of genes
RY:Number of genes
XUnit:1
Scale:
" if (!$isLinerMode);

print O "Y:Number of genes
RY:
XUnit:1
Scale:
" if ($isLinerMode);

for(@aXTickLabel){
    print O "$_\n";
}
print O ":End\n";
print O "Group:\n";
if($colum_proc ){
    print O "$colum_proc:Biological process\n";
}

if ($colum_comp){
    print O "$colum_comp:Cellular component\n";
}
if ($colum_func){
    print O "$colum_func:Molecular function\n";
}
print O ":End\n\n";

for (my $i=0; $i<$fileNumCount; $i++){
    print O "\nColor:$color_set[$i]\n";
    print O "YMark:r\n";
    print O "Start:0\n";
    print O "End:2\n";
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
}

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
my %opts=();
GetOptions(\%opts,"g:s","o:s","i:s","level:n","rc:s");

if (!$opts{g} or !$opts{o} or !$opts{i} or !$opts{level} or !$opts{rc}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input native go file, GENEID GO:######
        -g input obo_go file
        -o out file for plot
        -level #th level that you wanna retrieve
        -rc record functional categories in detail
----------------------------------------------------------------------\n";
    exit;
}


my $level=$opts{level};
die "EXIT: Level must be >=1\n" if (!$level);

#### global variables
my $start_time=time;
my $nGenes=0;
my @rect_width=("0.6","0.3","0.2","0.15","0.15","0.15");
my @moveper=("0.2","0.2","0.125","0.1","0.1","0.1");
my @color_set=("#9B97C8","#99330C","#0895F7","#F96611");
               #dark blue, purple,  Brown,   

my %hGO2Genes=();

open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
while(<IN>){
    next if (!/GO:\d+/);
    chomp;
    my ($gene,@aGOIDs)=split(/\t/,$_);
    $nGenes++;
    
    foreach my $goid (@aGOIDs){
        $hGO2Genes{$goid}{$gene}=1;
    }
}
close IN;
#### ontology parser
my %hNodesAndRelationship=();
my %hObsoleteTerms=();
&obo_parser($opts{g},\%hNodesAndRelationship,\%hObsoleteTerms);


#### searching for putative obsolete go items ..
foreach my $goid (keys %hGO2Genes){
    if (exists $hObsoleteTerms{$goid}){
        foreach my $alt_goid (@{$hObsoleteTerms{$goid}}){
            foreach my $gene_id (keys %{$hGO2Genes{$goid}}){
                $hGO2Genes{$alt_goid}{$gene_id}=1;
            }
        }
        delete $hGO2Genes{$goid};
    }
}

#### GO search by term
my %hGOResults=();
while (my ($goid,$refGeneHash)=each %hGO2Genes){
    my $allPath=&paths_to_root($goid);
    foreach my $local_path(@{$allPath}){
        if (defined $$local_path[$level-1]){
            my $current_node=$$local_path[$level-1]{acc};
            my $name_space=$hNodesAndRelationship{$current_node}{name_space};
            my $name=$hNodesAndRelationship{$current_node}{name};
            foreach my $gene_id (keys %{$refGeneHash}){
                $hGOResults{$name_space}{memmber}{$gene_id}=1;
                $hGOResults{$name_space}{subcat}{$name}{$gene_id}=1;
            }
        }
    }   
}

#### print to GRAPH format
open RECORD,">$opts{rc}" or die "Cannot create file: $opts{rc}!\n";

my @aYData=();#a hash contains pointers to arrays, its data structure is: $hYData{file_num}=@array of frequencies of eash x-element
my @aXTickLabel=();#array of labels for x-axis
my @colum_count=();#count of each three subcategories 

foreach my $top_level (sort keys %hGOResults){
    my $member_count=scalar keys %{$hGOResults{$top_level}{memmber}};
    my $subcat_count=scalar keys %{$hGOResults{$top_level}{subcat}};
    push @colum_count,$subcat_count;
    print RECORD $top_level,"\t",$member_count,"\t",$member_count/$nGenes*100,"\t",join (",", keys %{$hGOResults{$top_level}{memmber}}),"\n";
    foreach my $target_level (sort keys %{$hGOResults{$top_level}{subcat}}){
        my $subcat_member_count=scalar keys %{$hGOResults{$top_level}{subcat}{$target_level}};
        push @aXTickLabel,$target_level;
        push @aYData,$subcat_member_count/$nGenes*100;
        print RECORD "\t",$target_level,"\t",$subcat_member_count,"\t",$subcat_member_count/$nGenes*100,"\t",join(",",keys %{$hGOResults{$top_level}{subcat}{$target_level}}),"\n";
    }
}
close RECORD;
open O,">$opts{o}" or die "Cannot create file: $opts{o}!\n";

my $yStep=25;
my $yEnd=100;
my $yStart=0.01;
my @aRYTick=(int($nGenes/10000),int($nGenes/1000),int($nGenes/100),int($nGenes/10),$nGenes);
my ($colum_comp, $colum_func, $colum_proc)=@colum_count;
my $x_end=$colum_comp+$colum_func+$colum_proc;

my $svg_width=($x_end*45>=1200)?$x_end*45:1500;
my $svg_height=(300+$x_end*3)>400?(300+$x_end*3):400;
if ($x_end eq ''){
    warn"\n!ERROR:\n\tNo one class is generated,please check the input!\n";
    print"\n\tProgram abort!\n\n";
    exit;
}

&out($x_end,$rect_width[0],$moveper[1],0,$yStep,$yEnd,$yStart);

for(@aXTickLabel){
    print O "$_\n";
}
print O ":End\n";
print O "Group:\n";
if ($colum_comp){
    print O "$colum_comp:Cellular component\n";
}
if ($colum_func){
    print O "$colum_func:Molecular function\n";
}
if($colum_proc ){
    print O "$colum_proc:Biological process\n";
}
print O ":End\n\n";

print O "\nColor:$color_set[0]\n";
print O "YMark:r\n";
print O "Start:0\n";
print O "End:4\n";
print O "Step:1\n";
print O "Scale:\n";
for(@aRYTick){
    print O "$_\n";
}
print O ":End\n";
print O "Mark:\n";

for (my $j=0; $j<@aYData;$j++){
    print O "$j:$aYData[$j]\n";
}
close O;

my $used_time=time-$start_time;
print "time used $used_time\n";

#print "total items: \t", scalar (keys %hNodesAndRelationship),"\n";


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


sub out {
    my ($x_end,$width,$moveper,$linerMode,$ystep,$yend,$ystart)=@_;
print O "Type:Simple
Width:$svg_width
Height:$svg_height
";
print O "BothYAxis:0
MultiRY:0
#MultiY:1
" if ($linerMode);

print O "BothYAxis:1
MultiRY:1
#MultiY:1
" if (!$linerMode);

print O "ScaleLen:8
WholeScale:0.9
OffsetPer:$width
UnitPer:$width
MovePer:$moveper
XStep:1
";

print O "YStep:$ystep
RYStep:$ystep
XScalePos:0.5
XScaleRoate:75
XStart:0
YStart:$ystart
RYStart:$ystart
XEnd:$x_end
YEnd:$yend
RYEnd:$yend
";

print O "YNeedLog:10
" if (!$linerMode);

print O "YNeedLog:0
" if ($linerMode);

print O "MarkNoBorder:1
Note:
X:
";

print O "Y:Percent of genes
RY:Number of genes
XUnit:1
Scale:
" if (!$linerMode);

print O "Y:Number of genes
RY:
XUnit:1
Scale:
" if ($linerMode);
}
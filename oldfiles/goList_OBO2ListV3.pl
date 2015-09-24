#!/usr/bin/perl -w
#GO OBO file: (v1.2)
#Stat. the GO classification based on OBO file and sequenceName2GO_List file.
#Author:zhangbing.
#zhangbin@genomics.org.cn
#part_of relationship in OBO file is equal to is_a in current perl program.
#version:2007-03-21.
#version3
#added a function for alternative classification Level parameter.
if(@ARGV!=5){	
	print "Usage: $0 Gene_Ontology_OBO_file sequenceName2GoList_File outfile Classification_Level(>2) total_sequence_num\n";
#	print "\t(sequenceName2GoList_File format: sequenceName	GO_id\n)";
	die;
}
open(F,$ARGV[0])||die;
open(L,$ARGV[1])||die;
open(O,">$ARGV[2]")||die;
my $class=$ARGV[3];
print "running.............\n";

#######################################################
my @id;
my %name;
my %namespace;
my %parent;        #defined the parent hash between up-level and down-level of GO ID
my %alt;           #defined the alternative id hash.(redundant GO ID).
my $id;
my $anno;
my @topLevel;      #defined the three top classification Level.
my %son;           #defined the child hash between up-level and down-level of GO ID.
my %altLevel;
my %grade;

#defined the top three categories.
my $mf="GO:0003674";#molecular_function
my $bp="GO:0008150";#biological_process
my $cc="GO:0005575";#cellular_component

#read OBO file based on oneByOne records.
while(<F>)
{
	chomp;
	if(/\[Term\]/){
		$id='';
		$anno='';		
	}
	elsif(/^id: (\S+)/)
	{
		$id=$1;
		if($id eq $mf)
		{
			push @topLevel,$id;
		}
		elsif( $id eq $bp)
		{
			push @topLevel,$id;
		}
		elsif( $id eq $cc)
		{
			push @topLevel,$id;
		}
		else
		{
			push @id,$id;
		}		
	}
	elsif(/^alt_id: (\S+)/){
		push @{$alt{$id}},$1;
		$altLevel{$1}=$id;
		$name{$1}=$anno;
		$namespace{$1}=$namespace{$id};
	}
	elsif(/^name: (.+)/){
		$name{$id}=$1;
		$anno=$1;
	}
	elsif(/^namespace: (.+)/){
		$namespace{$id}=$1;
				
	}
	elsif(/^is_a: (\S+) \!/){
		push @{$parent{$id}},$1;
		push @{$son{$1}},$id;
	}
	elsif(/^relationship: part_of (\S+) \!/){
		push @{$parent{$id}},$1;
		push @{$son{$1}},$id;
	}
	elsif(/is_obsolete: true/)
	{
		$name{$id}="is_obsolete"."_$namespace{$id}";
		if($namespace{$id} eq 'molecular_function'){
			
			push @{$parent{$id}},$mf;			
			push @{$son{$mf}},$id;
		
		}
		elsif($namespace{$id} eq 'biological_process'){
			
			push @{$parent{$id}},$bp;			
			push @{$son{$bp}},$id;
		}
		elsif($namespace{$id} eq 'cellular_component'){
			
			push @{$parent{$id}},$cc;			
			push @{$son{$cc}},$id;
		}
	}
}
my $mfObs="is_obsolete_molecular_funtion";
push @{$son{$mf}},$mfObs;
$name{$mfObs}=$mfObs;
my $bpObs="is_obsolete_biological_process";
push @{$son{$bp}},$bpObs;
$name{$bpObs}=$bpObs;
my $ccObs="is_obsolete_cellular_component";
push @{$son{$cc}},$ccObs;
$name{$ccObs}=$ccObs;

######################################################################


#take count of the GO classification.
#GOListFile format: SampleName GO:0005575
my %totalNum;
my %sum;
my %obsolete;
my $totalSeq=0;
my %hashSeq;
open (T,">Log.info")||die;
while(<L>)        #read the sample2GOList file
{
	chomp;
	my $sampleName;
	my $goId;
    
	if(/(\S+)\s+(\S+)/)
	{
		$sampleName=$1;
		$goId=$2;
		print "$goId\n";
		#take count of total sequences.
		if(!exists $hashSeq{$sampleName})
		{
			$totalSeq++;
			$hashSeq{$sampleName}=1;
		}		
		##########
		
		if(exists $altLevel{$goId}) #change GO ID when goID is alternative ID.
		{
			my $temp=$altLevel{$goId};
			$goId=$temp;
		}
		
#		print $goId."\t"."Level:".$grade{$goId}."\n";
		print T $1."\t".$name{$goId}."\t|".$namespace{$goId}."\t"."\n";    #writing GO annotation of single sampleName into Log.info file.
		if($name{$goId}=~/is_obsolete/)#judging if the sample is is_obsolete.
		{
			&parentObsFound($goId,$sampleName);
		}
		else
		{
			&parentFound($goId,$sampleName);      #found  GO classification.
		}		
	}
	else
	{
		print "The format of GO list is not true.\n";
		$totalSeq++;
	}
	
}
close T;

#processing the obsolete GO annotation.
sub parentObsFound
{
	my($goId,$sampleName)=@_;
	my $parentId;
	foreach $parentId (@{$parent{$goId}})
	{
		my $sample_go=$sampleName.$parentId;
		if(exists $sum{$sample_go})
		{
			next;
		}
		elsif(!exists $totalNum{$parentId})
		{
			$totalNum{$parentId}=1;
		}
		else
		{
			$totalNum{$parentId}++;
		}
		$sum{$sample_go}=1;
		
	}
	if($name{$goId}=~/molecular_function/)
	{
		my $sample_go=$sampleName.$mfObs;
		if(exists $sum{$sample_go})
		{
		}
		elsif(!exists $totalNum{$mfObs})
		{
			$totalNum{$mfObs}=1;
			$sum{$sample_go}=1;
		}
		else
		{
			$totalNum{$mfObs}++;
			$sum{$sample_go}=1;
		}
	}
	elsif($name{$goId}=~/biological_process/)
	{
		my $sample_go=$sampleName.$bpObs;
		if(exists $sum{$sample_go})
		{
		}
		elsif(!exists $totalNum{$bpObs})
		{
			$totalNum{$bpObs}=1;
			$sum{$sample_go}=1;
		}
		else
		{
			$totalNum{$bpObs}++;
			$sum{$sample_go}=1;
		}
	}
	elsif($name{$goId}=~/cellular_component/)
	{
		my $sample_go=$sampleName.$ccObs;
		if(exists $sum{$sample_go})
		{
		}
		elsif(!exists $totalNum{$ccObs})
		{
			$totalNum{$ccObs}=1;
			$sum{$sample_go}=1;
		}
		else
		{
			$totalNum{$ccObs}++;
			$sum{$sample_go}=1;
		}
	}			
}
	

#found the first-($class)-level of GO classification.
sub parentFound
{
	my ($goId,$sampleName)=@_;
	my $sample_go=$sampleName.$goId;
	
	if(exists $sum{$sample_go})
	{
		return ();
	}
	if(!exists $totalNum{$goId})
	{
		$totalNum{$goId}=1;
		$sum{$sample_go}=1;
	}
	else
	{
		$totalNum{$goId}++;
		$sum{$sample_go}=1;
	}
	foreach my $up_goId (@{$parent{$goId}}){
		&parentFound($up_goId,$sampleName);
	}
}
######################################

$totalSeq=$ARGV[4];#redefined total number of sequences.



#writing GO results with tree-format.
foreach my $top1 (@topLevel)
{
	next if(!exists $totalNum{$top1});
	print O $name{$top1}."\t\|".$totalNum{$top1};
	my $percent=($totalNum{$top1}/$totalSeq)*100;
	printf O "(%.1f)\n",$percent;
	foreach my $sonId (@{$son{$top1}})
	{	
		my $k=2;
		&OutGrade($sonId,$k);
	}
	
	print O "------------------------------------------------\n";
}	

sub OutGrade
{
	my ($sonId,$k)=@_;
	if(!exists $totalNum{$sonId})
	{
		return ();
	}
	else
	{
		for(my $i=2;$i<=$k;$i++)
		{
			print O "\t";
		}
		print O $name{$sonId}."\t|".$totalNum{$sonId};
		my $percent=($totalNum{$sonId}/$totalSeq)*100;
		printf O "(%.1f)\n",$percent;
	}
	$k++;
	if($k<=$class)
	{
		my $j=$k;

		foreach my $son2Id (@{$son{$sonId}})
		{
			&OutGrade($son2Id,$j);
		}
	}
}
close O;
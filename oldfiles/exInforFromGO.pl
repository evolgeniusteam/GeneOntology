#!/usr/bin/perl

#此程序用来从GO分析的初始结果中提取出每个小分类中所含的UNIGENE数。
#输入文件为GO的初始分析结果，REF文件为GO的详细分类列表，保存在f:\00chen\smalltools\datasets目录下，分别为：
#component.id function.id process.id。在分析GO的结果时，选择相应的文件。只需输入文件名就可以了，程序会自动加上目录名。

use Getopt::Long;
use strict;

my $version=1.00;

my %opts;
GetOptions(\%opts,"i=s","o=s", "r=s","h");
if (!(defined $opts{i} and defined $opts{o} and defined $opts{r}) || defined $opts{h}) {			#necessary arguments
	&usage;
}

my $filein=$opts{'i'};
my $fileout=$opts{'o'};
my $ref=$opts{'r'};
$ref = "F:\\00chen\\smallTools\\datasets\\".$ref;

open IN,"<$filein";
my %comp;
while (my $line = <IN>) {
	chomp $line;
	my @arr = split (/\|/, $line);
	my $count = $arr[1];
	my @arr2 = split(/\)/, $arr[0]);
	my @arr3 = split (/-/, $arr2[0]);
	my $id = $arr3[0];
	$id =~ tr/\(//d;
	#print $id;
	if (exists $comp{$id}) {
		$comp{$id} += $count;
	}else {
		$comp{$id} = $count;
	}
}
while ((my $key , my $value) = each %comp) {
	print "$key => $value\n";
}
close IN;

my $ref_temp;
open REF, "<$ref" or die "cannot open ref file\n";
while (my $linein = <REF>) {
	if ($linein =~ /^\(/) {
		$ref_temp .= $linein;
	}
}
close REF;

my %ref_comp;
my @ref = split (/\n/, $ref_temp);
for (@ref) {
	my @arr9 = split (/:/, $_);
	$arr9[0] =~ tr/\(//d;
	my @temp = split(/\)/, $arr9[0]);
	$ref_comp {$temp[0]} = $temp[1];
}
while ((my $key , my $value) = each %ref_comp) {
	print "$key => $value\n";
}

open OUT,">$fileout";
while ((my $key, my $value) = each %comp) {
	print OUT "$ref_comp{$key}\t$value\n";
}
close OUT;


sub usage{
	print <<"USAGE";
Version $version
Usage:
	$0 -i <input file> -o <output file> -r <reference file>
options:
	-i input file
	-o output file
	-r reference file
		component.id 
		function.id 
		process.id
	-h help
USAGE
	exit(1);
}

#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s","d:s");

if (!$opts{i} or !$opts{o} or !$opts{d}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input  file
        -d dtabase
        -o out file
----------------------------------------------------------------------\n";
    exit;
    
}

my $MyDB=$opts{d};
###  Database access
my $MySQL   = "mysql --user=chenwh --password=mylonelyboy"; ### newstyle
##### database issues ...
my $dbh=DBI->connect("dbi:mysql:$MyDB:localhost","chenwh","mylonelyboy",{PrintError=>0,RaiseError=>1}) or die "Can't connect to mysql database: $DBI::errstr\n";

my %hTermName2ID=();
my $fetch_item=$dbh->prepare("SELECT name,acc from term");
$fetch_item->execute();
while(my ($name,$acc)=$fetch_item->fetchrow_array()){
    $hTermName2ID{$name}=$acc;
}

$hTermName2ID{'cellular component unknown'}='GO:0008372';
$hTermName2ID{'molecular function unknown'}='GO:0003674';
$hTermName2ID{'biological process unknown'}='GO:0008150';

open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
while(my $line=<IN>){
    next if ($line!~/\d+/);
    chomp $line;
    $line=~s/;/\t/g;
    $line=~ s/\s+$//;
    my ($nm_acc,$cluster,@aGOname)=split(/\t/,$line);
    if ((scalar @aGOname)<1){
        print "\t",$nm_acc,"\n";
        next;
    }
    print OUT $nm_acc;
    foreach my $goName (@aGOname){
        ($goName)=split(/\[/,$goName);
        next if (!$goName);
        $goName =~ s/^\s+//; #remove leading spaces
        $goName =~ s/\s+$//; #remove trailing spaces
        if (exists $hTermName2ID{$goName}){
            print OUT "\t",$hTermName2ID{$goName};
        }else{
            print $nm_acc,"\t",$goName,"\n";
        }
    }
    print OUT "\n";
}
close IN;
close OUT;


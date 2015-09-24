#!/usr/bin/perl -w
use strict;

use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s", "cutoff:f");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input file
        -o out file
      [optional]
        -cutoff = 0.05, default
----------------------------------------------------------------------\n";
    exit;
    
}

my $cutoff = defined $opts{cutoff} ? $opts{cutoff} : 0.05;

my $temp = 'tmp_'.int(rand(10000));
open TEMP, ">$temp" or die;
print TEMP "data <- read.table(file = \"$opts{i}\", header = TRUE, as.is = TRUE, sep = \"\\t\");
data[,\"dhyper\"] <- 1;
data[, \"note\"] <- '--';
pvalues <- c();
for(i in 1:dim(data)[1]){
    thisdata <- as.numeric(data[i, 4:7]);
    pval <- dhyper(thisdata[4], thisdata[3], thisdata[1]-thisdata[3], thisdata[2]);
    ratio <- thisdata[4] / thisdata[3] - thisdata[2] / thisdata[1];
    
    data[i,\"dhyper\"] <- pval;
    if(pval < $cutoff){
        if(ratio >= 0){
            data[i, \"note\"] <- \"enriched\";
        }else{
            data[i, \"note\"] <- \"depleted\";
        }        
    }
    
    pvalues[i] <- pval;
}

library(\"fdrtool\");
fdr <- fdrtool(pvalues, statistic = \"pvalue\", plot = FALSE);
data[,\"fdr\"] <- fdr\$qval;
write.table(data, file = \"$opts{o}\", sep = \"\\t\", quote=FALSE, row.names = FALSE);
\n";
close TEMP;
system("R < $temp --slave --vanilla");
unlink($temp);

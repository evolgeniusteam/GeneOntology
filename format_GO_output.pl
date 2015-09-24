use strict;
use Getopt::Long;
use Data::Dumper;
use Spreadsheet::WriteExcel;


my %opts=();
GetOptions(\%opts,"o:s","i:s","mark:s@");

if (!$opts{i} or !$opts{o} or !$opts{mark}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
      [note] : format input text file to Excel table, see example.txt for input example
        -i input file
        -o out file
        -mark input two 
----------------------------------------------------------------------\n";
    exit;
    
}

my @mark = @{$opts{mark}}; die "please input two marks!!!\n" if((scalar @mark) != 2);
my $outXlsFile = $opts{o} =~ /\.xls$/ ? $opts{o} : $opts{o}.'.xls';
my ($file1, $file2) = @mark;
my $sheetId = $file1."_vs_".$file2;
my $workbook = Spreadsheet::WriteExcel->new($outXlsFile);
my $worksheet = $workbook->add_worksheet($sheetId);

#################################################
## predefine some cell formats
#################################################
my $fTitle = $workbook->add_format(bold=>1, color=>'blue', size=>12);
my $fEnriched = $workbook->add_format(bg_color=>"red", pattern=>1);
my $fDepleted = $workbook->add_format(bg_color=>'green', pattern=>1);
my $fYellowBg = $workbook->add_format(bg_color=>'yellow');
# ----------------------------------------------

my $current_row = 0;
$worksheet->write($current_row, 0, "Note: depleted means genes in ".$file2." are depleted comparing to ".$file1);

$current_row += 2;
my @aOsTitles = ("aspect", "golevel", "goterm", $file1."Total", $file2."Total", $file1."InThisCat", $file2."InThisCat", "note", "p_fisher", "fdr");
for(my $col9875 = 0; $col9875 < @aOsTitles ; $col9875 ++){
    $worksheet->write($current_row, $col9875, $aOsTitles[$col9875], $fTitle);
}
$current_row ++;

open IN, $opts{i} or die;
$/ = "\n";
while(<IN>){
    chomp;
    next if(!/\d/);
    my $col2365 = 0;
    my ($aspect, $golevel, $goterm, $total1, $total2, $go1, $go2, $note, $p_fisher, $fdr) = split("\t", $_);
    $worksheet->write($current_row, $col2365, $aspect); $col2365 ++;
    $worksheet->write($current_row, $col2365, $golevel); $col2365 ++;
    $worksheet->write($current_row, $col2365, $goterm); $col2365 ++;
    $worksheet->write($current_row, $col2365, $total1); $col2365 ++;
    $worksheet->write($current_row, $col2365, $total2); $col2365 ++;
    $worksheet->write($current_row, $col2365, $go1); $col2365 ++;
    $worksheet->write($current_row, $col2365, $go2); $col2365 ++;
    
    if($note eq 'enriched'){
        $worksheet->write($current_row, $col2365, $note, $fEnriched); $col2365 ++;
    }elsif($note eq 'depleted'){
        $worksheet->write($current_row, $col2365, $note, $fDepleted); $col2365 ++;
    }else{
        $worksheet->write($current_row, $col2365, $note); $col2365 ++;
    }
    
    $worksheet->write($current_row, $col2365, $p_fisher); $col2365 ++;
    
    if($fdr <= 0.05){
        $worksheet->write($current_row, $col2365, $fdr, $fYellowBg); $col2365 ++;
    }else{
        $worksheet->write($current_row, $col2365, $fdr); $col2365 ++;
    }
    $current_row ++;
}
close IN;
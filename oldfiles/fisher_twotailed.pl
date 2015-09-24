use Text::NSP::Measures::2D::Fisher2::twotailed;

## n11 n12 | n1p
## n21 n22 | n2p
## --------|----
## np1 np2 | npp

my $n11 = 26603 - 851;
my $n12 =  851;
my $n21 = 9695 - 295;
my $n22 = 295;

my $n1p = $n11 + $n12;
my $n2p = $n21 + $n22;
my $np1 = $n11 + $n21;
my $np2 = $n12 + $n22;

my $npp = $n11 + $n12 + $n21 + $n22;

my $twotailed = calculateStatistic( n11=>$n11,
                                      n1p=>$n1p,
                                      np1=>$np1,
                                      npp=>$npp);

print "twotailed\t", $twotailed, "\n";

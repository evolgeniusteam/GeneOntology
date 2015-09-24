use strict;

# There are n "good" and m "bad" balls in an urn.
# Pick N of them. The probability of i or more successful selections:
# $n, $m, $N, $i
print hypergeom(851, 26603 - 851, 9695, 295),"\n";
## Total    UniTotal    TotalInThisGO   UniInThisGO
## 26603    9695	851	        295

##################################################################
# -------->
# Hyper Geometric distribution
sub hypergeom {
    my ($n, $m, $N, $i) = @_;
    # There are n "good" and m "bad" balls in an urn.
    # Pick N of them. The probability of i or more successful selections:
    # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
    my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
    my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
    return exp($loghyp1 - $loghyp2);
}

sub logfact {
    return gammln(shift(@_) + 1.0);
}

sub gammln {
    my $xx = shift;
    my @cof = (76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.12086509738661e-2, -0.5395239384953e-5);
    my $y = my $x = $xx;
    my $tmp = $x + 5.5;
    $tmp -= ($x + .5) * log($tmp);
    my $ser = 1.000000000190015;
    for my $j (0..5) {
        $ser += $cof[$j]/++$y;
    }
    return log(2.5066282746310005*$ser/$x) - $tmp;
}
# Hyper Geometric distribution
# <--------
##################################################################
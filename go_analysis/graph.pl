#BE CAREFUL THIS SCRIPTS ONLY PLOT INPUT FILE IN UNIX FORMAT UNDER LINUX CORRECTLY
use File::Basename;
use lib dirname (__FILE__);

print "\nBE CAREFUL THIS SCRIPTS ONLY PLOT INPUT FILE IN UNIX FORMAT UNDER LINUX CORRECTLY\n";
print "IF YOU RECEICED \" SOME THING WRONG\" MESSAGE ON THE SCREEN,\n";
print "PLEASE CHECK YOUR INPUT FILE FIRST!!!\n";

my ($in,$out)=@ARGV;
if ($out!~/.svg$/){
	$out.='.svg';
}
&graphing($in,$out);

sub graphing {
	
	my ($input,$output)=@_;
	#print $input,"\t",$output,"\n";
	my @zhk=();
	if ($input ne ''){
		push (@zhk,$input);
	}
	if ($output ne ''){
		push (@zhk,$output);
	}

#!/usr/local/bin/perl
#Author:Li Shengting
#E-mail:lishengting@genomics.org.cn
#Program Date:2002-7-10 11:31
#Last Update:2002-11-26 11:40
#Describe:���������͵ķֲ�ͼ
my $ver=2.2;#���Ը�������ֺ�
$ver=3.00;#����״ͼ
$ver=3.20;#�޸�һЩbug
$ver=3.40;#��ǿ�?�����Լ����������ʾ
$ver=3.60;#���������С��bug���߸�ΪԲ�ˣ�
		  #������ע�͹���
$ver=3.62;#�޸�������
$ver=3.63;#���Ӻ������ʾ
$ver=3.64;#�޸������ΪС��ʱ�Ĵ���
$ver=3.65;#��ߺ������ʾ��ȷ��
$ver=3.70;#���ӵ�ͼ����
		  #��������ָ����ʾ
		  #��ע��λ��ѡ��;
		  #����������
$ver=3.71;#�޸������ʾ
$ver=3.72;#�����������Html�ַ����ʾ����
$ver=3.73;#����powerpointɫ��
$ver=3.74;#�������ƫ������
$ver=3.75;#2002-9-20 13:55 ���Ӿ�������
$ver=3.76;#2002-9-20 13:58 ����ǿ�����Ͳ���
$ver=3.77;#2002-9-24 10:34 ������Ч���ֵ���������
$ver=3.80;#2002-9-28 10:41 �����Ҷ������ʾ
$ver=3.81;#2002-9-30 19:01 ����autofit
$ver=3.82;#2002-10-1 11:28 �ſ�����ʽҪ��
$ver=3.83;#2002-10-11 16:57 ����ִ�Сд;x y�������;�ɸı��߿�
$ver=3.84;#2002-10-12 22:52 ����ע��
$ver=3.85;#2002-10-18 11:16 ����һ��bug
$ver=3.86;#2002-10-20 14:40 ���Ӿֲ�����
$ver=3.90;#2002-10-20 16:46 y���ǿ��ʾ����
$ver=3.91;#2002-10-21 0:59 �����Զ�ȡlog����
$ver=3.92;#2002-10-22 13:46 ����bug
$ver=4.00;#2002-10-23 20:23 �������������ʾ
$ver=4.01;#2002-10-24 22:30 ���Ӵ�ֱ����ͼ
$ver=4.02;#2002-10-25 0:45 �޸�һ�Զ��bug(�޸�ɢ�б�ṹ)
$ver=4.03;#2002-10-26 23:56 �޸�YNeedLog!=0,YLog=0ʱbug
$ver=4.04;#2002-10-27 23:55 ������ɫ����
$ver=4.05;#2002-10-28 14:13 ����ָ���½�
$ver=4.06;#2002-10-28 14:37 �����ͼ�ηֿ��趨 YAxis,YMark
$ver=4.07;#2002-10-30 16:48 ���ӱ߿����
$ver=4.08;#2002-11-2 13:05 ����bug
$ver=4.10;#2002-11-5 11:50 ����͸������
$ver=4.11;#2002-11-5 20:07 ��������ע��λ��
$ver=4.12;#2002-11-7 4:10 ���������ͼ�ֿ�����
$ver=4.13;#2002-11-7 14:51 ����β���ո�
$ver=4.14;#2002-11-8 15:17 ������������
$ver=4.15;#2002-11-8 22:41 ���ӽ�׳��
$ver=4.16;#2002-11-9 1:07 �������������
$ver=4.17;#2002-11-9 4:53 ���ӷ����ͱ�ʶ
$ver=4.18;#2002-11-12 11:22 ����bug
$ver=4.19;#2002-11-18 14:24 ����bug
$ver=4.20;#2002-11-26 11:40 �����ͼ��ȱ��,ColorStep��Ϊ����
use strict;
use Getopt::Long;
use Svg::File;
use Svg::Graphics;

my %opts;
GetOptions(\%opts,"p!","t:s","f!");
die "Usage: $0 <list_file> <svg_file> [-p] [-t type] [-f]\n" if (@zhk<2);
#Constant
sub error;
my $XOFFSET=10;
my $YOFFSET=10;
my $CHRW=6;
my $CHRH=12;
my $XSPACE=2;
my $YSPACE=2;
my $SCALELEN=5;
my $MARK_SCALE=3/4;
my $MARK2_SCALE=3/4;
my $EXP_SCALE=0.3;
my $MARK_POINT_SIZE=5;
my $PI=3.1415926;

my $svg = Svg::File->new($zhk[1],\*OUT);
$svg->open("public","encoding","iso-8859-1");
my $g = $svg->beginGraphics();

#Globe Variable
#print map {"$_<-->$opts{$_}\n"} sort keys %opts;
my $pp=$opts{p};
my $type=$opts{t};
my $autofit=$opts{f};
my (%param,%rect,$ch,$vbW,$vbH,$xZero,$yZero,$ryZero);
my ($xDiv,$yDiv,$ryDiv,$x,$y,$ry,$colWidth,$xDDigits,$yDDigits,$ryDDigits,
$mark,$part,$offset,$unitPer,$tmp,$noConnect,$yMlen,$ryMlen,
$maxXScale,$maxYScale,$maxRYScale,$step,$realCh,
@mark,@tmp,@note,@mark2,@xScale,@yScale,@ryScale,@y,@ry,%raw,@group,@gxy);
my %d=(
	path=>{
		'stroke'=>"#000000",
		'fill'=>"none",
		'stroke-width'=>3
	},
	font=>{
		'fill'=>"#000000",
		'font-size'=>46,
		#'font-weight'=>'bold',
		'font-family'=>"Arial"	#modify 2002-12-03 zhk
		#'font-family'=>"ArialNarrow-bold"
	},
	rect=>{
		'fill'=>"none",
		'stroke'=>"black",
		'stroke-width'=>3
	},
	circle=>{
		'fill'=>"none",
		'stroke'=>"black"
	},
	polygon=>{
		'fill'=>"black",
		'stroke'=>"black"
	},
);

my @pShape=(
	'circle','rect',
	'circle','rect',
	'circle','rect',
	'circle','rect',
);

my @type=(
	'Rect',
	'Double',
	'Line',
	'Point',
	'Bar',
);

%param=();
@note=();
my @keys=(
"Type","Fill","Width","Height","WholeScale","BothYAxis","ScaleLen",
"MarkPos","MarkStyle","MarkNoBorder","MarkScale",
"Mark2Pos","Mark2Border","Mark2Scale",
"HaveMore","HaveLess","RightAngle","Part","OffsetPer","UnitPer","MovePer",
"FontSize","FontFamily","FontBold",
"XScalePos","XScaleLinePos","XUnit","XScaleRoate","XCut",
"YHasLow","YCut",
"XStart","XEnd","XStep","XLog","X","XDispStep","XDDigits",
"XDiv","XScaleDiv","XZeroPos","XZeroVal","XExp","XNeedLog",
"YStart","YEnd","YStep","YLog","Y","YDispStep","YDDigits",
"YDiv","YNum","YScaleDiv","YZeroPos","YZeroVal","YExp","YNeedLog",
"RYStart","RYEnd","RYStep","RYLog","RY","RYDispStep","RYDDigits",
"RYDiv","RYNum","RYScaleDiv","RYZeroPos","RYZeroVal","YExp","RYNeedLog",
"MultiY","MultiRY","Transparence",
"Note","AvailDigit","PointSize","LineWidth","VerticalLine",
"LineDash","NoLine","NoConnect",
);

open(F,$zhk[0]) || die "Can't open $zhk[0]!\n";
while (<F>) {
	last if ($_!~/\S/);
	if (/^\s*Note2:/i) {
		while (<F>) {
			last if ($_=~/^:End$/i);
			push(@note,$_);
		}
	}
	if (/^\s*Group:/i) {
		while (<F>) {
			last if ($_=~/^:End$/i);
			$_=~/^\s*(\d+)\s*:(.+)/;
			$#group++;
			$group[$#group]{len}=$1;
			$group[$#group]{mark}=$2;
		}
	}
	if (/^\s*Mark2:/i) {
		while (<F>) {
			last if ($_=~/^:End$/i);
			push(@mark2,$_);
		}
	}
	if ($_=~/^\s*Scale:/i || $_=~/^\s*XScale:/i) {
		while (<F>) {
			last if ($_=~/^:End$/i);
			chomp;
			push(@xScale,$_);
			$maxXScale=$_ if (length($maxXScale)<length($_));
		}
	}
	if (/^\s*YScale:/i) {
		while (<F>) {
			last if ($_=~/^:End$/i);
			chomp;
			push(@{$yScale[0]},$_);
			$maxYScale=$_ if (length($maxYScale)<length($_));
		}
	}
	if (/^\s*RYScale:/i) {
		while (<F>) {
			last if ($_=~/^:End$/i);
			chomp;
			push(@{$ryScale[0]},$_);
			$maxRYScale=$_ if (length($maxRYScale)<length($_));
		}
	}
	if (/(\S+?):(.+)/) {
		$param{lc($1)}=$2;
	}
}
foreach (@keys) {
	if (exists($param{lc($_)})) {
		$param{$_}=$param{lc($_)};
	}
}
my @stict_keys=("XStart","YStart","XEnd","YEnd","XStep","YStep");
foreach (@stict_keys) {
	if (!exists($param{$_})) {
		error "$_ must be defined!";
	}
}
#####################################################################################################################
#	SET param
#####################################################################################################################
%raw=%param;
error "XStart can't equal with XEnd!" if ($param{XStart}==$param{XEnd});
error "XStep must bigger than 0!" if ($param{XStep}<=0);
error "YStart can't equal with YEnd!" if ($param{YStart}==$param{YEnd});
error "YStep must bigger than 0!" if ($param{YStep}<=0);
if ($param{BothYAxis}) {
	error "RYStart can't equal with RYEnd!" if ($param{RYStart}==$param{RYEnd});
	error "RYStep must bigger than 0!" if ($param{RYStep}<=0);
}
$param{XLog}=10 if ($param{XLog}==1);
$param{YLog}=10 if ($param{YLog}==1);
$param{XExp}=10 if ($param{XExp}==1);
$param{YExp}=10 if ($param{YExp}==1);
$param{XExp}= ($param{XExp} ? $param{XExp} : 1);
$param{YExp}= ($param{YExp} ? $param{YExp} : 1);
$param{XNeedLog}=10  if ($param{XNeedLog}==1);
$param{YNeedLog}=10  if ($param{YNeedLog}==1);
$param{XLog}=$param{XNeedLog} if ($param{XNeedLog} && !exists($param{XLog}));
$param{YLog}=$param{YNeedLog} if ($param{YNeedLog} && !exists($param{YLog}));
error("XExp without XLog!\tError may occur!",1) if ($param{XExp}!=1 && !$param{XLog});
error("YExp without YLog!\tError may occur!",1) if ($param{YExp}!=1 && !$param{YLog});
if ($param{BothYAxis}) {
	@stict_keys=("RYStart","RYEnd","RYStep");
	foreach (@stict_keys) {
		if (!exists($param{$_})) {
			if ($param{MultiRY}) {
			}else{
				error "$_ must be defined!";
			}
		}
	}
	$param{RYLog}=10 if ($param{RYLog}==1);
	$param{RYExp}=10 if ($param{RYExp}==1);
	$param{RYExp}= ($param{RYExp} ? $param{RYExp} : 1);
	$param{RYNeedLog}=10  if ($param{RYNeedLog}==1);
	$param{RYLog}=$param{RYNeedLog} if ($param{RYNeedLog} && !exists($param{RYLog}));
	error("RYExp without RYLog!\tError may occur!",1) if ($param{RYExp}!=1 && !$param{RYLog});
}
if ($param{XDiv}) {
	if ($param{XLog} && !$param{XNeedLog}) {
		$param{XStart}-=log($param{XDiv})/log($param{XLog});
		$param{XEnd}-=log($param{XDiv})/log($param{XLog});
	}else{
		$param{XStart}/=$param{XDiv};
		$param{XEnd}/=$param{XDiv};
		$param{XStep}/=$param{XDiv};
		$param{XUnit}/=$param{XDiv};
	}
}
if ($param{XNeedLog}) {
	$step = rint(($param{XEnd}-$param{XStart})/$param{XStep});
	error("XStart and XEnd must bigger than 0!") if ($param{XStart}<=0 || $param{XEnd}<=0);
	$param{XStart}=log($param{XStart})/log($param{XNeedLog});
	$param{XEnd}=log($param{XEnd})/log($param{XNeedLog});
	$param{XStep}=availDigit(($param{XEnd}-$param{XStart})/$step,$param{AvailDigit},1,$param{XLog});
	$param{XUnit}=log($param{XUnit})/log($param{XNeedLog}) if ($param{XUnit});
}

if ($param{XStart}!=0 && !exists($param{XZeroPos}) && !exists($param{XZeroVal})) {
	$param{XZeroVal}=$param{XStart};
}
if ($param{YDiv}) {
	if ($param{YLog} && !$param{YNeedLog}) {
		$param{YStart}-=log($param{YDiv})/log($param{YLog});
		$param{YEnd}-=log($param{YDiv})/log($param{YLog});
	}else{
		$param{YStart}/=$param{YDiv};
		$param{YEnd}/=$param{YDiv};
		$param{YStep}/=$param{YDiv};
	}
}
if ($param{YNeedLog}) {
	#error("YStart,YEnd and YStep must bigger than 0!") if ($param{YStart}<=0 || $param{YEnd}<=0 || $param{YStep}<=0);
	$step = rint(($param{YEnd}-$param{YStart})/$param{YStep});
	$param{YStart}=log($param{YStart})/log($param{YNeedLog});
	$param{YEnd}=log($param{YEnd})/log($param{YNeedLog});
	$param{YStep}=availDigit(($param{YEnd}-$param{YStart})/$step,$param{AvailDigit},1,$param{YLog});
}
if ($param{YStart}!=0 && !exists($param{YZeroPos}) && !exists($param{YZeroVal})) {
	$param{YZeroVal}=$param{YStart};
}
if ($param{BothYAxis}) {
	if ($param{RYDiv}) {
		if ($param{RYLog} && !$param{RYNeedLog}) {
			$param{RYStart}-=log($param{RYDiv})/log($param{RYLog});
			$param{RYEnd}-=log($param{RYDiv})/log($param{RYLog});
		}else{
			$param{RYStart}/=$param{RYDiv};
			$param{RYEnd}/=$param{RYDiv};
			$param{RYStep}/=$param{RYDiv};
		}
	}
	if ($param{RYNeedLog}) {
		if ($param{RYStart}<=0 || $param{RYEnd}<=0 || $param{YStep}<=0) {
			if ($param{MultiRY}) {
			}else{
				error ("RYStart,RYEnd and RYStep must bigger than 0!");
			}
		}else{
			$step = rint(($param{RYEnd}-$param{RYStart})/$param{RYStep});
			$param{RYStart}=log($param{RYStart})/log($param{RYNeedLog});
			$param{RYEnd}=log($param{RYEnd})/log($param{RYNeedLog});
			$param{RYStep}=availDigit(($param{RYEnd}-$param{RYStart})/$step,$param{AvailDigit},1,$param{RYLog});
		}
	}
	if ($param{RYStart}!=0 && !exists($param{RYZeroPos}) && !exists($param{RYZeroVal})) {
		$param{RYZeroVal}=$param{RYStart};
	}
}
if ($param{FontBold}) {
	$d{font}{'font-weight'}='bold';
}
$param{Type}=$type if ($type ne '');
$param{Type}=~/\S+/;
$type=$&;
$type=$type[0] if ($type eq '' || $type eq 'Rect' || $type eq 'Rectangle');
$param{XZeroPos}=fx($param{XZeroPos}) if (exists($param{XZeroPos}));
$param{YZeroPos}=fy($param{YZeroPos}) if (exists($param{XZeroPos}));
$param{MarkPos}= ($param{MarkPos} ? $param{MarkPos} : 'right');
$param{MarkStyle}= ($param{MarkStyle} ? $param{MarkStyle} : 'v');
$param{MarkNoBorder}= ($param{MarkNoBorder} ? $param{MarkNoBorder} : 0);
$param{Mark2Pos}= ($param{Mark2Pos} ? $param{Mark2Pos} : 'left');
$param{Mark2Border}= ($param{Mark2Border} ? $param{Mark2Border} : 0);
$param{PointSize}= ($param{PointSize} ? $param{PointSize} : 5);
$param{XScaleDiv}= ($param{XScaleDiv} ? $param{XScaleDiv} : 1);
$param{YScaleDiv}= ($param{YScaleDiv} ? $param{YScaleDiv} : 1);
$param{XDispStep}= ($param{XDispStep} ? $param{XDispStep} : 1);
$param{YDispStep}= ($param{YDispStep} ? $param{YDispStep} : 1);
$param{WholeScale}= ($param{WholeScale} ? $param{WholeScale} : 1);
$param{Y}=~ s/([<&>])/sprintf("&#%d;", ord($1))/ge; 
$param{X}=~ s/([<&>])/sprintf("&#%d;", ord($1))/ge; 
$param{Note}=~ s/([<&>])/sprintf("&#%d;", ord($1))/ge;
if ($param{BothYAxis}) {
	$param{RY}=~ s/([<>])/sprintf("&#%d;", ord($1))/ge;
	$param{RYZeroPos}=fy($param{RYZeroPos});
	$param{RYScaleDiv}= ($param{RYScaleDiv} ? $param{RYScaleDiv} : 1);
	$param{RYDispStep}= ($param{RYDispStep} ? $param{RYDispStep} : 1);
}
#####################################################################################################################
#	SET
#####################################################################################################################
if ($pp) {
	$d{font}{fill}='#F0F000';
}
$MARK_SCALE=$param{MarkScale} if ($param{MarkScale} ne '');
$MARK2_SCALE=$param{Mark2Scale} if ($param{Mark2Scale} ne '');
$SCALELEN=$param{ScaleLen} if ($param{ScaleLen});
$noConnect=(($param{NoLine} || $param{NoConnect})? 1 : 0);
$d{font}{'font-size'}=sprintf("%.2f",$param{FontSize}) if ($param{FontSize} != 0);
if ($param{FontFamily} ne '') {
	$param{FontFamily}=~/^\s*(.+?)\s*$/;
	$d{font}{'font-family'}=$1;
	#print "\"$d{font}{'font-family'}\"\n";
}
	$g->setFontSize($d{font}{'font-size'});
	$g->setFontFamily($d{font}{'font-family'});
	$g->setFontColor($d{font}{fill});
	$g->setCharSpacing(0);
	$g->setWordSpacing(0);
$param{XEnd}+=$param{XUnit}*uint($param{XScaleLinePos},1);
$xDiv=($param{XEnd}-$param{XStart})/($param{XStep} ? $param{XStep}:$param{XUnit});
$yDiv=($param{YEnd}-$param{YStart})/($param{YStep} ? $param{YStep}:1);
if ($param{BothYAxis}) {
	$ryDiv=($param{RYEnd}-$param{RYStart})/($param{RYStep} ? $param{RYStep}:1);
}
$ch=$CHRH*$d{font}{'font-size'}/10;
$rect{top}=$YOFFSET+$ch+$YSPACE*7;
$rect{width}=int($param{Width});
$rect{height}=int($param{Height});

$yMlen=$ryMlen=0;
if ($param{YNum} && $param{YExp}==1) {
	$yMlen=txtWidth($param{YNum});
	#print "$yMlen\n";
}elsif ($param{YExp}!=1) {
	$yMlen=txtWidth($param{YExp})+txtWidth($param{YNum} ? $param{YNum} : '9')*5/9;
}elsif ($param{YLog}) {
	$yMlen=txtWidth(int($param{YLog}**$param{YEnd}));
}elsif ($param{YNeedLog}) {
	$yMlen=(txtWidth(int($param{YEnd}))>txtWidth($param{YStep}) ? txtWidth(int($param{YEnd})) : txtWidth($param{YStep}));
}else{
	$yMlen=(txtWidth($param{YEnd})>txtWidth($param{YStep}) ? txtWidth($param{YEnd}) : txtWidth($param{YStep}));
	$yMlen=txtWidth($maxYScale) if ($yMlen<txtWidth($maxYScale));
	#print "($param{YEnd}>$param{YStep})\t$yMlen\n";
}
if ($param{BothYAxis}) {
	if ($param{RYNum} && $param{RYExp}==1) {
		$ryMlen=txtWidth($param{RYNum});
	}elsif ($param{RYExp}!=1) {
		$ryMlen=txtWidth($param{RYExp})+txtWidth($param{RYNum} ? $param{RYNum} : '9')*5/9;
	}elsif ($param{RYLog}) {
		$ryMlen=txtWidth(int($param{RYLog}**$param{RYEnd}));
	}else{
		$ryMlen=(txtWidth($param{RYEnd})>txtWidth($param{RYStep}) ? txtWidth($param{RYEnd}) : txtWidth($param{RYStep}));
		$ryMlen=txtWidth($maxRYScale) if ($ryMlen<txtWidth($maxRYScale));
	}
}
if ($autofit) {
	$rect{left}=$XOFFSET+$ch+$XSPACE;
	$rect{left}+=$yMlen;
	$rect{left}+=max((txtWidth($param{Note}) < ($rect{width}+$rect{left}*2) ? 0 : ((txtWidth($param{Note})-$rect{width})/2-$rect{left})),
				 (txtWidth($param{X}) < ($rect{width}+$rect{left}*2) ? 0 : ((txtWidth($param{X})-$rect{width})/2-$rect{left})));
}else{
	$rect{left}=$XOFFSET*5+$ch+$XSPACE*20;
}
if ($autofit) {
	$rect{right}=$XOFFSET*6;
	if ($param{BothYAxis}) {
		$rect{right}+=$ryMlen;
	}
	$rect{right}+=max((txtWidth($param{Note}) < ($rect{width}+$rect{right}*2) ? 0 : ((txtWidth($param{Note})-$rect{width})/2-$rect{right})),
				  (txtWidth($param{X}) < ($rect{width}+$rect{right}*2) ? 0 : ((txtWidth($param{X})-$rect{width})/2-$rect{right})));
	$rect{right}+=txtWidth($param{XEnd})/2+($param{HaveMore} ? txtWidth(' ') : 0);
}else{
	$rect{right}=$XOFFSET*10;
	if ($param{BothYAxis}) {
		$rect{right}+=$rect{left};
	}
}
$realCh=$ch*8/10;
$vbW=$rect{width}+$rect{left}+$rect{right};
$maxXScale=$param{XEnd} if ($maxXScale eq '');
$vbH=$rect{height}+$rect{top}+$ch+$ch*cos($param{XScaleRoate}/180*$PI)+txtWidth($maxXScale)*sin($param{XScaleRoate}/180*$PI);#+$YOFFSET;
#if ($param{XExp}!=1) {
#	$vbH+=$ch*$EXP_SCALE;
#}
$colWidth=$rect{width}/($param{XEnd}-$param{XStart})*$param{XUnit};		#��λ���
$xZero=$rect{width}*$param{XZeroPos}+$colWidth*$param{XScaleLinePos};
$yZero=$rect{height}*$param{YZeroPos};
$xDDigits=10**$param{XDDigits};
$yDDigits=10**$param{YDDigits};
if ($param{BothYAxis}) {
	$ryZero=$rect{height}*$param{RYZeroPos};
	$ryDDigits=10**$param{RYDDigits};
}

#####################################################################################################################
#	��ͼ
#####################################################################################################################

$g->b("svg","viewBox","0 0 ".($vbW*$param{WholeScale})." ".(($vbH+$ch*@note)*$param{WholeScale}),"width",$vbW*$param{WholeScale},"height",($vbH+($ch+$YSPACE)*@note)*$param{WholeScale});
$g->svgPrint("\n<Author>Li Shengting</Author>");
$g->svgPrint("\n<E-mail>lishengting\@genomics.org.cn</E-mail>");
$g->svgPrint("\n<Version>$ver</Version>");
my $drawer=getlogin()."@".(`hostname`);
chomp $drawer;
$g->svgPrint("\n<Drawer>$drawer</Drawer>");
$g->svgPrint("\n<Date>".(localtime())."</Date>");
$g->b("g","transform","scale($param{WholeScale})",'xml:space',"preserve");

$g->b("g","transform","translate(".$rect{left}.",".$rect{top}.")");

$offset=$param{OffsetPer} >0 ? $param{OffsetPer} :0;
$unitPer=$param{UnitPer} >0 ? $param{UnitPer} :1;

if ($param{Part}>0) {
	$part=$param{Part};
	$offset=$unitPer/$part if (!$offset);
}elsif ($type =~ /$type[0]/) {
	$part=1;
}elsif ($type =~ /$type[1]/) {
	$part=2;
	$offset=0.3;
	$unitPer=0.8
}else{
	$part=1;
}
if ($pp) {
	$d{rect}{fill}='#DDDDDD';
	$d{rect}{stroke}='none';
	$g->d("rect","xval",0,"yval",0,"width",$rect{width},"height",$rect{height},"style",style($d{rect}));
}
#$g->b('clipPath', 'id', 'RectClip');
	#$g->d("rect","xval",0,"yval",0,"width",$rect{width},"height",$rect{height});
#$g->e();  

#print "$part $offset\n";
my $i=0;
my @shape=@pShape;
while (!eof(F)) {
	while (<F>) {
		if (/\S/) {
			seek(F,-length,1);
			last;
		}
	}
	draw(($i*$offset+$param{MovePer})*$colWidth,$colWidth*($unitPer-$offset*($part-1)));
	$i++;
}
#####################################################################################################################
#	�����Ͻ�ע��
#####################################################################################################################
#print "Ok1!\n";
if ($pp) {
	$d{font}{fill}='#000000';
}
if (@mark2) {
	@tmp=sort {txtWidth($a)<=>txtWidth($b)} @mark2;
	my ($xx,$yy,$i,$ww);
	$tmp=txtWidth($tmp[$#tmp]);
	$ww=$tmp+$XSPACE*6;
	if ($param{Mark2Pos}=~/r/i) {
		$xx = $rect{width}-$SCALELEN-$ww*$MARK2_SCALE-$XSPACE;
	}else{
		$xx = $SCALELEN+$XSPACE;
	}
	if ($param{Mark2Pos}=~/b/i) {
		$yy = $rect{height} - (@tmp*($realCh+$YSPACE)+$YSPACE*6)*$MARK2_SCALE - $SCALELEN - $YSPACE;
	}else{
		$yy = $SCALELEN+$YSPACE;
	}
	$g->b("g","transform","translate(".$xx.",".$yy.") scale($MARK2_SCALE)");
	if ($param{Mark2Border}) {
		$d{rect}{fill}='none';
		$d{rect}{stroke}='black';
		$d{rect}{'stroke-width'}=3;
		$g->d("rect","xval",0,"yval",0,"width",$ww,"height",@tmp*($realCh+$YSPACE)+$YSPACE*6,"style",style($d{rect}));
	}

	$x=$XSPACE*3;
	for ($i=0;$i<@mark2;$i++) {
		$g->d("txtLM",$mark2[$i],"xval",txtWidth(' '),"yval",($i+0.5)*($realCh),"style",style($d{font}));
	}
	$g->e();
}

@tmp=sort {txtWidth($a->{name})<=>txtWidth($b->{name})} @mark;
if ($tmp[$#tmp]{name} ne '') {
	shift @tmp until ($tmp[0] ne '');
	my ($xx,$yy,$i,$j,$ww,$chNum);
	for ($i=$#mark;$i>=0;$i--) {
		if (!$mark[$i]{name}) {
			for ($j=$i;$j<$#mark;$j++) {
				$mark[$j]=$mark[$j+1];
			}
			$#mark--;
		}
	}
	$chNum=5;
	$tmp=txtWidth($tmp[$#tmp]{name});
	if ($param{MarkStyle}=~/^[hH]/) {
		$ww=0;
		foreach (@mark) {
			$ww +=txtWidth($_->{name})+txtWidth(' ' x ($chNum+2));
		}
		$ww+=$XSPACE*6;
	}else{
		$ww=$tmp+txtWidth(' ' x ($chNum+1))+$XSPACE*6;
	}
	if ($param{MarkPos}=~/r/i) {
		$xx = $rect{width}-$ww*$MARK_SCALE-$SCALELEN-$XSPACE;
	}else{
		$xx = $SCALELEN+$XSPACE;
	}
	if ($param{MarkPos}=~/b/i) {
		if ($param{MarkStyle}=~/^[hH]/) {
			$yy = $rect{height} - ($realCh+$YSPACE+$YSPACE*6)*$MARK_SCALE - $SCALELEN - $YSPACE;
		}else{
			$yy = $rect{height} - (@mark*($realCh+$YSPACE)+$YSPACE*6)*$MARK_SCALE - $SCALELEN - $YSPACE;
		}
	}else{
		$yy = $SCALELEN+$YSPACE;
	}
	$g->b("g","transform","translate(".$xx.",".$yy.") scale($MARK_SCALE)");
	if (!$param{MarkNoBorder}) {
		$d{rect}{fill}='none';
		$d{rect}{stroke}='black';
		$d{rect}{'stroke-width'}=3;
		#$y=$realCh+$YSPACE;
		if ($param{MarkStyle}=~/^[hH]/) {
			$g->d("rect","xval",0,"yval",0,"width",$ww,"height",$realCh+$YSPACE*6,"style",style($d{rect}));
		}else{
			$g->d("rect","xval",0,"yval",0,"width",$ww,"height",@mark*($realCh+$YSPACE)+$YSPACE*6,"style",style($d{rect}));
		}
	}

	$x=$XSPACE*3;
	for ($i=0;$i<@mark;$i++) {
		next if ($mark[$i]{name} eq '');
		if ($param{MarkStyle}=~/^[hH]/) {
			$y=$realCh/2;
			#print "$x\t$y\n";
			$g->d("txtLM",$mark[$i]{name},"xval",txtWidth(' ' x ($chNum+1))+$x,"yval",$y,"style",style($d{font}));
			$y+=$YSPACE*3;
			if ($type =~ /$type[2]/) { #Line
				$d{path}{stroke}=$mark[$i]{color};
				$d{path}{'stroke-width'}=3/$MARK_SCALE;
				$g->d("line","x1",$x+txtWidth(' ')/2,"y1",$y,"x2",$x+txtWidth(' ' x $chNum),"y2",$y,"style",style($d{path}));
			}elsif ($type =~ /$type[3]/) { #Point
				$d{path}{stroke}=$mark[$i]{color};
				$d{path}{'stroke-width'}=3/$MARK_SCALE;
				if (!$noConnect) {
					$g->d("line","x1",$x+txtWidth(' ')/2,"y1",$y,"x2",$x+txtWidth(' ' x $chNum),"y2",$y,"style",style($d{path}));
				}
				$d{$pShape[$i]}{fill}=$mark[$i]{color};
				$d{$pShape[$i]}{stroke}='none';
				if ($pShape[$i] eq 'circle') {
					$g->d("circle","cx",$x+txtWidth(' ' x $chNum)/2,'cy',$y,'r',$MARK_POINT_SIZE,"style",style($d{circle}));
				}else{
					$g->d("rect","xval",$x+txtWidth(' ' x $chNum)/2-$MARK_POINT_SIZE,
								 "yval",$y-$MARK_POINT_SIZE,
								 "width",$MARK_POINT_SIZE*2,
								 "height",$MARK_POINT_SIZE*2,
								 "style",style($d{rect}));
				}
			}else{
				$d{rect}{fill}=$mark[$i]{color};
				$g->d("rect","xval",$x+txtWidth(' ' x $chNum)/2-$realCh*3/8,"yval",$y-$realCh*3/8,"width",$realCh*3/4,"height",$realCh*3/4,"style",style($d{rect}));
			}
			$x+=txtWidth(' ' x ($chNum+2))+txtWidth($mark[$i]{name});
		}else{
			$y=($i+0.5)*$realCh;
			$g->d("txtLM",$mark[$i]{name},"xval",txtWidth(' ' x ($chNum+1)),"yval",$y,"style",style($d{font}));
			$y+=$YSPACE*3;
			if ($type =~ /$type[2]/ && !$param{Fill}) { #Line
				$d{path}{stroke}=$mark[$i]{color};
				$d{path}{'stroke-width'}=3/$MARK_SCALE;
				$g->d("line","x1",$x+txtWidth(' ')/2,"y1",$y,"x2",txtWidth(' ' x $chNum),"y2",$y,"style",style($d{path}));
			}elsif ($type =~ /$type[3]/) { #Point
				$d{path}{stroke}=$mark[$i]{color};
				$d{path}{'stroke-width'}=3/$MARK_SCALE;
				if (!$noConnect) {
					$g->d("line","x1",$x+txtWidth(' ')/2,"y1",$y,"x2",txtWidth(' ' x $chNum),"y2",$y,"style",style($d{path}));
				}
				$d{$pShape[$i]}{fill}=$mark[$i]{color};
				$d{$pShape[$i]}{stroke}='none';
				if ($pShape[$i] eq 'circle') {
					$g->d("circle","cx",$x+txtWidth(' ' x $chNum)/2,'cy',$y,'r',$MARK_POINT_SIZE,"style",style($d{circle}));
				}else{
					$g->d("rect","xval",$x+txtWidth(' ' x $chNum)/2-$MARK_POINT_SIZE,
								 "yval",$y-$MARK_POINT_SIZE,
								 "width",$MARK_POINT_SIZE*2,
								 "height",$MARK_POINT_SIZE*2,
								 "style",style($d{rect}));
				}
			}else{
				$d{rect}{fill}=$mark[$i]{color};
				$g->d("rect","xval",$x+txtWidth(' ' x $chNum)/2-$realCh*3/8,"yval",$y-$realCh*3/8,"width",$realCh*3/4,"height",$realCh*3/4,"style",style($d{rect}));
			}
		}
	}
	$g->e();
}
if ($pp) {
	$d{font}{fill}='#F0F000';
}
#####################################################################################################################
#	���������
#####################################################################################################################
#print "Ok2!\n";
$d{path}{stroke}='#000000';
$d{path}{fill}='none';
$d{path}{'stroke-width'}=3;

my ($xMark,$yMark,$ryMark,$num,$sNum);
my ($tmp1,$tmp2,$mvOff,$yZeroVal,$yStep,$ryZeroVal,$ryStep);
#############
#	Y
#############
$x= - $XSPACE*2;
$yZeroVal=$param{YZeroVal};
$yStep=$param{YStep};
$sNum=@y;
$num=0;
while (1) {
	if ($param{MultiY}) {
		$yStep=$y[0]{Step};
		$yZeroVal=$y[0]{ZeroVal};
		$yDiv=($y[0]{End}-$yZeroVal)/($yStep ? $yStep:1);
		error ("MultiY usage wrong!") if (!$yDiv);
	}
	for (my $i=0;$i<=4*$yDiv*$param{YScaleDiv};$i++) {
		foreach (($i/$param{YScaleDiv},-$i/$param{YScaleDiv})) {
			$y = $yZero+$_*($rect{height}/$yDiv);
			#print "$y\n";
			if (!($i % $param{YScaleDiv})) {	#�����
				$yMark=$_*$yStep+$yZeroVal;
				if ($param{YLog}) {
					$yMark=availDigit($param{YLog}**$yMark,$param{AvailDigit},1,$param{YLog});
				}
				if (!($i % $param{YDispStep})) {
					if ($param{YExp}!=1) {
						#print "$yMark\t$param{YExp}\n";
						if ($yMark) {
							$tmp=log($yMark)/log($param{YExp});
						}else{
							$tmp=0;
						}
						$yMark=rint($tmp,1/$yDDigits);# if (exists($param{YDDigits}));
						if ($param{YLog}) {
							$y=$yZero
								+((log($param{YExp}**$yMark)/log($param{YLog}))-$yZeroVal)/$yStep
								*($rect{height}/$yDiv);
						}else{
							$y=$yZero
								+($param{YExp}**$yMark-$yZeroVal)/$yStep
								*($rect{height}/$yDiv);
						}
						$y=cut($y);
						next if ($y<0);
						next if ($y>$rect{height});
						$y=$rect{height}-$y;
						#$yMark=$param{YExp}."<tspan dy=-$ch transform=\"scale($EXP_SCALE)\" >".$yMark."</tspan>";
						$d{font}{'font-size'}*=5/9;
						$g->setFontSize($d{font}{'font-size'});
						$mvOff=txtWidth($yMark);
						$g->d("txtRB",($yMark),
							  "xval",$x,"yval",$y,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
						#$y+=$ch*$EXP_SCALE;
						$d{font}{'font-size'}/=5/9;
						$g->setFontSize($d{font}{'font-size'});
						#print ";;;$y+$mvOff\n";
						$g->d("txtRM",$param{YExp},
							"xval",$x-$mvOff,"yval",$y,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
					}else{
						if ($param{YLog}) {
							$yMark=rint($yMark,1/$yDDigits) if (exists($param{YDDigits}));
							next if (!$yMark);
							$y=$yZero
								+((log($yMark)/log($param{YLog}))-$yZeroVal)/$yStep
								*($rect{height}/$yDiv);
						}elsif (exists($param{YDDigits})) {
							$yMark=rint($yMark,1/$yDDigits);
							$y=$yZero
								+($yMark-$yZeroVal)/$yStep
								*($rect{height}/$yDiv);
						}
						$y=cut($y);
						next if ($y<0);
						next if ($y>$rect{height});
						$y=$rect{height}-$y;
						if (@yScale && @{$yScale[0]}) {
							$yMark=$yScale[0][$i / $param{YScaleDiv}];
						}
						if ($param{MultiY}) {
							$tmp=$d{font}{fill};
							$d{font}{fill}=$y[0]{Color};
							$g->d("txtRM",($yMark),
								  "xval",$x,"yval",$y+($num-$sNum/2+0.5)*($realCh)-$YSPACE*3,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
							$d{font}{fill}=$tmp;
						}else{
							$g->d("txtRM",($yMark),
								  "xval",$x,"yval",$y-$YSPACE*3,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
						}
					}
				}else{
					$y=$rect{height}-$y;
				}
				$g->d("line","x1",0,"y1",fy($y),"x2",$SCALELEN,"y2",fy($y),"style",style($d{path}));
				if (!$param{BothYAxis}) {
					$g->d("line","x1",$rect{width},"y1",fy($y),"x2",$rect{width}-$SCALELEN,"y2",fy($y),"style",style($d{path}));
				}
			}else{								#С���
				if ($param{YLog}) {
					$tmp1=uint($_-1)*$yStep+$yZeroVal;
					$tmp1=availDigit($param{YLog}**$tmp1,$param{AvailDigit},1,$param{YLog});
					$tmp2=uint($_)*$yStep+$yZeroVal;
					$tmp2=availDigit($param{YLog}**$tmp2,$param{AvailDigit},1,$param{YLog});
					$tmp=($tmp2-$tmp1)/$param{YScaleDiv};
					$yMark=$tmp1+$tmp*($i % $param{YScaleDiv});
					if ($param{YExp}!=1) {
						$tmp1=rint(log($tmp1)/log($param{YExp}),1/$yDDigits);
						$tmp2=rint(log($tmp2)/log($param{YExp}),1/$yDDigits);
						$tmp=($param{YExp}**$tmp2-$param{YExp}**$tmp1)/$param{YScaleDiv};
						$yMark=$param{YExp}**$tmp1+$tmp*($i % $param{YScaleDiv});
					}
					#$tmp=$_*$yStep+$yZeroVal;
					#$y-=($tmp-log($yMark)/log(10))/$yStep*($rect{width}/$yDiv); #��λ�� ($rect{width}/$yDiv)/$yStep
					$y=$yZero
						+((log($yMark)/log($param{YLog}))-$yZeroVal)/$yStep
						*($rect{height}/$yDiv);
				}else{
					if ($param{YExp}!=1) {
						$tmp1=uint($_-1)*$yStep+$yZeroVal;
						$tmp2=uint($_)*$yStep+$yZeroVal;
						$tmp1=rint(log($tmp1)/log($param{YExp}),1/$yDDigits);
						$tmp2=rint(log($tmp2)/log($param{YExp}),1/$yDDigits);
						$tmp=($param{YExp}**$tmp2-$param{YExp}**$tmp1)/$param{YScaleDiv};
						$yMark=$param{YExp}**$tmp1+$tmp*($i % $param{YScaleDiv});
					}
				}
				$y=cut($y);
				next if ($y<0);
				next if ($y>$rect{height});
				$y=$rect{height}-$y;
				$g->d("line","x1",0,"y1",fy($y),"x2",$SCALELEN/2,"y2",fy($y),"style",style($d{path}));
				if (!$param{BothYAxis}) {
					$g->d("line","x1",$rect{width},"y1",fy($y),"x2",$rect{width}-$SCALELEN/2,"y2",fy($y),"style",style($d{path}));
				}
			}
			last if ($_==0);
		}
	}
	shift (@y);
	shift (@yScale);
	last if (@y==0);
	$num++;
}
#############
#	RY
#############
$x= $rect{width} + $XSPACE*2;
$ryZeroVal=$param{RYZeroVal};
$ryStep=$param{RYStep};
$sNum=@ry;
$num=0;
if ($param{BothYAxis}) {
	while (1) {
		if ($param{MultiRY}) {
			$ryStep=$ry[0]{Step};
			$ryZeroVal=$ry[0]{ZeroVal};
			$ryDiv=($ry[0]{End}-$ry[0]{Start})/($ryStep ? $ryStep:1);
			#print"$ryStep\t$ryZeroVal\t$ryDiv\n";
			error ("MultiRY usage wrong!") if (!$ryDiv);
		}
		for (my $i=0;$i<=4*$ryDiv*$param{RYScaleDiv};$i++) {
			foreach (($i/$param{RYScaleDiv},-$i/$param{RYScaleDiv})) {
				$ry = $ryZero+$_*($rect{height}/$ryDiv);
				#print "$ry\n";
				if (!($i % $param{RYScaleDiv})) {	#�����
					$ryMark=$_*$ryStep+$ryZeroVal;
					if ($param{RYLog}) {
						$ryMark=availDigit($param{RYLog}**$ryMark,$param{AvailDigit},1,$param{RYLog});
					}
					#print "$ryMark\n";
					if (!($i % $param{RYDispStep})) {
						if ($param{RYExp}!=1) {
							if ($ryMark) {
								$tmp=log($ryMark)/log($param{RYExp});
							}else{
								$tmp=0;
							}
							$ryMark=rint($tmp,1/$ryDDigits);# if (exists($param{RYDDigits}));
							if ($param{RYLog}) {
								$ry=$ryZero
									+((log($param{RYExp}**$ryMark)/log($param{RYLog}))-$ryZeroVal)/$ryStep
									*($rect{height}/$ryDiv);
							}else{
								$ry=$ryZero
									+($param{RYExp}**$ryMark-$ryZeroVal)/$ryStep
									*($rect{height}/$ryDiv);
							}
							$ry=cut($ry);
							next if ($ry<0);
							next if ($ry>$rect{height});
							$ry=$rect{height}-$ry;
							#$ryMark=$param{RYExp}."<tspan dy=-$ch transform=\"scale($EXP_SCALE)\" >".$ryMark."</tspan>";
							$mvOff=txtWidth($param{RYExp});
							$d{font}{'font-size'}*=5/9;
							$g->setFontSize($d{font}{'font-size'});
							$g->d("txtLB",($ryMark),
								  "xval",$x+$mvOff,"yval",$ry,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
							#$ry+=$ch*$EXP_SCALE;
							$d{font}{'font-size'}/=5/9;
							$g->setFontSize($d{font}{'font-size'});
							#print ";;;$ry+$mvOff\n";
							$g->d("txtLM",$param{RYExp},
								"xval",$x,"yval",$ry,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
						}else{
							if ($param{RYLog}) {
								$ryMark=rint($ryMark,1/$ryDDigits) if (exists($param{RYDDigits}));
								next if (!$ryMark);
								$ry=$ryZero
									+((log($ryMark)/log($param{RYLog}))-$ryZeroVal)/$ryStep
									*($rect{height}/$ryDiv);
							}elsif (exists($param{RYDDigits})) {
								$ryMark=rint($ryMark,1/$ryDDigits);
								$ry=$ryZero
									+($ryMark-$ryZeroVal)/$ryStep
									*($rect{height}/$ryDiv);
							}
							$ry=cut($ry);
							next if ($ry<0);
							next if ($ry>$rect{height});
							$ry=$rect{height}-$ry;
							if (@ryScale && @{$ryScale[0]}) {
								$ryMark=$ryScale[0][$i / $param{RYScaleDiv}];
							}
							if ($param{MultiRY}) {
								$tmp=$d{font}{fill};
								$d{font}{fill}=$ry[0]{Color};
								$g->d("txtLM",($ryMark),
									  "xval",$x,"yval",$ry+($num-$sNum/2+0.5)*($realCh)-$YSPACE*3,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
								$d{font}{fill}=$tmp;
							}else{
								$g->d("txtLM",($ryMark),
									  "xval",$x,"yval",$ry-$YSPACE*3,"style",style($d{font}));#,"onclick","alert('y:$tmp')");
							}
						}
					}else{
						$ry=$rect{height}-$ry;
					}
					$g->d("line","x1",$rect{width},"y1",fy($ry),"x2",$rect{width}-$SCALELEN,"y2",fy($ry),"style",style($d{path}));
				}else{								#С���
					if ($param{RYLog}) {
						$tmp1=uint($_-1)*$ryStep+$ryZeroVal;
						$tmp1=availDigit($param{RYLog}**$tmp1,$param{AvailDigit},1,$param{RYLog});
						$tmp2=uint($_)*$ryStep+$ryZeroVal;
						$tmp2=availDigit($param{RYLog}**$tmp2,$param{AvailDigit},1,$param{RYLog});
						$tmp=($tmp2-$tmp1)/$param{RYScaleDiv};
						$ryMark=$tmp1+$tmp*($i % $param{RYScaleDiv});
						if ($param{RYExp}!=1) {
							$tmp1=rint(log($tmp1)/log($param{RYExp}),1/$ryDDigits);
							$tmp2=rint(log($tmp2)/log($param{RYExp}),1/$ryDDigits);
							$tmp=($param{RYExp}**$tmp2-$param{RYExp}**$tmp1)/$param{RYScaleDiv};
							$ryMark=$param{RYExp}**$tmp1+$tmp*($i % $param{RYScaleDiv});
						}
						#$tmp=$_*$ryStep+$ryZeroVal;
						#$ry-=($tmp-log($ryMark)/log(10))/$ryStep*($rect{width}/$ryDiv); #��λ�� ($rect{width}/$ryDiv)/$ryStep
						$ry=$ryZero
							+((log($ryMark)/log($param{RYLog}))-$ryZeroVal)/$ryStep
							*($rect{height}/$ryDiv);
					}else{
						if ($param{RYExp}!=1) {
							$tmp1=uint($_-1)*$ryStep+$ryZeroVal;
							$tmp2=uint($_)*$ryStep+$ryZeroVal;
							$tmp1=rint(log($tmp1)/log($param{RYExp}),1/$ryDDigits);
							$tmp2=rint(log($tmp2)/log($param{RYExp}),1/$ryDDigits);
							$tmp=($param{RYExp}**$tmp2-$param{RYExp}**$tmp1)/$param{RYScaleDiv};
							$ryMark=$param{RYExp}**$tmp1+$tmp*($i % $param{RYScaleDiv});
						}
					}
					$ry=cut($ry);
					next if ($ry<0);
					next if ($ry>$rect{height});
					$ry=$rect{height}-$ry;
					$g->d("line","x1",$rect{width},"y1",fy($ry),"x2",$rect{width}-$SCALELEN/2,"y2",fy($ry),"style",style($d{path}));
				}
				last if ($_==0);
			}
		}
		shift (@ry);
		shift (@ryScale);
		last if (@ry==0);
		$num++;
	}
}

#############
#	X
#############
$y=$rect{height} + $YSPACE*2;
for (my $i=0;$i<=4*$xDiv*$param{XScaleDiv};$i++) {
	foreach (($i/$param{XScaleDiv},-$i/$param{XScaleDiv})) {
		$x = $xZero+$_*($rect{width}/$xDiv);
		#print "$_\t$x\n";
		if (!($i % $param{XScaleDiv})) {	#�����
			$xMark=$_*$param{XStep}+$param{XZeroVal};
			if ($param{XLog}) {
				$xMark=availDigit($param{XLog}**$xMark,$param{AvailDigit},1,$param{XLog});
			}
			if (!($i % $param{XDispStep})) {
				if ($param{XExp}!=1) {
					if ($xMark) {
						$tmp=log($xMark)/log($param{XExp});
					}else{
						$tmp=0;
					}
					$xMark=rint($tmp,1/$xDDigits);# if (exists($param{XDDigits}));
					if ($param{XLog}) {
						$x=$xZero
							+((log($param{XExp}**$xMark)/log($param{XLog}))-$param{XZeroVal})/$param{XStep}
							*($rect{width}/$xDiv);
					}else{
						$x=$xZero
							+($param{XExp}**$xMark-$param{XZeroVal})/$param{XStep}
							*($rect{width}/$xDiv);
					}
					$x=cut($x);
					next if ($x<0);
					next if ($x>$rect{width}-$param{MovePer}*$colWidth);
					#print "$xMark\n";
					#$xMark=$param{XExp}."<tspan dy=-$ch transform=\"scale($EXP_SCALE)\" >".$xMark."</tspan>";
					$mvOff=txtWidth('10');
					$d{font}{'font-size'}*=5/9;
					$g->setFontSize($d{font}{'font-size'});
					$mvOff-=txtWidth($xMark);
					$mvOff/=2;
					$g->d("txtLT",($xMark)
								 .((!$param{XZeroPos} && $param{HaveLess} && ($xZero+($_-1)*($rect{width}/$xDiv)) < 0) ? '-' : '')
								 .((!$param{XZeroPos} && $param{HaveMore} &&  ($xZero+($_+1)*($rect{width}/$xDiv)) > $rect{width}) ? '+' : ''),
						"xval",$x+$mvOff+$colWidth*$param{XScalePos},"yval",$y,"style",style($d{font}));#,"onclick","alert('x:$tmp')");
					#$y+=$ch*$EXP_SCALE;
					$d{font}{'font-size'}/=5/9;
					$g->setFontSize($d{font}{'font-size'});
					#print ";;;$x+$mvOff\n";
					$g->d("txtRT",$param{XExp},
						"xval",$x+$mvOff+$colWidth*$param{XScalePos},"yval",$y,"style",style($d{font}));#,"onclick","alert('x:$tmp')");
				}else{
					#print "$xMark\t";
					if ($param{XLog}) {
						$xMark=rint($xMark,1/$xDDigits) if (exists($param{XDDigits}));
						next if (!$xMark);
						$x=$xZero
							+((log($xMark)/log($param{XLog}))-$param{XZeroVal})/$param{XStep}
							*($rect{width}/$xDiv);
					}elsif (exists($param{XDDigits})) {
						$xMark=rint($xMark,1/$xDDigits);
						$x=$xZero
							+($xMark-$param{XZeroVal})/$param{XStep}
							*($rect{width}/$xDiv);
					}
					$x=cut($x);
					next if ($x<0);
					next if ($x>$rect{width}-$param{MovePer}*$colWidth);
					if (@xScale) {
						$xMark=$xScale[$i / $param{XScaleDiv}];
					}
					#print "$xMark\n";
					if ($param{XScaleRoate}) {
						$g->d("txtRM",($xMark)
									 .((!$param{XZeroPos} && $param{HaveLess} && ($xZero+($_-1)*($rect{width}/$xDiv)) < 0) ? '-' : '')
									 .((!$param{XZeroPos} && $param{HaveMore} &&  ($xZero+($_+1)*($rect{width}/$xDiv)) > $rect{width}) ? '+' : ''),
							"xval",0,"yval",0,"transform","translate(".($x+$colWidth*$param{XScalePos}).",$y) rotate(-$param{XScaleRoate})","style",style($d{font}));#,"onclick","alert('x:$tmp')");
						if (@group) {
							#print int($i / $param{XScaleDiv})."\t$xMark\n";
							$gxy[$i / $param{XScaleDiv}]{x}=($x+$colWidth*$param{XScalePos})-txtWidth("  ".$xMark)*cos($param{XScaleRoate}/180*$PI);
							$gxy[$i / $param{XScaleDiv}]{y}=$y+txtWidth("  ".$xMark)*sin($param{XScaleRoate}/180*$PI);
						}
					}else{
						$g->d("txtCT",($xMark)
									 .((!$param{XZeroPos} && $param{HaveLess} && ($xZero+($_-1)*($rect{width}/$xDiv)) < 0) ? '-' : '')
									 .((!$param{XZeroPos} && $param{HaveMore} &&  ($xZero+($_+1)*($rect{width}/$xDiv)) > $rect{width}) ? '+' : ''),
							"xval",$x+$colWidth*$param{XScalePos},"yval",$y,"style",style($d{font}));#,"onclick","alert('x:$tmp')");
						if (@group) {
							$gxy[$i / $param{XScaleDiv}]{x}=($x+$colWidth*$param{XScalePos});
							$gxy[$i / $param{XScaleDiv}]{y}=$y+$ch+$YSPACE*2;
						}
					}
				}
			}
			$g->d("line","x1",$x,"y1",$rect{height},"x2",$x,"y2",$rect{height}-$SCALELEN,"style",style($d{path}));#,"onclick","alert('x:$tmp')");
			$g->d("line","x1",$x,"y1",0,"x2",$x,"y2",$SCALELEN,"style",style($d{path}));#,"onclick","alert('x:$tmp')");
		}else{								#С���
			if ($param{XLog}) {
				$tmp1=uint($_-1)*$param{XStep}+$param{XZeroVal};
				$tmp1=availDigit($param{XLog}**$tmp1,$param{AvailDigit},1,$param{XLog});
				$tmp2=uint($_)*$param{XStep}+$param{XZeroVal};
				$tmp2=availDigit($param{XLog}**$tmp2,$param{AvailDigit},1,$param{XLog});
				$tmp=($tmp2-$tmp1)/$param{XScaleDiv};
				$xMark=$tmp1+$tmp*($i % $param{XScaleDiv});
				if ($param{XExp}!=1) {
					$tmp1=rint(log($tmp1)/log($param{XExp}),1/$xDDigits);
					$tmp2=rint(log($tmp2)/log($param{XExp}),1/$xDDigits);
					$tmp=($param{XExp}**$tmp2-$param{XExp}**$tmp1)/$param{XScaleDiv};
					$xMark=$param{XExp}**$tmp1+$tmp*($i % $param{XScaleDiv});
				}
				#$tmp=$_*$param{XStep}+$param{XZeroVal};
				#$x-=($tmp-log($xMark)/log(10))/$param{XStep}*($rect{width}/$xDiv); #��λ�� ($rect{width}/$xDiv)/$param{XStep}
				$x=$xZero
					+((log($xMark)/log($param{XLog}))-$param{XZeroVal})/$param{XStep}
					*($rect{width}/$xDiv);
			}else{
				if ($param{XExp}!=1 && $_>0) {
					$tmp1=uint($_-1)*$param{XStep}+$param{XZeroVal};
					$tmp2=uint($_)*$param{XStep}+$param{XZeroVal};
					$tmp1=rint(log($tmp1)/log($param{XExp}),1/$xDDigits);
					$tmp2=rint(log($tmp2)/log($param{XExp}),1/$xDDigits);
					$tmp=($param{XExp}**$tmp2-$param{XExp}**$tmp1)/$param{XScaleDiv};
					$xMark=$param{XExp}**$tmp1+$tmp*($i % $param{XScaleDiv});
				}
			}
			$x=cut($x);
			next if ($x<0);
			next if ($x>$rect{width});
			$g->d("line","x1",$x,"y1",$rect{height},"x2",$x,"y2",$rect{height}-$SCALELEN/2,"style",style($d{path}));#,"onclick","alert('x:$tmp')");
			$g->d("line","x1",$x,"y1",0,"x2",$x,"y2",$SCALELEN/2,"style",style($d{path}));#,"onclick","alert('x:$tmp')");
		}
		last if ($_==0);
	}
}

if (@group) {
	my (@tmp,$angle,$maxY,$x,$y,$lastx,$lasty,$i);
	$angle=$param{XScaleRoate} ? $param{XScaleRoate} : 90;
	@tmp=sort {$a->{y}<=>$b->{y}} @gxy;
	$maxY=$tmp[$#tmp]{y};
	if ($pp) {
		$d{path}{stroke}='#F0F000';
	}
	$d{rect}{'stroke-width'}=5;
	#print map {"=>$_->{x}\t$_->{y}<=\n"} @gxy;
	foreach (@group) {
		$lastx=$gxy[0]{x}-($maxY-$gxy[0]{y})/sin($param{XScaleRoate}/180*$PI)*cos($param{XScaleRoate}/180*$PI);
		$lasty=$maxY;
		$g->d("line","x1",fx($gxy[0]{x}),"y1",fy($gxy[0]{y}),"x2",fx($lastx),"y2",fy($lasty),"style",style($d{path}));
		for ($i=0;$i<$_->{len}-1;$i++) {
			last if ($#gxy<=0);
			shift(@gxy);
		}
		$x=$gxy[0]{x}-($maxY-$gxy[0]{y})/sin($param{XScaleRoate}/180*$PI)*cos($param{XScaleRoate}/180*$PI);
		$y=$maxY;
		$g->d("line","x1",fx($gxy[0]{x}),"y1",fy($gxy[0]{y}),"x2",fx($x),"y2",fy($y),"style",style($d{path}));
		$g->d("line","x1",fx($lastx),"y1",fy($lasty),"x2",fx($x),"y2",fy($y),"style",style($d{path}));
		$g->d("txtCT",$_->{mark},
			  "xval",fx(($lastx+$x)/2),
			  "yval",fy($maxY+$YSPACE),
			  "style",style($d{font}),
			 );
		shift(@gxy) if ($#gxy>0);
	}
	$d{rect}{'stroke-width'}=3;
}

$d{rect}{fill}='none';
$d{rect}{stroke}='black';
$d{rect}{'stroke-width'}=3;
$g->d("rect","xval",0,"yval",0,"width",$rect{width},"height",$rect{height},"style",style($d{rect}));
$g->e();

if (@note) {
	$g->b("g","transform","translate(".$XOFFSET.",".($vbH+$YOFFSET).") scale(0.95,1)");
	$d{font}{'font-family'}='ArialNarrow-Bold';
	for (my $i=0;$i<@note;$i++) {
		$g->d('txtLT',$note[$i],'xval',0,'yval',$i*($ch+$YSPACE),"style",style($d{font}));
	}
	$g->e();
	$d{font}{'font-family'}='ArialNarrow-Bold';
}

$x=$rect{left}-$XSPACE*8;
$x-=$yMlen;
$y=$rect{height}/2+$rect{top};
$g->d("txtCB",$param{Y},"xval",0,"yval",0,"transform",
	"translate(".$x.",".$y.") rotate(-90)", "style",style($d{font}));
$g->d("txtCB",$param{X},"xval",($rect{width}/2+$rect{left}),"yval",$vbH-$YOFFSET*1.5,"style",style($d{font}));
$g->d("txtCT",$param{Note},"xval",($rect{width}/2+$rect{left}),"yval",$YOFFSET, "style",style($d{font}));
if ($param{BothYAxis}) {
	$g->d("txtCT",$param{RY},"xval",0,"yval",0,"transform",
		"translate("
		.($rect{left}+$rect{width}+$ryMlen+$XSPACE*1)
		.",".($rect{height}/2+$rect{top})
		.") rotate(-90)", "style",style($d{font}));
}
#print "($x,$y)\t$yMlen\t$ryMlen\n";
$g->e();
$g->e();
$svg->close($g);
close(F);

#####################################################################################################################
#	��
#####################################################################################################################

sub draw{
	my($wOff,$w)=@_;
	my($color,$h,$lastx,$lasty,$lst,$x,$y,$ry,$y2,$tmpx,$tmpy,
	   $ltmpx,$ltmpy,$i,$tmp,$tmp1,$tmp2,$tmp3,$zeroy,$dotNum,$maxDotNum,$maxDot,
	   $shape,$ryAxis,$ryMark,$tmpScale,$maxLenScale,
	   @other,@tmp,%draw,%xy,@x,%y,%ry,@ys,%dots,$dot,$overflow,
	   $drawType,$noconnect);
	my ($hex1,$hex2,$hex3);
	my @lkeys=(
		"Type","Color","Mark","YAxis","YMark","YHasLow",
		"Start","End","Step","ZeroVal","ZeroPos",
		"ColorStep","LightColor","SmoothColor","Transparence",
		"NoLine","NoConnect","LineDash","LineWidth",		
		#Line
		"Fill",
		#Rect
		"NoFill",
		#Point
	);

	$tmpScale=[];
	while (<F>) {
		if ($_=~/^\s*[\+-]?\d/ || $_!~/\S/) {
			seek(F,-length,1);
			last;
		}
		if (/^\s*Scale:/i) {
			while (<F>) {
				last if ($_=~/^:End$/i);
				chomp;
				push(@{$tmpScale},$_);
				$maxLenScale=$_ if (length($maxLenScale)<length($_));
			}
		}
		if (/(\S+?):(.+)/) {
			$draw{lc($1)}=$2;
		}
	}
	foreach (@lkeys) {
		if (exists($draw{lc($_)})) {
			$draw{$_}=$draw{lc($_)};
		}
	}

	$draw{Transparence}=$param{Transparence} if (!exists($draw{Transparence}) && $param{Transparence});
	$draw{LineDash}=$param{LineDash} if (!exists($draw{LineDash}) && $param{LineDash});
	$drawType=$draw{Type} ? $draw{Type} : $type ;
	$color=$draw{Color} ? $draw{Color} : '#000000' ;
	$noconnect=$draw{NoLine} || $draw{NoConnect} || $noConnect;
	$#mark++;
	$mark[$#mark]{name}=$draw{Mark};
	$mark[$#mark]{color}=$color;
	if (!exists($draw{YMark})) {
		$draw{YMark} = $draw{YAxis};
	}
	$ryAxis=($draw{YAxis}=~/r/i) ? 1 : 0;
	$ryMark=($draw{YMark}=~/r/i) ? 1 : 0;

	if ($param{MultiY} && !$ryMark) {
		$y{Start}=exists($draw{Start}) ? $draw{Start} : $raw{YStart};
		$y{End}=exists($draw{End}) ? $draw{End} : $raw{YEnd};
		$y{Step}=exists($draw{Step}) ? $draw{Step} : $raw{YStep};
		$y{Color}=$color;
		error ("YStart can't equal with YEnd!") if ($y{Start}==$y{End});
		error ("YStep must bigger than 0!") if ($y{Step}<=0);
		if ($param{YDiv}) {
			if ($param{YLog} && !$param{YNeedLog}) {
				$y{Start}-=log($param{YDiv})/log($param{YLog});
				$y{End}-=log($param{YDiv})/log($param{YLog});
			}else{
				$y{Start}/=$param{YDiv};
				$y{End}/=$param{YDiv};
				$y{Step}/=$param{YDiv};
			}	
		}
		if ($param{YNeedLog}) {
			error ("YStart,YEnd and YStep must bigger than 0!") if ($y{Start}<=0 || $y{End}<=0 || $y{Step}<=0);
			$step = rint(($y{End}-$y{Start})/$y{Step});
			$y{Start}=log($y{Start})/log($param{YNeedLog});
			$y{End}=log($y{End})/log($param{YNeedLog});
			$y{Step}=availDigit(($y{End}-$y{Start})/$step,$param{AvailDigit},1,$param{YLog});
		}
		if ($y{Start}!=0 && !exists($draw{ZeroPos}) && !exists($raw{YZeroPos}) 
						 && !exists($draw{ZeroVal}) && !exists($raw{YZeroVal})) {
			$y{ZeroVal}=$y{Start};
		}else{
			$y{ZeroVal}=exists($draw{ZeroVal}) ? $draw{ZeroVal} : $raw{YZeroVal};
		}
		$y{ZeroPos}=exists($draw{ZeroPos}) ? $draw{ZeroPos} : $raw{YZeroPos};
		push(@y,\%y);
		push(@yScale,$tmpScale) ;#if (@{$tmpScale});
		$yMlen=txtWidth($maxLenScale) if ($yMlen<txtWidth($maxLenScale));
	}
	if ($param{MultiRY} && $ryMark) {
		$ry{Start}=exists($draw{Start}) ? $draw{Start} : $raw{RYStart};
		$ry{End}=exists($draw{End}) ? $draw{End} : $raw{RYEnd};
		$ry{Step}=exists($draw{Step}) ? $draw{Step} : $raw{RYStep};
		$ry{Color}=$color;
		error ("RYStart can't equal with RYEnd!") if ($ry{Start}==$ry{End});
		error ("RYStep must bigger than 0!") if ($ry{Step}<=0);
		if ($param{RYDiv}) {
			if ($param{RYLog} && !$param{RYNeedLog}) {
				$ry{Start}-=log($param{RYDiv})/log($param{RYLog});
				$ry{End}-=log($param{RYDiv})/log($param{RYLog});
			}else{
				$ry{Start}/=$param{RYDiv};
				$ry{End}/=$param{RYDiv};
				$ry{Step}/=$param{RYDiv};
			}
		}
		if ($param{RYNeedLog}) {
			error ("RYStart,RYEnd and RYStep must bigger than 0!") if ($ry{Start}<=0 || $ry{End}<=0 || $ry{Step}<=0);
			$step = rint(($ry{End}-$ry{Start})/$ry{Step});
			$ry{Start}=log($ry{Start})/log($param{RYNeedLog});
			$ry{End}=log($ry{End})/log($param{RYNeedLog});
			$ry{Step}=availDigit(($ry{End}-$ry{Start})/$step,$param{AvailDigit},1,$param{RYLog});
			#print "$step\t$ry{Start}\t$ry{End}\n";
		}
		if ($ry{Start}!=0 && !exists($draw{ZeroPos}) && !exists($raw{RYZeroPos})
						  && !exists($draw{ZeroVal}) && !exists($raw{RYZeroVal})) {
			$ry{ZeroVal}=$ry{Start};
		}else{
			$ry{ZeroVal}=exists($draw{ZeroVal}) ? $draw{ZeroVal} : $raw{RYZeroVal};
		}
		$ry{ZeroPos}=exists($draw{ZeroPos}) ? $draw{ZeroPos} : $raw{RYZeroPos};
		push(@ry,\%ry);
		push(@ryScale,$tmpScale) ;#if (@{$tmpScale});
		$ryMlen=txtWidth($maxLenScale) if ($ryMlen<txtWidth($maxLenScale));
	}

	%xy=();
	$tmpx=$tmpy=0;
	while (<F>) {
		last if ($_!~/\S/);
		$_=~s/\s//g;
		if (/([^:]+):([^:]+):?(\S*)/) {
			$x=$1;
			$y=$2;
			#print "$x,$y\n";
			@other=split(/:/,$3);
			$y2=@other[0];
			#print @other;
			#print "\n";
			if ($param{XDiv}) {
				if ($param{XLog} && !$param{XNeedLog}) {
					$x-=log($param{XDiv})/log($param{XLog});
				}else{
					$x/=$param{XDiv};
				}
			}
			if ($param{XNeedLog}) {
				$x=log($x)/log($param{XNeedLog});
			}
			#print "$x:$y\n";
			if ($ryAxis) {
				if ($param{RYDiv}) {
					if ($param{RYLog} && !$param{RYNeedLog}) {
						$y-=log($param{RYDiv})/log($param{RYLog});
						$y2-=log($param{RYDiv})/log($param{RYLog}) if ($y2);
					}else{
						$y/=$param{RYDiv};
						$y2/=$param{RYDiv};
					}
				}
				if ($param{RYNeedLog}) {
					$y=log($y)/log($param{RYNeedLog}) if ($y);
					$y2=log($y2)/log($param{RYNeedLog}) if ($y2);
				}
				$zeroy=$param{RYZeroVal};
			}else{
				if ($param{YDiv}) {
					if ($param{YLog} && !$param{YNeedLog}) {
						$y-=log($param{YDiv})/log($param{YLog});
						$y2-=log($param{YDiv})/log($param{YLog}) if ($y2);
					}else{
						$y/=$param{YDiv};
						$y2/=$param{YDiv};
					}
				}
				if ($param{YNeedLog}) {
					$y=log($y)/log($param{YNeedLog}) if ($y);
					$y2=log($y2)/log($param{YNeedLog}) if ($y2);
				}
				$zeroy=$param{YZeroVal};
			}
			$x=cut($x);
			$y=cut($y);
			if ($draw{YHasLow} || $param{YHasLow}) {
				$zeroy=cut($y2);
			}
			#print "$x\t$y\t$param{XZeroPos}\n";
			#��(С)��$param{XEnd}($param{XStart})����ݴ���
			$tmp=$param{XUnit}*uint($param{XScaleLinePos},1);
			$tmp1=cut($param{XEnd})-($param{XLog} ? 0 : ($tmp?$tmp:$param{XUnit}));
			#print "$x>$tmp1\n" if (cut($x)>cut($tmp1));
			if ($param{XZeroPos}) {
				$xy{$x}{$y}{n}++;
				$xy{$x}{$y}{low}=$zeroy;
			}elsif (cut($x)>cut($tmp1)) {
				#print "$x\t$param{XEnd}\t$tmp?$tmp:$param{XUnit}\t";
				#print cut($param{XEnd})."\n";
				$tmpy+=$y;
				if ($x>cut($param{XEnd}) && !$param{XCut}) {
					$param{HaveMore}=1;
				}
				#print "$x++++++$y\n";
			}elsif (cut($x)<cut($param{XStart})) {
				$tmpx+=$y;
				if (!$param{XCut}) {
					$param{HaveLess}=1;
				}
			}else{
				$xy{$x}{$y}{n}++;
				$xy{$x}{$y}{low}=$zeroy;
				#print "$x\t$y\n";
			}
		}
	}


	foreach $x (keys %xy) {
		foreach $y (keys %{$xy{$x}}) {
			$maxDotNum=$xy{$x}{$y}{n} if ($maxDotNum < $xy{$x}{$y}{n});
			$dots{$xy{$x}{$y}{n}}=1;
		}
	}
	$maxDot=0;
	foreach (keys %dots) {
		$maxDot++;
	}
	if ($draw{SmoothColor} eq '') {
		$maxDot=$maxDotNum;
	}
	
	if ($draw{ColorStep} eq 'auto') {
		$draw{LightColor}='#FFFFFF';
	}
	if ($draw{LightColor} ne '' && $maxDot>1) {
		$hex1=(hex(substr($draw{LightColor},1,2))-hex(substr($color,1,2)));
		$hex2=(hex(substr($draw{LightColor},3,2))-hex(substr($color,3,2)));
		$hex3=(hex(substr($draw{LightColor},5,2))-hex(substr($color,5,2)));
		$hex1=$hex1/($maxDot-1);
		$hex2=$hex2/($maxDot-1);
		$hex3=$hex3/($maxDot-1);
		$draw{ColorStep}='auto';
		#print "$maxDot\t$hex1\t$hex2\t$hex3\n";
	}

	if ($draw{ColorStep} ne '' && $draw{ColorStep} ne 'auto') {
		$overflow=0;
		$tmp=(hex(substr($color,1))+hex(substr($draw{ColorStep},1))*($maxDot-1));
		if ($tmp>=hex('FFFFFF')) {
			error("Dot color too light!",1);
			$overflow=1;
		}
	}
	
	if ($param{HaveMore}) {
		$tmp=cut($param{XEnd})-($tmp?$tmp:$param{XUnit});
		@ys=sort {$a<=>$b} keys %{$xy{$tmp}};
		#print "$tmp\n";
		delete($xy{$tmp});
		$xy{$tmp}{$tmpy+$ys[$#ys]}{n}++;
		$xy{$tmp}{$tmpy+$ys[$#ys]}{low}=$zeroy;
		#print "$tmp[$#tmp]\t$tmpy\n";
	}
	if ($param{HaveLess}) {
		$tmp=cut($param{XStart});
		@ys=sort {$a<=>$b} keys %{$xy{$tmp}};
		delete($xy{$param{XStart}});
		$xy{$tmp}{$tmpx+$ys[0]}{n}++;
		$xy{$tmp}{$tmpx+$ys[0]}{low}=$zeroy;
	}

	$d{path}{stroke}=$color;
	$d{path}{fill}='none';
	$d{path}{'stroke-width'}=$draw{LineWidth} ? $draw{LineWidth} : ($param{LineWidth} ? $param{LineWidth} : 3);
	$d{path}{'stroke-linecap'}='round';
	$d{path}{'stroke-dasharray'}=$draw{LineDash} if ($draw{LineDash} ne '');
	if ($drawType =~ /$type[2]/) {	#Line
		$d{polygon}{fill}=$color;
		$d{polygon}{stroke}=$color;
	}elsif ($drawType =~ /$type[3]/) {	#Point
		$shape=shift(@shape);
		$d{$shape}{fill}=$color;
		$d{$shape}{stroke}='none';
	}else{
		if ($draw{NoFill}) {
			$d{rect}{fill}='none';
			$d{rect}{stroke}=$color;
			$d{rect}{'stroke-width'}=$param{LineWidth} ? $param{LineWidth} :1;
		}else{
			$d{rect}{fill}=$color;
			$d{rect}{stroke}='none';
		}
	}
	$lst=1;
	$tmpx=$tmpy=0;
	if ($draw{Transparence}) {
		$g->b("g","opacity",1-$draw{Transparence});
	}
	@x=sort {$a<=>$b} keys %xy;
	$dot=0;
	for ($dotNum=1;$dotNum<=$maxDotNum;$dotNum++) {
		next if (!$dots{$dotNum});
		if ($draw{SmoothColor}) {
			$dot++;
		}else{
			$dot=$dotNum;
		}
		#print "$dotNum\t$dot\n";
		foreach (@x) {
			$x = $rect{width}*($_-$param{XZeroVal})/($param{XEnd}-$param{XStart})+$wOff+$rect{width}*$param{XZeroPos};
			#print "$_\t$x\n";
			#@ys=sort {$xy{$_}{$a}{n}<=>$xy{$_}{$b}{n}} keys %{$xy{$_}};
			foreach $i (keys %{$xy{$_}}) {
				next if ($xy{$_}{$i}{n}!=$dotNum);
				#print "$_\t$i\n";
				#print "$x = $rect{width}*$_/($param{XEnd}-$param{XStart})+$rect{width}*$param{XZeroPos}+$wOff\n";
				if ($ryAxis) {
					if ($param{MultiRY} && $ryMark) {
						$tmp = $rect{height}*($i-$ry{ZeroVal})/($ry{End}-$ry{Start})+$rect{height}*$ry{ZeroPos};
						$h = $rect{height}*(($i-$xy{$_}{$i}{low}))/($ry{End}-$ry{Start})+$rect{height}*$ry{ZeroPos};
					}else{
						$tmp = $rect{height}*($i-$param{RYZeroVal})/($param{RYEnd}-$param{RYStart})+$rect{height}*$param{RYZeroPos};
						$h = $rect{height}*(($i-$xy{$_}{$i}{low}))/($param{RYEnd}-$param{RYStart})+$rect{height}*$param{RYZeroPos};
					}
				}else{
					if ($param{MultiY} && !$ryMark) {
						$tmp = $rect{height}*($i-$y{ZeroVal})/($y{End}-$y{Start})+$rect{height}*$y{ZeroPos};
						$h = $rect{height}*(($i-$xy{$_}{$i}{low}))/($y{End}-$y{Start})+$rect{height}*$y{ZeroPos};
					}else{
						$tmp = $rect{height}*($i-$param{YZeroVal})/($param{YEnd}-$param{YStart})+$rect{height}*$param{YZeroPos};
						$h = $rect{height}*(($i-$xy{$_}{$i}{low}))/($param{YEnd}-$param{YStart})+$rect{height}*$param{YZeroPos};
					}
					#print "\n$h = $rect{height}*(($i-$xy{$_}{$i}{low})-$param{RYZeroVal})/($param{RYEnd}-$param{RYStart})+$rect{height}*$param{RYZeroPos}\n";
				}
				#print "$tmp\t$h\n= $rect{height}*($i-$param{YZeroVal})/($param{YEnd}-$param{YStart})+$rect{height}*$param{YZeroPos}\n";
				#next if ($h==0);
				$y = $rect{height}-$tmp;#-$d{rect}{'stroke-width'}/2;
				if ($param{YCut}) {
					next if ($y<0 || $y>$rect{height} || $h<0 || $h>$rect{height});
				}
				$tmpx= ($param{XLog} ? $param{XLog}**$_ : $_);
				if ($ryAxis) {
					$tmpy= ($param{RYLog} ? $param{RYLog}**$i : $i);
				}else{
					$tmpy= ($param{YLog} ? $param{YLog}**$i : $i);
				}
				if ($drawType =~ /$type[2]/) { #Line
					if ($lst) {
						$lst=0;
					}else{
						if ($param{RightAngle}) {	#����
								$g->d("line","x1",fx($lastx),"y1",fy($lasty),"x2",fx($x),"y2",fy($lasty),"style",style($d{path}),"onclick","alert('x:$ltmpx\\ny:$ltmpy')");
							if (!$noconnect) {
								$g->d("line","x1",fx($x),"y1",fy($lasty),"x2",fx($x),"y2",fy($y),"style",style($d{path}));
							}
						}elsif ($draw{Fill} || $param{Fill}) {
							$g->d("polygon","points",fx($lastx)." ".fy($lasty).",".fx($x)." ".fy($y).",".fx($x)." ".fy($rect{height}).",".fx($lastx)." ".fy($rect{height}),
								"style",style($d{polygon}));#,"onclick","alert('".$msg."');");
						}else{
							if (!$noconnect) {
								$g->d("line","x1",fx($lastx),"y1",fy($lasty),"x2",fx($x),"y2",fy($y),"style",style($d{path}));
							}
						}
					}
					#print "$x\t$y\n";
					if (!$w || $drawType =~ /$type[4]/) {
						$lastx=$x;
						#$y+=$d{rect}{'stroke-width'}/2;
					}else{
						$g->d("line","x1",fx($x),"y1",fy($y),"x2",fx($x+$w),"y2",fy($y),"style",style($d{path}),"onclick","alert('x:$tmpx\\ny:$tmpy')");
						$lastx=$x+$w;
					}
					$lasty=$y;
					$ltmpx=$tmpx;
					$ltmpy=$tmpy;
					#print "------>$lastx\t$lasty\n";
				}elsif ($drawType =~ /$type[3]/) { #Point
					if ($lst) {
						$lst=0;
					}elsif ($param{VerticalLine}) {
						#print "$lastx\t$x\n";
						if ($lastx==$x) {
							$g->d("line","x1",fx($lastx),"y1",fy($lasty),"x2",fx($x),"y2",fy($y),"style",style($d{path})); #����
						}
					}elsif (!$noconnect) {
						$g->d("line","x1",fx($lastx),"y1",fy($lasty),"x2",fx($x),"y2",fy($y),"style",style($d{path})); #����
					}

					if ($draw{ColorStep} ne '') {
						if ($draw{ColorStep} eq 'auto') {
							$tmp1=int($hex1*($maxDot-$dot)+0.5);
							$tmp2=int($hex2*($maxDot-$dot)+0.5);
							$tmp3=int($hex3*($maxDot-$dot)+0.5);
							$tmp=(hex(substr($color,1))+($tmp1*hex('10000')+$tmp2*hex('100')+$tmp3));
						}else{
							$tmp=(hex(substr($color,1))+hex(substr($draw{ColorStep},1))*($maxDot-$dot));
						}
						#print substr($color,1)."\n";
						#print sprintf("%X",$tmp>=0 ? $tmp : 0);
						$d{$shape}{'fill'}='#'.sprintf("%06X",$tmp>=0 ? ($tmp<=hex('FFFFFF')?$tmp:hex('FFFFFF')) : 0);
					}
					if ($shape eq 'circle') {
						$g->d("circle","cx",fx($x),'cy',fy($y),'r',$param{PointSize},"style",style($d{circle}),"onclick","alert('x:$tmpx\\ny:$tmpy".($maxDot>1 ? "\\ndot:$dotNum" : '')."')");
					}else{
						$g->d("rect","xval",fx($x)-$param{PointSize},
									 "yval",fy($y)-$param{PointSize},
									 "width",$param{PointSize}*2,
									 "height",$param{PointSize}*2,
									 "style",style($d{rect}),
									 #"opacity",1-$param{Transparence},
									 "onclick","alert('x:$tmpx\\ny:$tmpy".($maxDotNum>1 ? "\\ndot:$dotNum" : '')."')",
							 );
					}
					$lastx=$x;
					$lasty=$y;
					$ltmpx=$tmpx;
					$ltmpy=$tmpy;
					#print "------>$lastx\t$lasty\n";
				}else{	#Rectangle
					#print "$x\t$w\n";
					#print "$x,$y,$w,$h\n";
					$g->d("rect","xval",fx($x),"yval",fy($y),"width",fx($w),"height",fy($h),"style",style($d{rect}),"onclick","alert('x:$tmpx\\ny:$tmpy')");
				}
			}
		}
	}
	if ($draw{Transparence}) {
		$g->e();
	}
	delete($d{path}{'stroke-dasharray'});
}

sub error{
	my($str,$ret)=@_;
	warn "\n|~|~|~|~|~|~|~|~|~|~|~|~|~|~|\n";
	warn "$str\n";
	warn "|_|_|_|_|_|_|_|_|_|_|_|_|_|_|\n";
	die "\n" if (!$ret);
	print "\n";
	return $str;
}

sub style{
	my($style)=@_;
	my $str="";
	for (keys %$style) {
		$str.="$_:$style->{$_};";
	}
	chop($str);
	return $str;
}

sub availDigit{
	my ($num,$n,$rounding,$exp)=@_;
	$exp=10 if (!$exp);
	#print "$num\n";
	my $log10=log($num)/log($exp);
	$n=1 if (!$n);
	$log10 = $log10 >0 ? int($log10) : int($log10)==$log10 ? $log10 : int($log10-1);
	return int($num*$exp**($n-1)/$exp**$log10+0.5*$rounding)*$exp**$log10/$exp**($n-1);
}

sub uint{
	my($num,$div)=@_;
	$div = 1 if (!$div);
	my $tmp=$num/$div;
	if ($tmp!=oint($tmp)) {
		$tmp=oint($tmp+1);
	}
	$tmp*=$div;
	return $tmp;
}

sub rint{
	my($num,$dot)=@_;
	if ($dot) {
		$num/=$dot;
		return $num>=0 ? int($num+0.5)*$dot : int($num-0.5)*$dot;
	}else{
		return $num>=0 ? int($num+0.5) : int($num-0.5);
	}
}

sub oint{
	my($num)=@_;
	if ($num>0) {
		return int($num);
	}else{
		return int($num)==$num ? $num : int($num)-1;
	}
}

sub fx{
	my ($num)=@_;
	#return uint($num,1);
	return sprintf("%f",$num);
}

sub fy{
	my ($num)=@_;
	#return uint($num,1);
	return sprintf("%f",$num);
}

sub max{
	my($m1,$m2)=@_;
	$m1 > $m2 ? $m1 : $m2;
}

sub min{
	my($m1,$m2)=@_;
	$m1 < $m2 ? $m1 : $m2;
}

sub cut{
	my ($num)=@_;
	return int($num*100000000+0.5)/100000000;
}

sub txtWidth{
	my($str)=@_;
	return $g->textWidth($d{font}{'font-size'},$g->{charspacing},$g->{wordspacing},$g->{hscaling},$str);
}



}
1;

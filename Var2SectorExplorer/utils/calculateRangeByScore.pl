use strict;
use warnings;
# score file format: "","pbmc.Celltype","pbmc.Groups","pbmc.Lightharvest1"
#"WT2_TATCTGTGTAGGCAAC","Epidermal cell","WT",3.05009022860266


my $scorefile=shift;
my $window =shift;

my ($n,$w,$m)=(0,0,0);
my $sum=0;my $wc=0;
open IN, $scorefile;
<IN>;
print "AverageScore($window)\twildRatio\tmutantRatio\n";
while(<IN>){
	chomp;
	$_=~s/\"//g;
	my @t=split(/\t/,$_);
	$n++;
	if($n % $window == 0){
		my $ave=sprintf("%.3f",$sum/$window);
		my $wratio=sprintf("%.3f",$w/($w+$m));
		my $mratio=1-$wratio;
		print "Bin$wc\t$ave\t$wratio\t$mratio\n";
		$sum=0; $m=0; $w=0;
		$wc++;
	}
	else{
		$sum+=$t[-1];
		if($t[4] eq "WT"){
			$w++;
		}
		if($t[4] eq "mutant"){
			$m++;
		}
	}
}
my $number=$n  % $window;
my $ave=sprintf("%.3f",$sum/$number);
my $wratio=sprintf("%.3f",$w/($w+$m));
my $mratio=1-$wratio;
print "Bin$wc\t$ave\t$wratio\t$mratio\n";

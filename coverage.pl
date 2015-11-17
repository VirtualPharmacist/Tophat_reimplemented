use diagnostics;
open FILE2, "./$ARGV[0]/coverage.pileup";

open FILE1, "window.bed";

open FILE3, ">./$ARGV[0]/histogram.txt";


$window = $ARGV[1];


while($line1 = <FILE1>){
	$line1 =~ s/\n//;
	$line1 =~ s/\R//;
	($chr,$start,$_) = split(/\t/,$line1);
	$hash1{$chr."\t".$start} = $line1;
	$hash2{$line1} = 0;
	}
while($line2 =<FILE2>){
	$line2 =~ s/\n//;
	$line2 =~ s/\R//;
	($c,$posi,$_,$count)	= split(/\t/,$line2);
	$tmp = int($posi / $window);
	$tmp = $tmp * $window + 1;
	
	$hash2{$hash1{$c."\t".$tmp}} += $count if(exists($hash1{$c."\t".$tmp}));
	}

foreach $key (sort keys %hash2){
	$value = $hash2{$key} / $window;
	print FILE3 "$key\t$value\n";
	}	


	
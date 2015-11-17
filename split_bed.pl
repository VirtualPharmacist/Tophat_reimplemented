
open FILE1, "chromosome_length.txt";

open FILE2, ">window.bed";

$window = $ARGV[0];

while($line1 = <FILE1>){
	$line1 =~ s/\n//;
	$line1 =~ s/\R//;
	$i = 0;
	($chr,$length) = split(/\t/,$line1);
	$label = $length / $window;
	$label = int $label;
	while($i < $label){
		$start = $i * $window + 1;
		$end = $window * ($i + 1);
		print FILE2 "$chr\t$start\t$end\n";
		$i++;
		}
		if($i == $label){
			$start = $i * $window + 1;
			$end = $length;
			print FILE2 "$chr\t$start\t$end\n";
			}
		
	}

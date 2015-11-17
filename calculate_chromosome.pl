open FILE1, "/1_disk/public_resources/hg19.fa";

open FILE2, ">chromosome_length.txt";

while($line1 = <FILE1>){
	$line1 =~ s/\n//;
	$line1 =~ s/\R//;
	if($line1 =~ m/>/){
		$line1 =~ s/>//;
		$chr = $line1;
		$hash{$chr} = 0;
		next;
		
		}
		$hash{$chr} = $hash{$chr}	+ length($line1);
	}
	
foreach $key (sort keys %hash){
	print FILE2 "$key\t$hash{$key}\n";
	}
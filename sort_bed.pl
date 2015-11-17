use diagnostics;
open FILE1, "./$ARGV[0]/histogram.txt";

open FILE2, ">./$ARGV[0]/histogram_sort.txt";




$tmp = "";
while($line1 = <FILE1>){
	$line1 =~ s/\n//;
	$line1 =~ s/\R//;
	($chr,$start,$end,$count) = split(/\t/,$line1);
	if($tmp ne $chr){
		if($tmp ne ""){
			foreach $key (sort {$a <=> $b} (keys %hash)){
				$tmp =~ s/chr/hs/;
				print FILE2 "$tmp\t$hash{$key}\n";
				}

			$tmp = $chr;
			undef(%hash);						
			}
		if($tmp eq ""){
			$tmp = $chr;
			undef(%hash);	
			}

					
		}
	$hash{$start} = "$start\t$end\t$count";		
	if(eof(FILE1)){
		foreach $key (sort {$a <=> $b} (keys %hash)){
			$tmp =~ s/chr/hs/;
			print FILE2 "$tmp\t$hash{$key}\n";
			}			
			
		}		
	
	}
	



	
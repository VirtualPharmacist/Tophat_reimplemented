use diagnostics;
open FILE1, "./$ARGV[0]/final_output.sam";

open FILE2, ">./$ARGV[0]/read_statistics";

$all = 0;
while($line1 = <FILE1>){
	$line1 =~ s/\n//;
	$line1 =~ s/\R//;
	next if($line1 =~ m/^@/);
	@a =split(/\t/,$line1);
	
	
	$chrom = $a[2];
	$posi = $a[3];
	next if($a[2] eq '*');
	
	$len = length($a[9]);
	print FILE2 "$chrom\t$posi\t$len\n";
	$all += $len;
	}
open FILE3, ">./$ARGV[0]/all_read";
$all = $all / 100;
print FILE3 "$all\n";
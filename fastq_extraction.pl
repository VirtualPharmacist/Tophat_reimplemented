use diagnostics;
open FILE1, "$ARGV[0]";
print "$ARGV[0]\t$ARGV[1]unalign1.fq\n";

open FILE21, ">$ARGV[1]unalign1.fq";
open FILE22, ">$ARGV[1]unalign2.fq";
open FILE23, ">$ARGV[1]unalign3.fq";
open FILE24, ">$ARGV[1]unalign4.fq";


$tmp = 0;	
while($line1 = <FILE1>){
	$tmp++;
	$line1 =~ s/\n//;
	$line1 =~ s/\R//;
	@a = split(/\t/,$line1);
	$reads = $a[9];
	$qual = $a[10];	
	$r = $a[0];
	$r1 = substr($reads,0,25);
	$r2 = substr($reads,25,25);
	$r3 = substr($reads,50,25);
	$r4 = substr($reads,75,25);			
	$h1 = '@'."$r|$tmp|1";
	$h2 = '@'."$r|$tmp|2";
	$h3 = '@'."$r|$tmp|3";
	$h4 = '@'."$r|$tmp|4";
	$q1 = substr($qual,0,25);
	$q2 = substr($qual,25,25);
	$q3 = substr($qual,50,25);
	$q4 = substr($qual,75,25);		
	print FILE21 "$h1\n$r1\n".'+'."\n$q1\n";
	print FILE22 "$h2\n$r2\n".'+'."\n$q2\n";
	print FILE23 "$h3\n$r3\n".'+'."\n$q3\n";
	print FILE24 "$h4\n$r4\n".'+'."\n$q4\n";
	}

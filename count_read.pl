use diagnostics;
open FILE1, "./$ARGV[0]/final_output.sam";
open FILE2, ">./$ARGV[0]/all_read.txt";

$count = 0;
while($line1 = <FILE1>){
	$line1 =~ s/\R//;
	$line1 =~ s/\n//;
	$line1 =~ s/\"//g;
	next if($line1 =~ m/^@/);
	($read,$f,$chr,$_) = split(/\t/,$line1);
	next if($chr eq '*');
	$count++ if($chr ne '*');
	
	}

print FILE2 "$count\n"
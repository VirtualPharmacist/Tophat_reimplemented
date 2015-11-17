use diagnostics;
open FILE1, "./hg19_exon_unique2.gtf";


%len = ();
############read two hash table
while($line1 = <FILE1>){
	$line1 =~ s/\R//;
	$line1 =~ s/\n//;
	$line1 =~ s/\"//g;
	$line1 =~ s/gene_id//;
	my($char, $start,$end,$gene) = split(/\t/,$line1);
	#print "$gene\n";	
	$result{$gene} = 0;
	$len{$gene} = $end - $start if(!exists($len{$gene}));
	$len{$gene} = $len{$gene} + $end -$start if(exists($len{$gene}));
	for(my $i = $start;$i <= $end;$i++){
		#$count{$char.':'.$i} = 0;
		$name{$char.':'.$i} = $gene;
		}
	
	}

open FILE2, "./$ARGV[0]/final.sam";

while($line2 = <FILE2>){
	$line2 =~ s/\R//;
	$line2 =~ s/\n//;
	next if($line2 =~ m/^@/);
	
	@bb = split(/\t/,$line2);
	$char = $bb[2];
	$posi = $bb[3];
	$i = $posi;
	next if($char eq "*");
	
	if(exists($name{$char.':'.$i})){
		$result{$name{$char.':'.$i}} = 1 if(!exists($result{$name{$char.':'.$i}}));		
		$result{$name{$char.':'.$i}} = $result{$name{$char.':'.$i}} + 1 if(exists($result{$name{$char.':'.$i}}));
			
		}		
		

	$all_read++;
	}



open FILE3, ">./$ARGV[0]/out_gtf_new_new.gtf";

foreach $key (keys %result){
	$rpkm = 1000000000*$result{$key}/($all_read*$len{$key});
	
	
	print FILE3 "$key\t$rpkm\n";
	
	}

	
	
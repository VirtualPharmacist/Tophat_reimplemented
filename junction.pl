use diagnostics;

######load whole genome coordinate;

open FILE, '/1_disk/public_resources/hg19.fa';
@test = <FILE>;	
close(FILE);

######################3

################build the chromosome index
open FILE2, '/1_disk/rhcheng/RNA_project/tophat_v1/coordinate.txt';
$num =0;	
while($line = <FILE2>){		
	$line =~ s/\R//;
	$line =~ s/\n//;
	
	@array = split(/\t/,$line);
	$hash{$array[0]} = $num+1;
	#print "$array[1]\n";
	$num += $array[1];
	
	}
close(FILE2);
########################

############# get the chromosome sequence at specific position
sub cut{
	my ($seed,$chrom) = @_;
	#print "$hash{$chrom}\n";
	my $breath = 50;
	#print "1:$seed\t2:$breath\t3:$chrom\t4:$hash{$chrom}\n";
	my $l1 = int($seed / $breath)+$hash{$chrom};
	my $l2 = $seed % $breath;
	my $t1 = $test[$l1-1];
	$t1 =~ s/\n//;
	$t1 =~ s/\R//;
	
	my $t2 = $test[$l1];
	$t2 =~ s/\n//;
	$t2 =~ s/\R//;
	
	my $t3 = $test[$l1+1];
	$t3 =~ s/\n//;
	$t3 =~ s/\R//;
	
	
	my $f = $t1.$t2.$t3;
	
	my $str = uc(substr($f,$l2+49,25));
	
	
	return $str;		
	}
##############################
$intron_length = 150000;
$leg = "AS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:25\tYT:Z:UU";
$mid = "*\t0\t0";
$mapq = 42;
sub check{
	my($str1,$str2,$tmp_seg) = @_;
	#$leg = "AS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:25\tYT:Z:UU";
	#$head1 = "$h1.1\t16";
	#$head2 = "$h1.2\t16";
	#$mid = "*\t0\t0";
	#$mapq = 42;
	
	
	#print "$str1\t$str2\t$tmp_seg\t$chrom\n";
	$match =0;
	$mismatch = 0;
	$threshold = 1;
	$limit = 5;
	for(my $i=0;$i<length($str1);$i++){
		$s1 = substr($str1,$i,1);
		$s2 = substr($tmp_seg,$i,1);
		if($s1 eq $s2){
			$match++;
			next;
			
			}
		if($s1 ne $s2 && $mismatch<$threshold && $match <= $limit){
			$mismatch++;
			$match++;
			next;			
			}
		if($s1 ne $s2 && $mismatch == $threshold){
			last;
			
			}
								
		}
	if($match != 0){
		$final_seg1 = substr($tmp_seg,0,$match);
		#$final_qual = 'I' x $match;
		#$des = $match.'M';
		#$posi = $p1;
		#print FILE3 "$head1\t$chrom\t$posi\t$mapq\t$des\t$mid\t$final_seg\t$final_qual\t$leg\n";
		
		}	
	###############这个是反着找
	$match =0;
	
	
	$tmp_seg = reverse($tmp_seg);
	$str2_1	= reverse($str2);	
	for(my $i=0;$i<length($str2_1);$i++){
		$s1 = substr($str2_1,$i,1);
		$s2 = substr($tmp_seg,$i,1);
		if($s1 eq $s2){
			$match++;
			next;
			
			}
		if($s1 ne $s2 && $mismatch<$threshold && $match <= $limit){
			$mismatch++;
			$match++;
			next;			
			}
		if($s1 ne $s2 && $mismatch == $threshold){
			last;
			
			}
								
		}
	if($match != 0){
		$final_seg2 = reverse(substr($tmp_seg,0,$match));
		#$final_qual = 'I' x $match;
		#$des = $match.'M';
		#$posi = $p2 - $match;
		#print FILE3 "$head2\t$chrom\t$posi\t$mapq\t$des\t$mid\t$final_seg\t$final_qual\t$leg\n";
		
		}
	if($final_seg2 && $final_seg1){
		return ($final_seg1,$final_seg2);
		
		}		
	}



open F1, "$ARGV[0]/unalign1.sam";
open F2, "$ARGV[0]/unalign2.sam";
open F3, "$ARGV[0]/unalign3.sam";
open F4, "$ARGV[0]/unalign4.sam";
open FILE3, ">$ARGV[0]/junction_output_new_test_mis.sam";

while($line1 = <F1>){
	$line2 = <F2>;
	$line3 = <F3>;
	$line4 = <F4>;
	next if($line1 =~ m/^@/);
	$line1 =~ s/\R//;
	$line1 =~ s/\n//;
	$line2 =~ s/\R//;
	$line2 =~ s/\n//;
	$line3 =~ s/\R//;
	$line3 =~ s/\n//;
	$line4 =~ s/\R//;
	$line4 =~ s/\n//;			
	@a1 = split(/\t/,$line1);		#$a1[2] $a1[3] position $a1[9] read
	@a2 = split(/\t/,$line2);
	@a3 = split(/\t/,$line3);
	@a4 = split(/\t/,$line4);
	
	#print FILE3 "$line1\n" if($a1[2] ne "*");
	#print FILE3 "$line2\n" if($a2[2] ne "*");
	#print FILE3 "$line3\n" if($a3[2] ne "*");
	#print FILE3 "$line4\n" if($a4[2] ne "*");
	#print "$chrom_a[0]\t$chrom_a[1]\t$chrom_a[2]\t$chrom_a[3]\n";
	
		if($a1[2] ne '*'&& $a2[2] eq '*' && $a1[2] eq $a3[2] && $a4[2] eq $a3[2] && abs($a3[3]-$a1[3]) <= $intron_length){	#seg1,3,4 有map，并且在同一个染色体的的同一个位置


			if($a3[3] - $a4[3] == -25){	#在正链的情况,暂时先不考虑 insertion deletion 
				
				$seed1 = $a1[3]+25;
				$seed2 = $a3[3]-25;
				#print "$seed1\t$seed2\n";
				$str1 = cut($seed1,$a1[2]);
				$str2 = cut($seed2,$a1[2]);				
				$ss1 = "";
				$ss2 = "";

				($ss1,$ss2) = check($str1,$str2,$a2[9]);
				if($ss1 && $ss2){

					$num_n = $a3[3] - $a1[3] - length($a1[9]) - length($ss1) - length($ss2);
					$cigar = length($a1[9].$ss1)."M".$num_n."N".length($ss2.$a3[9].$a4[9])."M";
					$flag = 0;
					$tt_seg = $a1[9].$ss1.$ss2.$a3[9].$a4[9];
					$tt_qual = 'I' x length($tt_seg);
					print FILE3 "$a1[0]\t$flag\t$a1[2]\t$a1[3]\t$mapq\t$cigar\t$mid\t$tt_seg\t$tt_qual\t$leg\n" if($num_n>0);
					
					}
				}

			if($a3[3] - $a4[3] == 25){	#在负链的情况,暂时先不考虑 insertion deletion 
				
				$seed1 =  $a3[3]+25;
				$seed2 = $a1[3]-25;
				$str1 = cut($seed1,$a1[2]);
				$str2 = cut($seed2,$a1[2]);
				$ss1 = "";
				$ss2 = "";									
				$reverse_seg = $a2[9];
				$reverse_seg =~ tr/AGCT/TCGA/;
				$reverse_seg = reverse($reverse_seg);
				
				($ss1,$ss2) = check($str1,$str2,$reverse_seg);
				if($ss1 && $ss2){
					
					$num_n = $a1[3] - $a3[3] - length($a3[9]) - length($ss1) - length($ss2);
					$cigar = length($a4[9].$a3[9].$ss1)."M".$num_n."N".length($ss2.$a1[9])."M";
					$flag = 16;
					$tt_seg = $a4[9].$a3[9].$ss1.$ss2.$a1[9];
					$tt_qual = 'I' x length($tt_seg);
					print FILE3 "$a4[0]\t$flag\t$a4[2]\t$a4[3]\t$mapq\t$cigar\t$mid\t$tt_seg\t$tt_qual\t$leg\n" if($num_n>0);
					
					}				
				}



			}

		if($a2[2] ne '*'&& $a3[2] eq '*' && $a2[2] eq $a4[2] && $a1[2] eq $a2[2] && abs($a2[3]-$a4[3]) <= $intron_length){	#seg1,2,4 有map，并且在同一个染色体的的同一个位置
			#print "$posi_a[1]\n";
			if($a1[3] - $a2[3] == -25){	#在正链的情况,暂时先不考虑 insertion deletion 
				
				$seed1 =  $a2[3]+25;
				$seed2 = $a4[3]-25;
				$str1 = cut($seed1,$a2[2]);
				$str2 = cut($seed2,$a2[2]);					
				$ss1 = "";
				$ss2 = "";				

				
				($ss1,$ss2) = check($str1,$str2,$a3[9]);
				if($ss1 && $ss2){

					$num_n = $a4[3] - $a2[3] - length($a2[9]) - length($ss1) - length($ss2);
					
					$cigar = length($a1[9].$a2[9].$ss1)."M".$num_n."N".length($ss2.$a4[9])."M";
					$flag = 0;
					$tt_seg = $a1[9].$a2[9].$ss1.$ss2.$a4[9];
					$tt_qual = 'I' x length($tt_seg);
					print FILE3 "$a1[0]\t$flag\t$a1[2]\t$a1[3]\t$mapq\t$cigar\t$mid\t$tt_seg\t$tt_qual\t$leg\n" if($num_n>0);
					
					}
				}
			if($a1[3] - $a2[3] == 25){	#在负链的情况,暂时先不考虑 insertion deletion 
				
				$seed1 =  $a4[3]+25;
				$seed2 = $a2[3]-25;
				$str1 = cut($seed1,$a2[2]);
				$str2 = cut($seed2,$a2[2]);					
				$ss1 = "";
				$ss2 = "";
				$reverse_seg = $a3[9];
				$reverse_seg =~ tr/AGCT/TCGA/;
				$reverse_seg = reverse($reverse_seg);


				($ss1,$ss2) = check($str1,$str2,$reverse_seg);
				if($ss1 && $ss2){
					
					$num_n = $a2[3] - $a4[3] - length($a4[9]) - length($ss1) - length($ss2);
					
					$cigar = length($a4[9].$ss1)."M".$num_n."N".length($ss2.$a2[9].$a1[9])."M";
					$flag = 16;
					$tt_seg = $a4[9].$ss1.$ss2.$a2[9].$a1[9];
					$tt_qual = 'I' x length($tt_seg);
					print FILE3 "$a4[0]\t$flag\t$a4[2]\t$a4[3]\t$mapq\t$cigar\t$mid\t$tt_seg\t$tt_qual\t$leg\n" if($num_n>0);
					
					}				
				}
	
			#print "$str1\t$str2\n";	
			#print "$seg_a[2]\n";
			
			}
			
		if($a1[2] eq $a2[2] && $a2[2] eq $a3[2] && $a3[2] eq $a4[2] && $a1[1] == $a2[1] && $a2[1] == $a3[1] && $a3[1] == $a4[1]){	#正好切断的状态
			my $d12 = $a1[3] - $a2[3];
			my $d13 = $a1[3] - $a3[3];
			my $d14 = $a1[3] - $a4[3];
			$num_n = 0;
	
			
			if($d12 < 0 && $a1[1] == 0){	#正链的情况

				$tt_seg = $a1[9].$a2[9].$a3[9].$a4[9];
				$tt_qual = 'I' x length($tt_seg);
				
				if($d12 < -25 && $d12 >= $intron_length && $d13 - $d14 == -25 && $d12 - $d13 == -25){ #25 + 75 正链
					
					$num_n = $a2[3] - $a1[3] - 25;
					
					$cigar = length($a1[9])."M".$num_n."N".length($a2[9].$a3[9].$a4[9])."M";
				
					
					
					}
								
				if($d12 == -25 && $d13 < -50 && $d13 >= $intron_length && $d13 - $d14 == -25){ #50 + 50 正链

					$num_n = $a3[3] - $a2[3] - 25;
					
					$cigar = length($a1[9].$a2[9])."M".$num_n."N".length($a3[9].$a4[9])."M";
					
					
					}
	
				if($d12 == -25 && $d13 == -50 && $d14 >= $intron_length && $d14 < -75){ #75 + 25 正链
					
					$num_n = $a4[3] - $a3[3] - 25;
					
					$cigar = length($a1[9].$a2[9].$a3[9])."M".$num_n."N".length($a4[9])."M";
					
					}
									
				print FILE3 "$a1[0]\t$a1[1]\t$a1[2]\t$a1[3]\t$mapq\t$cigar\t$mid\t$tt_seg\t$tt_qual\t$leg\n" if($num_n>0);
				
				}

			if($d12 > 0 && $a1[1] == 16){	#负链的情况

				$tt_seg = $a4[9].$a3[9].$a2[9].$a1[9];
				$tt_qual = 'I' x length($tt_seg);
				
				if($d12 > 25 && $d12 <= $intron_length && $d13 - $d14 == 25 && $d12 - $d13 == 25){ #25 + 75 负链 
					
					$num_n = $a1[3] - $a2[3] - 25;
					
					$cigar = length($a4[9].$a3[9].$a2[9])."M".$num_n."N".length($a1[9])."M";
				
					
					
					}
								
				if($d12 == 25 && $d13 > 50 && $d13 <= $intron_length && $d13 - $d14 == 25){ #50 + 50 负链

					$num_n = $a2[3] - $a3[3] - 25;
					
					$cigar = length($a4[9].$a3[9])."M".$num_n."N".length($a2[9].$a1[9])."M";
					
					
					}
	
				if($d12 == 25 && $d13 == 50 && $d14 <= $intron_length && $d14 > 75){ #75 + 25 负链
					
					$num_n = $a3[3] - $a4[3] - 25;
					
					$cigar = length($a4[9])."M".$num_n."N".length($a3[9].$a2[9].$a1[9])."M";
					
					}
									
				print FILE3 "$a4[0]\t$a4[1]\t$a4[2]\t$a4[3]\t$mapq\t$cigar\t$mid\t$tt_seg\t$tt_qual\t$leg\n" if($num_n>0);
				
				}

						
			}	
	}



=pod
foreach $val (sort keys %hash){
	print "$val\t$hash{$val}\n";
	}

=cut
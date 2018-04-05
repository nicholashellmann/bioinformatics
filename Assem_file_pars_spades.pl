#!/usr/bin/perl

system "ls assemb_spades_k\*.txt > list_assem_out";# captures all the different outputs from the assemblethon program
open (INFILE, "list_assem_out");
open (OUTFILE, ">table_assem_out.txt") or die;# creates a output file to print the table we want to analyze in excel to

while ($line_assem_out = <INFILE>){

	$line_assem_out =~ s/\s+$//;# removes the space that follows eac of the output files from the assemblethon program
	push @out_assem, $line_assem_out;# adds the 'corrected' file on to the end of the array each time it iterates through the while loop 
}


foreach $out_file (@out_assem){# values in the array are the files
	open (INFILE, "$out_file") or die;# so the files need to be opened to get at the data
	$out_file =~ m/assemb_spades_k(\d+).txt/;# in order to define the keys in the hash before I assign a vaule to the hash table I need to make 	
	$kmer = $1;# a key match to the kmer value by matching the number part of out_assem_50X_.txt
	while ($inside_out_file = <INFILE>){# need to loop through each line of out_assem_50X_.txt file that I just opened
		if ($inside_out_file =~ m/Total size of scaffolds\s+(\d+)/){# the () will store the number following the "Total size of scaffolds" as a scalar called $1			
			$genome_hash{$kmer} = $1;# assigns a value to the keys, kmer, of the hash. I don't have to change the scalar created from the match because it is in a loop
		}
		if ($inside_out_file =~ m/N50 contig length\s+(\d+)/){		
			$N50_hash{$kmer} = $1;
		}
		if ($inside_out_file =~ m/Number of scaffolds\s+(\d+)/){			
			$Num_scaf_hash{$kmer} = $1;
		}	
	}
	close INFILE;# frees up the file handle to be used again
}

foreach $kmer (sort {$a <=> $b} keys %genome_hash){# because $kmer is the key for all the hashes I can loop through the keys for all the hashes with the $kmer key
	print OUTFILE "$kmer\t$genome_hash{$kmer}\t";# prints the table to the output file I created at the begining
	print OUTFILE "$N50_hash{$kmer}\t";
	print OUTFILE "$Num_scaf_hash{$kmer}\n";
}


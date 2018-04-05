#!/usr/bin/perl

open(INFILE, "prodigal_ref_K12_genes.fasta") or die;# using ARGV[0] allows us to test multiple gene files 

while ($line = <INFILE>){
	if ($line =~ m/^>/){# will match the first line of a gene in the FASTA format
		$count_genes++;}# adds 1 to the count when the if statement is sastified for each line of the INFILE
}

for ($ID = 1; $ID <= 32; $ID++){

		system "perl assemblathon.pl $ID.scaffolds.fasta > $ID.assemblathon.txt";
	
	open (INFILE1, "$ID.assemblathon.txt") or die;
	while ($line_assem = <INFILE1>){
		if ($line_assem =~ m/Total size of scaffolds\s+(\d+)/){# will match the assembly size
			$assem_size_hash{$ID} = $1;}
		if ($line_assem =~ m/Number of contigs\s+(\d+)/){			
			$Num_contigs_hash{$ID} = $1;}
		if ($line_assem =~ m/Number of scaffolds\s+(\d+)/){			
			$Num_scaf_hash{$ID} = $1;}
		if ($line_assem =~ m/N50 contig length\s+(\d+)/){		
			$N50_contig_hash{$ID} = $1;}
		if ($line_assem =~ m/N50 scaffold length\s+(\d+)/){		
			$N50_scaffold_hash{$ID} = $1;}
	}	
	
	system "prodigal -i $ID.contigs.fasta -d $ID.prodigal_gene.txt";
	open (INFILE2, "$ID.prodigal_gene.txt")or die;
	while ($line_prodigal = <INFILE2>){
			if ($line_prodigal =~ m/^>/){# will match the first line of a gene in the FASTA format
				$Num_pred_genes{$ID}++;}# adds 1 to the count when the if statement is sastified for each line of the INFILE
			if ($line_prodigal =~ m/partial=00/){# will match each line that has a complete gene not truncated	
				$Num_complete_genes{$ID}++;}
	}
#everything up till here works, as in prints to the table

	system "nucmer $ID.contigs.fasta K12_MG1655_ref_genome.fna";# performs the alignment
	system "mv out.delta $ID.out.delta";# need to do this so out.delta that is created by default is not overwritten each time
	system "delta-filter -1 $ID.out.delta > $ID.filtered.delta";# restricts output to non-overlapping hits
	system "show-snps $ID.filtered.delta > $ID.filtered.snps";# creates a table of SNPs and indels in aligned regions
	
	open (INFILE3, "$ID.filtered.snps" ) or die;
	<INFILE3>; 
	<INFILE3>;
	<INFILE3>;
	<INFILE3>;# skips the first few lines of the table

	$indel_count{$ID} = 0;
	$snps_count{$ID} = 0;

	while ($lines_file_snps = <INFILE3>){
		@array_snps = split /\s+/, $lines_file_snps;
		if ($array_snps[2] =~ m/\./){# matches the first column of the sub columns searching for '.', the symbol for an indel
			$indel_count{$ID}++}
		elsif ($array_snps[3] =~ m/\./){# matches the second column of the sub columns searching for '.', the symbol for an indel
			$indel_count{$ID}++}
		else {$snps_count{$ID}++} # if the line doesn't have a '.' it will be counted as a SNP
	}

	system "show-coords $ID.filtered.delta > $ID.filtered.coords";# create table of aligned regions	
	%contig_hash = (); # initializes for each genome	
	
	open (INFILE4, "$ID.contigs.fasta") or die;
	while ($line_contigs = <INFILE4>){		
		if ($line_contigs =~ m/^>(\w+)/){
			$contig_ID = $1;
			$contig_hash{$contig_ID} = 0;}# the keys of the hash will be every contig in every contig file, the values will all be initially 0
	}

	open (INFILE5, "$ID.filtered.coords") or die;
	<INFILE5>; 
	<INFILE5>;
	<INFILE5>;
	<INFILE5>;# skips the first few lines of the table

	while ($lines_unaligned = <INFILE5>){
		@array_prct = split /\s+/, $lines_unaligned;# splits the coords file into an array
		$contig_hash{$array_prct[13]}++;# adds one to the value of $contig_hash when the key matches array position 12 in the coords file
	}
	
	$cannot_align{$ID} = 0;
	$misassem{$ID} = 0;
	
	foreach $contig_ID (keys %contig_hash){
		if ($contig_hash{$contig_ID} == 0){# when no match is found between the coords file and the keys of the hash then the contig hasn't been matched to the reference genome
			$cannot_align{$ID}++;}# add one to the count for the number of contigs not aligning to the ref genome
		if ($contig_hash{$contig_ID} > 1){# if the coords file has more than one match of a single contig to the keys of the hash then that means that there are misassemblies
			$misassem{$ID} = $misassem{$ID} + ($contig_hash{$contig_ID} - 1);}# the first match between the coords file and the keys of the hash for a single contig is a correct alignment but every match after that is a misassembly
	}
}

open (OUTFILE, ">table_assign_4.txt") or die;			
foreach $ID (sort {$a <=> $b} keys %assem_size_hash){# because $ID is the key for all the hashes I can loop through the keys for all the hashes with the $kmer key
	print OUTFILE "$ID\t$assem_size_hash{$ID}\t";# prints the table to the output file I created at the begining
	print OUTFILE (($assem_size_hash{$ID} / 4641652) * 100), "%\t";# 4641652 is the size of the MG1655 genome
	print OUTFILE "$Num_contigs_hash{$ID}\t";
	print OUTFILE "$Num_scaf_hash{$ID}\t";
	print OUTFILE "$N50_contig_hash{$ID}\t";
	print OUTFILE "$N50_scaffold_hash{$ID}\t";
	print OUTFILE "$Num_pred_genes{$ID}\t";
	print OUTFILE (($Num_pred_genes{$ID}/ $count_genes) * 100), "%\t";
	print OUTFILE (($Num_pred_genes{$ID}/ $Num_complete_genes{$ID}) *100), "%\t";	
	print OUTFILE "$snps_count{$ID}\t";
	print OUTFILE "$indel_count{$ID}\t";
	print OUTFILE "$cannot_align{$ID}\t";
	print OUTFILE "$misassem{$ID}\t\n";	
}

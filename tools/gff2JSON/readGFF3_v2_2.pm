#!/usr/bin/env perl
#
#	Read GFF file
#
# AUTHOR:Gemy George Kaithakottil (gemy.kaithakottil@tgac.ac.uk || gemygk@gmail.com)

#package gff2ensembl;
package readGFF3_v2;
use strict;
use warnings;


use base 'Exporter';
our @EXPORT = qw(readGFF);


my $usage = "
Script to read GFF file

	Usage: perl $0 <file.gff>

OUTPUT to STDOUT

NOTE:
   Just make sure the - strand (col 7) \"increase\" from one line to the following
   e.g.   mRNA	10131	15124
          exon	10131	10348
          exon	12173	12310
          exon	15037	15124
\n";

my $gff = "";

my $gene_id="";
my $mrna_id="";

my @gene_array=();
my @mrna_array=();
my @exon_array=();

my $gene_line="";
my $mrna_line="";
my @exon_line=();

my $prev_gene_line="";
my $prev_mrna_line="";
#my $prev_exon_line="";

my %prev_parent_id_hash=();


my $temp="";			# for mRNA id
my $prev_parent_id="";	# for mRNA parent
my $boolean=0; 
my $gene_boo=0;

my @exon_info=();
my @cds_info=();
#my $exon_count=1;


sub readGFF{
#print "readGFF\n";

	$gff = $_[0];

	# Check the input GFF file for Errors
	# GFF3 file
	open(CHECK,$gff) or die "$!";
	while (<CHECK>) { # Read lines from file(s) specified on command line. Store in $_.
	    s/#.*//; # Remove comments from $_.
	    next unless /\S/; # \S matches non-whitespace.  If not found in $_, skip to next line.
		chomp;
		my @f = split (/\t/); # split $_ at tabs separating fields.
		if( $f[8] !~ /ID/ ){  # check if col 9 has got "NO" ID 
			print "ERROR:\tNo ID field for line\n\t\"$_\"\n";
			print "ERROR: Please make sure all features has an unique ID tag\n";
			exit;
		}
	}
	close(CHECK);

	# Now that file has passed the check, proceeding with it!!!

	# GFF3 file
	open(FILE,$gff) or die "$!";
	while (<FILE>) { # Read lines from file(s) specified on command line. Store in $_.
	    s/#.*//; # Remove comments from $_.
	    next unless /\S/; # \S matches non-whitespace.  If not found in $_, skip to next line.
		chomp;
		my @f = split (/\t/); # split $_ at tabs separating fields.

		# To take GFF with and without gene feature
		if($f[2] eq "gene") {
#print "gene\n";
			$gene_boo=1; # There is gene feature for the GFF

			# Gene ID
			($gene_id) = $f[8] =~ /ID\s*=\s*([^;]+)/;
			$gene_id =~ s/\s+$//;

			# Gene Name - not bothered at the moment
			#$f[8] =~ /Name\s*=\s*([^;]+)/;
			#my $gene_name = $1;
			#$gene_name =~ s/\s+$//;

			# Gene Note
			my $gene_note="";
			if( $f[8] =~ /Note\s*=\s*([^;]+)/ ) {
				#$f[8] =~ /Note\s*=\s*([^;]+)/;
				$gene_note = $1;
				$gene_note =~ s/\s+$//;
			}
			# Gene biotype:
			my $gene_biotype="";
			if ( $f[8] =~ /biotype\s*=\s*([^;]+)/ ) {
				$gene_biotype = $1;
				$gene_biotype =~ s/\s+$//;
			}
			#print "$f[2]\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tID=$gene_id;Name=$gene_id;biotype=$gene_biotype\n";

	#	#	load_gff::load_gene(@gene_array);  #-- Uncomment this when loading to db

			($gene_id eq "") ? ($gene_id = "NULL") : ($gene_id = $gene_id);
			($gene_note eq "") ? ($gene_note = "NULL") : ($gene_note = $gene_note);
			($gene_biotype eq "") ? ($gene_biotype = "NULL") : ($gene_biotype = $gene_biotype);

			$gene_line=join (" ",$f[2],$f[0],$f[1],$f[2],$f[3],$f[4],$f[5],$f[6],$f[7],$gene_id,$gene_note,$gene_biotype);

#print "$gene_line\n";
		}
		elsif($f[2] eq "mRNA" || $f[2] eq "transcript") {
#print "mRNA\n";	
		# mRNA ID
			$f[8] =~ /ID\s*=\s*([^;]+)/;
			$mrna_id = $1;
			$mrna_id =~ s/\s+$//;

			# Only if there is gene we need to look for parent
			my $parent_id="";
			if($gene_boo) { 
				# mRNA parent
				$f[8] =~ /Parent\s*=\s*([^;]+)/;
				$parent_id = $1;
				$parent_id =~ s/\s+$//;
			}
			# mRNA Note
			my $mrna_note="";
			if( $f[8] =~ /Note\s*=\s*([^;]+)/ ) {
				#$f[8] =~ /Note\s*=\s*([^;]+)/;
				$mrna_note = $1;
				$mrna_note =~ s/\s+$//;
			}
			# mRNA Name
			my $mrna_name=undef;
			if( $f[8] =~ /Name\s*=\s*([^;]+)/ ) {
				$mrna_name = $1;
				$mrna_name =~ s/\s+$//;
			}
			# mRNA biotype:
			my $mrna_biotype=undef;
			if ( $f[8] =~ /biotype\s*=\s*([^;]+)/ ) {
				$mrna_biotype = $1;
				$mrna_biotype =~ s/\s+$//;
			}

			#print "$f[2]\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tID=$mrna_id;Parent=$parent_id;Name=$mrna_name;biotype=$mrna_biotype\n";

		#	#print "@mrna_array\n";

			($mrna_id eq "") ? ($mrna_id = "NULL") : ($mrna_id = $mrna_id);
			($mrna_note eq "") ? ($mrna_note = "NULL") : ($mrna_note = $mrna_note);
			($mrna_biotype eq "") ? ($mrna_biotype = "NULL") : ($mrna_biotype = $mrna_biotype);
			($parent_id eq "") ? ($parent_id = "NULL") : ($parent_id = $parent_id);

			$mrna_line = join (" ",$f[2],$f[0],$f[1],$f[2],$f[3],$f[4],$f[5],$f[6],$f[7],$mrna_id,$mrna_note,$mrna_biotype,$parent_id);

	# This is where we get the exon and CDS linking
			if($mrna_id ne $temp) {
				if ($boolean) {
					unless (exists $prev_parent_id_hash{$prev_parent_id}) {
						$prev_parent_id_hash{$prev_parent_id}=$prev_gene_line;
						#print "$prev_gene_line\n";	# MAIN
						my @load_gene=();
						push (@load_gene,split(/" "/,$prev_gene_line));
						#print "@load_gene\n";

						# LOADING GENE
						load_gff::load_gene(@load_gene);
					}
					#print "$prev_mrna_line\n";		# MAIN
					my @load_mrna=();
					push (@load_mrna,split(/" "/,$prev_mrna_line));
					#print "@load_mrna\n";

					# LOADING mRNA
					load_gff::load_transcript(@load_mrna);

					# Get the first and last line of the CDS attribute - assuming incremental and sorted gff3 file
					my @cds_firststLine = split(/\t/,$cds_info[0]);
					my @cds_lastLine = split(/\t/,$cds_info[$#cds_info]);
 					my @load_cds=();
					foreach my $exn_line (@exon_line) {
						#print "$exn_line\n";		# MAIN
						my @load_exon=();
						push (@load_exon,split(/" "/,$exn_line));
						#print "@load_exon\n";

						# LOADING EXON
						load_gff::load_exon(@load_exon);

						my @exn_line = split(/\s/,$exn_line);
# 2nd is the index of @cds_firststLine is the CDS start
						# 4,5 is the index of @exn_line are the exon coords
						if ($cds_firststLine[2] > $exn_line[4] && $cds_firststLine[2] <= $exn_line[5]) {
							#print join (" ", "utr",$exn_line[1],$exn_line[2],"CDS",$exn_line[4],$cds_firststLine[2],$exn_line[6],$exn_line[7],$exn_line[8],$cds_firststLine[0],$exn_line[9]) ."\n";	# MAIN
#							my @load_cds=();
						push (@load_cds,$cds_firststLine[1]." ".($cds_firststLine[2]-$exn_line[4])." ".$exn_line[9]);							
						#print "@load_cds\n";

							# LOADING UTR
						#	load_gff::load_utr(@load_cds);
						}
						# 3rd is the index of @cds_lastLine is the CDS stop
						# 4,5 is the index of @exn_line are the exon coords
						if ($cds_lastLine[3] >= $exn_line[4] && $cds_lastLine[3] < $exn_line[5]) {
							#print join (" ", "utr",$exn_line[1],$exn_line[2],"CDS",$cds_lastLine[3],$exn_line[5],$exn_line[6],$exn_line[7],$exn_line[8],$cds_lastLine[0],$exn_line[9]) ."\n";	# MAIN
						#	my @load_cds=();
							push (@load_cds,($exn_line[5]-$cds_lastLine[3])." ".$exn_line[9]);

							#print "@load_cds\n";

							# LOADING UTR
							load_gff::load_utr(@load_cds);
						}
					}
				}
				$boolean=1;
				$temp=$mrna_id;  				# Main get the transcript ID
				$prev_parent_id=$parent_id;		# Gene ID of mRNA
				$prev_gene_line=$gene_line;
				$prev_mrna_line=$mrna_line;
				@exon_line=();
				@cds_info=();

				@gene_array=();
				@mrna_array=();
				@exon_array=();
			}
		}
		elsif($f[2] eq "exon") {
			# exon ID
			my ($exon_id) = $f[8] =~ /ID\s*=\s*([^;]+)/;
			$exon_id =~ s/\s+$//;
			
			# exon parent
			my ($parent_id) = $f[8] =~ /Parent\s*=\s*([^;]+)/;
			$parent_id =~ s/\s+$//;

			#print "$f[2]\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tID=$exon_id;Parent=$parent_id\n";

			# Store the exon information to an array to loop through
			my $exon_lines="$exon_id\t$parent_id\t$f[3]\t$f[4]";
			push(@exon_info,$exon_lines);

			($exon_id eq "") ? ($exon_id = "NULL") : ($exon_id = $exon_id);
			($parent_id eq "") ? ($parent_id = "NULL") : ($parent_id = $parent_id);

			my $exon_line = join (" ",$f[2],$f[0],$f[1],$f[2],$f[3],$f[4],$f[5],$f[6],$f[7],$exon_id,$parent_id);
			push (@exon_line,$exon_line);

		}
		elsif($f[2] eq "CDS") {
			# cds ID
			my ($cds_id) = $f[8] =~ /ID\s*=\s*([^;]+)/;
			$cds_id =~ s/\s+$//;
		
			# cds parent
			my ($parent_id) = $f[8] =~ /Parent\s*=\s*([^;]+)/;
			$parent_id =~ s/\s+$//;

			# Store the cds information to an array to loop through
			my $cds_lines="$cds_id\t$parent_id\t$f[3]\t$f[4]";
			push(@cds_info,$cds_lines);
		}
		if (eof(FILE)) {
			unless (exists $prev_parent_id_hash{$prev_parent_id}) {
				$prev_parent_id_hash{$prev_parent_id}=$prev_gene_line;
				#print "$prev_gene_line\n";	# MAIN
				my @load_gene=();
				push (@load_gene,split(/" "/,$prev_gene_line));
				#print "@load_gene\n";

				# LOADING GENE
				load_gff::load_gene(@load_gene);
			}
			#print "$prev_mrna_line\n";		# MAIN
			my @load_mrna=();
			push (@load_mrna,split(/" "/,$prev_mrna_line));
			#print "@load_mrna\n";

			# LOADING mRNA
			load_gff::load_transcript(@load_mrna);

			# Get the first and last line of the CDS attribute - assuming incremental and sorted gff3 file
			my @cds_firststLine = split(/\t/,$cds_info[0]);
			my @cds_lastLine = split(/\t/,$cds_info[$#cds_info]);

			my @load_cds=();
			
			foreach my $exn_line (@exon_line) {
				#print "$exn_line\n";		# MAIN
				my @load_exon=();
				push (@load_exon,split(/" "/,$exn_line));
				#print "@load_exon\n";

				# LOADING EXON
				load_gff::load_exon(@load_exon);

				my @exn_line = split(/\s/,$exn_line);
				# 2nd is the index of @cds_firststLine is the CDS start
				# 4,5 is the index of @exn_line are the exon coords
				if ($cds_firststLine[2] > $exn_line[4] && $cds_firststLine[2] <= $exn_line[5]) {
					#print join (" ", "utr",$exn_line[1],$exn_line[2],"CDS",$exn_line[4],$cds_firststLine[2],$exn_line[6],$exn_line[7],$exn_line[8],$cds_firststLine[0],$exn_line[9]) ."\n";	# MAIN
					push (@load_cds,$cds_firststLine[1]." ".($cds_firststLine[2]-$exn_line[4])." ".$exn_line[9]);							
					#push (@load_cds,"utr",$exn_line[1],$exn_line[2],"CDS",$exn_line[4],$cds_firststLine[2],$exn_line[6],$exn_line[7],$exn_line[8],$cds_firststLine[0],$exn_line[9]);
					#print "@load_cds\n";

					# LOADING UTR
					#load_gff::load_utr(@load_cds);
				}
				# 3rd is the index of @cds_lastLine is the CDS stop
				# 4,5 is the index of @exn_line are the exon coords
				if ($cds_lastLine[3] >= $exn_line[4] && $cds_lastLine[3] < $exn_line[5]) {
					#print join (" ", "utr",$exn_line[1],$exn_line[2],"CDS",$cds_lastLine[3],$exn_line[5],$exn_line[6],$exn_line[7],$exn_line[8],$cds_lastLine[0],$exn_line[9]) ."\n";	# MAIN
					
					push (@load_cds,($exn_line[5]-$cds_lastLine[3])." ".$exn_line[9]);
					#push (@load_cds,"utr",$exn_line[1],$exn_line[2],"CDS",$cds_lastLine[3],$exn_line[5],$exn_line[6],$exn_line[7],$exn_line[8],$cds_lastLine[0],$exn_line[9]);
					#print "@load_cds\n";

					# LOADING UTR
					load_gff::load_utr(@load_cds);
				}
			}
		}
	}
	close(FILE);
}
1;

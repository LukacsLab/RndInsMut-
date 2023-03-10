#!/usar/bin/perl
# Programmer: Elton Vasconcelos (20/Aug/2019)
# Edited by: Julia Gala de Pablo (04/May/2022)
# Usage: perl getClusters4VL.pl VL-allSamples.merged.bedLabCPM.tsv [CPM_CUTOFF] [window_size] [number_of_samples] >outfile.tsv
# Script that finds clusters of reads in different samples in the same genomic locus (considering a 100 kb window and 1 CPM per sample, as well as hal of the samples present on the locus).
# This script was first run within /nobackup/fbsev/LeedsOmics/VLukacs-pumpPriming/03-mergeBed/ directory on arc3 server.
# Please read cmd-line.txt file under that directory in order to see how the input file "VL-allSamples.merged.bedLabCPM.tsv" was generated.

use List::MoreUtils qw(uniq);
open (INFILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
my $line = <INFILE>;
chomp($line);
my(@array, @coord_start, @coord_end, @samples, @polarity, @rawreads, @unique_samples, @cpm, $chr); 
my $c = 1;

while ($line ne "") {
	@array = split(/\t/, $line);
	$chr = $array[0];
	while ($array[0] eq $chr && $line ne "") {
		if ($array[6] >= $ARGV[1] && $array[5]>=5) {
			push(@coord_start, $array[1]);
			push(@coord_end, $array[2]);
			push(@samples, $array[3]);
      push(@polarity, $array[4]);
      push(@rawreads, $array[5]);
			push(@cpm, $array[6]);
			$line = <INFILE>;
			chomp($line);
			@array = split(/\t/, $line);
			while ($array[1] <= $coord_start[0] + $ARGV[2] && $array[0] eq $chr && $line ne "") {
  			push(@coord_start, $array[1]);
  			push(@coord_end, $array[2]);
  			push(@samples, $array[3]);
        push(@polarity, $array[4]);
        push(@rawreads, $array[5]);
  			push(@cpm, $array[6]);
				$line = <INFILE>;
				chomp($line);
				@array = split(/\t/, $line);
			}
      for ($i=0;$i<$ARGV[3]*20;$i++){
        if ($cpm[$i]<$ARGV[1]){
          $samples[$i]=""}}          
			@unique_samples = uniq @samples;
			if (@unique_samples >= $ARGV[3]) {
				print("\# CLUSTER_$c \n");
				for ($i=0;$i<$ARGV[3]*20;$i++) {
          if ($cpm[$i] >= $ARGV[1] && $rawreads[$i] >= 5) {
					  print("$chr\t$coord_start[$i]\t$coord_end[$i]\t$samples[$i]\t$cpm[$i]\t$rawreads[$i]\t$polarity[$i]\n");
				  }}
				@coord_start = ();
				@coord_end = ();
				@samples = ();
        @polarity = ();
        @rawreads = ();
				@cpm = ();
				$c++;}		
			else {
				@coord_start = ();
				@coord_end = ();
				@samples = ();
        @polarity = ();
        @rawreads = ();
				@cpm = ();
			}
		}
		else {
			$line = <INFILE>;
			chomp($line);
			@array = split(/\t/, $line);
		}
	}
}
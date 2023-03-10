#$ -cwd
#$ -V
#$ -l h_rt=3:00:00,h_vmem=2G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-768G
# 2) Merging positions from each individual bed file, considering a maximum distance of 800nt between aligned reads, and counting the occurrences 
# based on 6th column (strand), but could actually have been based on any column, since the operation is just for counting.�
# Changing the code to merge by polarity (forward/reverse) and to report polarity in the output file
# This code additionally obtains the list of start and end positions, and calculates the median of those to calculate an accurate start and end.
for i in `ls ../02-bwa2020/*sorted.bed | uniq`;\
do /nobackup/fbsev/bioinformatics-tools/bedtools2/bin/mergeBed -d 600 -s -c 2,3,6,6,4 -o collapse,collapse,distinct,count,collapse -i $i >`echo $i | sed 's/..\/02-bwa2020\///g' | sed 's/sorted/merged/g'`;\
done
# Remove the long name of the insertion, keeping only the polarity
# Remove those lines that correspond to read 2 - this ensures they are only counted once and keep the polarity of the read of interest
# Substitutes the sample names with the read (1 or 2), only one per line (they are all the same). # sed 's/M[^,]*\//R/g' keeps all the Read numbers for a certain merge, but they are all the same. Careful! The names start with M, but there is also a chrM!
# Add the clone name to the end of the line
for i in `ls *merged.bed | sed 's/.merged.bed//g' | uniq`; 
do less `echo $i`.merged.bed | sed 's/M[^\t]*\//R/g' > `echo $i`_1.merged.bed;\ # This is different for CHO: here the locations start with a K
less `echo $i`_1.merged.bed | sed '/R2/d' > `echo $i`_2.merged.bed;\
less `echo $i`_2.merged.bed | sed 's/\tR1//g' > `echo $i`_1.merged.bed;\
less `echo $i`_1.merged.bed | sed "s/$/\t$i/g" > `echo $i`_3.merged.bed;\
rm `echo $i`_1.merged.bed;\
rm `echo $i`_2.merged.bed;\
# We use the mode of the starting position of R1s to calculate the start position, rather than the output from merge.
cut -f 5 `echo $i`_3.merged.bed > `echo $i`_pos_end.tsv
while read -r line; do echo $line | grep -Po '([0-9]*)' | sed 's\,\\g' | uniq --count | sort -n | tail -n1 | awk '{ print $2}' ; done < `echo $i`_pos_end.tsv > `echo $i`_mode_end.tsv
cut -f 4 `echo $i`_3.merged.bed > `echo $i`_pos_start.tsv
while read -r line; do echo $line | grep -Po '([0-9]*)' | sed 's\,\\g' | uniq --count | sort -n | tail -n1 | awk '{ print $2}' ; done < `echo $i`_pos_start.tsv > `echo $i`_mode_start.tsv
cut -f 1 `echo $i`_3.merged.bed > w
cut -f 1 `echo $i`_mode_start.tsv > x1
cut -f 1 `echo $i`_mode_end.tsv > x2
cut -f 6,7,8 `echo $i`_3.merged.bed > y
paste w x1 x2 y > `echo $i`.merged.bed
rm w x1 x2 y
rm `echo $i`_3.merged.bed
rm `echo $i`_pos_start.tsv
rm `echo $i`_pos_end.tsv
rm `echo $i`_mode_start.tsv
rm `echo $i`_mode_end.tsv
done
# 5) Concatenating all libraries/clones in a single merged.bedLab.tsv file (cmd 1012), and then flipping columns 4 (library name) and 5 # now we need to flip columns 6 with 5 instead:
# (counts) and getting rid of intermediate unnecessary files (cmds 1013-1021).
cat *bed | sort -k1,1 -k2,2n >CHO6-11LIB-Tol2.merged.bedLab.tsv; # Concatenating them all in the same file
cut -f 1,2,3 CHO6-11LIB-Tol2.merged.bedLab.tsv >x;
cut -f 6 CHO6-11LIB-Tol2.merged.bedLab.tsv >y;
cut -f 4,5 CHO6-11LIB-Tol2.merged.bedLab.tsv >z;
paste x y z >CHO6-11LIB-Tol2_2.merged.bedLab.tsv;
rm x y z
rm CHO6-11LIB-Tol2.merged.bedLab.tsv;
mv CHO6-11LIB-Tol2_2.merged.bedLab.tsv CHO6-11LIB-Tol2.merged.bedLab.tsv;
# 6) Counting the total number of mapped reads in each library
wc -l ../02-bwa2020/*sorted.bed | sed 's/^ *//g'| sed 's/..\/02-bwa2020\///g' | sed 's/.sorted.bed//' | sed 's/ /\t/g' >totalAlignedReads.tab;
grep -v 'total$' totalAlignedReads.tab | sort -k1,1nr | sed 's/.sorted.bed//g' >totalAlignedReads-sorted.tab;
# 7) Getting normalized CPM values
perl -e 'open(FILE, "CHO6-11LIB-Tol2.merged.bedLab.tsv");
open(FILE2, "totalAlignedReads-sorted.tab");
while(<FILE2>) {chomp($_); @array2 = split(/\t/, $_); $hash{$array2[1]} = $array2[0];} while(<FILE>) {chomp($_);
@array = split(/\t/, $_);
if($hash{$array[3]} ne "") {$cpm = ($array[5] / $hash{$array[3]}) * 1000000; $cpm = int($cpm); print("$_\t$cpm\n");
} else {print STDERR ("$_\t***Error***\n");
}}' >CHO6-11LIB-Tol2_CPM.tsv
# 8) Getting clusters
# perl getClusters4VL.pl VL-allSamples.merged.bedLabCPM.tsv [CPM_cutoff] [window_size] [number_of_samples] >outfile.tsv
perl getClusters_CPM.pl CHO6-11LIB-Tol2_CPM.tsv 1000 100000 1 >CHO6-11LIB-Tol2_CPM_clusters.tsv
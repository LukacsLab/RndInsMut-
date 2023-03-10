#$ -cwd
#$ -V
#$ -l h_rt=10:00:00,h_vmem=4G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-128G
# To clean any non-paired alignments
for i in `ls ./*.sam | sed -r 's/.\///g'| sed -r 's/.sam//g'| uniq`; \
do /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools fixmate -rc `echo $i`.sam `echo $i`-fixmates.sam
echo 'Fixmates was run in file' $i
# Running fixmate to fix some mess-ups by BWA MEM
grep -Po '^(.*?)\t(?:.*?)XS:i:(?:[^0]{1,6})' `echo $i`-fixmates.sam | sed -r 's/\t.*//g' | uniq > `echo $i`-locs.txt;\
echo 'Locations obtained is XS flag>0 for file' $i
# Obtaining reads where the XS tag is > 0 (multiple possible alignments) and saving them to a file.
grep -v -f `echo $i`-locs.txt `echo $i`-fixmates.sam > `echo $i`-unique.sam;\
#
## Alternative method, but very slow!
#cp `echo $i`-fixmates.sam `echo $i`-unique.sam
#while IFS= read -r line; do
#	sed -i "/^$line/d" `echo $i`-unique.sam
#done < `echo $i`-locs.txt
#echo 'Locations with XS flag>0 removed from file' $i
#
# Remove both pairs for those locations
/nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view -F12 -q 40 -h -o `echo $i`-paired.sam `echo $i`-unique.sam;\
echo 'Paired reads obtained for file' $i
done
# Quality control to obtain number of positive and negative reads
echo -e "Clone Name\tSam file\tFixmates\tUnique\tPaired" > ClusterStats.txt
for i in `ls ./*.sam | sed -r 's/.\///g'| sed -r 's/-paired//g'| sed -r 's/-unique//g' | sed -r 's/-fixmates//g' | sed -r 's/.sam//g'| uniq`; \
do \
printf `echo $i` >> ClusterStats.txt; printf "\t" >> ClusterStats.txt;\
sed -i "\${s/$/\t$( /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view `echo $i`.sam | wc -l  )/}" ClusterStats.txt;\
sed -i "\${s/$/\t$( /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view `echo $i`-fixmates.sam | wc -l  )/}" ClusterStats.txt;\
sed -i "\${s/$/\t$( /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view `echo $i`-unique.sam | wc -l  )/}" ClusterStats.txt;\
sed -i "\${s/$/\t$( /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view `echo $i`-paired.sam | wc -l  ) \n/}" ClusterStats.txt;\
rm `echo $i`-unique.sam
rm `echo $i`-fixmates.sam
rm `echo $i`-locs.txt
done
echo 'Cluster stats file produced'
# To convert sam to bam
for i in `ls ./*paired.sam | sed -r 's/.\///g'| sed -r 's/-paired//g'| sed -r 's/.sam//g'| uniq`; \
do /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view -S -b `echo $i`-paired.sam  >  `echo $i`.bam;\
/nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools sort `echo $i`.bam -o `echo $i`.sorted.bam;\
/nobackup/fbsev/bioinformatics-tools/bedtools2/bin/bamToBed -i `echo $i`.sorted.bam >`echo $i`.sorted.bed;\
echo 'Sorted bed file obtained for file' $i
done
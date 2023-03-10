#$ -cwd
#$ -V
#$ -l h_rt=3:00:00,h_vmem=2G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-128G
# To clean any non-paired alignments
for i in `ls ./*.sam | sed -r 's/.\///g'| sed -r 's/.sam//g'| uniq`; \
do /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view -F12 -q 40 -h -o `echo $i`-paired.sam `echo $i`.sam;\
done
# Quality control to obtain number of positive and negative reads
for i in `ls ./*.sam | sed -r 's/.\///g'| sed -r 's/-paired//g'| sed -r 's/.sam//g'| uniq`; \
do \
printf `echo $i` >> ClusterStats.txt; printf "\t" >> ClusterStats.txt;\
sed -i "\${s/$/\t$( /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view `echo $i`.sam | wc -l  )/}" ClusterStats.txt;\
sed -i "\${s/$/\t$( /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view `echo $i`-paired.sam | wc -l  )/}" ClusterStats.txt;\
sed -i "\${s/$/\t$( /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view -F0x10 `echo $i`-paired.sam | wc -l  ) \n/}" ClusterStats.txt;\
done
# To convert sam to bam
for i in `ls ./*paired.sam | sed -r 's/.\///g'| sed -r 's/.sam//g'| uniq`; \
do /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools view -S -b `echo $i`.sam  >  `echo $i`.bam;\
done
# The code below is for sorting
for i in `ls ./*.bam | sed -r 's/.\///g'| uniq`; \
do \
/nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools sort `echo $i` -o `echo $i | sed 's/-paired.bam/.sorted.bam/g'`;\
done
# Now to convert bam to bed
for i in `ls ./*.sorted.bam | sed -r 's/.\///g'| uniq`; \
do /nobackup/fbsev/bioinformatics-tools/bedtools2/bin/bamToBed -i $i >`echo $i | sed 's/bam$/bed/g'`; \
done
#$ -cwd
#$ -V
#$ -l h_rt=3:00:00,h_vmem=2G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-128G
# To clean any non-paired alignments
for i in `ls /nobackup/phyjga/CHO6-11LIB/02-bwa/*paired.bam | sed -r 's/\/nobackup\/phyjga\/CHO6-11LIB\/02-bwa\///g'| sed -r 's/.bam//g'| uniq`; \
do /nobackup/fbsev/bioinformatics-tools/samtools-1.9/samtools index > `echo $i`.bai;\
done
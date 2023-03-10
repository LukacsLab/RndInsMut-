#$ -cwd
#$ -V
#$ -l h_rt=3:00:00,h_vmem=2G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-768G
for i in `ls ../01-trimmomatics/*_paired* | sed 's/_R[12]-trimmed_paired.fq.gz//g' | sed -r 's/..\/01-trimmomatics\///g'| uniq`; \
do /nobackup/fbsev/bioinformatics-tools/bwa-0.7.17/bwa mem -A 45 -t 20 -k 40 -w 4 -c 1 -O 100 -o `echo $i`.sam \
    GCF_003668045_3.fa \
	../01-trimmomatics/`echo $i`_R1-trimmed_paired.fq.gz	../01-trimmomatics/`echo $i`_R2-trimmed_paired.fq.gz; \
done
# Code explanation
# $ bwa mem index_prefix [input_reads.fastq|input_reads_pair_1.fastq input_reads_pair_2.fastq] [options]
# index_prefix is the index for the reference genome generated from bwa index (in this case, /nobackup/phyjga/02-bwa/chroms/hg38)
# input_reads.fastq¦ input_reads_pair_1.fastq, input_reads_pair_2.fastq are the input files of sequencing data that can be single-end or paired-end respectively. 
# We need to run here the output files for trimmomatics.
# (from the help of BWA-MEM)
# -t is the number of threads (20)
# -k minimum seed alignment
# -A sets a threshold for the mapping quality score
# -U sets a penalty for unpaired alignments (I set it very high)
# -o needs to be before the output file name
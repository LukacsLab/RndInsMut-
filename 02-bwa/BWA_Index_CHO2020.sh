#$ -cwd
#$ -V
#$ -l h_rt=12:00:00,h_vmem=2G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-768G
/nobackup/fbsev/bioinformatics-tools/bwa-0.7.17/bwa index -a bwtsw GCF_003668045_3.fa
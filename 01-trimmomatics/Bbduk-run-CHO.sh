#$ -cwd
#$ -V
#$ -l h_rt=1:00:00,h_vmem=2G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-768G
for i in `ls ../*fastq.gz | sed 's/_L001_R[12]_001.fastq.gz//g' | sed -r 's/..\///g' | uniq`; \
do \
/nobackup/phyjga/Tools/bbmap/bbduk.sh in1=../`echo $i`_L001_R1_001.fastq in2=../`echo $i`_L001_R2_001.fastq out1=`echo $i`_L001_R1_001-1.fastq out2=`echo $i`_L001_R2_001-1.fastq ref=Bbduk_adaptors_Left_R1.fa ktrim=l k=20 hdist=1;
/nobackup/phyjga/Tools/bbmap/bbduk.sh in1=`echo $i`_L001_R1_001-1.fastq in2=`echo $i`_L001_R2_001-1.fastq out1=`echo $i`_L001_R1_001-2.fastq out2=`echo $i`_L001_R2_001-2.fastq ftl=1
/nobackup/phyjga/Tools/bbmap/bbduk.sh in1=`echo $i`_L001_R1_001-2.fastq in2=`echo $i`_L001_R2_001-2.fastq out1=`echo $i`_L001_R1_001-3.fastq out2=`echo $i`_L001_R2_001-3.fastq ref=Bbduk_adaptors_Right_R1.fa ktrim=r k=15 hdist=0 mink=10 rcomp=f minlen=15;
/nobackup/phyjga/Tools/bbmap/bbduk.sh in1=`echo $i`_L001_R1_001-3.fastq in2=`echo $i`_L001_R2_001-3.fastq out1=`echo $i`_L001_R1_001.fastq out2=`echo $i`_L001_R2_001.fastq ref=Bbduk_adaptors_Right_R2.fa ktrim=r k=15 hdist=0 rcomp=f
done
rm *-1.fastq
rm *-2.fastq
rm *-3.fastq


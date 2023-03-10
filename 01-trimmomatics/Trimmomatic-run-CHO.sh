#$ -cwd
#$ -V
#$ -l h_rt=1:00:00,h_vmem=2G
#$ -l nodes=1,ppn=24
#$ -m abe
#$ -M phyjga@leeds.ac.uk
#$ -l node_type=24core-768G
#for i in `ls ../SingleClones/*fastq.gz | sed 's/_L001_R[12]_001.fastq.gz//g' | sed -r 's/..\/SingleClones\///g' | uniq`; do java -jar /nobackup/fbsev/bioinformatics-tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 ../SingleClones/`echo $i`_L001_R1_001.fastq.gz ../SingleClones/`echo $i`_L001_R2_001.fastq.gz `echo $i`_R1-trimmed_paired.fq.gz `echo $i`_R1-trimmed_unpaired.fq.gz `echo $i`_R2-trimmed_paired.fq.gz `echo $i`_R2-trimmed_unpaired.fq.gz ILLUMINACLIP:contaminant_list.fasta:2:30:7 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:30; done
for i in `ls ../*fastq.gz | sed 's/_L001_R[12]_001.fastq.gz//g' | sed -r 's/..\///g' | uniq`; \
do java -jar /nobackup/fbsev/bioinformatics-tools/Trimmomatic-0.36/trimmomatic-0.36.jar \
PE -threads 24 -phred33 `echo $i`_L001_R1_001.fastq `echo $i`_L001_R2_001.fastq \
`echo $i`_R1-trimmed_paired.fq.gz \
`echo $i`_R1-trimmed_unpaired.fq.gz \
`echo $i`_R2-trimmed_paired.fq.gz \
`echo $i`_R2-trimmed_unpaired.fq.gz \
TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:30; \
done
#
# Previous settings
# ILLUMINACLIP:TruSeq3-PE-Adaptors.fa:2:30:7 LEADING:20 HEADCROP:26 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:30; \
# I removed LEADING to avoid trimming the adaptor by mistake. I want to make sure I can trust the start of the aligment.
# I removed HEADCROPT to avoid trimming twice (adaptor and also following 26 bp). This will help me trust the start position.
# I changed the TRAILING values from 20 down to 3.
#
# ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clipthreshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
# our data is Pair End, so we choose PE and we will require two input files (forward and reverse reads) and 4 output files (forward paired, unpaired, reverse paired, unpaired)
# seed mismatches: specifies the maximum mismatch count which will still allow a full match to be performed
# palindrome clipthreshold:  specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
# simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
# minAdapterLength - default as 8, but can be reduced to remove smaller adapter cuts for the reverse reads.
#
# SLIDINGWINDOW:<windowSize>:<requiredQuality>
#
#  We choose 24 threads (cores) to run our data
#  the quality score character encoding is chosen between phred33 or the standard phred64. More about it here: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
#  the saving location is set to ../SingleClones/`echo $i`_L001_R1_001.fastq.gz, where `echo $i` is the clone name that is running at the moment from the previously trimmed list of clones.
# input file forward corresponds to read 1
# input file reverse corresponds to read 2
# output file paired forward
# output file unpaired forward
# output file paired reverse
# output file unpaired reverse
# LEADING:20 \ # removed low quality (<20) bases from the beginning
# TRAILING:20 \ # removes low quality (<20) bases from the end
# SLIDINGWINDOW:4:15 \ # Scans the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 
# MINLEN:30; \ # removes reads less than 36 bases long
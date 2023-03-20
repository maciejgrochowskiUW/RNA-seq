cd ~
cd workdir/
gunzip *.gz
source activate fastq
#quality control with fastqc, -o stands for output#
for fastqfiles in *.fastq
do
	fastqc $fastqfiles -o fastqc/
done

#this script prepares your mapped files for analysis in IGV software#
###############
cd ~
cd workdir/mapped
source activate fastq
###############
#convert your .sam files to .bam#
###############
for files in s_*.sam
do
	base=${files%.sam}
	short=${base#s_}
	samtools view -Sb $files -o ${short}.bam
done
###############
#and index it#
#-b generates .bai file#
###############
for files in *.bam
do
	samtools index -b $files
done
###############
#now u can analyse your mapping in IGV#
###############
source activate igv
igv

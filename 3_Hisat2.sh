#this script edits raw fastq files and once those are ready,#
#it aligns them to reference genome#
#in also marks PCR duplicates in your alignment file# 
#before alignment u have to prepare a reference genome,#
#download fasta (*fa.gz) file that contains DNA sequence#
#for S. pombe https://fungi.ensembl.org/Schizosaccharomyces_pombe/Info/Index#
#once you download it, copy it to Documents/workdir/reference#
#################
cd ~
cd Documents/workdir/
source activate fastq
#################
#move UMIs to headers, here it moves 8 bases, to change it UMI1:<put number of bases to be moved here>#
#################
for files in *.fastq
do
	je clip F=$files 'RL=<UMI1:8><SAMPLE1:x>' 'OL=1:U1:S1' OUTPUT_DIR=edit
	gunzip edit/out_1.txt.gz
	mv edit/out_1.txt edit/c_$files
	gzip $files
done
ls -l edit/
#################
#for each read remove first nucleotide that comes from ligation step (-f 2)#
#and last one that usually has odd base distribution (-t 1)#
#################
cd edit/
for files in c_*.fastq
do
	fastx_trimmer -f 2 -i $files -o t_$files
	rm $files
	fastx_trimmer -t 1 -i t_$files -o tt_$files
	rm t_$files
done
ls -l
#################
#now index your reference genome#
#################
cd ..
cd reference/
gunzip *.gz
hisat2-build *.fa reference
#################
#this divides your file into several smaller once named reference?.ht2#
#your files are ready to be aligned#
#################
cd ..
cd edit/
for file in *_R1.fastq
do 
	base=${file%_R1*}
	x=${base}_R1.fastq y=${base}_R2.fastq;
	echo $x $y
	short=${base#*_c_}
	echo $short
	hisat2 -q -x ../reference/reference -1 $x -2 $y --rna-strandness FR -S ../mapped/${short}.sam --summary-file  ../mapped/hisat2/${short}.txt --new-summary -p 16
	rm $x
	rm $y
done
#################
#-q indicates that input file is in .fastq format,#
#-x ../prefix idicates which files will be used as reference,#
#-S is used to define output file#
#introducing "base" is neccessery to pair reads,#
#"%_" shorten name by cutting off everything that stands after %#
#I also introduced "short" to make .sam file look like <genename.sam>#
#"#*_c_" shorten name by cutting off everything that stands pre "_c_" including this motive#
#################
#now u want to mark PCR duplicates and check out output file with statistics#
#at first sort .sam files with samtools, then mark duplicates with je#
#################
cd ../mapped
for file in *.sam
do
	samtools sort $file -o s_$file -@ 8
	rm $file
done
for file in s_*.sam
do
	base=${file#s_}
	red='\033[0;31m'
	white='\033[0m'
	printf "${red}$base${white}"
	short=${base%.sam}
	je markdupes I=$file O=m_$file M=markdupes/dups_report_${short}.txt MM=0 
	rm $file 
done

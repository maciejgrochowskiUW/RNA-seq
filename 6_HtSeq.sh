#this script will create a file, that shows how many reads are mapped#
#to any particular feature annotated in reference genome.#
#in order to make htseq analysis u need to prepare your .sam file#
###############
cd ~
cd workdir/mapped
source activate fastq
###############
#first sort your reads by names#
###############
for files in m_s_*.sam
do
	base=${files#m_s_}
	samtools sort -n $files -o sbn_$base
	rm $files
done
###############
#By default feature=exon (use -t <featurename> option to change it).#
#Script searches for feature name in third column of .gff3 file,#
#(REMEMBER to change code if your reference has different extension than .gff3)#
#then jumps to 8th column and by deafult looks for gene_id,#
#which might be not present there. In such a case use -i <8thcolumnstartingword> option.#
#Depending on the method of library preparation Read1 might be equal#
#to feature sequence or complementar to it. By default its treated as equal,#
#here we change it with -s <reverse> option.#
###############
for files in sbn_*.sam
do
	base=${files%.sam}
	short=${base#sbn_}
	red='\033[0;31m'
	white='\033[0m'
	printf "${red}$short${white}" 
	htseq-count -t gene -i ID -s reverse $files ../reference2/*.gff3 > htseq/counttable_gene_${short}.txt
done
###############
#this command below prints the sum of every score#
#located in second collumn of the .txt file, exept last 5 lanes#
cd htseq
for file in *.txt
do
	red='\033[0;31m'
	white='\033[0m'
	printf "${red}$file${white}" 
	head -n -5 $file | cut -f 2 | awk '{total += $0} END{print "sum="total}'
done

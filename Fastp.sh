#fastp merges R1 with R2

for files in *R1.fastq
do
base=${files%_R1*}
echo $base
x=${base}_R1.fastq
y=${base}_R2.fastq
echo $x $y
#short=${base#_c_}
#echo $short
fastp -i $x -o ${short}.fastq [-I $y] [-m][--discard_unmerged][--overlap_len_require 10][-j fastp/${short}.json][-h fastp/${short}.html][-A]

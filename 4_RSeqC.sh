#this script allows to check theoritical distances between paired reads#
#download annotated reference genome#
#it should have .ggf3 extension or similar, if it differs change code below#
#move file to Documents/workdir/reference2#
##################
source activate fastq
cd ~
cd workdir/reference2
gff2bed < *.gff3 > reference.gff3.bed
##################
#use RSeqC - inner_distance, it will create 4 files, one of which is run with R#
#-i is input, -r is reference, -o is output#
##################
cd ../mapped
for files in m_s_*.sam
do
	base=${files#*m_s_}
	inner_distance.py -i $files -r ../reference2/file_bed* -o rseqc/innerdistance_$base
done
##################

for files in m_s_*.sam
do
	base=${files#*m_s_}
	junction_annotation.py -i $files -o rseqc/junctionannotation_$base -r ../reference2/file_bed*
done
##################
#junction_saturation tells u how                     
for files in m_s_*.sam
do
	base=${files#*m_s_}
	junction_saturation.py -i $files -o rseqc/junctionsaturation_$base -r ../reference2/file_bed*
done
##################
cd ../
geneBody_coverage.py -r reference2/file_bed* -i mapped/  -o rseqc/genebodycoverage_$base
cd ~

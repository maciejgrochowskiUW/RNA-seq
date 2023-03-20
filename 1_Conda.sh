#download miniconda https://conda.io/miniconda.html and run command#
bash Downloads/Miniconda3-latest-Linux-x86_64.sh
#restart terminal (this is not tested) and test instalation#
exec bash
read -p "Press [Enter] key to start backup..."
conda list
#now create environment for conda packages and install them#
conda create --name fastq
source activate fastq
for packages in seqkit fastqc je-suite HISAT2 samtools bedops rseqc htseq cutadapt
do
	conda install -c bioconda $packages
done
conda install -c r r
conda install -c biobuilds fastx-toolkit
source deactivate fastq
#now create another environment for igv (it can't be installed to fastq env)#
conda create --name igv
source activate igv
conda install -c bioconda igv
source deactivate igv

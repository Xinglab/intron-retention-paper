# intron retention project

The scripts used in the manuscript entitled 'Tracking pre-mRNA maturation across subcellular compartments identifies developmental gene regulation through intron retention and nuclear anchoring' are included in the present folders for the purpose of reproducibility and transparency of the results from the manuscript.

The scripts are broken down by the order of analysis.

Please download the datasets from https://doi.org/10.5281/zenodo.4540589 and run the analysis.

The code is provided without a warranty. Codes from parts of the analysis are rewritten to make them easy to understand thus may cause errors when running them. If you encounter any problems or questions, please report them to the author (zcpan1016@gmail.com).

## download data

* Download from Zenodo: https://doi.org/10.5281/zenodo.4540589
  + `cd data`
  + `curl -L -o Intron.feature.annotation https://zenodo.org/record/4540589/files/Intron.feature.annotation?download=1`
  + `curl -L -o intron_FI_combined_results.txt https://zenodo.org/record/4540589/files/intron_FI_combined_results.txt?download=1`
  + `curl -L -o Intron_transcript.txt https://zenodo.org/record/4540589/files/Intron_transcript.txt?download=1`
  + `curl -L -o mm10.fa.gz https://zenodo.org/record/4540589/files/mm10.fa.gz?download=1`
  + `curl -L -o Mus_musculus.GRCm38.91.chr.gtf.gz https://zenodo.org/record/4540589/files/Mus_musculus.GRCm38.91.chr.gtf.gz?download=1`
* Download other data:
  + `curl -L -O ftp://ftp.ensembl.org/pub/release-91/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz`
* Unzip
  + `gunzip -c ./mm10.fa.gz > mm10.fa`
  + `gunzip -c ./Mus_musculus.GRCm38.91.chr.gtf.gz > ./Mus_musculus.GRCm38.91.chr.gtf`
  + `gunzip -c ./Mus_musculus.GRCm38.cdna.all.fa.gz > ./Mus_musculus.GRCm38.cdna.all.fa`
  + `cd ..`

## conda environment

* `conda create --prefix ./conda_env_py2`
* `conda activate ./conda_env_py2`
* `conda install -c conda-forge -c bioconda python=2 r-base=4 sra-tools star kallisto samtools=1.11 pysam=0.16.0 numpy=1.15.4 tqdm keras scikit-learn h5py r-rtsne r-ggplot2 r-ggthemes r-scales r-ggpubr`
* `R`
* `> repos <- "http://cran.us.r-project.org"`
* `> install.packages("circlize", repos=repos)`
* `> quit()`
* `conda deactivate`
* `conda create --prefix ./conda_env_py3`
* `conda activate ./conda_env_py3`
* `conda install -c conda-forge -c bioconda python=3 r-base=4 numpy scipy pyclustering r-ggplot2 r-ggthemes r-scales r-ggpubr`
* `conda deactivate`

## `quantification_intronRetention_1`

Download fastq files

* `conda activate ./conda_env_py2`
* `cd quantification_intronRetention_1`
* `prefetch SRR12883492`
* `fastq-dump --split-files ./SRR12883492/SRR12883492.sra`

Create indices

* `STAR --runMode genomeGenerate --runThreadN 4 --genomeDir ./star_index --genomeFastaFiles ../data/mm10.fa --sjdbGTFfile ../data/Mus_musculus.GRCm38.91.chr.gtf`
* `kallisto index ../data/Mus_musculus.GRCm38.cdna.all.fa --i Mus_musculus.GRCm38.cdna.all.fa.idx`

Setup for run

* `mkdir -p results/mESC/Chr`
* `mkdir -p data/mESC/fastq`
* `mv ./SRR12883492_1.fastq data/mESC/fastq/Chr_1.fastq`
* `mv ./SRR12883492_2.fastq data/mESC/fastq/Chr_2.fastq`
* edit `PARENT_DIR` variable at top of `RNAseq_Analysis_Pipeline.py`
* `chmod +x SIRI/bin/SIRI`

Run

* 8 threads 64 GB: `python RNAseq_Analysis_Pipeline.py`
* `conda deactivate`
* `cd ..`

## `xmeans_cluster_U_intron_2`

Output: cluster.png should look like what is shown in the .ipynb

* `conda activate ./conda_env_py3`
* `cd xmeans_cluster_U_intron_2`
* `python run_scripts.py`
* `conda deactivate`
* `cd ..`

## `ptc_analysis_3`

Output: PTC_figure.png should look like what is shown in the .ipynb

* `conda activate ./conda_env_py2`
* `cd ptc_analysis_3`
* `python run_scripts.py`
* `conda deactivate`
* `cd ..`

## `deep_learning_4`

Output: ./results/performance.png should look like what is shown in the .ipynb

* `conda activate ./conda_env_py2`
* `cd deep_learning_4`
* 8 threads 32 GB: `python run_scripts.py`
* `conda deactivate`
* `cd ..`

## `tsne_analysis_5`

Output: tsne_plot.pdf should look like the one uploaded to the repo

* `conda activate ./conda_env_py2`
* `cd tsne_analysis_5`
* `Rscript tsne_plot.R`
* `conda deactivate`
* `cd ..`

## `circos_plot_6`

Output: circos.pdf should look like the one uploaded to the repo

* `conda activate ./conda_env_py2`
* `cd circos_plot_6`
* `Rscript circos_plot.R`
* `conda deactivate`
* `cd ..`

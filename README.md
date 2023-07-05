# CYTO-SV-ML
#### CYTO-SV-ML pipeline design:
![CYTO-SV-ML PIPELINE WORKFLOW](workflow.png)

### CYTO-SV-ML snakemake workflow
![CYTO-SV-ML Snakemake Snakemake workflow](cyto-sv-ml_Snakemake_workflow.png)

<!--ts-->
   * [Installation](#Installation)
   * [Pipeline Usage](#Pipeline-usage)
   * [Download SV related Resource](#Download-SV-related-Resource)

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Installation
============
### Install conda environment
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
# conda update -n base -c defaults conda # in case that conda version is old
conda install -n base -c conda-forge mamba
```

### Install parliament and chromoseq Docker
```
# install docker in the linux system
sudo yum install -y yum-utils
sudo yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
sudo yum install docker-ce docker-ce-cli containerd.io
sudo systemctl start docker
# download parliament and chromoseq Docker
sudo docker pull docker.io/dnanexus/parliament2:latest # Parliament: https://github.com/dnanexus/parliament2
sudo docker pull  docker.io/zatawada/docker-basespace_chromoseq_v2:master # ChromoSeq: https://github.com/genome/docker-basespace_chromoseq
```

### Install CYTO-SV-ML Snakemake pipeline
```
git clone https://github.com/tzhang-nmdp/CYTO-SV-ML.git
cd CYTO-SV-ML
mamba env create -n py27 -f py27.yaml
mamba env create -n py39 -f py39.yaml
mamba env create -n cyto-sv-ml -f cyto-sv-ml.yaml
conda activate cyto-sv-ml
pip install --upgrade snakemake # in case that snakemake version is old
# note: python-graphviz might be conflicted with the pre-existing packages in your environment. If it happened, please remove it from requirement file and install it separately.
```

### Download hg38 reference genome 
```
cd CYTO-SV-ML/reference
mkdir CYTO-SV-ML/reference/hg38
# hg38 from Broad Institute Google Cloud
wget https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
mv Homo_sapiens_assembly38.fasta reference/hg38/hs38.fasta
mv Homo_sapiens_assembly38.fasta.fai reference/hg38/hs38.fasta.fai
# VEP GRCh38 from ensembl
curl -O https://ftp.ensembl.org/pub/release-90/variation/VEP/homo_sapiens_vep_90_GRCh38.tar.gz
tar xzf homo_sapiens_vep_90_GRCh38.tar.gz
cd CYTO-SV-ML/reference/homo_sapiens/90_GRCh38
```

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Pipeline Usage
============
### 1. Run CYTO-SV-ML Snakemake preprocess pipeline for each sample
Please change the config.yaml according to your own environment settings:                                                    
{your_work_dir} for input/output dir   
{your_work_dir}/in/${sample}.cram # input for chromseq pipeline <br/>
{your_work_dir}/in/${sample}.bam # input for parliament2 pipeline <br/>
{cyto_sv_ml_dir} for cyto_sv_ml dir and software subdir and reference subdir
```
conda activate cyto-sv-ml
snakemake --core ${number_of_cores} -s cyto-sv-ml-preprocess.smk --use-conda --config sample=${sample} gender=${gender}
```
### 2. Run CYTO-SV-ML Snakemake modeling pipeline for whole cohort
```
snakemake --core ${number_of_cores} -s cyto-sv-ml-automl.smk --config cohort_name=${cohort_name}
snakemake --core ${number_of_cores} -s cyto-sv-ml-automl.smk --config cohort_name=${cohort_name} --report ${your_work_dir}/out/${cohort_name}_report.html
```
![CYTO-SV-ML snakemake report](cyto-sv-ml_snakemake_report.png)

### 3. Run CYTO-SV-ML Snakemake interface pipeline for Web-portal application
```
snakemake --core ${number_of_cores} -s cyto-sv-ml-interface.smk  --config cohort_name=${cohort_name} k=${kfolds}
# to run the docker image in the local machine and open user interface with "http://localhost:8000/"
sudo docker run -d -p 8000:80 cyto-sv-ml-app:${sample_all}
```
The analysis summary of 494 MDS cohort using CYTO-SV-ML pipeline 
### Online Web-portal
```
http://cyto-sv-ml.b12x.org/
```
![CYTO-SV-ML Rshiny web-portal](cyto-sv-ml_web-portal.png)

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
Download SV related Resource
=============================

## SV database Download (The official websites contain detailed information)
```
gnomAD: https://gnomad.broadinstitute.org/downloads
1000g: https://www.internationalgenome.org/phase-3-structural-variant-dataset
CytoAtlas: https://github.com/genome/docker-basespace_chromoseq/blob/master/workflow_files/chromoseq_translocations.bedpe
COSMIC: https://cancer.sanger.ac.uk/cosmic/download
```

## SV Tool Download
```
SURVIVOR: https://github.com/fritzsedlazeck/SURVIVOR
SVTyper: https://github.com/hall-lab/svtyper
SeqComplex: https://github.com/caballero/SeqComplex
Komplexity: https://github.com/eclarke/komplexity
```

## ML Pipeline Download
```
AUTOML: https://github.com/mljar/mljar-supervised
```

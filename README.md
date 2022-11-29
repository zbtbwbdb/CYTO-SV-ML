# CYTO-SV-ML
#### CYTO-SV-ML PIPELINE WORKFLOW:
![CYTO-SV-ML PIPELINE WORKFLOW](workflow.png)

# Run CYTO-SV-ML Snakemake pipeline
### CYTO-SV-ML Snakemake pipeline workflow
![CYTO-SV-ML Snakemake Snakemake workflow](cyto-sv-ml_Snakemake_workflow.svg)
### Install CYTO-SV-ML Snakemake pipeline
```
git clone https://github.com/tzhang-nmdp/CYTO-SV-ML.git
cd CYTO-SV-ML
conda env create cyto-sv-ml -f cyto-sv-ml.yaml
```

### Run CYTO-SV-ML Snakemake preprocess pipeline
change the config.yaml to your own environment settings
```
conda activate cyto-sv-ml
snakemake --core ${number_of_cores} -s cyto-sv-ml-preprocess.smk --config sample=${sample} gender=${gender}
```
### Run CYTO-SV-ML Snakemake AUTO-ML pipeline
change the config.yaml for your own settings
```
snakemake --core ${number_of_cores} -s cyto-sv-ml-modeling.smk
```

# Run Shiny Web-portal 
The summarized visualization of CYTO-SV-ML pipeline 
### Install Shiny Web-portal
```
git clone https://github.com/tzhang-nmdp/CYTO-SV-ML.git
cd CYTO-SV-ML
docker build -t CYTO-SV-ML:main .
```

### Start Shiny Web-portal
```
docker run -d -p 8000:80  CYTO-SV-ML:main
```

### open Shiny Web-portal
```
http://127.0.0.1:8000/ # in a web browser 
```

### Online Web-portal
http://cyto-sv-ml.b12x.org/


# Related Resource Download
## SV Calling Pipeline Download
```
Parliament: docker pull docker.io/dnanexus/parliament2:latest #https://github.com/dnanexus/parliament2
ChromoSeq: docker pull  docker.io/zatawada/docker-basespace_chromoseq_v2:master #https://github.com/genome/docker-basespace_chromoseq
```

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

# CYTO-SV-ML
#### CYTO-SV-ML PIPELINE WORKFLOW:
![CYTO-SV-ML PIPELINE WORKFLOW](Workflow.png)

# Run Shiny Application 

## Install R Packages

CYTO-SV-ML Shiny Application requires the following R packages to be installed.
```
RUN install2.r remotes
RUN R -e 'install.packages("lineupjs",dependencies=TRUE, repos="https://cran.r-project.org/")'
#RUN R -e 'remotes::install_github("tzhang-nmdp/lineup_htmlwidget")'
RUN install2.r ggforce plotly DT ggplot2 stringr shinydashboard shinyjs shinyWidgets dplyr Hmisc
```

# Local Run

### Change directory to `cyto-sv-ml`
```
cd cyto-sv-ml
```

### Start Shiny Application from R Console

```
shiny::runApp()
```

# Online Web-portal
http://cyto-sv-ml.b12x.org/

# Related Resource Download
## SV Calling Pipeline Download
```
Parliament: docker pull dongwonlee/parliament2-sing:v0.12 #https://github.com/dnanexus/parliament2
ChromoSeq: https://github.com/genome/docker-basespace_chromoseq
```

## SV database Download
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

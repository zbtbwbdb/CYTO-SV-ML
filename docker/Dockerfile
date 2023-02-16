FROM conda/miniconda3
MAINTAINER "Tao Zhang tzhang@nmdp.org"

RUN conda install -n base -c conda-forge mamba

# create conda env
COPY . /CYTO-SV-ML/
WORKDIR CYTO-SV-ML

RUN mamba env create cyto-sv-ml -f cyto-sv-ml.yaml
RUN pip install --upgrade snakemake

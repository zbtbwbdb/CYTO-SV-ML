#!/bin/bash
main_dir=$1
cyto_sv_ml_dir=$2
sample_all=$3
cyto_band_file=$4
k=$5

echo ${main_dir} ${cyto_sv_ml_dir} ${sample_all}
echo "# prepare the data file for cyto-sv-ml application" && date   
python ${cyto_sv_ml_dir}/Pipeline_script/interface_data.py -t ${main_dir}'/out/'${sample_all}'/cyto_sv_ml/'${sample_all}'_trs_'${k} -n ${main_dir}'/out/'${sample_all}'/cyto_sv_ml/'${sample_all}'_nontrs_'${k} -i ${main_dir}'/out/'${sample_all}'/'${sample_all}
Rscript ${cyto_sv_ml_dir}/Pipeline_script/interface_data.R -i ${main_dir}'/out/'${sample_all}'/'${sample_all} -o ${cyto_sv_ml_dir}"/docker/cyto-sv-ml/data/cyto_sv_ml.RData" -c ${cyto_band_file}
cd ${cyto_sv_ml_dir}"/docker"
echo "# build R-shiny based user interface in docker container" && date   
sudo docker build -t cyto-sv-ml-app:${sample_all} .
# to run the docker image in the local machine and open user interface with "http://localhost:8000/"
# sudo docker run -d -p 8000:80 cyto-sv-ml-app:${sample_all}

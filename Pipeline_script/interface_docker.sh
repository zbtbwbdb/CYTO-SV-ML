#!/bin/bash
main_dir=$1
cyto_sv_ml_dir=$2
sample_all=$3

echo ${main_dir} ${cyto_sv_ml_dir} ${sample_all}
echo "# build R-shiny based user interface in docker container " && date   
python ${cyto_sv_ml_dir}/Pipeline_script/interface_data.py -t ${main_dir}'/out/cyto_sv_ml/'${sample_all}'_trs_'${k}'_trs_EXP' -n ${main_dir}'/out/cyto_sv_ml/'${sample_all}'_nontrs_'${k}'_nontrs_EXP' -i ${main_dir}'/out/'${sample_all}'.sv.all.combine_all'
Rscript ${cyto_sv_ml_dir}/Pipeline_script/interface_data.R 
cd ${cyto_sv_ml_dir}/cyto-sv-ml/
sudo docker build -t cyto-sv-ml-app:${sample_all} .
# to run the docker image in the local machine and open user interface with "http://localhost:8000/"
# sudo docker run -d -p 8000:80 cyto-sv-ml-app:${sample_all}

#!/bin/bash
main_dir=$1
sample_all=$2

echo ${sample_all}
echo "# build R-shiny based user interface in docker container " && date   
python ${cyto_sv_ml_dir}/Pipeline_script/interface_docker.py -t OUTPUT_DIR+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_trs_EXP -n  OUTPUT_DIR+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_nontrs_EXP -i outdir+'/'+str(cohort_name)+'.sv.all.combine_all'
cd ${cyto_sv_ml_dir}/cyto-sv-ml/
sudo docker build -t cyto-sv-ml-app:${cohort_name} .
# to run the docker in the local machine and open user interface with "http://localhost:8000/"
# sudo docker run -d -p 8000:80 cyto-sv-ml-app:${cohort_name}

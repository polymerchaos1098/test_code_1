#!/bin/bash

### PARAMETERS

########################################
##### 1
sample_set_name=AML_new_samples

# "tumor_xxxxx"
tumor_list=(IID_H121248_T03_02_WG01 IID_H133128_T02_02_WG01 IID_H133566_T02_02_WG01 IID_H133567_T02_02_WG01 IID_H134807_T02_02_WG01) 

# "normal_xxxxx"
normal_list=(IID_H121248_N02_02_WG01 IID_H133128_N01_02_WG01 IID_H133566_N01_02_WG01 IID_H133567_N01_02_WG01 IID_H134807_N01_02_WG01) 

# "batch_xxxxx", using the patient ID
batch_list=(IID_H121248_T03 IID_H133128_T02 IID_H133566_T02 IID_H133567_T02 IID_H134807_T02) 

########################################
########################################
########################################

# save folder
target_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/results/bcbio-GRCh37/$sample_set_name
if [ ! -d $target_folder ]; then mkdir $target_folder; fi

# tumor vs normal pair analysis
# https://stackoverflow.com/questions/27787536/how-to-pass-a-variable-containing-slashes-to-sed
# "normal_fastq_xxxxx"
normal_fastq_folder=/projectsp/foran/Team_Chan/project142_AML_6_10_2019
# replace with slash "/", add "\" before it
normal_fastq_folder=$(echo "$normal_fastq_folder" | sed 's/\//\\\//g')

# "tumor_fastq_xxxxx"
tumor_fastq_folder=/projectsp/foran/Team_Chan/project142_AML_6_10_2019
# replace with slash "/", add "\" before it
tumor_fastq_folder=$(echo "$tumor_fastq_folder" | sed 's/\//\\\//g')

#echo $tumor_fastq_folder

code_date="2019-07-16" # yyyy-mm-dd, "code_date_xxxxx"
code_name=$sample_set_name # "code_name_xxxxx"

# template_file
template_file=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/code/bcbio/series_job_run/bam/bcbio-template-bam.yaml

# generate yaml for each tumor sample
# and create corresponding folders
#for i in ${tumor_list[@]}; do 
#for i in ${#tumor_list[@]}; do 

i=0

while [ $i -lt ${#tumor_list[@]} ]; do
	echo $i
	tumor_sample=${tumor_list[$i]}
	normal_sample=${normal_list[$i]}
	batch_id=${batch_list[$i]}
	
	final_save_dir=$target_folder/$batch_id; if [ ! -d $final_save_dir ]; then mkdir $final_save_dir; fi
	# "final_save_dir_xxxxx"
	final_save_dir=$target_folder/$batch_id/final; if [ ! -d $final_save_dir ]; then mkdir $final_save_dir; fi
	final_save_dir=$(echo "$final_save_dir" | sed 's/\//\\\//g')

	# save this yaml file in below folder
	config_dir=$target_folder/$batch_id/config; if [ ! -d $config_dir ]; then mkdir $config_dir; fi

	# replace in template
	cp $template_file $config_dir/$batch_id-SV-dectect.yaml
	
	# replace with slash "/", add "\" before it
	sed -i "s/normal_fastq_xxxxx/$normal_fastq_folder/g" $config_dir/$batch_id-SV-dectect.yaml
	sed -i "s/tumor_fastq_xxxxx/$tumor_fastq_folder/g" $config_dir/$batch_id-SV-dectect.yaml
	
	sed -i "s/normal_xxxxx/$normal_sample/g" $config_dir/$batch_id-SV-dectect.yaml
	sed -i "s/tumor_xxxxx/$tumor_sample/g" $config_dir/$batch_id-SV-dectect.yaml
	
	sed -i "s/batch_xxxxx/$batch_id/g" $config_dir/$batch_id-SV-dectect.yaml
	sed -i "s/code_date_xxxxx/$code_date/g" $config_dir/$batch_id-SV-dectect.yaml
	sed -i "s/code_name_xxxxx/$code_name/g" $config_dir/$batch_id-SV-dectect.yaml
	sed -i "s/final_save_dir_xxxxx/$final_save_dir/g" $config_dir/$batch_id-SV-dectect.yaml

	i=$(($i+1))
	
done


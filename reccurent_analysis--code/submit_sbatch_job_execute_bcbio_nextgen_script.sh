#!/bin/bash

### PARAMETERS

# save log files
cd /projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/code/bcbio/series_job_run/bam/log

########################################
##### 1
sample_set_name=AML_new_samples

tumor_list=(IID_H121248_T03 IID_H133128_T02 IID_H133566_T02 IID_H133567_T02 IID_H134807_T02) 

tumor_list=(IID_H133128_T02) 

########################################
########################################
########################################
#
process_number=12


# parallel running. not using. don't know how to use in the perceval.
parallelization_type="ipython"
ipython_parallelization_scheduler="pbspro"

#
work_folder=/scratch/chanc3/bcbio-result-intermediate-work-folder/$sample_set_name
if [ ! -d $work_folder ]; then mkdir $work_folder; fi

#
config_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/results/bcbio-GRCh37/$sample_set_name

run_dir=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/code/bcbio/series_job_run/bam

for i in ${tumor_list[@]}; do 

	intermediate_folder=$work_folder/$i; if [ ! -d $intermediate_folder ]; then mkdir $intermediate_folder; fi
	intermediate_folder=$intermediate_folder/work; if [ ! -d $intermediate_folder ]; then mkdir $intermediate_folder; fi
	
	job_yaml_file=$config_folder/$i/config/$i-SV-dectect.yaml
	
	### command: --partition perceval   --cpus-per-task $process_number  --mem=124000
	#	--partition perceval --nodes 1 --cpus-per-task 12 --mem-per-cpu=124000  \
    # 5-00:00:00
	sbatch \
	--partition perceval --nodes 1 --cpus-per-task 1 --ntasks 1 --mem-per-cpu=124000  \
	--time 7-00:00:00 \
    --export=ALL --job-name="$i" --error="$i".err --output="$i".out \
	$run_dir/execute_bcbio_nextgen_script.sh $intermediate_folder $job_yaml_file $process_number $parallelization_type $ipython_parallelization_scheduler 
		
	#sbatch \
	#--partition mem_p --nodes 1 --cpus-per-task 1 --ntasks 1 --mem-per-cpu=1240000  \
	
	# pass args:
	# $1: intermediate_folder in /scratch/chanc3; 
	# $2: job_yaml_file; 
	# $3: total number of processes to use, -n; 
	# $4: type of parallelization to use, -t; 
	# $5: -s, scheduler for ipython parallelization (lsf, sge, slurm, torque, pbspro); 
	
done



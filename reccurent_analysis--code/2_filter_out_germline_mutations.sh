#!/usr/bin/bash

# 1. put germline & somatic mutations in differrent folders

###
# file name
file_name_list=(IID_H121248_T03 IID_H133128_T02 IID_H133566_T02 IID_H133567_T02 IID_H134807_T02 I-H-121248-T2-1-D1-1 I-H-121248-T3-1-D1-1 I-H-133128-T1-1-D1-1 I-H-133128-T2-1-D1-1 I-H-133128-T3-1-D1-1 I-H-133566-T1-1-D1-1 I-H-133566-T1-2-D1-1 I-H-133566-T2-1-D1-1 I-H-133567-T1-1-D1-1 I-H-133567-T2-1-D1-1 I-H-133567-T3-1-D1-1 I-H-134807-T1-1-D1-1 I-H-134807-T2-1-D1-1) 


######
data_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/reccurent_analysis/data/mutations_only

######
save_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/reccurent_analysis/data/mutations_only

######
summary_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/reccurent_analysis/results

echo -e "sample_ID"'\t'"germline_number"'\t'"somatic_number" > $summary_folder/germline_somatic_mutations_summary.txt

#
if [ ! -d "$save_folder/somatic" ]; then mkdir $save_folder/somatic; fi
if [ ! -d "$save_folder/germline" ]; then mkdir $save_folder/germline; fi

######
i=0

while [ $i -lt ${#file_name_list[@]} ]; do

#while [ $i -lt 1 ]; do

	echo $i

	# 
	file_name=${file_name_list[$i]}

	cd $data_folder
	
	# save
	grep -P "\trs" data-$file_name-ensemble-annotated.vcf > $save_folder/germline/germline-$file_name-ensemble-annotated.vcf
	
	grep -P -v "\trs" data-$file_name-ensemble-annotated.vcf > $save_folder/somatic/somatic-$file_name-ensemble-annotated.vcf
	
	# summary
	germline_number=$(grep -P "\trs" data-$file_name-ensemble-annotated.vcf | wc -l) 
	
	somatic_number=$(grep -P -v "\trs" data-$file_name-ensemble-annotated.vcf | wc -l) 

	echo -e $file_name'\t'$germline_number'\t'$somatic_number >> $summary_folder/germline_somatic_mutations_summary.txt

	i=$(($i+1))

done






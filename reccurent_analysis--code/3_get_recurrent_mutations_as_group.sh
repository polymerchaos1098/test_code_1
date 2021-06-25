#!/usr/bin/bash

# 1. get mutations for given grouped file name

###
# file name
grouped_sample_list=(I-H-121248-T2-1-D1-1 I-H-121248-T3-1-D1-1) 
group_ID="I-H-121248"
mutation_type="somatic"

######
data_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/reccurent_analysis/data/mutations_only/$mutation_type

######
save_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/reccurent_analysis/results/recurrent_mutations/$mutation_type

if [ ! -d "$save_folder" ]; then mkdir $save_folder; fi
if [ ! -d "$save_folder/$group_ID" ]; then mkdir $save_folder/$group_ID; fi

######
summary_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/reccurent_analysis/results

#echo -e "sample_ID"'\t'"germline_number"'\t'"somatic_number" > $summary_folder/germline_somatic_mutations_summary.txt

#
cd $data_folder

# 1. get concatenated ID as 1,2,4,5 column combined
i=0

while [ $i -lt ${#grouped_sample_list[@]} ]; do

	echo $i

	# 
	file_name=${grouped_sample_list[$i]}

	# get concatenated ID
	cat $mutation_type-$file_name-ensemble-annotated.vcf | awk '{print  $1"-"$2"-"$4"-"$5"\t"$0}' > $save_folder/$group_ID/with_ID_$mutation_type-$file_name-ensemble-annotated.vcf

	i=$(($i+1))

done


###
# 2. get recurrent mutations
cd $save_folder/$group_ID

# the 1st file
file_name=${grouped_sample_list[0]}
cp with_ID_$mutation_type-$file_name-ensemble-annotated.vcf recurrent-with_ID_$mutation_type-$group_ID-ensemble-annotated.vcf

i=1

while [ $i -lt ${#grouped_sample_list[@]} ]; do

	echo $i

	# 
	file_name=${grouped_sample_list[$i]}
	
#https://www.unix.com/unix-for-dummies-questions-and-answers/230809-intersect-two-columns-two-separate-files.html
# the 
	awk 'NR==FNR{A[$1];next}$1 in A' recurrent-with_ID_$mutation_type-$group_ID-ensemble-annotated.vcf with_ID_$mutation_type-$file_name-ensemble-annotated.vcf > tmp-recurrent-with_ID_$mutation_type-$group_ID-ensemble-annotated.vcf 

	# rm tmp file
	cp tmp-recurrent-with_ID_$mutation_type-$group_ID-ensemble-annotated.vcf recurrent-with_ID_$mutation_type-$group_ID-ensemble-annotated.vcf
	
	rm tmp-recurrent-with_ID_$mutation_type-$group_ID-ensemble-annotated.vcf 
	
	#
	i=$(($i+1))

done

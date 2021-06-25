#!/bin/bash

data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/ALL---ORIGINAL---data/Manta-original"
#result_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV"

mkdir -p $result_folder

# "save_list"
save_list=(IID_H134807_T02 IID_H133566_T02 IID_H133567_T02  IID_H133128_T02 IID_H121248_T03 I-H-113385-T1-1-D1-1_vs_I-H-113385-N1-1-D1-1  I-H-113386-T1-1-D1-1_vs_I-H-113386-N1-1-D1-1 I-H-113387-T1-1-D1-1_vs_I-H-113387-N1-1-D1-1 I-H-113388-T1-1-D1-1_vs_I-H-113388-N1-1-D1-1 I-H-113389-T1-1-D1-1_vs_I-H-113389-N1-1-D1-1 I-H-113390-T1-1-D1-1_vs_I-H-113390-N1-1-D1-1 I-H-112855-T1-1-D1-1_vs_I-H-112855-N1-1-D1-1) 

# "file_list"
#file_list=(134807 133566 133567 133128 121248 113385  113386 113387 113388 113389 113390 112855)


###
i=0

while [ $i -lt ${#save_list[@]} ]; do

	echo "---------"$i"----------"

	# 1. names
	#file_name=${file_list[$i]}
	save_name=${save_list[$i]}
    echo $save_name
	
	cd $data_folder/"$save_name"/results/variants
    cp somaticSV.vcf.gz backup---somaticSV.vcf.gz
	zcat somaticSV.vcf.gz > somaticSV.vcf
	
	i=$(($i+1))

done



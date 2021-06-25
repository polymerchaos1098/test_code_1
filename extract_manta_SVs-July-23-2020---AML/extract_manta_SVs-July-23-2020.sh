#!/bin/bash


# 1. copy manta somaticSV files into 1 folder
# 2. extract header & mutations data separatively with columns

###
# "save_list"
save_list=(IID_H134807_T02 IID_H133566_T02 IID_H133567_T02  IID_H133128_T02 IID_H121248_T03 I-H-113385-T1-1-D1-1_vs_I-H-113385-N1-1-D1-1  I-H-113386-T1-1-D1-1_vs_I-H-113386-N1-1-D1-1 I-H-113387-T1-1-D1-1_vs_I-H-113387-N1-1-D1-1 I-H-113388-T1-1-D1-1_vs_I-H-113388-N1-1-D1-1 I-H-113389-T1-1-D1-1_vs_I-H-113389-N1-1-D1-1 I-H-113390-T1-1-D1-1_vs_I-H-113390-N1-1-D1-1 I-H-112855-T1-1-D1-1_vs_I-H-112855-N1-1-D1-1) 

# "file_list"
sample_list=(134807 133566 133567 133128 121248 113385  113386 113387 113388 113389 113390 112855)

######
vcf_folder=/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/ALL---ORIGINAL---data/Manta-original

######
save_folder=/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV
analysis_folder=/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary

mkdir -p $save_folder; mkdir -p $analysis_folder

###### step 1: copy files
######
i=0

while [ $i -lt ${#save_list[@]} ]; do

	echo $i

	# 1. copy SomaticSV files into 1 folder
	sample_name=${sample_list[$i]}
	save_name=${save_list[$i]}

	# copy 
	cp $vcf_folder/$save_name/results/variants/somaticSV.vcf.gz $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf.gz
	cp $vcf_folder/$save_name/results/variants/somaticSV.vcf $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf

	i=$(($i+1))

done

###### step 2: 1). filter out "PASS" candidates;
######         2). remove the header and combine all samples' "PASS" candidates into 1 file;
mkdir $analysis_folder/only-PASS-candidates
mkdir $analysis_folder/only-PASS-candidates-for-snpEff

# record the somatic SV #
file_1=$analysis_folder/somaticSV_number_statistics.txt
echo "Sample Total_somaticSV SomaticSV_with_PASS" > $file_1

# combine all "PASS" candidates into a file
file_2=$analysis_folder/ALL_PASS_somaticSV_candidates.txt
#touch $file_2
echo "Sample CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOR" > $file_2

# for snpEff
#file_3=$analysis_folder/ALL_PASS_somaticSV_candidates-for-snpEff.vcf
#touch $file_3

###
i=0

while [ $i -lt ${#sample_list[@]} ]; do

	echo $i

	# 1. names
	sample_name=${sample_list[$i]}

	# copy 
    total_number=$(cat $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf | grep -v "#" | wc -l)
    pass_candidate_number=$(cat $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf | grep -v "#" | grep "PASS" | wc -l)
	
	comment_number=$(cat $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf | grep "#" | wc -l)

	# 
	cat $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf | grep "#CHROM" > $analysis_folder/only-PASS-candidates/only-PASS-candidates---"$sample_name"-tumor-WGS---somaticSV.vcf

	cat $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf | grep -v "#" | grep "PASS" >> $analysis_folder/only-PASS-candidates/only-PASS-candidates---"$sample_name"-tumor-WGS---somaticSV.vcf
	
	# for snpEff
	cat $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf | grep "#" > $analysis_folder/only-PASS-candidates-for-snpEff/only-PASS-candidates---"$sample_name"-tumor-WGS---somaticSV.vcf

	cat $save_folder/"$sample_name"-tumor-WGS---somaticSV.vcf | grep -v "#" | grep "PASS" >> $analysis_folder/only-PASS-candidates-for-snpEff/only-PASS-candidates---"$sample_name"-tumor-WGS---somaticSV.vcf
	
	#
	echo "$sample_name $total_number $pass_candidate_number" >> $file_1

	# line by line
	j=1
	while read -r line; do
	    #echo $j
	    if [ $j -gt 1 ]; then
	        echo $sample_name" "$line >> $file_2
			#echo $line >> $file_3
	    fi
	    j=$(($j+1))
	done < $analysis_folder/only-PASS-candidates/only-PASS-candidates---"$sample_name"-tumor-WGS---somaticSV.vcf
	
	
	i=$(($i+1))

done

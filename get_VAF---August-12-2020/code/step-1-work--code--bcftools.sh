#
save_folder=/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/get_VAF---August-12-2020/result/data
mkdir -p $save_folder

ying_vcf_folder=/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/somatic_mutation---August-12-2020/result-filter/filter_test_1/all_files

########################################
# "tumor"
tumor_list=(IID_H121248_T03 IID_H133128_T02 IID_H133566_T02 IID_H133567_T02 IID_H134807_T02) 

label_list=(121248 133128 133566 133567 134807) 

i=0

while [ $i -lt ${#tumor_list[@]} ]; do
	
	echo $i
	
    #47123.somatic.filtered.snpEff.PASS.vcf.gz
    tumor_sample=${tumor_list[$i]}
	#echo $tumor_sample
    file_name=diploid--"$tumor_sample"-mutect2-annotated.vcf 
	#echo $file_name
	
	# 
	save_name=${label_list[$i]}_VAF.txt; 
	#echo $save_name
	
	#
	cat $ying_vcf_folder/$file_name | grep "#CHROM" | awk '{print $1, $2, $(NF-1), $NF}' > $save_folder/$save_name

	# only use "biallelic" SNP: very few have "triallelic" SNP; separate with ","
	cat $ying_vcf_folder/$file_name | bcftools query -f '%CHROM %POS[ %AF]\n' | grep -v "," >> $save_folder/$save_name
	
	#
	i=$(($i+1))

done



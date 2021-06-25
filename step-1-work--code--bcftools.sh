#
save_folder=/projectsp/foran/Team_Chan/Hua_files/Human_LFS_VAF-9-5-2019/result/data

vcf_folder=/projectsp/foran/Team_Chan/Hua_files/Human_LFS_VAF-9-5-2019/vcf

########################################
##### 1
sample_set_name=LFS

# "tumor"
tumor_list=(RGH7005504-T) 

normal_list=(RGH7005504-N) 

label_list=(RGH7005504-T) 

#
i=0

while [ $i -lt ${#tumor_list[@]} ]; do
	
	echo $i
	
    #47123.somatic.filtered.snpEff.PASS.vcf.gz
    tumor_sample=${tumor_list[$i]}
	normal_sample=${normal_list[$i]}

	#echo $tumor_sample
    file_name="$tumor_sample"--"$normal_sample".snv.union.v5B.annotated.vcf; 
	#echo $file_name
	
	# 
	save_name=${label_list[$i]}_VAF.txt; 
	#echo $save_name
	
	#
	#zcat $vcf_folder/$file_name | grep "#CHROM" | awk '{print $1, $2, $(NF-1), $NF}' > $save_folder/$save_name
	cat $vcf_folder/$file_name | grep "#CHROM" | awk '{print $1, $2, $(NF-1), $NF}' > $save_folder/$save_name

	# remove germline
	
	cat $vcf_folder/$file_name | grep -P -v "\trs" > $vcf_folder/remove-germline-$file_name
	
	cat $vcf_folder/$file_name | grep -P "\trs" > $vcf_folder/only-germline-$file_name

	# only use "biallelic" SNP: very few have "triallelic" SNP; separate with ","
	#zcat $vcf_folder/$file_name | bcftools query -f '%CHROM %POS[ %AF]\n' | grep -v "," >> $save_folder/$save_name
	#cat $vcf_folder/remove-germline-$file_name | bcftools query -f '%CHROM %POS[ %mutect_AD]\n' | grep -v "," >> $save_folder/$save_name
	cat $vcf_folder/remove-germline-$file_name | bcftools query -f '%CHROM %POS[ %mutect_AD]\n' >> $save_folder/$save_name
	cat $save_folder/$save_name | grep -v '\.' > $save_folder/remove-no-mutect_AD-$save_name

	#
	i=$(($i+1))

done


#
mutation_number=$(cat $save_folder/remove-no-mutect_AD-$save_name | wc -l)
echo "mutation #:" $mutation_number

echo -e "#CHROM\tPOS\tNormal_RGH7005504-N\tTumor_RGH7005504-T" > $save_folder/final-all-$save_name
echo -e "#CHROM\tPOS\tNormal_RGH7005504-N\tTumor_RGH7005504-T" > $save_folder/final-only-AF-$save_name

j=2 # file line begins with 1

# each mutation have AD, but have no "AF"
# search AD and AF
# if having no AD, skip
# https://stackoverflow.com/questions/8009664/how-to-split-a-delimited-string-into-an-array-in-awk
while [ $j -le $mutation_number ]; do
#while [ $j -le 10 ]; do	

	chr=$(awk -v var="$j" 'NR==var' $save_folder/remove-no-mutect_AD-$save_name | awk '{printf $1}')
	pos=$(awk -v var="$j" 'NR==var' $save_folder/remove-no-mutect_AD-$save_name | awk '{printf $2}')

	normal_column=$(awk -v var="$j" 'NR==var' $save_folder/remove-no-mutect_AD-$save_name | awk '{printf $3}')
	tumor_column=$(awk -v var="$j" 'NR==var' $save_folder/remove-no-mutect_AD-$save_name | awk '{printf $4}')

	ref_1=$(echo $normal_column | awk '{
	n = split($0, t, ",");print t[1]}')
	
	alt_1=$(echo $normal_column | awk '{
	n = split($0, t, ",");print t[2]}')
	
	ref_2=$(echo $tumor_column | awk '{
	n = split($0, t, ",");print t[1]}')
	
	alt_2=$(echo $tumor_column | awk '{
	n = split($0, t, ",");print t[2]}')	
	
	echo $j #$normal_column $tumor_column $ref_1 $alt_1 $ref_2 $alt_2 $(echo "scale=4; $alt_1/($ref_1+$alt_1)" | bc -l) $(echo "scale=4; $alt_2/($ref_2+$alt_2)" | bc -l)
		
	echo -e $chr"\t"$pos"\t"$normal_column"\t"$tumor_column"\t"$ref_1"\t"$alt_1"\t"$ref_2"\t"$alt_2"\t"$(echo "scale=4; $alt_1/($ref_1+$alt_1)" | bc -l)"\t"$(echo "scale=4; $alt_2/($ref_2+$alt_2)" | bc -l) >> $save_folder/final-all-$save_name
	
	echo -e $chr"\t"$pos"\t"$(echo "scale=4; $alt_1/($ref_1+$alt_1)" | bc -l)"\t"$(echo "scale=4; $alt_2/($ref_2+$alt_2)" | bc -l)  >> $save_folder/final-only-AF-$save_name

	j=$(($j+1))	

done

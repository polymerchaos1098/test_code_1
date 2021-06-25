#!/usr/bin/bash

#tumor_list=(GPK0267-2000 GPK0266-2000 GPK0266-2001 GPK0266-2002 GPK0266-2003 GPK0266-2004 GPK0266-2005 GPK0266-2006 GPK0266-2007 GPK0265-2000 GPK0264-2000 GPK0264-2001 GPK0263-2004 GPK0263-2003 GPK0263-2002 GPK0263-2001 GPK0263-2000 31252T 32545T 34854T 33175T 38783T 33014T 36266T 32537T 34059T 47123T 33059T 36585T) 

#
#tumor_list=(31252T 32545T 34854T 33175T 38783T 33014T 36266T 32537T 34059T 47123T 33059T 36585T) 
#tumor_list=(47123T)
tumor_list=(GPK0267-2000 GPK0266-2000 GPK0266-2001 GPK0266-2002 GPK0266-2003 GPK0266-2004 GPK0266-2005 GPK0266-2006 GPK0266-2007 GPK0265-2000 GPK0264-2000 GPK0264-2001 GPK0263-2004 GPK0263-2003 GPK0263-2002 GPK0263-2001 GPK0263-2000) 

### 
#filter_test_1:
#i). vcf "FILTER" column as "PASS"; ii). TUMOR Depth, REF + ALT > 10 & REF + ALT < 200; iii). NORMAL Depth, REF + ALT > 8 & REF + ALT < 200; iv). Tumor AF > 0.2; 5). Normal AF < 0.05


### vcf_folder
#vcf_folder=/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/mutation_signature/test-new-condition-June-30-2020---oncocytomas/result-filter
vcf_folder=/scratch/hk618/test-new-condition-June-30-2020---oncocytomas/result-filter

#save_folder=/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/mutation_signature/test-new-condition-June-30-2020---oncocytomas/result---step-4-mutation-analysis/all_vcf_from_folder---result-filter
save_folder=/scratch/hk618/test-new-condition-June-30-2020---oncocytomas/result---step-4-mutation-analysis/all_vcf_from_folder---result-filter
#mkdir $save_folder
mkdir -p $save_folder/tmp

########################################
#cd /scratch/hk618/test-new-condition-June-30-2020---oncocytomas/updated-step-4-mutation-analysis
cd /projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/mutation_signature/test-new-condition-June-30-2020---oncocytomas/updated-step-4-mutation-analysis

#
echo -e "Sample"'\t'"CHROM"'\t'"POS"'\t'"ID"'\t'"REF"'\t'"ALT"'\t'"QUAL"'\t'"FILTER"'\t'"INFO"'\t'"FORMAT"'\t'"NORMAL"'\t'"TUMOR" > $save_folder/../oncocytomas_all_PASS_vcf_in_a_file.txt

#
i=0

while [ $i -lt ${#tumor_list[@]} ]; do
	
	echo $i

	#
	tumor_sample=${tumor_list[$i]}
	
	# filter_condition_1
	filter_condition="filter_test_1"

	cp $vcf_folder/$filter_condition/$tumor_sample/PASS_"$tumor_sample".somatic.filtered.snpEff.PASS.vcf $save_folder/filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.vcf
    
	cat $vcf_folder/$filter_condition/$tumor_sample/PASS_"$tumor_sample".somatic.filtered.snpEff.PASS.vcf | grep -v "#" > $save_folder/tmp/filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.txt
	
	# rows
	cat $save_folder/tmp/filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.txt | wc -l
	row_num=$(cat $save_folder/tmp/filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.txt | wc -l)

	# write the sample name
	touch $save_folder/tmp/sample_name---filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.txt
	j=0
	
	while [ $j -lt $row_num ]; do
	    echo $tumor_sample >> $save_folder/tmp/sample_name---filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.txt
		j=$(($j+1))
	done	
	
	#paste $result_folder/tmp/sample_name-$data_file $result_folder/tmp/$data_file | column -s $'\t' -t > $result_folder/$data_file
	
	#awk '{print $0, FILENAME}' $result_folder/tmp/sample_name-$data_file $result_folder/tmp/$data_file > $result_folder/$data_file
	awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],"\t",$0}' $save_folder/tmp/sample_name---filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.txt $save_folder/tmp/filtered-out---"$tumor_sample".somatic.filtered.snpEff.PASS.txt >> $save_folder/../oncocytomas_all_PASS_vcf_in_a_file.txt
		
	#
	i=$(($i+1))

done


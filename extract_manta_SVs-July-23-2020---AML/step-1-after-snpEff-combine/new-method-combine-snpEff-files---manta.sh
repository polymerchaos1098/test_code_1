#!/bin/bash

# 
data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/snpEff-result-only-PASS-candidates"

result_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method"
mkdir $result_folder
mkdir $result_folder/tmp
mkdir $result_folder/final


# "tumor_xxxxx"
tumor_list=(121248 133128 133566 133567 134807 112855 113385 113386 113387 113388 113389 113390) 


###
i=0

while [ $i -lt ${#tumor_list[@]} ]; do

	echo "---------"$i"----------"

	# 1. names
	tumor_sample=${tumor_list[$i]}

	data_file=$(echo snpEff---only-PASS-candidates---"$tumor_sample"-tumor-WGS---somaticSV.vcf)

    #
    cat $data_folder/$data_file | grep -v "#" > $result_folder/tmp/$data_file
	
	# rows
	cat $result_folder/tmp/$data_file | wc -l
	row_num=$(cat $result_folder/tmp/$data_file | wc -l)

	# write the sample name
	touch $result_folder/tmp/sample_name-$data_file
	j=0
	
	while [ $j -lt $row_num ]; do
	    echo $tumor_sample >> $result_folder/tmp/sample_name-$data_file
		j=$(($j+1))
	done	
	
	#paste $result_folder/tmp/sample_name-$data_file $result_folder/tmp/$data_file | column -s $'\t' -t > $result_folder/$data_file
	
	#awk '{print $0, FILENAME}' $result_folder/tmp/sample_name-$data_file $result_folder/tmp/$data_file > $result_folder/$data_file
	awk 'NR==FNR{a[NR]=$0;next}{print a[FNR],"\t",$0}' $result_folder/tmp/sample_name-$data_file $result_folder/tmp/$data_file > $result_folder/final/$data_file
	
	i=$(($i+1))

done
 
#
echo -e "Sample"'\t'"CHROM"'\t'"POS"'\t'"ID"'\t'"REF"'\t'"ALT"'\t'"QUAL"'\t'"FILTER"'\t'"INFO"'\t'"FORMAT"'\t'"NORMAL"'\t'"TUMOR" > $result_folder/manta_SVs_with_PASS---snpEfff.txt
cat $result_folder/final/* >> $result_folder/manta_SVs_with_PASS---snpEfff.txt

cat $result_folder/manta_SVs_with_PASS---snpEfff.txt | wc -l

#
#echo -e "Sample"'\t'"CHROM"'\t'"POS"'\t'"ID"'\t'"REF"'\t'"ALT"'\t'"QUAL"'\t'"FILTER"'\t'"INFO"'\t'"FORMAT"'\t'"NORMAL"'\t'"TUMOR" > $result_folder/chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt

# not work
#cat $result_folder/manta_SVs_with_PASS---snpEfff.txt | grep -v "#" | grep -E "chr1	|chr2	|chr3	|chr4	|chr5	|chr6	|chr7	|chr8	|chr9	|chr10	|chr11	|chr12	|chr13	|chr14	|chr15	|chr16	|chr17	|chr18	|chr19	|chr20	|chr21	|chr22	|chrX	|chrY	|chrM	" >> $result_folder/chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt

# try this
cat $result_folder/manta_SVs_with_PASS---snpEfff.txt | grep -v "#" | grep -v "chrUn_" | grep -v "_random" | grep -v "_alt" >> $result_folder/chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt


cat $result_folder/chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | wc -l



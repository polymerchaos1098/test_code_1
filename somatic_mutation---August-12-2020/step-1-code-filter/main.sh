#!/usr/bin/bash

#
tumor_list=(IID_H121248_T03 IID_H133128_T02 IID_H133566_T02 IID_H133567_T02 IID_H134807_T02) 

### 
#filter_test_1:
#i). vcf "FILTER" column as "PASS"; ii). TUMOR Depth, REF + ALT > 10 & REF + ALT < 200; iii). NORMAL Depth, REF + ALT > 8 & REF + ALT < 200; iv). Tumor AF > 0.2; 5). Normal AF < 0.05

########################################
cd /projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/somatic_mutation---August-12-2020/step-1-code-filter

i=0

while [ $i -lt ${#tumor_list[@]} ]; do
	
	echo $i

	#
	tumor_sample=${tumor_list[$i]}
	
	# filter_condition_1
	filter_condition="filter_test_1"
	filter_1_value="PASS"
	filter_2_value=10
	filter_3_value=200
	filter_4_value=8
	filter_5_value=200
	filter_6_value=0.20
	filter_7_value=0.05

	./get_filters.sh $tumor_sample $filter_condition $filter_1_value $filter_2_value $filter_3_value $filter_4_value $filter_5_value $filter_6_value $filter_7_value

	#
	i=$(($i+1))

done


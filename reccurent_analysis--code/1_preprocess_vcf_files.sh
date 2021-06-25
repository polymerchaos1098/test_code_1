#!/usr/bin/bash

# 1. copy vcf.gz files into 1 folder and unzip
# 2. extract header & mutations data separatively with columns

###
# file name
file_name_list=(IID_H121248_T03 IID_H133128_T02 IID_H133566_T02 IID_H133567_T02 IID_H134807_T02 I-H-121248-T2-1-D1-1 I-H-121248-T3-1-D1-1 I-H-133128-T1-1-D1-1 I-H-133128-T2-1-D1-1 I-H-133128-T3-1-D1-1 I-H-133566-T1-1-D1-1 I-H-133566-T1-2-D1-1 I-H-133566-T2-1-D1-1 I-H-133567-T1-1-D1-1 I-H-133567-T2-1-D1-1 I-H-133567-T3-1-D1-1 I-H-134807-T1-1-D1-1 I-H-134807-T2-1-D1-1) 


# "tumor_xxxxx"
tumor_list=(IID_H121248_T03_02_WG01 IID_H133128_T02_02_WG01 IID_H133566_T02_02_WG01 IID_H133567_T02_02_WG01 IID_H134807_T02_02_WG01 I-H-121248-T2-1-D1-1_5ebd0be6ea371bff_W0036501F_251195_Blast_D_IGO_07224_G_1_S9_001 I-H-121248-T3-1-D1-1_120d8e4e637b0c98_W0038835F_M16-14754_IGO_07224_M_2_S66_001 I-H-133128-T1-1-D1-1_5a6ca23ca016b95f_W0036509F_F17-4614_Blast_D_IGO_07224_G_9_S17_001 I-H-133128-T2-1-D1-1_4ce4e6f9c5cb0f08_W0038838F_311472_IGO_07224_M_4_S68_001 I-H-133128-T3-1-D1-1_6ef3d3021bdfad0a_361878_IGO_07224_Y_8_S65_001 I-H-133566-T1-1-D1-1_5419b406cb8c400d_W0038080F_F17-2460_Blast_25P_D_IGO_07224_I_8_S24_001 I-H-133566-T1-2-D1-1_8af1036c4cba10c_W0038081F_F17-2460_Blast_25N_D_IGO_07224_I_9_S25_001 I-H-133566-T2-1-D1-1_1eb393fdfc605c26_282955_IGO_07224_O_3_S20_001 I-H-133567-T1-1-D1-1_54a91c3d0a9d8d61_W0038082F_F17-4344_Blast_D_IGO_07224_I_10_S16_001 I-H-133567-T2-1-D1-1_6340aeeec7f6df9c_298215_IGO_07224_O_13_S24_001 I-H-133567-T3-1-D1-1_7d03ad0c364288a3_M18-621_IGO_07224_Y_18_S58_001 I-H-134807-T1-1-D1-1_605dd86f38bd104a_W0042912F_344629_Blasts_IGO_07224_Q_4_S35_001 I-H-134807-T2-1-D1-1_5e8cb2ed5106605b_M17-23672_IGO_07224_Y_2_S59_001) 

# "normal_xxxxx"
normal_list=(IID_H121248_N02_02_WG01 IID_H133128_N01_02_WG01 IID_H133566_N01_02_WG01 IID_H133567_N01_02_WG01 IID_H134807_N01_02_WG01 I-H-121248-N2-1-D1-1_459f9a005d6a101f_W0038836F_127205_IGO_07224_M_6_S47_001 I-H-121248-N2-1-D1-1_459f9a005d6a101f_W0038836F_127205_IGO_07224_M_6_S47_001 I-H-133128-N1-1-D1-1_493564628c966c83_W0038837F_162449_IGO_07224_M_7_S48_001 I-H-133128-N1-1-D1-1_493564628c966c83_W0038837F_162449_IGO_07224_M_7_S48_001 I-H-133128-N1-1-D1-1_493564628c966c83_W0038837F_162449_IGO_07224_M_7_S48_001 I-H-133566-N1-1-D1-1_6fb9d144f5c6a7a5_162534_IGO_07224_O_11_S22_001 I-H-133566-N1-1-D1-1_6fb9d144f5c6a7a5_162534_IGO_07224_O_11_S22_001 I-H-133566-N1-1-D1-1_6fb9d144f5c6a7a5_162534_IGO_07224_O_11_S22_001 I-H-133567-N1-1-D1-1_39f443d212b3aab6_162533_IGO_07224_O_10_S21_001 I-H-133567-N1-1-D1-1_39f443d212b3aab6_162533_IGO_07224_O_10_S21_001 I-H-133567-N1-1-D1-1_39f443d212b3aab6_162533_IGO_07224_O_10_S21_001 I-H-134807-N1-1-D1-1_3a72ecbfc446f348_195915_IGO_07224_Y_10_S48_001 I-H-134807-N1-1-D1-1_3a72ecbfc446f348_195915_IGO_07224_Y_10_S48_001)

# this are not done yet
#tumor_list=(I-H-133566-T2-1-D1-1_1eb393fdfc605c26_282955_IGO_07224_O_3_S20_001); # normal: I-H-133566-N1-1-D1-1_6fb9d144f5c6a7a5_162534_IGO_07224_O_11_S22_001; sample name: I-H-133566-T2-1-D1-1

######
data_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/results/bcbio-GRCh37/AML_new_samples

######
save_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/reccurent_analysis/data

mkdir $save_folder/vcf_original
mkdir $save_folder/mutations_only
mkdir $save_folder/header_only

######
i=0

while [ $i -lt ${#tumor_list[@]} ]; do

#while [ $i -lt 1 ]; do

	echo $i

	# 1. copy vcf.gz files into 1 folder and unzip
	file_name=${file_name_list[$i]}
	tumor_sample=${tumor_list[$i]}
	normal_sample=${normal_list[$i]}

	vcf_folder=/projectsp/foran/Team_Chan/Hua_files/analysis_for--project142_AML_6_10_2019/results/bcbio-GRCh37/AML_new_samples/$file_name/final/2019-07-16_AML_new_samples
	
	# replace in template
	cd $vcf_folder
	cp $file_name-ensemble-annotated.vcf.gz $save_folder/vcf_original/$file_name-ensemble-annotated.vcf.gz
	
	cd $save_folder/vcf_original;
	gunzip $file_name-ensemble-annotated.vcf.gz
	
	# 2. extract header & mutations data separatively with columns
	grep "#" $file_name-ensemble-annotated.vcf | grep -v "#CHROM" > $save_folder/header_only/header-$file_name-ensemble-annotated.vcf

	grep "#CHROM" $file_name-ensemble-annotated.vcf > $save_folder/header_only/column_name-$file_name-ensemble-annotated.vcf
	
	grep -v "#" $file_name-ensemble-annotated.vcf > $save_folder/mutations_only/data-$file_name-ensemble-annotated.vcf
	
	#
	col_list=$(grep "#CHROM" $file_name-ensemble-annotated.vcf)
	#echo $col_list
	
	normal_location=$(echo $col_list | grep -bo $normal_sample | sed 's/:.*$//')
	#echo "$normal_location"
	tumor_location=$(echo $col_list | grep -bo $tumor_sample | sed 's/:.*$//')
	#echo "$tumor_location"

	#
	cd $save_folder
	
	echo $col_list >> checking_if_normal_tumor_columns_swapped.txt
	
	if [ $normal_location -gt $tumor_location ]; then 
		echo "normal_location" "tumor_location" "$normal_location" "$tumor_location"
		echo $file_name >> sample_list_with_swapped_normal_tumor_columns.txt	
	fi
	
	i=$(($i+1))

done



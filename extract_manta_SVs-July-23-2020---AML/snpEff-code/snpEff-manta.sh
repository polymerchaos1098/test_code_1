#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -p main
#SBATCH --mem=124000
#SBATCH --constraint=oarc 
#SBATCH --export=ALL
#SBATCH -D /projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/snpEff-code
#SBATCH -o snpeff.out
#SBATCH -e snpeff.err
# this may work, check

# 
snpEff_folder="/cache/home/hk618/tools/snpEff"

data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/only-PASS-candidates-for-snpEff"

result_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/snpEff-result-only-PASS-candidates"
mkdir $result_folder

# "tumor_xxxxx"
tumor_list=(121248 133128 133566 133567 134807 112855 113385 113386 113387 113388 113389 113390) 

###
i=0

while [ $i -lt ${#tumor_list[@]} ]; do

echo $i

# 1. names
tumor_sample=${tumor_list[$i]}

data_file=$(echo only-PASS-candidates---"$tumor_sample"-tumor-WGS---somaticSV.vcf)
save_name=$(echo snpEff---only-PASS-candidates---"$tumor_sample"-tumor-WGS---somaticSV.vcf)
#echo $data_file
#echo $save_name
#
java -Xmx16g -jar $snpEff_folder/snpEff.jar GRCh37.75 $data_folder/$data_file > $result_folder/$save_name

i=$(($i+1))

done



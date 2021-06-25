#!/usr/bin/bash

#filter 1: "PASS" in vcf "FILTER" column
#filter 2: Alt (mutations) reads > 6 in tumor sample.     Check AD (Allelic depths) term in tumor sample.
#filter 3: AF > 0.2 in tumor sample.                                         Check AF (Allele fractions) term in tumor sample.
#filter 4: AF < 0.01 in NORMAL sample.                                         Check AF (Allele fractions) term in NORMAL sample.

# genotype information
# https://www.internationalgenome.org/wiki/Analysis/vcf4.0/
###
#sample_name=47123T; #sample_name=33175T;#sample_name=GPK0263-2000

# filter condition
#filter_condition="filter_condition_1"

#
#filter_test_1:
#i). vcf "FILTER" column as "PASS"; ii). TUMOR Depth, REF + ALT > 10 & REF + ALT < 200; iii). NORMAL Depth, REF + ALT > 8 & REF + ALT < 200;  iv). Normal AF < 0.05

###
sample_name=$1
filter_condition=$2
filter_1_value=$3
filter_2_value=$4
filter_3_value=$5
filter_4_value=$6
filter_5_value=$7
filter_6_value=$8
filter_7_value=$8

echo $1 $2 $3 $4 $5 $6 $7 $8 $9

### original_vcf_folder
original_vcf_folder=/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/ALL---ORIGINAL---data/Mutect2-original

### save_folder
#save_folder=/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/mutation_signature/test-new-condition-June-30-2020---oncocytomas/result-filter
save_folder=/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/somatic_mutation---August-12-2020/result-filter
mkdir -p $save_folder
summary_folder=$save_folder

# mkdir for specific sample
if [ ! -d "$save_folder/$filter_condition" ]; then mkdir $save_folder/$filter_condition; fi
if [ ! -d "$save_folder/$filter_condition/$sample_name" ]; then mkdir $save_folder/$filter_condition/$sample_name; fi
if [ ! -d "$save_folder/$filter_condition/all_files" ]; then mkdir $save_folder/$filter_condition/all_files; fi

# word dir
cd $save_folder/$filter_condition/$sample_name

########

# copy vcf into work folder
original_vcf=$sample_name-mutect2-annotated.vcf 
cp $original_vcf_folder/$sample_name-mutect2-annotated.vcf $original_vcf

# notice:
# using diploid mutations. deconstructSigs can ONLY do diploid analysis
#https://github.com/samtools/bcftools/issues/118
# bcftools +setGT
#https://samtools.github.io/bcftools/howtos/filtering.html
#diploid 
#bcftools filter -e 'GT ="."' $original_vcf > diploid_$original_vcf
# only include ".|.|."； . does not work
# GT[0] is the normal GT; GT[1] is the tumor GT

# GT genotype, encoded as alleles values separated by either of ”/” or “|”, e.g. The allele values are 0 for the reference allele (what is in the reference sequence), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1 or 1|0 etc. For haploid calls, e.g. on Y, male X, mitochondrion, only one allele value should be given. All samples must have GT call information; if a call cannot be made for a sample at a given locus, ”.” must be specified for each missing allele in the GT field (for example ./. for a diploid). The meanings of the separators are:
#    / : genotype unphased
#    | : genotype phased
# Phased data are ordered along one chromosome and so from these data you know the haplotype. Unphased data are simply the genotypes without regard to which one of the pair of chromosomes holds that allele.

###########################
###########################
###########################
### without filter
####### header
header_file=header_$original_vcf
grep "#" $original_vcf > $header_file

# original output vcf with "PASS":
PASS_file=data_PASS_$original_vcf
cat $original_vcf | grep -v "#" | grep "PASS" > $PASS_file

###########################
###########################
###########################
### REMOVE "rs" germline mutations
PASS_file=remove_germline_data_PASS_$original_vcf

cat data_PASS_$original_vcf | grep -v -P '\trs' > $PASS_file 

###########################
###########################
###########################

PASS_file=PASS_$original_vcf
cat $header_file remove_germline_data_PASS_$original_vcf > $PASS_file

#cp $PASS_file $save_folder/$filter_condition/all_files/

# get diploid mutations
bcftools filter -i 'GT[1]="0/1" | GT[1]="1/0" | GT[1]="1/1"' PASS_$original_vcf > diploid--$original_vcf

# using the diploid file
original_vcf=diploid--$original_vcf

cp $original_vcf $save_folder/$filter_condition/all_files/

# test
#bcftools filter -i 'GT="0/1"' $original_vcf > test-2.vcf
#bcftools filter -i 'GT[1]="0/1" | GT[1]="1/0" | GT[1]="1/1"' $original_vcf |grep -v "#" | head -2

####### header
header_file=header_$original_vcf
grep "#" $original_vcf > $header_file

#######
# filter 1, vcf with "PASS":
filter_1_file=$filter_1_value--$original_vcf
#grep "PASS" $original_vcf > $filter_1_file
bcftools view -f $filter_1_value $original_vcf > $filter_1_file
cat $filter_1_file | wc -l

# filter 2:
filter_2_file=tumor_depth--between-$filter_2_value-$filter_3_value--$filter_1_file
bcftools filter -i '(AD[1:0]+AD[1:1]) > '$filter_2_value'  & (AD[1:0]+AD[1:1]) < '$filter_3_value'' $filter_1_file > $filter_2_file
cat $filter_2_file | wc -l

# filter 3: AF[x:y] meaning. x [0/1], 0 is the normal columan; 1 is the tumor info columan. y [0,1,2]: an info term may have more than 1 values
filter_3_file=normal_depth--between-$filter_4_value-$filter_5_value--$filter_2_file
bcftools filter -i '(AD[0:0]+AD[0:1]) > '$filter_4_value'  & (AD[0:0]+AD[0:1]) < '$filter_5_value'' $filter_2_file > $filter_3_file
cat $filter_3_file | wc -l

# filter 4:
filter_4_file=AF_tumor_greater--$filter_3_file
bcftools filter -i 'AF[1:0] > '$filter_6_value'' $filter_3_file > $filter_4_file
cat $filter_4_file | wc -l

# filter 5:
filter_5_file=AF_normal_less--$filter_4_file
bcftools filter -i 'AF[0:0] < '$filter_7_value'' $filter_4_file > $filter_5_file
cat $filter_5_file | wc -l

# filter "chrUn_":
filter_chrUn_file=remove_chrUn--$filter_5_file
cat $filter_5_file | grep "#" > $filter_chrUn_file
# NO "chr"
#cat $filter_5_file | grep -v "#" | grep -E "1	|2	|3	|4	|5	|6	|7	|8	|9	|10	|11	|12	|13	|14	|15	|16	|17	|18	|19	|20	|21	|22	|X	|Y	|M	" >> $filter_chrUn_file
cat $filter_5_file | grep -v "#" | awk '{if ($1=="1"||$1=="2"||$1=="3"||$1=="4"||$1=="5"||$1=="6"||$1=="7"||$1=="8"||$1=="9"||$1=="10"||$1=="11"||$1=="12"||$1=="13"||$1=="14"||$1=="15"||$1=="16"||$1=="17"||$1=="18"||$1=="19"||$1=="20"||$1=="21"||$1=="22"||$1=="X"||$1=="Y"||$1=="M") print $0;}' >> $filter_chrUn_file
cat $filter_chrUn_file | wc -l


cp $filter_chrUn_file $save_folder/$filter_condition/all_files/

######################################
# count row

# . without any filter
#cat $original_vcf | grep -v "#" | wc -l

# 1. filter 1
#cat $filter_1_file | grep -v "#" | wc -l

# 2. filter 1 & 2
#cat $filter_2_file | grep -v "#" | wc -l

# 3. filter 1 & 2 & 3
#cat $filter_3_file | grep -v "#" | wc -l

# 4. filter 1 & 2 & 3 & 4
#cat $filter_4_file | grep -v "#" | wc -l



######################################
# save count row

i0=$(cat $original_vcf | grep -v "#" | wc -l); echo $i0
i1=$(cat $filter_1_file | grep -v "#" | wc -l); echo $i1
i2=$(cat $filter_2_file | grep -v "#" | wc -l); echo $i2
i3=$(cat $filter_3_file | grep -v "#" | wc -l); echo $i3
i4=$(cat $filter_4_file | grep -v "#" | wc -l); echo $i4
i5=$(cat $filter_5_file | grep -v "#" | wc -l); echo $i5
i6=$(cat $filter_chrUn_file | grep -v "#" | wc -l); echo $i6

echo -e $sample_name'\t'$i0'\t'$i1'\t'$i2'\t'$i3'\t'$i4'\t'$i5'\t'$i6 >> $summary_folder/$filter_condition--summary.txt

echo -e "ID"'\t'"No_Filter"'\t'"Filter_1"'\t'"Filter_2"'\t'"Filter_3"'\t'"Filter_4"'\t'"Filter_5"'\t'"Filter_chrUn" >> $summary_folder/$filter_condition--column_and_filter_information.txt

echo -e "--"'\t'"No_Filter"'\t'"mutations labeled as PASS"'\t'"mutations with tumor_depth--between- "$filter_2_value" - "$filter_3_value'\t'"mutations normal_depth--between- "$filter_4_value" - "$filter_5_value'\t'"mutations with AF_tumor_greater "$filter_6_value'\t'"mutations with AF_normal_less "$filter_7_value >> $summary_folder/$filter_condition--column_and_filter_information.txt

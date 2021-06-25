#!/bin/bash
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -p main
#SBATCH --mem=124000
#SBATCH --constraint=oarc 
#SBATCH --export=ALL
#SBATCH -D /projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/extract_manta_SVs-June-11-2020/step-1-after-snpEff-combine
#SBATCH -o combine.out
#SBATCH -e combine.err
# this may work, check

# 
data_folder="/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/results/extract_manta_SVs-June-11-2020/somaticSV-analysis-and-summary/snpEff-result-only-PASS-candidates"

result_folder="/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/results/extract_manta_SVs-June-11-2020/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff"
mkdir $result_folder
mkdir $result_folder/tmp

#
save_name="manta_with_SOMATICSCORE_greater_40.txt"
touch $result_folder/$save_name


# "tumor_xxxxx"
tumor_list=(GPK0267-2000 GPK0266-2000 GPK0266-2001 GPK0266-2002 GPK0266-2003 GPK0266-2004 GPK0266-2005 GPK0266-2006 GPK0266-2007 GPK0265-2000 GPK0264-2000 GPK0264-2001 GPK0263-2004 GPK0263-2003 GPK0263-2002 GPK0263-2001 GPK0263-2000 31252T 32545T 34854T 33175T 38783T 33014T 36266T 32537T 34059T 47123T 33059T 36585T) 

# "normal_xxxxx"
normal_list=(GPK0267-2001 GPK0266-2008 GPK0266-2008 GPK0266-2008 GPK0266-2008 GPK0266-2008 GPK0266-2008 GPK0266-2008 GPK0266-2008 GPK0265-2001 GPK0264-2003 GPK0264-2003 GPK0263-2005 GPK0263-2005 GPK0263-2005 GPK0263-2005 GPK0263-2005 31252N 32545N 34854N 33175N 38783N 33014N 36266N 32537N 34059N 47123N 33059N 36585N) 


###
i=0

while [ $i -lt ${#tumor_list[@]} ]; do

	echo "---------"$i"----------"

	# 1. names
	tumor_sample=${tumor_list[$i]}
	normal_sample=${normal_list[$i]}

	data_file=$(echo snpEff---only-PASS-candidates---"$tumor_sample"_vs_"$normal_sample"---somaticSV.vcf)

    #
    cat $data_folder/$data_file | grep -v "#" > $result_folder/tmp/$data_file
	
	# line by line
	j=1
	while read -r line; do
	    echo $j
	    if [ $j -gt 0 ]; then
	        echo $tumor_sample" "$line >> $result_folder/$save_name
			#echo $line >> $file_3
	    fi
	    j=$(($j+1))
	done < $result_folder/tmp/$data_file
	
	
	i=$(($i+1))

done



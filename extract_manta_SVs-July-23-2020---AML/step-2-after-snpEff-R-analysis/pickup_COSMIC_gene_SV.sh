cd /projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/results/extract_manta_SVs-June-11-2020/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method

filename="COSMIC---SOMATICSCORE-greater-44---extracted-information---Thyroid---ALL_PASS_somaticSV_candidates.txt"
touch $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33059T" | grep "chr1" | grep "205625650" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33059T" | grep "chr6" | grep "82791982" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33059T" | grep "chr6" | grep "82791991" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33059T" | grep "chr9" | grep "65437801" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33059T" | grep "chr9" | grep "65438087" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33059T" | grep "chr17" | grep "7241460" >> $filename

cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33175T" | grep "chr10" | grep "43115990" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33175T" | grep "chr10" | grep "43116194" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "33175T" | grep "chr14" | grep "21378612" >> $filename

cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "34854T" | grep "chr8" | grep "108203503" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "34854T" | grep "chr10" | grep "43115256" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "34854T" | grep "chr10" | grep "43115259" >> $filename
cat chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff.txt | grep "34854T" | grep "chr19" | grep "34172761" >> $filename

cat $filename | wc -l



Sample	CHROM	POS
33059T 	 chr1	205625650
33059T 	 chr6	82791982
33059T 	 chr6	82791991
33059T 	 chr9	65437801
33059T 	 chr9	65438087
33059T 	 chr17	7241460
33175T 	 chr10	43115990
33175T 	 chr10	43116194
33175T 	 chr14	21378612
34854T 	 chr8	108203503
34854T 	 chr10	43115256
34854T 	 chr10	43115259
34854T 	 chr19	34172761

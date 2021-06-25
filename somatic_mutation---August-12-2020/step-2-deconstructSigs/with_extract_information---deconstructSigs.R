#####
#####  NOTICE
#####

# must install the github package
# using install.packages()
# the function vcf.to.sigs.input
# has a bug and produce the error as below

# 
#library(devtools)
#install_github("raerose01/deconstructSigs")

### error 1
#> sigs.input = vcf.to.sigs.input(vcf_file)
#Error in .geno2geno(lst) : only diploid variants are supported

### error 2
#> sigs.input = vcf.to.sigs.input(vcf_file)
#Error in data.frame(sample = sample, chr = chr[alt1], pos = pos[alt1],  :
#  arguments imply differing number of rows: 1, 0


#https://chuansongme.com/n/2648933752117
#
#http://www.bio-info-trainee.com/2518.html

library("deconstructSigs")
library("BSgenome.Hsapiens.UCSC.hg38")
vcf.to.sigs.input

# work dir
data_dir_folder = "/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/mutation_signature/test-new-condition-June-30-2020---oncocytomas/result-filter/"

#
save_folder = "/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/mutation_signature/test-new-condition-June-30-2020---oncocytomas/result-deconstructSigs/"

#
#tumor_list = c("GPK0267-2000", "GPK0266-2000", "GPK0266-2001", "GPK0266-2002", "GPK0266-2003", "GPK0266-2004", "GPK0266-2005", "GPK0266-2006", "GPK0266-2007", "GPK0265-2000", "GPK0264-2000", "GPK0264-2001", "GPK0263-2004", "GPK0263-2003", "GPK0263-2002", "GPK0263-2001", "GPK0263-2000", "31252T", "32545T", "34854T", "33175T", "38783T", "33014T", "36266T", "32537T", "34059T", "47123T", "33059T", "36585T") 
#tumor_list = c("GPK0267-2000")
tumor_list = c("GPK0267-2000", "GPK0266-2000", "GPK0266-2001", "GPK0266-2002", "GPK0266-2003", "GPK0266-2004", "GPK0266-2005", "GPK0266-2006", "GPK0266-2007", "GPK0265-2000", "GPK0264-2000", "GPK0264-2001", "GPK0263-2004", "GPK0263-2003", "GPK0263-2002", "GPK0263-2001", "GPK0263-2000") 


filter_condition = c("filter_condition_1", "filter_condition_2", "filter_condition_3")

filter_condition = c("filter_test_1")

for (j in 1:length(filter_condition)) {
    dir.create(paste0(save_folder, "/", filter_condition[j]), showWarnings = FALSE)
	file_path = paste0(save_folder, "/", filter_condition[j], "-summary", ".txt")
    title = c("sample_ID,sample_name,mutation_signature"); write(title,file=file_path,append=TRUE)

}


for (i in 1:length(tumor_list)) {

	sample_name = tumor_list[i]
	
	for (j in 1:length(filter_condition)) {
		data_folder = paste0(data_dir_folder, "/", filter_condition[j], "/", sample_name)
		setwd(data_folder)
		
		#
		vcf_file = paste0("remove_chrUn--AF_normal_less--AF_tumor_greater--normal_depth--between-8-200--tumor_depth--between-10-200--PASS--diploid--", sample_name, ".somatic.filtered.snpEff.PASS.vcf")
		
		# this one
		sigs.input = vcf.to.sigs.input(vcf_file, bsg=BSgenome.Hsapiens.UCSC.hg38)
		sum(sigs.input)		
		#sample_1 = whichSignatures(tumor.ref = sigs.input, signatures.ref = signatures.nature2013, sample.id = "47123T", contexts.needed = TRUE,tri.counts.method = 'default')
		
		# Plot example
		plot_example <- whichSignatures(tumor.ref = sigs.input, 
		#signatures.ref = signatures.nature2013, 
		signatures.ref = signatures.cosmic,
		sample.id = sample_name, 
		contexts.needed = TRUE,
		tri.counts.method = 'default')
		
		# Plot output
		#png(filename=paste0(paste0(save_folder, "/", filter_condition[j]),"/", sample_name, "--", filter_condition[j], ".png"), width = 800, height = 600)
		pdf(paste0(paste0(save_folder, "/", filter_condition[j]),"/", sample_name, "--", filter_condition[j], ".pdf"), width = 8, height = 8)

		plotSignatures(plot_example, sub = paste0("analyzed with ", sum(sigs.input), " diploid mutations"))
		dev.off()
		
		
		###  Extract the signature information
		weights <- plot_example[["weights"]]
		tmp <- which(weights != 0)
		c <- paste(colnames(weights)[tmp[1]], " : ", round(weights[tmp[1]],
		3), sep = "")
		
		if (length(tmp) > 1) {
			for (i in tmp[2:length(tmp)]) {
		        c <- paste(c, " & ", colnames(weights)[i], " : ",
		        round(weights[i], 3), sep = "")
		    }
		}
		c
		
		write(paste0(sample_name, ",", sample_name, ",", c),file=file_path,append=TRUE)
	
	}

}

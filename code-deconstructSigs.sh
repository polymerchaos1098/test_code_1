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
library("BSgenome.Hsapiens.NCBI.GRCh38")
library("BSgenome.Hsapiens.UCSC.hg19")

vcf.to.sigs.input

# work dir
data_dir_folder = "/projectsp/foran/Team_Chan/Hua_files/Human_LFS_VAF-9-5-2019/mutation_signature/result-UAC/data/signature_data/deconstructSigs/"

sample_name = "group_1"

#data_folder = paste0(data_dir_folder, "/", sample_name)
data_folder = data_dir_folder

setwd(data_folder)

#
vcf_file = paste0(sample_name, "-data.vcf")

# this one
sigs.input = vcf.to.sigs.input(vcf_file, bsg=BSgenome.Hsapiens.UCSC.hg19)
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
plotSignatures(plot_example, sub = paste0("analyzed with ", sum(sigs.input), " diploid mutations"))



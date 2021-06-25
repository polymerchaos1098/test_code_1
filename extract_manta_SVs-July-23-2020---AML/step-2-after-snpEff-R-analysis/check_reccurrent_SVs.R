library(plyr)

data_folder="/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/results/extract_manta_SVs-June-11-2020/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method/"

save_folder="/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/results/extract_manta_SVs-June-11-2020/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method/"

# data
somatic_data = read.csv(paste0(data_folder, "SOMATICSCORE-greater-44---extracted-information---Thyroid---ALL_PASS_somaticSV_candidates", ".txt"), header=TRUE, sep=" ")
dim(somatic_data)


colnames(somatic_data)

somatic_data = data.frame(id=paste(trimws(as.character(somatic_data$CHROM)),"-",  trimws(as.character(somatic_data$POS)), sep=""), somatic_data)
somatic_data = somatic_data[, -10]
head(somatic_data); dim(somatic_data); 

a=count(somatic_data$id)
b = a[order(a$freq, decreasing=TRUE),]
c=b[b$freq>1,]


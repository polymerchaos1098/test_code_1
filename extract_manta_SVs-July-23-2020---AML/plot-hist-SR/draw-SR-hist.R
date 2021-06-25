#
library(stringr)

data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/"

save_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/"

# data
somatic_data = read.csv(paste0(data_folder, "ALL_PASS_somaticSV_candidates", ".txt"), header=FALSE, sep=" ")
somatic_data = data.frame(somatic_data, REF_SR=-1, ALT_SR=-1, Ratio_ALT_divided_REF=-1, SOMATICSCORE=-1)
head(somatic_data); dim(somatic_data); 

ref_list = c()
alt_list = c()
ratio_list =c()
score_list =c()

#
i = 1
i = 3
i = 439

for (i in 1:nrow(somatic_data)) {
    tmp = as.character(somatic_data$V12[i])
	tmp = unlist(strsplit(tmp, ":"))[2]

	ref = as.numeric(unlist(strsplit(tmp, ","))[1])
	alt = as.numeric(unlist(strsplit(tmp, ","))[2])
	ref; alt
	
	ref_list = c(ref_list, ref)
	alt_list = c(alt_list, alt)
	ratio_list = c(ratio_list, round(alt/ref,3))
	
	somatic_data$REF_SR[i] = ref
	somatic_data$ALT_SR[i] = alt
	somatic_data$Ratio_ALT_divided_REF[i] = round(alt/ref,3)
	
	# SOMATICSCORE
	tmp = as.character(somatic_data$V9[i])
	tmp = substr(tmp, str_locate(tmp, ";SOMATICSCORE=")[2] + 1, nchar(tmp)) 
	tmp
	
	# check if it is in the end
	pos = str_locate(tmp, ";")
	
	if (is.na(pos[1])) {
	    somatic_data$SOMATICSCORE[i] = tmp
	} else {
		somatic_data$SOMATICSCORE[i] = substr(tmp, 1, pos[1] - 1) 
	}
	score_list = c(score_list, somatic_data$SOMATICSCORE[i])
	
}

write.csv(somatic_data, paste0(data_folder, "AML---ALL_PASS_somaticSV_candidates", ".csv"))

#
length(ratio_list)
ratio_list = na.omit(ratio_list)
length(ratio_list)


#https://blog.csdn.net/weixin_40628687/article/details/79254791

# ref
png(filename=paste0(data_folder, "SR-REF-READ---AML-manta", ".png"), width=1000, height=600)
par(mar=c(5,5,5,5))
hist(ref_list, 
	 main=paste0("REF READ from: SR, ", length(ref_list), " candidates"), col.main="red", cex.main=2.2,
     ylab="Frequency", cex.lab=1.9,cex.axis=1.5,
     xlab="REF READ from SR", 
     border="blue", 
     col="green"
	)
dev.off()



# alt
png(filename=paste0(data_folder, "SR-ALT-READ---AML-manta", ".png"), width=1000, height=600)
par(mar=c(5,5,5,5))
hist(alt_list, 
	 main=paste0("ALT READ from: SR, ", length(alt_list), " candidates"), col.main="red", cex.main=2.2,
     ylab="Frequency", cex.lab=1.9,cex.axis=1.5,
     xlab="ALT READ from SR", 
     border="blue", 
     col="blue",
	 xlim=c(0,100),
	 breaks=seq(0,225,5)
	)
dev.off()


# ratio
png(filename=paste0(data_folder, "ALT-REF-ratio---AML-manta", ".png"), width=1000, height=600)
par(mar=c(5,5,5,5))
hist(ratio_list[ratio_list<1.0], 
	 main=paste0("ALT/REF ratio from: SR, ", length(ratio_list[ratio_list<1.0]), " candidates, some >= 1 not shown here"), col.main="red", cex.main=2.2,
     ylab="Frequency", cex.lab=1.9,cex.axis=1.5,
     xlab="ALT/REF ratio from SR", 
     border="blue", 
     col="red",
	 xlim=c(0,1),
	)
dev.off()


# score
score_list = as.numeric(na.omit(score_list)); min(score_list); max(score_list)
length(score_list); score_list = score_list[score_list<=120]
png(filename=paste0(data_folder, "SOMATICSCORE---AML-manta", ".png"), width=800, height=600)
par(mar=c(5,5,2,5))
hist(score_list, 
     main="",
	 #main=paste0("SOMATICSCORE, ", length(score_list), " candidates, interval = 2"), col.main="red", cex.main=2.2,
     ylab="Frequency", cex.lab=1.9,cex.axis=1.5,
     xlab="SOMATICSCORE", 
     border="blue", 
     col="hotpink1",
	 xlim=c(0,120),
	 breaks=seq(0,120,2),
	 xaxt ="n",
	 yaxt ="n"
	)
axis(2,seq(0, 1200, 100), seq(0, 1200, 100), las = 1, font=2.4)
axis(1,seq(0, 120, 5), seq(0, 120, 5), las = 1, font=2.4)
box()
dev.off()












hist(ratio_list)
hist(ratio_list[ratio_list<1.0])
min(ratio_list)

########################################
# tumor
tumor_list=c("GPK0267-2000", "GPK0266-2000", "GPK0266-2001", "GPK0266-2002", "GPK0266-2003", "GPK0266-2004", "GPK0266-2005", "GPK0266-2006", "GPK0266-2007", "GPK0265-2000", "GPK0264-2000", "GPK0264-2001", "GPK0263-2004", "GPK0263-2003", "GPK0263-2002", "GPK0263-2001", "GPK0263-2000", "31252T", "32545T", "34854T", "33175T", "38783T", "33014T", "36266T", "32537T", "34059T", "47123T", "33059T", "36585T")


########################################
# get Sequenza purity estimation value:
sequenza_purity = read.csv("/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/get_VAF/code/best_fitting_tumor_purity-with_correct_X_Y.csv", header = TRUE)
head(sequenza_purity); dim(sequenza_purity)
	dev.off()


#
i = 1

for (i in 1:length(tumor_list)) {
	print(i)
	tmp = read.csv(paste0(data_folder, tumor_list[i], "_VAF.txt"), header=FALSE, sep=" ")
	head(tmp); dim(tmp)
	
	# 1 row is the names
	idx = which(tumor_list[i] == tmp[1,])
	
	#
	target_data = tmp[-1,]
	head(target_data); dim(target_data); head(as.numeric(as.character(trimws(target_data[,idx],"both"))))

	a = as.numeric(as.character(trimws(target_data[,idx],"both")))
	
	# get sequenza purity
	p=0
	p = which(tumor_list[i]==sequenza_purity$ID); p
	sequenza_purity$cellularity[p]
	
	###
	#dev.new(width=10, height=6)

	png(filename=paste0("/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/get_VAF/result/figure/", tumor_list[i], "_", "column_is_", idx, ".png"), width=1000, height=600)

	par(mar=c(5,5,5,5))
	
	hist(a, 
	 main=paste0("Histogram of VAF from: ", tumor_list[i], "\n", sequenza_purity$ID[p], " tumor purity (Sequenza) is: ", sequenza_purity$cellularity[p]), col.main="red", cex.main=2.2,
     ylab="Frequency", cex.lab=1.9,cex.axis=1.5,
     xlab="Variant Allele Frequency (VAF)", 
     border="blue", 
     col="green",
     xlim=c(0,1),
     ylim=c(0,300),
     las=1, 
	 breaks=seq(0,1,0.005)
	)
	
	dev.off()
 
}




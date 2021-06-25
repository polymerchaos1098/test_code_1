#
data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/get_VAF---August-12-2020/result/data/"

save_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/get_VAF---August-12-2020/result/figure/"
dir.create(save_folder, showWarnings = FALSE)

########################################
# tumor
tumor_list=c("121248", "133128", "133566", "133567", "134807") 

#
i = 1

for (i in 1:length(tumor_list)) {
	print(i)
	tmp = read.csv(paste0(data_folder, tumor_list[i], "_VAF.txt"), header=TRUE, sep=" ")
	head(tmp); dim(tmp)
	
	# 1 row is the names
	#idx = which(tumor_list[i] == tmp[1,])
	
	#
	#target_data = tmp[-1,]
	#head(target_data); dim(target_data); head(as.numeric(as.character(trimws(target_data[,idx],"both"))))

	#a = as.numeric(as.character(trimws(target_data[,idx],"both")))

	# col 4 is the tumor
	idx = 4
	a =  as.numeric(as.character(trimws(tmp[,idx],"both")))
	
	###
	#dev.new(width=10, height=6)

	png(filename=paste0(save_folder, tumor_list[i], "_", "column_is_", idx, ".png"), width=1000, height=600)

	par(mar=c(5,5,5,5))
	
	hist(a, 
	 main=paste0("Sample-", tumor_list[i], ": 1). used ", length(a), " mutations marked as \"PASS\" & \n", "2). only used  biallelic mutations (removed few multiallelic mut.)"), col.main="blue", cex.main=2.2,
     ylab="Frequency", cex.lab=1.9,cex.axis=1.5,
     xlab="Variant Allele Frequency (VAF), biallelic mutations", 
     border="blue", 
     col="green",
     xlim=c(0,1),
     ylim=c(0,75),
     las=1, 
	 breaks=seq(0,1,0.005)
	)
	
	dev.off()
 
}




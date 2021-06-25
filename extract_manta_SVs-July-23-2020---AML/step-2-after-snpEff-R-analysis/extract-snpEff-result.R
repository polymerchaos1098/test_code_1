#
### COSMIC gene: cosmic-grch37-v84-cancer_gene_census.csv
### too long, separate them into 2 list and concatenate them
COSMIC_gene_list_1 = c("A1CF", "ABI1", "ABL1", "ABL2", "ACKR3", "ACSL3", "ACSL6", "ACVR1", "ACVR2A", "AFF1", "AFF3", "AFF4", "AKAP9", "AKT1", "AKT2", "AKT3", "ALDH2", "ALK", "AMER1", "ANK1", "APC", "APOBEC3B", "AR", "ARAF", "ARHGAP26", "ARHGAP5", "ARHGEF10", "ARHGEF10L", "ARHGEF12", "ARID1A", "ARID1B", "ARID2", "ARNT", "ASPSCR1", "ASXL1", "ASXL2", "ATF1", "ATIC", "ATM", "ATP1A1", "ATP2B3", "ATR", "ATRX", "AXIN1", "AXIN2", "B2M", "BAP1", "BARD1", "BAX", "BAZ1A", "BCL10", "BCL11A", "BCL11B", "BCL2", "BCL2L12", "BCL3", "BCL6", "BCL7A", "BCL9", "BCL9L", "BCLAF1", "BCOR", "BCORL1", "BCR", "BIRC3", "BIRC6", "BLM", "BMP5", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRD3", "BRD4", "BRIP1", "BTG1", "BTK", "BUB1B", "C15ORF65", "C2ORF44", "CACNA1D", "CALR", "CAMTA1", "CANT1", "CARD11", "CARS", "CASC5", "CASP3", "CASP8", "CASP9", "CBFA2T3", "CBFB", "CBL", "CBLB", "CBLC", "CCDC6", "CCNB1IP1", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCR4", "CCR7", "CD209", "CD274", "CD28", "CD74", "CD79A", "CD79B", "CDC73", "CDH1", "CDH10", "CDH11", "CDH17", "CDK12", "CDK4", "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2C", "CDX2", "CEBPA", "CEP89", "CHCHD7", "CHD2", "CHD4", "CHEK2", "CHIC2", "CHST11", "CIC", "CIITA", "CLIP1", "CLP1", "CLTC", "CLTCL1", "CNBD1", "CNBP", "CNOT3", "CNTNAP2", "CNTRL", "COL1A1", "COL2A1", "COL3A1", "COX6C", "CPEB3", "CREB1", "CREB3L1", "CREB3L2", "CREBBP", "CRLF2", "CRNKL1", "CRTC1", "CRTC3", "CSF1R", "CSF3R", "CSMD3", "CTCF", "CTNNA2", "CTNNB1", "CTNND1", "CTNND2", "CUL3", "CUX1", "CXCR4", "CYLD", "CYP2C8", "CYSLTR2", "DAXX", "DCAF12L2", "DCC", "DCTN1", "DDB2", "DDIT3", "DDR2", "DDX10", "DDX3X", "DDX5", "DDX6", "DEK", "DGCR8", "DICER1", "DNAJB1", "DNM2", "DNMT3A", "DROSHA", "DUX4L1", "EBF1", "ECT2L", "EED", "EGFR", "EIF1AX", "EIF3E", "EIF4A2", "ELF3", "ELF4", "ELK4", "ELL", "ELN", "EML4", "EP300", "EPAS1", "EPHA3", "EPHA7", "EPS15", "ERBB2", "ERBB3", "ERBB4", "ERC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERG", "ESR1", "ETNK1", "ETV1", "ETV4", "ETV5", "ETV6", "EWSR1", "EXT1", "EXT2", "EZH2", "EZR", "FAM131B", "FAM135B", "FAM46C", "FAM47C", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FAS", "FAT1", "FAT3", "FAT4", "FBLN2", "FBXO11", "FBXW7", "FCGR2B", "FCRL4", "FEN1", "FES", "FEV", "FGFR1", "FGFR1OP", "FGFR2", "FGFR3", "FGFR4", "FH", "FHIT", "FIP1L1", "FKBP9", "FLCN", "FLI1", "FLNA", "FLT3", "FLT4", "FNBP1", "FOXA1", "FOXL2", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FOXR1", "FSTL3", "FUBP1", "FUS", "GAS7", "GATA1", "GATA2", "GATA3", "GLI1", "GMPS", "GNA11", "GNAQ", "GNAS", "GOLGA5", "GOPC", "GPC3", "GPC5", "GPHN", "GRIN2A", "GRM3", "H3F3A", "H3F3B", "HERPUD1", "HEY1", "HIF1A", "HIP1", "HIST1H3B", "HIST1H4I", "HLA-A", "HLF", "HMGA1", "HMGA2", "HMGN2P46", "HNF1A", "HNRNPA2B1", "HOOK3", "HOXA11", "HOXA13", "HOXA9", "HOXC11", "HOXC13", "HOXD11", "HOXD13", "HRAS", "HSP90AA1", "HSP90AB1", "ID3", "IDH1", "IDH2", "IGF2BP2", "IGH", "IGK", "IGL", "IKBKB", "IKZF1", "IL2", "IL21R", "IL6ST", "IL7R", "IRF4", "IRS4", "ISX", "ITGAV", "ITK", "JAK1", "JAK2", "JAK3", "JAZF1", "JUN", "KAT6A", "KAT6B", "KAT7", "KCNJ5", "KDM5A", "KDM5C", "KDM6A", "KDR", "KDSR", "KEAP1", "KIAA1549", "KIAA1598", "KIF5B", "KIT", "KLF4", "KLF6", "KLK2", "KMT2A", "KMT2C", "KMT2D", "KNSTRN", "KRAS", "KTN1", "LARP4B", "LASP1", "LCK", "LCP1", "LEF1", "LEPROTL1", "LHFP", "LIFR", "LMNA", "LMO1", "LMO2", "LPP", "LRIG3", "LRP1B", "LSM14A", "LYL1", "LZTR1", "MAF", "MAFB", "MALAT1", "MALT1", "MAML2", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAPK1", "MAX", "MB21D2", "MDM2", "MDM4", "MDS2", "MECOM", "MED12")

COSMIC_gene_list_2 = c("MEN1", "MET", "MGMT", "MITF", "MKL1", "MLF1", "MLH1", "MLLT1", "MLLT10", "MLLT11", "MLLT3", "MLLT4", "MLLT6", "MN1", "MNX1", "MPL", "MSH2", "MSH6", "MSI2", "MSN", "MTCP1", "MTOR", "MUC1", "MUC16", "MUC4", "MUTYH", "MYB", "MYC", "MYCL", "MYCN", "MYD88", "MYH11", "MYH9", "MYO5A", "MYOD1", "N4BP2", "NAB2", "NACA", "NBEA", "NBN", "NCKIPSD", "NCOA1", "NCOA2", "NCOA4", "NCOR1", "NCOR2", "NDRG1", "NF1", "NF2", "NFATC2", "NFE2L2", "NFIB", "NFKB2", "NFKBIE", "NIN", "NKX2-1", "NONO", "NOTCH1", "NOTCH2", "NPM1", "NR4A3", "NRAS", "NRG1", "NSD1", "NT5C2", "NTHL1", "NTRK1", "NTRK3", "NUMA1", "NUP214", "NUP98", "NUTM1", "NUTM2A", "NUTM2B", "OLIG2", "OMD", "P2RY8", "PABPC1", "PAFAH1B2", "PALB2", "PAX3", "PAX5", "PAX7", "PAX8", "PBRM1", "PBX1", "PCBP1", "PCM1", "PDCD1LG2", "PDE4DIP", "PDGFB", "PDGFRA", "PDGFRB", "PER1", "PHF6", "PHOX2B", "PICALM", "PIK3CA", "PIK3CB", "PIK3R1", "PIM1", "PLAG1", "PLCG1", "PML", "PMS1", "PMS2", "POLD1", "POLE", "POLG", "POLQ", "POT1", "POU2AF1", "POU5F1", "PPARG", "PPFIBP1", "PPM1D", "PPP2R1A", "PPP6C", "PRCC", "PRDM1", "PRDM16", "PRDM2", "PREX2", "PRF1", "PRKACA", "PRKAR1A", "PRKCB", "PRPF40B", "PRRX1", "PSIP1", "PTCH1", "PTEN", "PTK6", "PTPN11", "PTPN13", "PTPN6", "PTPRB", "PTPRC", "PTPRD", "PTPRK", "PTPRT", "PWWP2A", "QKI", "RABEP1", "RAC1", "RAD17", "RAD21", "RAD51B", "RAF1", "RALGDS", "RANBP2", "RAP1GDS1", "RARA", "RB1", "RBM10", "RBM15", "RECQL4", "REL", "RET", "RFWD3", "RGPD3", "RGS7", "RHOA", "RHOH", "RMI2", "RNF213", "RNF43", "ROBO2", "ROS1", "RPL10", "RPL22", "RPL5", "RPN1", "RSPO2", "RSPO3", "RUNDC2A", "RUNX1", "RUNX1T1", "S100A7", "SALL4", "SBDS", "SDC4", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SEPT5", "SEPT6", "SEPT9", "SET", "SETBP1", "SETD1B", "SETD2", "SF3B1", "SFPQ", "SFRP4", "SGK1", "SH2B3", "SH3GL1", "SIRPA", "SIX1", "SIX2", "SKI", "SLC34A2", "SLC45A3", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1", "SMARCD1", "SMARCE1", "SMC1A", "SMO", "SND1", "SOCS1", "SOX2", "SOX21", "SPECC1", "SPEN", "SPOP", "SRC", "SRGAP3", "SRSF2", "SRSF3", "SS18", "SS18L1", "SSX1", "SSX2", "SSX4", "STAG1", "STAG2", "STAT3", "STAT5B", "STAT6", "STIL", "STK11", "STRN", "SUFU", "SUZ12", "SYK", "TAF15", "TAL1", "TAL2", "TBL1XR1", "TBX3", "TCEA1", "TCF12", "TCF3", "TCF7L2", "TCL1A", "TEC", "TERT", "TET1", "TET2", "TFE3", "TFEB", "TFG", "TFPT", "TFRC", "TGFBR2", "THRAP3", "TLX1", "TLX3", "TMEM127", "TMPRSS2", "TNC", "TNFAIP3", "TNFRSF14", "TNFRSF17", "TOP1", "TP53", "TP63", "TPM3", "TPM4", "TPR", "TRA", "TRAF7", "TRB", "TRD", "TRIM24", "TRIM27", "TRIM33", "TRIP11", "TRRAP", "TSC1", "TSC2", "TSHR", "U2AF1", "UBR5", "USP44", "USP6", "USP8", "VAV1", "VHL", "VTI1A", "WAS", "WHSC1", "WHSC1L1", "WIF1", "WNK2", "WRN", "WT1", "WWTR1", "XPA", "XPC", "XPO1", "YWHAE", "ZBTB16", "ZCCHC8", "ZEB1", "ZFHX3", "ZMYM3", "ZNF198", "ZNF278", "ZNF331", "ZNF384", "ZNF429", "ZNF479", "ZNF521", "ZNRF3", "ZRSR2")

COSMIC_gene_list = c(COSMIC_gene_list_1, COSMIC_gene_list_2)
length(COSMIC_gene_list); head(COSMIC_gene_list); tail(COSMIC_gene_list)

######
######
######
library(stringr)
library(plyr)

data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method/"

save_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method/"

# data
somatic_data = read.csv(paste0(data_folder, "chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff", ".txt"), header=TRUE, sep="\t")
head(somatic_data); dim(somatic_data); 

somatic_data = data.frame(somatic_data, REF_SR=-1, ALT_SR=-1, Ratio_ALT_divided_REF=-1, SOMATICSCORE=-1, SVTYPE="-", Involved_genes="-")
somatic_data$SVTYPE = as.character(somatic_data$SVTYPE)
somatic_data$Involved_genes = as.character(somatic_data$Involved_genes)

head(somatic_data); dim(somatic_data); 

ref_list = c()
alt_list = c()
ratio_list =c()
score_list =c()

#
i = 1
i = 3
i = 439
i = 529

for (i in 1:nrow(somatic_data)) {
	print(i)
    tmp = as.character(somatic_data$TUMOR[i])
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
	tmp = as.character(somatic_data$INFO[i])
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
	
	
	# ;SVTYPE=
	tmp = as.character(somatic_data$INFO[i])
	tmp = substr(tmp, str_locate(tmp, "SVTYPE=")[2] + 1, nchar(tmp)) 
	if (!is.na(tmp)) {
		tmp
	
		# check if it is in the end
		pos = str_locate(tmp, ";")
	
		if (is.na(pos[1])) {
			somatic_data$SVTYPE[i] = tmp
		} else {
			somatic_data$SVTYPE[i] = substr(tmp, 1, pos[1] - 1) 
		}
	}
	
	# ;ANN=
	tmp = as.character(somatic_data$INFO[i])
	tmp = substr(tmp, str_locate(tmp, ";ANN=")[2] + 1, nchar(tmp)) 
	
	###########
	# few candidates have "LOF="
	# remove them, cause it has different format than "ANN"
	
	lof_pos = str_locate(tmp, "LOF=")
	
	if (!is.na(lof_pos[1])) {
		tmp = substr(tmp, 1, lof_pos[1]-1) 
	}
	###########
	
	if (!is.na(tmp)) {
	   tmp
	
	   all_ANN_terms = unlist(strsplit(tmp, ","))
	
		# extract "gene"
		gene = c()
		j = 20
		for (j in 1:length(all_ANN_terms)) {
			one_ANN = all_ANN_terms[j]
			# 4th is gene
			a = unlist(strsplit(one_ANN, "\\|"))[4]
			if (a!="") gene = c(gene, a)
		}
		gene = paste(unique(gene), collapse = '; ')
	
		# check if it is in the end
		pos = str_locate(tmp, ";")
	
		if (!is.null(gene)) somatic_data$Involved_genes[i] = str_replace_all(gene, "&", "; ")
	}
	
}


write.table(somatic_data, paste0(save_folder, "ALL---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

# filter out SOMATICSCORE>38
dim(somatic_data)
res = subset(somatic_data, SOMATICSCORE>38)
dim(res)

write.table(res, paste0(save_folder, "SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

# ANN is too huge, not save
res_1 = res[, -9]
write.table(res_1, paste0(save_folder, "no-ANN---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

write.csv(res_1, paste0(save_folder, "no-ANN---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)


# get each sample hist for type
sample_list = as.character(unique(res_1$Sample))
sample_list

need = matrix("", nrow = length(sample_list)+1, ncol = 5) 
i = 1
table(res_1$SVTYPE)

# 4 categories: BND DEL DUP INV
categories = c("BND", "DEL", "DUP" ,"INV", "INS")

for (i in 1:length(sample_list)) {
	
	for (j in 1:length(categories)) {
		tmp = subset(res_1, as.character(Sample) == sample_list[i] & as.character(SVTYPE)==categories[j])
		need[i, j] = nrow(tmp)
    }
	
}

for (j in 1:length(categories)) {
	tmp = subset(res_1, as.character(SVTYPE)==categories[j])
	need[length(sample_list)+1, j] = nrow(tmp)
}

colnames(need) = categories
rownames(need) = c(sample_list, "ALL")

need 

write.csv(need, paste0(save_folder, "SV-categories---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = TRUE)




# check reccurrent SVs
dim(somatic_data)
res = subset(somatic_data, SOMATICSCORE>38)
dim(res)

# test
#test = somatic_data
#test = data.frame(uid=paste(trimws(as.character(test$CHROM)),"-",  trimws(as.character(test$POS)), sep=""), test)
#test1 = count(test$uid)
#test1


#
colnames(res)

res = data.frame(uid=paste(trimws(as.character(res$CHROM)),"-",  trimws(as.character(res$POS)), sep=""), res)
#res = res[, -10]
head(res); dim(res); 

a=count(res$uid)
b = a[order(a$freq, decreasing=TRUE),]
c=b[b$freq>1,]
c


# reccurent subset
x = subset(res, is.element(uid, c$x))
dim(x)

#
write.csv(x, paste0(save_folder, "reccurrent---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)

####################
####################
####################
# get COSMIC genes


# some mutations have several involved genes, here, separate them into individidual rows
gene_statistics = c()

i = 1
for (i in 1:nrow(x)) {

    if (i%%10==1) print(i)
	
    tmp = as.character(x$Involved_genes[i])
	all_Involved_genes = trimws(unlist(strsplit(tmp, ";")))
	all_Involved_genes
	for (j in 1:length(all_Involved_genes)) {
		one_row = data.frame(sample_gene=paste0(trimws(as.character(x$Sample[i])), "-", all_Involved_genes[j]), gene=all_Involved_genes[j], x[i, ])
		gene_statistics = rbind(gene_statistics, one_row)
	}
	
}

head(gene_statistics); dim(gene_statistics)

#
cosmic_x = subset(gene_statistics, is.element(toupper(gene), toupper(COSMIC_gene_list)))
dim(cosmic_x)

write.csv(cosmic_x, paste0(save_folder, "COSMIC---reccurrent---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)


### comparison inside a group
### used two table above: need & x
group_list = c("GPK0263", "GPK0264", "GPK0265", "GPK0266", "GPK0267")
group_number = c(5, 2, 1, 8, 1)
sample_list = c("GPK0263-2000", "GPK0263-2001", "GPK0263-2002", "GPK0263-2003", "GPK0263-2004", "GPK0264-2000", "GPK0264-2001", "GPK0265-2000", "GPK0266-2000", "GPK0266-2001", "GPK0266-2002", "GPK0266-2003", "GPK0266-2004", "GPK0266-2005", "GPK0266-2006", "GPK0266-2007", "GPK0267-2000")

i = 3
for (i in 1:length(group_list)) {

	target_table = matrix(-1, nrow = group_number[i], ncol = group_number[i])
	target_sample_set = sample_list[grep(group_list[i], sample_list)]
	colnames(target_table) = target_sample_set
	rownames(target_table) = target_sample_set
	
	# add the total number in diagonal line
	# check the need table
	j = 1
	for (j in 1:nrow(target_table)) {
		idx = which(target_sample_set[j]==trimws(rownames(need)))
		target_table[j, j] = sum(as.numeric(need[idx,]))		
	}
	
	# check pair recurrent variant
	j = 1
	if (nrow(target_table) > 1) {
	for (j in 1:(nrow(target_table)-1)) {
        s1 = subset(x, trimws(as.character(Sample))==rownames(target_table)[j])	
		dim(s1)
		for (k in (j+1):nrow(target_table)) {
           t1 = subset(x, trimws(as.character(Sample))==rownames(target_table)[k])	
		   dim(t1)
		   a = intersect(s1$uid, t1$uid)
		   target_table[j, k] = length(a)
		   # get information		   
		}
	}
	} # if
	
	target_table[target_table==-1] = "-"
	write.csv(target_table, paste0(save_folder, group_list[i], "---inside_group_statistic---SOMATICSCORE-greater-38---AML", ".csv"), row.names = TRUE)
	
}


####################
####################
####################
# check the COSMIC gene
head(res);dim(res)


# some mutations have several involved genes, here, separate them into individidual rows
gene_statistics = c()

i = 1
for (i in 1:nrow(res)) {

    if (i%%10==1) print(i)
	
    tmp = as.character(res$Involved_genes[i])
	all_Involved_genes = trimws(unlist(strsplit(tmp, ";")))
	all_Involved_genes
	for (j in 1:length(all_Involved_genes)) {
		one_row = data.frame(sample_gene=paste0(trimws(as.character(res$Sample[i])), "-", all_Involved_genes[j]), gene=all_Involved_genes[j], res[i, ])
		gene_statistics = rbind(gene_statistics, one_row)
	}
	
}

head(gene_statistics); dim(gene_statistics)

cosmic_candidates = subset(gene_statistics, is.element(toupper(gene), toupper(COSMIC_gene_list)))
dim(cosmic_candidates); dim(gene_statistics);

# uid, some mutations include many COSMIC genes
res_cosmic = cosmic_candidates[!duplicated(cosmic_candidates[ , c("uid")]),]
dim(res_cosmic);
res_cosmic = data.frame(all_COSMIC_GENES="", res_cosmic)
res_cosmic$all_COSMIC_GENES = as.character(res_cosmic$all_COSMIC_GENES)
### collect all COSMIC gene in a column
i = 1
for (i in 1:nrow(res_cosmic)) {
	tmp = subset(cosmic_candidates, trimws(as.character(uid))==trimws(as.character(res_cosmic$uid[i])))
	res_cosmic$all_COSMIC_GENES[i] = paste(trimws(as.character(tmp$gene)), collapse="; ")
	print(dim(tmp))
}

head(res_cosmic);dim(res_cosmic);

#write.csv(res_cosmic, paste0(save_folder, "COSMIC---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)
write.table(res_cosmic, paste0(save_folder, "COSMIC---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".txt"), row.names = FALSE, quote=FALSE, sep = " ")
write.table(res_cosmic[, -13], paste0(save_folder, "no-ANN---COSMIC---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".txt"), row.names = FALSE, quote=FALSE, sep = " ")

write.csv(res_cosmic[, -13], paste0(save_folder, "no-ANN---COSMIC---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)

# for COSMIC
need = matrix("", nrow = length(sample_list)+1, ncol = 5) 

# 5 categories: BND DEL DUP INV INS
categories = c("BND", "DEL", "DUP" ,"INV", "INS")

i = 2; j = 4

for (i in 1:length(sample_list)) {
	
	for (j in 1:length(categories)) {
		tmp = subset(res_cosmic, as.character(trimws(Sample)) == sample_list[i] & as.character(trimws(SVTYPE))==categories[j])
		dim(tmp)
		need[i, j] = nrow(tmp)
    }
	
}

for (j in 1:length(categories)) {
	tmp = subset(res_cosmic, as.character(trimws(SVTYPE))==categories[j])
	need[length(sample_list)+1, j] = nrow(tmp)
}

colnames(need) = categories
rownames(need) = c(sample_list, "ALL")

need 

write.csv(need, paste0(save_folder, "COSMIC---SV-categories---SOMATICSCORE-greater-38---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = TRUE)

 











#
length(ratio_list)
ratio_list = na.omit(ratio_list)
length(ratio_list)


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
score_list = as.numeric(score_list); min(score_list); max(score_list)
png(filename=paste0(data_folder, "SOMATICSCORE---AML-manta", ".png"), width=1000, height=600)
par(mar=c(5,5,5,5))
hist(score_list, 
	 main=paste0("SOMATICSCORE, ", length(score_list), " candidates"), col.main="red", cex.main=2.2,
     ylab="Frequency", cex.lab=1.9,cex.axis=1.5,
     xlab="SOMATICSCORE", 
     border="blue", 
     col="hotpink1",
	 xlim=c(0,200),
	 breaks=seq(0,200,2)
	)
dev.off()




hist(ratio_list)
hist(ratio_list[ratio_list<1.0])
min(ratio_list)

########################################
# tumor
tumor_list=c("GPK0267-2000", "GPK0266-2000", "GPK0266-2001", "GPK0266-2002", "GPK0266-2003", "GPK0266-2004", "GPK0266-2005", "GPK0266-2006", "GPK0266-2007", "GPK0265-2000", "GPK0264-2000", "GPK0264-2001", "GPK0263-2004", "GPK0263-2003", "GPK0263-2002", "GPK0263-2001", "GPK0263-2000", "31252T", "32545T", "34854T", "33175T", "38783T", "33014T", "36266T", "32537T", "34059T", "47123T", "33059T", "36585T")


########################################
# get Sequenza purity estimation value:
sequenza_purity = read.csv("/projectsp/foran/Team_Chan/AML_SV_Hua_Job/AML/codes/get_VAF/code/best_fitting_tumor_purity-with_correct_X_Y.csv", header = TRUE)
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

	png(filename=paste0("/projectsp/foran/Team_Chan/AML_SV_Hua_Job/AML/codes/get_VAF/result/figure/", tumor_list[i], "_", "column_is_", idx, ".png"), width=1000, height=600)

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




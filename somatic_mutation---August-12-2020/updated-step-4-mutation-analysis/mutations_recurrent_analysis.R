#
library(stringr)
library(plyr)

data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/somatic_mutation---August-12-2020/result---step-4-mutation-analysis/"

save_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/somatic_mutation---August-12-2020/result---step-4-mutation-analysis/"

# data
somatic_data = read.csv(paste0(data_folder, "AML_all_filtered_out_vcf_in_a_file", ".txt"), header=TRUE, sep="\t")
head(somatic_data); dim(somatic_data); 

somatic_data = data.frame(somatic_data, Involved_genes="-", HGVS_c="-", HGVS_p="-", AF=0, Annot="-", Annot_Impact="-")
somatic_data$Involved_genes = as.character(somatic_data$Involved_genes)
somatic_data$HGVS_c = as.character(somatic_data$HGVS_c)
somatic_data$HGVS_p = as.character(somatic_data$HGVS_p)
somatic_data$AF = as.character(somatic_data$AF)
somatic_data$Annot = as.character(somatic_data$Annot)
somatic_data$Annot_Impact = as.character(somatic_data$Annot_Impact)

head(somatic_data); dim(somatic_data); 

#
i = 1
i = 3
i = 439
i = 529

for (i in 1:nrow(somatic_data)) {
	print(i)

	# 1. Involved_genes, hgvs_c, hgvs_p
	# ;ANN=
	tmp = as.character(somatic_data$INFO[i])
	tmp = substr(tmp, str_locate(tmp, ";ANN=")[2] + 1, nchar(tmp)) 
	
	if (!is.na(tmp)) {
	   tmp
	
	   all_ANN_terms = unlist(strsplit(tmp, ","))
	
		# extract "gene"
		gene = c()
		hgvs_c = c()
		hgvs_p = c()
		Annot = c()
		Annot_Impact = c()

		j = 9
		for (j in 1:length(all_ANN_terms)) {
			one_ANN = all_ANN_terms[j]
			
			# 2nd is Annotation
			a = unlist(strsplit(one_ANN, "\\|"))[2]
			if (!is.na(a)) if (a!="") Annot = c(Annot, a)

			# 3rd is Annotation_Impact
			a = unlist(strsplit(one_ANN, "\\|"))[3]
			if (!is.na(a)) if (a!="") Annot_Impact = c(Annot_Impact, a)

			# 4th is gene
			a = unlist(strsplit(one_ANN, "\\|"))[4]
			if (!is.na(a)) if (a!="") gene = c(gene, a)
			
			# 10th is hgvs_c
			a = unlist(strsplit(one_ANN, "\\|"))[10]
			if (!is.na(a)) if (a!="") hgvs_c = c(hgvs_c, a)
			
			# 11th is hgvs_p
			a = unlist(strsplit(one_ANN, "\\|"))[11]
			if (!is.na(a)) if (a!="") hgvs_p = c(hgvs_p, a)
		}

		#
		if (!is.null(Annot)) {
			Annot = paste(unique(Annot), collapse = '; ')
		    somatic_data$Annot[i] = Annot			
		}
		
		if (!is.null(Annot_Impact)) {
			Annot_Impact = paste(unique(Annot_Impact), collapse = '; ')
		    somatic_data$Annot_Impact[i] = Annot_Impact			
		}
			
		if (!is.null(gene)) {
			gene = paste(unique(gene), collapse = '; ')
		    somatic_data$Involved_genes[i] = gene			
		}
		
		if (!is.null(hgvs_c)) {
			hgvs_c = paste(unique(hgvs_c), collapse = '; ')
		    somatic_data$HGVS_c[i] = hgvs_c			
		}
		
		if (!is.null(hgvs_p)) {
			hgvs_p = paste(unique(hgvs_p), collapse = '; ')
		    somatic_data$HGVS_p[i] = hgvs_p			
		}
		
	}
	
	# 2. AF tumor
	tmp = as.character(somatic_data$FORMAT[i])
	all_FORMAT_terms = unlist(strsplit(tmp, ":"))
	pos = which("AF"==all_FORMAT_terms)
	
	tmp = as.character(somatic_data$TUMOR[i])
	if (length(pos)!=0) somatic_data$AF[i] = unlist(strsplit(tmp, ":"))[pos]
	
}

head(somatic_data); dim(somatic_data)

somatic_data = data.frame(uid=paste(trimws(as.character(somatic_data$CHROM)),"-",  trimws(as.character(somatic_data$POS)), sep=""), somatic_data)
head(somatic_data); dim(somatic_data); 

write.table(somatic_data, paste0(save_folder, "AML---ALL_filtered_out_mutations", ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

write.csv(somatic_data, paste0(save_folder, "AML---ALL_filtered_out_mutations", ".csv"), row.names = FALSE)

###############
# distribution of mutations in chr 1-22, X, Y, M
# no "chr" in the beginning
chr_list = paste("", c(c(1:22), "X", "Y", "M"), sep="")
chr_list
sample_list = trimws(as.character(unique(somatic_data$Sample)))
sample_list = sort(sample_list)

res = matrix("", nrow=length(chr_list), ncol=length(sample_list))

for (i in 1:length(sample_list)) {
    for (j in 1:length(chr_list)) {
		tmp = subset(somatic_data, trimws(as.character(Sample))==sample_list[i]&trimws(as.character(CHROM))==chr_list[j])
		res[j, i] = nrow(tmp)		
	}
}

colnames(res) = sample_list
rownames(res) = chr_list
head(res); dim(res); 

sum(as.numeric(res))

write.csv(res, paste0(save_folder, "AML---mutation_distribution", ".csv"), row.names = TRUE)

###############
# recurrent mutations
uid_statistics = count(somatic_data$uid); x=uid_statistics[uid_statistics$freq==1,]; dim(x); head(x); tail(x)
head(uid_statistics); dim(uid_statistics); 
uid_statistics = uid_statistics[uid_statistics$freq>1,]

uid_statistics = uid_statistics[order(uid_statistics$freq, decreasing = FALSE),]
head(uid_statistics); dim(uid_statistics); 
colnames(uid_statistics)[1] = "recurrent_mutation_position"

need = count(uid_statistics$freq)
colnames(need)[1] = "#_of_samples_in_this_recurrent"
need

write.csv(uid_statistics, paste0(save_folder, "AML---recurrent_mutation_statistics", ".csv"), row.names = FALSE)
write.csv(need, paste0(save_folder, "AML---recurrent_mutation_distribution", ".csv"), row.names = FALSE)

# candidates
recurrent_candidates = subset(somatic_data, is.element(trimws(as.character(uid)), trimws(as.character(uid_statistics[,1]))))
head(recurrent_candidates); dim(recurrent_candidates)
write.csv(recurrent_candidates, paste0(save_folder, "AML---recurrent_mutation_ALL_CANDIDATES", ".csv"), row.names = FALSE)

###############
###############
###############
# get recurrent table for 2 samples

rec_table = matrix(0, nrow=length(sample_list), ncol=length(sample_list))
rec_table_gene = matrix("-", nrow=length(sample_list), ncol=length(sample_list))

colnames(rec_table) = sample_list
rownames(rec_table) = sample_list
colnames(rec_table_gene) = sample_list
rownames(rec_table_gene) = sample_list

need_uid = uid_statistics[uid_statistics$freq==2,]
head(need_uid); dim(need_uid); 

i = 1
for (i in 1:length(need_uid[,1])) {
	tmp = subset(recurrent_candidates, trimws(as.character(uid))==trimws(as.character(need_uid[i, 1])))
	rec_samples = sort(trimws(as.character(tmp$Sample)))
	
	row_idx = which(rec_samples[1]==rownames(rec_table)); row_idx
	col_idx = which(rec_samples[2]==colnames(rec_table)); col_idx
	rec_table[row_idx, col_idx] = rec_table[row_idx, col_idx] + 1
	rec_table_gene[row_idx, col_idx] = paste0(rec_table_gene[row_idx, col_idx], "; ", tmp$Involved_genes[1])
}

write.csv(rec_table_gene, paste0(save_folder, "AML---GENE---Table_with_mutation_statistics_of_2_recurrent_samples", ".csv"), row.names = TRUE)

# unique id
# later should remove the # with more than 2 reccurent samples
i = 1
for (i in 1:length(sample_list)) {
	tmp = subset(somatic_data, trimws(as.character(Sample))==trimws(as.character(sample_list[i])))
	row_idx = which(sample_list[i]==rownames(rec_table)); row_idx
	rec_table[row_idx, row_idx] = nrow(tmp) - sum(rec_table[row_idx, ])
}

# later should remove the # with more than 2 reccurent samples
write.csv(rec_table, paste0(save_folder, "AML---Table_with_mutation_statistics_of_2_recurrent_samples", ".csv"), row.names = TRUE)

# candidates with recurrent samples > 2
specific_uid = uid_statistics[uid_statistics$freq>2,]
specific_uid

tmp = subset(somatic_data, is.element(trimws(as.character(uid)), trimws(as.character(specific_uid$recurrent_mutation_position))))

write.csv(tmp, paste0(save_folder, "AML---mutations_with_MORE_THAN_3_recurrent_samples", ".csv"), row.names = TRUE)

# not do this
# later should remove the # with more than 2 reccurent samples
#a = count(tmp$Sample)
#i = 1
#for (i in 1:length(a[,1])) {
#	row_idx = which(trimws(as.character(a[i, 1]))==rownames(rec_table)); row_idx
#    rec_table[row_idx, row_idx] = rec_table[row_idx, row_idx] - a[i, 2]
#}

#write.csv(rec_table, paste0(save_folder, "AML---Table_with_mutation_statistics_of_2_recurrent_samples", ".csv"), row.names = TRUE)


######
head(somatic_data); dim(somatic_data); 

# gene checking
gene_list=c("CDKN2A", "CDKN2B", "CCNE1", "CD274", "PDCD1LG2", "JAK2", "RAC1", "KDR", "KIT", "PDGFRA", "TERT", "PIK3CA", "EIF1AX", "BRAF", "PPM1D", "CHEK2", "TP53", "RET", "MAPK1", "PTEN", "NRAS", "HRAS", "KRAS", "PPM1D", "ARID1B", "MLL", "BDP1", "TG", "ZFHX3", "ATM", "RB1", "TSHR", "EZH1", "MEN1", "CDH4", "SPOP", "MLL3", "APC", "PPARG", "NTRK1", "NTRK3", "ALK", "LTK", "MET", "FGFR2", "THADA", "PAX8", "AKT2", "AKT1", "NF1", "NF2", "AGK", "FAM114A2", "LTK", "PTC1", "PTC3", "PIK3A", "MMR", "TBX3", "LEAP1", "RBM10", "GNAS", "ETV6")

#
save_folder_1 = paste0(save_folder, "/", "specific_gene_data")
dir.create(save_folder_1, showWarnings = FALSE)

gene_candidates = subset(somatic_data, is.element(trimws(as.character(Involved_genes)), trimws(as.character(gene_list))))
head(gene_candidates); dim(gene_candidates); 

write.csv(gene_candidates, paste0(save_folder, "AML---mutations_in_specific_genes", ".csv"), row.names = TRUE)


###############################################################
###############################################################
###############################################################
# some mutations have several involved genes, here, separate them into individidual rows
# the columns of "Involved_genes", "Annot", "Annot_Impact" have multiple terms, hard to do statistics anaysis
# separate them first as normalized data


# 1. separate "Involved_genes"
work_data = somatic_data
try_statistics = c()

i = 1
for (i in 1:nrow(work_data)) {

    if (i%%1000==1) print(i)
	
    tmp = as.character(work_data$Involved_genes[i])
	all_Involved_genes = trimws(unlist(strsplit(tmp, ";")))
	all_Involved_genes
	
	# here only consider, "all_Involved_genes" is not null, otherwise, need to check
	for (j in 1:length(all_Involved_genes)) {
		one_row = data.frame(sample_gene=paste0(trimws(as.character(work_data$Sample[i])), "-", all_Involved_genes[j]), gene=all_Involved_genes[j], work_data[i, ])
		try_statistics = rbind(try_statistics, one_row)
	}
	
}

normalized_data_statistics = try_statistics
head(normalized_data_statistics); dim(normalized_data_statistics)


# 2. separate "Annot"
work_data = normalized_data_statistics
try_statistics = c()

i = 1
for (i in 1:nrow(work_data)) {

    if (i%%1000==1) print(i)
	
    tmp = as.character(work_data$Annot[i])
	all_Annot = trimws(unlist(strsplit(tmp, ";")))
	all_Annot
	
	# here only consider, "Annot" is not null, otherwise, need to check
	for (j in 1:length(all_Annot)) {
		one_row = data.frame(separate_Annot=all_Annot[j], work_data[i, ])
		try_statistics = rbind(try_statistics, one_row)
	}
	
}

normalized_data_statistics = try_statistics
head(normalized_data_statistics); dim(normalized_data_statistics)


# 3. separate "Annot_Impact"
work_data = normalized_data_statistics
try_statistics = c()

i = 1
for (i in 1:nrow(work_data)) {

    if (i%%1000==1) print(i)
	
    tmp = as.character(work_data$Annot_Impact[i])
	all_Annot_Impact = trimws(unlist(strsplit(tmp, ";")))
	all_Annot_Impact
	
	# here only consider, "Annot_Impact" is not null, otherwise, need to check
	for (j in 1:length(all_Annot_Impact)) {
		one_row = data.frame(separate_Annot_Impact=all_Annot_Impact[j], work_data[i, ])
		try_statistics = rbind(try_statistics, one_row)
	}
	
}

normalized_data_statistics = try_statistics
head(normalized_data_statistics); dim(normalized_data_statistics)

write.csv(normalized_data_statistics, paste0(save_folder, "AML---ALL-candidates---with-normalized-format", ".csv"), row.names = TRUE)
write.table(normalized_data_statistics, paste0(save_folder, "AML---ALL-candidates---with-normalized-format", ".txt"), 
quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)




######################################################
######################################################
######################################################
### COSMIC gene: cosmic-grch37-v84-cancer_gene_census.csv
# occurrence frequency of annotated genes 

### COSMIC gene: cosmic-grch37-v84-cancer_gene_census.csv
### too long, separate them into 2 list and concatenate them
COSMIC_gene_list_1 = c("A1CF", "ABI1", "ABL1", "ABL2", "ACKR3", "ACSL3", "ACSL6", "ACVR1", "ACVR2A", "AFF1", "AFF3", "AFF4", "AKAP9", "AKT1", "AKT2", "AKT3", "ALDH2", "ALK", "AMER1", "ANK1", "APC", "APOBEC3B", "AR", "ARAF", "ARHGAP26", "ARHGAP5", "ARHGEF10", "ARHGEF10L", "ARHGEF12", "ARID1A", "ARID1B", "ARID2", "ARNT", "ASPSCR1", "ASXL1", "ASXL2", "ATF1", "ATIC", "ATM", "ATP1A1", "ATP2B3", "ATR", "ATRX", "AXIN1", "AXIN2", "B2M", "BAP1", "BARD1", "BAX", "BAZ1A", "BCL10", "BCL11A", "BCL11B", "BCL2", "BCL2L12", "BCL3", "BCL6", "BCL7A", "BCL9", "BCL9L", "BCLAF1", "BCOR", "BCORL1", "BCR", "BIRC3", "BIRC6", "BLM", "BMP5", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRD3", "BRD4", "BRIP1", "BTG1", "BTK", "BUB1B", "C15ORF65", "C2ORF44", "CACNA1D", "CALR", "CAMTA1", "CANT1", "CARD11", "CARS", "CASC5", "CASP3", "CASP8", "CASP9", "CBFA2T3", "CBFB", "CBL", "CBLB", "CBLC", "CCDC6", "CCNB1IP1", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCR4", "CCR7", "CD209", "CD274", "CD28", "CD74", "CD79A", "CD79B", "CDC73", "CDH1", "CDH10", "CDH11", "CDH17", "CDK12", "CDK4", "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2C", "CDX2", "CEBPA", "CEP89", "CHCHD7", "CHD2", "CHD4", "CHEK2", "CHIC2", "CHST11", "CIC", "CIITA", "CLIP1", "CLP1", "CLTC", "CLTCL1", "CNBD1", "CNBP", "CNOT3", "CNTNAP2", "CNTRL", "COL1A1", "COL2A1", "COL3A1", "COX6C", "CPEB3", "CREB1", "CREB3L1", "CREB3L2", "CREBBP", "CRLF2", "CRNKL1", "CRTC1", "CRTC3", "CSF1R", "CSF3R", "CSMD3", "CTCF", "CTNNA2", "CTNNB1", "CTNND1", "CTNND2", "CUL3", "CUX1", "CXCR4", "CYLD", "CYP2C8", "CYSLTR2", "DAXX", "DCAF12L2", "DCC", "DCTN1", "DDB2", "DDIT3", "DDR2", "DDX10", "DDX3X", "DDX5", "DDX6", "DEK", "DGCR8", "DICER1", "DNAJB1", "DNM2", "DNMT3A", "DROSHA", "DUX4L1", "EBF1", "ECT2L", "EED", "EGFR", "EIF1AX", "EIF3E", "EIF4A2", "ELF3", "ELF4", "ELK4", "ELL", "ELN", "EML4", "EP300", "EPAS1", "EPHA3", "EPHA7", "EPS15", "ERBB2", "ERBB3", "ERBB4", "ERC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERG", "ESR1", "ETNK1", "ETV1", "ETV4", "ETV5", "ETV6", "EWSR1", "EXT1", "EXT2", "EZH2", "EZR", "FAM131B", "FAM135B", "FAM46C", "FAM47C", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FAS", "FAT1", "FAT3", "FAT4", "FBLN2", "FBXO11", "FBXW7", "FCGR2B", "FCRL4", "FEN1", "FES", "FEV", "FGFR1", "FGFR1OP", "FGFR2", "FGFR3", "FGFR4", "FH", "FHIT", "FIP1L1", "FKBP9", "FLCN", "FLI1", "FLNA", "FLT3", "FLT4", "FNBP1", "FOXA1", "FOXL2", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FOXR1", "FSTL3", "FUBP1", "FUS", "GAS7", "GATA1", "GATA2", "GATA3", "GLI1", "GMPS", "GNA11", "GNAQ", "GNAS", "GOLGA5", "GOPC", "GPC3", "GPC5", "GPHN", "GRIN2A", "GRM3", "H3F3A", "H3F3B", "HERPUD1", "HEY1", "HIF1A", "HIP1", "HIST1H3B", "HIST1H4I", "HLA-A", "HLF", "HMGA1", "HMGA2", "HMGN2P46", "HNF1A", "HNRNPA2B1", "HOOK3", "HOXA11", "HOXA13", "HOXA9", "HOXC11", "HOXC13", "HOXD11", "HOXD13", "HRAS", "HSP90AA1", "HSP90AB1", "ID3", "IDH1", "IDH2", "IGF2BP2", "IGH", "IGK", "IGL", "IKBKB", "IKZF1", "IL2", "IL21R", "IL6ST", "IL7R", "IRF4", "IRS4", "ISX", "ITGAV", "ITK", "JAK1", "JAK2", "JAK3", "JAZF1", "JUN", "KAT6A", "KAT6B", "KAT7", "KCNJ5", "KDM5A", "KDM5C", "KDM6A", "KDR", "KDSR", "KEAP1", "KIAA1549", "KIAA1598", "KIF5B", "KIT", "KLF4", "KLF6", "KLK2", "KMT2A", "KMT2C", "KMT2D", "KNSTRN", "KRAS", "KTN1", "LARP4B", "LASP1", "LCK", "LCP1", "LEF1", "LEPROTL1", "LHFP", "LIFR", "LMNA", "LMO1", "LMO2", "LPP", "LRIG3", "LRP1B", "LSM14A", "LYL1", "LZTR1", "MAF", "MAFB", "MALAT1", "MALT1", "MAML2", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAPK1", "MAX", "MB21D2", "MDM2", "MDM4", "MDS2", "MECOM", "MED12")

COSMIC_gene_list_2 = c("MEN1", "MET", "MGMT", "MITF", "MKL1", "MLF1", "MLH1", "MLLT1", "MLLT10", "MLLT11", "MLLT3", "MLLT4", "MLLT6", "MN1", "MNX1", "MPL", "MSH2", "MSH6", "MSI2", "MSN", "MTCP1", "MTOR", "MUC1", "MUC16", "MUC4", "MUTYH", "MYB", "MYC", "MYCL", "MYCN", "MYD88", "MYH11", "MYH9", "MYO5A", "MYOD1", "N4BP2", "NAB2", "NACA", "NBEA", "NBN", "NCKIPSD", "NCOA1", "NCOA2", "NCOA4", "NCOR1", "NCOR2", "NDRG1", "NF1", "NF2", "NFATC2", "NFE2L2", "NFIB", "NFKB2", "NFKBIE", "NIN", "NKX2-1", "NONO", "NOTCH1", "NOTCH2", "NPM1", "NR4A3", "NRAS", "NRG1", "NSD1", "NT5C2", "NTHL1", "NTRK1", "NTRK3", "NUMA1", "NUP214", "NUP98", "NUTM1", "NUTM2A", "NUTM2B", "OLIG2", "OMD", "P2RY8", "PABPC1", "PAFAH1B2", "PALB2", "PAX3", "PAX5", "PAX7", "PAX8", "PBRM1", "PBX1", "PCBP1", "PCM1", "PDCD1LG2", "PDE4DIP", "PDGFB", "PDGFRA", "PDGFRB", "PER1", "PHF6", "PHOX2B", "PICALM", "PIK3CA", "PIK3CB", "PIK3R1", "PIM1", "PLAG1", "PLCG1", "PML", "PMS1", "PMS2", "POLD1", "POLE", "POLG", "POLQ", "POT1", "POU2AF1", "POU5F1", "PPARG", "PPFIBP1", "PPM1D", "PPP2R1A", "PPP6C", "PRCC", "PRDM1", "PRDM16", "PRDM2", "PREX2", "PRF1", "PRKACA", "PRKAR1A", "PRKCB", "PRPF40B", "PRRX1", "PSIP1", "PTCH1", "PTEN", "PTK6", "PTPN11", "PTPN13", "PTPN6", "PTPRB", "PTPRC", "PTPRD", "PTPRK", "PTPRT", "PWWP2A", "QKI", "RABEP1", "RAC1", "RAD17", "RAD21", "RAD51B", "RAF1", "RALGDS", "RANBP2", "RAP1GDS1", "RARA", "RB1", "RBM10", "RBM15", "RECQL4", "REL", "RET", "RFWD3", "RGPD3", "RGS7", "RHOA", "RHOH", "RMI2", "RNF213", "RNF43", "ROBO2", "ROS1", "RPL10", "RPL22", "RPL5", "RPN1", "RSPO2", "RSPO3", "RUNDC2A", "RUNX1", "RUNX1T1", "S100A7", "SALL4", "SBDS", "SDC4", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SEPT5", "SEPT6", "SEPT9", "SET", "SETBP1", "SETD1B", "SETD2", "SF3B1", "SFPQ", "SFRP4", "SGK1", "SH2B3", "SH3GL1", "SIRPA", "SIX1", "SIX2", "SKI", "SLC34A2", "SLC45A3", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1", "SMARCD1", "SMARCE1", "SMC1A", "SMO", "SND1", "SOCS1", "SOX2", "SOX21", "SPECC1", "SPEN", "SPOP", "SRC", "SRGAP3", "SRSF2", "SRSF3", "SS18", "SS18L1", "SSX1", "SSX2", "SSX4", "STAG1", "STAG2", "STAT3", "STAT5B", "STAT6", "STIL", "STK11", "STRN", "SUFU", "SUZ12", "SYK", "TAF15", "TAL1", "TAL2", "TBL1XR1", "TBX3", "TCEA1", "TCF12", "TCF3", "TCF7L2", "TCL1A", "TEC", "TERT", "TET1", "TET2", "TFE3", "TFEB", "TFG", "TFPT", "TFRC", "TGFBR2", "THRAP3", "TLX1", "TLX3", "TMEM127", "TMPRSS2", "TNC", "TNFAIP3", "TNFRSF14", "TNFRSF17", "TOP1", "TP53", "TP63", "TPM3", "TPM4", "TPR", "TRA", "TRAF7", "TRB", "TRD", "TRIM24", "TRIM27", "TRIM33", "TRIP11", "TRRAP", "TSC1", "TSC2", "TSHR", "U2AF1", "UBR5", "USP44", "USP6", "USP8", "VAV1", "VHL", "VTI1A", "WAS", "WHSC1", "WHSC1L1", "WIF1", "WNK2", "WRN", "WT1", "WWTR1", "XPA", "XPC", "XPO1", "YWHAE", "ZBTB16", "ZCCHC8", "ZEB1", "ZFHX3", "ZMYM3", "ZNF198", "ZNF278", "ZNF331", "ZNF384", "ZNF429", "ZNF479", "ZNF521", "ZNRF3", "ZRSR2")

COSMIC_gene_list = c(COSMIC_gene_list_1, COSMIC_gene_list_2)
length(COSMIC_gene_list); head(COSMIC_gene_list); tail(COSMIC_gene_list)

### for all mutations
# statistic for mutations number in each gene
# too many
#x = count(normalized_data_statistics$gene)
#head(x); dim(x)

###  COSMIC genes
cosmic_candidates = subset(normalized_data_statistics, is.element(toupper(gene), toupper(COSMIC_gene_list)))
dim(cosmic_candidates); dim(normalized_data_statistics);

# normalized_data_statistics may count the same mutations multiple times
# so, check back the original "somatic_data" data
cosmic_uid_set = trimws(as.character(cosmic_candidates$uid))

cosmic_candidates_original = subset(somatic_data, is.element(uid, cosmic_uid_set))
head(cosmic_candidates_original); dim(cosmic_candidates_original)

# save these
write.csv(cosmic_candidates_original, paste0(save_folder, "AML---COSMIC_GENE_mutations", ".csv"), row.names = FALSE)
write.table(cosmic_candidates_original, paste0(save_folder, "AML---COSMIC_GENE_mutations", ".txt"), 
quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

# statisticcally analyze the occurence frequency for COSMMIC gene using the data "cosmic_candidates"
# filter out with COSMIC genes
# the columns of "Involved_genes", "Annot", "Annot_Impact" have multiple terms, hard to do statistics anaysis
# separate them first as normalized data
# this process repeat the same mutations many time. 
# remove duplicated with  "uid"        "Sample"
n1 = which("uid"==colnames(cosmic_candidates))
n2 = which("Sample"==colnames(cosmic_candidates))
df = cosmic_candidates
df = df[!duplicated(df[c(n1,n2)]),]

x = count(df$gene)
x = x[order(x$fre, decreasing = TRUE), ]
head(x); dim(x)

write.csv(x, paste0(save_folder, "AML---statistics---COSMIC_GENE_mutations", ".csv"), row.names = FALSE)
write.table(x, paste0(save_folder, "AML---statistics---COSMIC_GENE_mutations", ".txt"), 
quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)


######################################################
######################################################
### RECURRENT mutations
# check "somatic_data/cosmic_candidates_original" to remove unique
all_uid = as.character(cosmic_candidates_original$uid)
length(all_uid)

#
recurrent_uid = all_uid[duplicated(all_uid)]
length(recurrent_uid); length(unique(recurrent_uid)) 

# 
unique_uid_for_all = all_uid[!duplicated(all_uid)]
unique_uid = setdiff(unique_uid_for_all, unique(recurrent_uid))
length(unique_uid)

# save
recurrent_mutation_set = subset(cosmic_candidates_original, is.element(uid, unique(recurrent_uid)))
head(recurrent_mutation_set); dim(recurrent_mutation_set)

dim(cosmic_candidates_original)

unique_mutation_set = subset(cosmic_candidates_original, is.element(uid, unique(unique_uid)))
head(unique_mutation_set); dim(unique_mutation_set)


# save these
write.csv(recurrent_mutation_set, paste0(save_folder, "RECURRENT---AML---COSMIC_GENE_mutations", ".csv"), row.names = FALSE)
write.table(recurrent_mutation_set, paste0(save_folder, "RECURRENT---AML---COSMIC_GENE_mutations", ".txt"), 
quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

write.csv(unique_mutation_set, paste0(save_folder, "UNIQUE---AML---COSMIC_GENE_mutations", ".csv"), row.names = FALSE)
write.table(unique_mutation_set, paste0(save_folder, "UNIQUE---AML---COSMIC_GENE_mutations", ".txt"), 
quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)


# statisticcally analyze the occurence frequency for COSMMIC gene using the data "recurrent_mutation_set"
# filter out with COSMIC genes
recurrent_cosmic_candidates = subset(cosmic_candidates, is.element(uid, unique(recurrent_uid)))
head(recurrent_cosmic_candidates); dim(recurrent_cosmic_candidates)

# the columns of "Involved_genes", "Annot", "Annot_Impact" have multiple terms, hard to do statistics anaysis
# separate them first as normalized data
# this process repeat the same mutations many time. 
# remove duplicated with  "uid"        "Sample"
n1 = which("uid"==colnames(recurrent_cosmic_candidates))
n2 = which("Sample"==colnames(recurrent_cosmic_candidates))
df = recurrent_cosmic_candidates
df = df[!duplicated(df[c(n1,n2)]),]

x = count(df$gene)
x = x[order(x$fre, decreasing = TRUE), ]
head(x); dim(x)

write.csv(x, paste0(save_folder, "RECURRENT---AML---statistics---COSMIC_GENE_mutations", ".csv"), row.names = FALSE)
write.table(x, paste0(save_folder, "RECURRENT---AML---statistics---COSMIC_GENE_mutations", ".txt"), 
quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)


#########
#########  the annotation as “missense_variant” or the annotation impact as “High”.
# cosmic_candidates
head(cosmic_candidates); dim(cosmic_candidates)

pickup_missense_high = subset(cosmic_candidates, separate_Annot_Impact == "HIGH" | separate_Annot == "missense_variant")
dim(pickup_missense_high)

# this process repeat the same mutations many time. 
# remove duplicated with  "uid"        "Sample"
n1 = which("uid"==colnames(pickup_missense_high))
n2 = which("Sample"==colnames(pickup_missense_high))
pickup_missense_high = pickup_missense_high[!duplicated(pickup_missense_high[c(n1,n2)]),]
dim(pickup_missense_high)

write.csv(pickup_missense_high, paste0(save_folder, "missense_variant_or_HIGH---AML---COSMIC_GENE_mutations", ".csv"), row.names = FALSE)
write.table(pickup_missense_high, paste0(save_folder, "missense_variant_or_HIGH---AML---COSMIC_GENE_mutations", ".txt"), 
quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)














#########

#########

# not filter with "missense"
a = count(somatic_data$Involved_genes); head(a)
b = subset(a, is.element(trimws(as.character(x)), trimws(as.character(gene_list))))

gene_table =  data.frame(b, Information="")
gene_table$Information = as.character(gene_table$Information)
colnames(gene_table)[1:2] = c("Gene", "Freq")
head(gene_table); dim(gene_table); 

i = 1
for (i in 1:nrow(gene_table)) {
	
	tmp = subset(somatic_data, trimws(as.character(Involved_genes))==gene_table$Gene[i])
	mut_info = c()
	for (j in 1:nrow(tmp)) {
	    x = paste0(trimws(as.character(tmp$Sample[j])), ", ", trimws(as.character(tmp$CHROM[j])), ":", tmp$POS[j], ", REF=", tmp$REF[j], ", ALT=", tmp$ALT[j], ", HGVS.c=", tmp$HGVS_c[j], ", HGVS.p=", tmp$HGVS_p[j], ", AF=", tmp$AF[j])
	    mut_info = c(mut_info, x)
	}
	gene_table$Information[i] = paste(mut_info, collapse="; ")
	
}

write.csv(gene_table, paste0(save_folder, "AML---STATISTICS_for_mutations_in_specific_genes", ".csv"), row.names = FALSE)


# missense
a = count(somatic_data$Involved_genes); head(a)
b = subset(a, is.element(trimws(as.character(x)), trimws(as.character(gene_list))))

# Annotation_Term="-", Annotation_Impact="-"
gene_table =  data.frame(b, Information="")
gene_table$Information = as.character(gene_table$Information)
colnames(gene_table)[1:2] = c("Gene", "Freq")
gene_table$Freq = as.character(gene_table$Freq)

head(gene_table); dim(gene_table); 

i = 1

for (i in 1:nrow(gene_table)) {
	
	tmp = subset(somatic_data, trimws(as.character(Involved_genes))==gene_table$Gene[i])
	mut_info = c()
	
	# count 
	n = 0
	
	for (j in 1:nrow(tmp)) {
		Annot = tmp$Annot[j]
	    Annot_Impact = tmp$Annot_Impact[j]
		
		if (Annot=="missense_variant") {
	        x = paste0(trimws(as.character(tmp$Sample[j])), ", ", trimws(as.character(tmp$CHROM[j])), ":", tmp$POS[j], ", REF=", tmp$REF[j], ", ALT=", tmp$ALT[j], ", HGVS.c=", tmp$HGVS_c[j], ", HGVS.p=", tmp$HGVS_p[j], ", AF=", tmp$AF[j], ", Annotation=", tmp$Annot[j], ", Annotation_Impact=", tmp$Annot_Impact[j])
	        mut_info = c(mut_info, x)
			n = n + 1
		}	 
	}
	
	gene_table$Information[i] = paste(mut_info, collapse="; ")
	gene_table$Freq[i] = n

}


write.csv(gene_table, paste0(save_folder, "missense-AML---STATISTICS_for_mutations_in_specific_genes", ".csv"), row.names = FALSE)


###############
###############
###############
# statistics for recurreng genes


############################################
###########    ----filter out "missense"------
############################################
missense_set = c()
i = 1
for (i in 1:nrow(cosmic_candidates)) {
   tmp = trimws(as.character(cosmic_candidates$Annot[i]))
   pos = str_locate(tmp, "missense_variant")
   if (!is.na(pos[1])) missense_set = rbind(missense_set, cosmic_candidates[i,])
}
head(missense_set); dim(missense_set)

#
write.csv(missense_set, paste0(save_folder, "missense---ALL-candidates----AML---in-COSMIC_GENEs", ".csv"), row.names = TRUE)







############################################
###########    follow ups with NOT filter out "missense"
############################################

# statistics for reccurent genes, where there are mutations in different samples
# one sample may have more than 1 mutations in 1 gene
# remove such row
tmp = cosmic_candidates[!duplicated(cosmic_candidates[ , c("sample_gene")]),]
dim(tmp)

#tmp = as.character(unique(cosmic_candidates$sample_gene))
#length(cosmic_candidates$sample_gene); length(tmp)

y = count(tmp$gene)
y = y[order(y$fre, decreasing = TRUE), ]
head(y); dim(y)

y = data.frame(y, Sample=" ")
y$Sample = as.character(y$Sample)


i = 1
for (i in 1:nrow(y)) {
   target = subset(tmp, trimws(as.character(gene))==trimws(as.character(y[i,1])))
   sample_list = sort(target$Sample)
   y$Sample[i] = paste(trimws(as.character(sample_list)), collapse="; ")
}

write.csv(y, paste0(save_folder, "AML---recurrent_genes---COSMIC_mutations", ".csv"), row.names = FALSE)


### filter out the 2 recurrent mutations with COSMIC genes

rec_table = matrix(0, nrow=length(sample_list), ncol=length(sample_list))
rec_table_gene = matrix("-", nrow=length(sample_list), ncol=length(sample_list))

colnames(rec_table) = sample_list
rownames(rec_table) = sample_list
colnames(rec_table_gene) = sample_list
rownames(rec_table_gene) = sample_list

need_uid = uid_statistics[uid_statistics$freq==2,]
head(need_uid); dim(need_uid); 

a = c()

i = 1
for (i in 1:length(need_uid[,1])) {
	tmp = subset(recurrent_candidates, trimws(as.character(uid))==trimws(as.character(need_uid[i, 1])))
	rec_samples = sort(trimws(as.character(tmp$Sample)))
	
	row_idx = which(rec_samples[1]==rownames(rec_table)); row_idx
	col_idx = which(rec_samples[2]==colnames(rec_table)); col_idx
	rec_table[row_idx, col_idx] = rec_table[row_idx, col_idx] + 1
	
	# COSMIC gene
	pos = which(trimws(as.character(tmp$Involved_genes[1]))==COSMIC_gene_list)
	if (length(pos)!=0) {
	    rec_table_gene[row_idx, col_idx] = paste0(rec_table_gene[row_idx, col_idx], "; ", tmp$Involved_genes[1])
		a = rbind(a, tmp)
	}
}

write.csv(rec_table_gene, paste0(save_folder, "AML---COSMIC_GENE---Table_with_mutation_statistics_of_2_recurrent_samples", ".csv"), row.names = TRUE)


write.csv(a, paste0(save_folder, "condidtates--for---AML---COSMIC_GENE---Table_with_mutation_statistics_of_2_recurrent_samples", ".csv"), row.names = TRUE)































#
a = count(somatic_data$Involved_genes); head(a)
b = subset(a, is.element(trimws(as.character(x)), trimws(as.character(gene_list))))

gene_table =  data.frame(b, Samples="")
gene_table$Samples = as.character(gene_table$Samples)
colnames(gene_table)[1:2] = c("Gene", "Freq")
head(gene_table); dim(gene_table); 

i = 1
for (i in 1:nrow(gene_table)) {
	
	tmp = subset(somatic_data, trimws(as.character(Involved_genes))==gene_table$Gene[i])
	gene_table$Samples[i] = print(paste(trimws(as.character(tmp$Sample)), collapse=";"))
	
}


write.csv(gene_table, paste0(save_folder, "AML---STATISTICS_for_mutations_in_specific_genes", ".csv"), row.names = TRUE)





#
head(somatic_data); dim(somatic_data)
i = 1
for (i in 1:length(gene_list)) {
    res = c()
    for (j in 1:nrow(somatic_data)) {
		tmp = trimws(as.character(somatic_data$Involved_genes))
		get_flag = str_locate(tmp, gene_list[i])
		if (!is.na(get_flag[1])) {
		    res = rbind(res, somatic_data[j, ])
		}		
    }
	
}













# filter out SOMATICSCORE>40
dim(somatic_data)
res = subset(somatic_data, SOMATICSCORE>40)
dim(res)

write.table(res, paste0(save_folder, "SOMATICSCORE-greater-40---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

# ANN is too huge, not save
res_1 = res[, -9]
write.table(res_1, paste0(save_folder, "no-ANN---SOMATICSCORE-greater-40---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)

write.csv(res_1, paste0(save_folder, "no-ANN---SOMATICSCORE-greater-40---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)


# get each sample hist for type
sample_list = as.character(unique(res_1$Sample))
sample_list

need = matrix("", nrow = length(sample_list)+1, ncol = 4) 
i = 1
table(res_1$SVTYPE)

# 4 categories: BND DEL DUP INV
categories = c("BND", "DEL", "DUP" ,"INV")

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

write.csv(need, paste0(save_folder, "SV-categories---SOMATICSCORE-greater-40---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = TRUE)




# check reccurrent SVs
dim(somatic_data)
res = subset(somatic_data, SOMATICSCORE>40)
dim(res)

#
colnames(res)

res = data.frame(uid=paste(trimws(as.character(res$CHROM)),"-",  trimws(as.character(res$POS)), sep=""), res)
res = res[, -10]
head(res); dim(res); 

a=count(res$uid)
b = a[order(a$freq, decreasing=TRUE),]
c=b[b$freq>1,]

# reccurent subset
x = subset(res, is.element(uid, c$x))


write.csv(x, paste0(save_folder, "reccurrent---SOMATICSCORE-greater-40---extracted-information---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)







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




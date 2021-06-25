#
library(stringr)
library(plyr)

#data_folder="/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/extract_manta_SVs-July-6-2020---oncocytomas/result/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method/"
data_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/"

#save_folder="/projectsp/foran/Team_Chan/oncocytomas_SV_Hua_Job/oncocytomas/codes/extract_manta_SVs-July-6-2020---oncocytomas/result/somaticSV-analysis-and-summary/follow_up_analysis_with_snpEff-new-method/"
save_folder="/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/extract_manta_SVs-July-23-2020---AML/result/somaticSV-analysis-and-summary/"

# data
# somatic_data = read.csv(paste0(data_folder, "chr-1-22-X-Y-M---manta_SVs_with_PASS---snpEfff", ".txt"), header=TRUE, sep="\t")

# data with only "PASS"
somatic_data = read.csv(paste0(data_folder, "ALL_PASS_somaticSV_candidates", ".txt"), header=TRUE, sep=" ")
head(somatic_data); dim(somatic_data); 
unique(somatic_data$CHROM)

# filter out candidates with "chr1-22, X, Y, M"
chr_list = c(paste("", c(1:22), sep=""), "X", "Y", "M")
chr_list

somatic_data = subset(somatic_data, is.element(as.character(CHROM), chr_list))
head(somatic_data); dim(somatic_data); 
unique(somatic_data$CHROM)

#####
somatic_data = data.frame(somatic_data, BreakPoint_1=-1, BreakPoint_2=-1, BreakPoint_Gene_1="-", BreakPoint_Gene_2="-", SOMATICSCORE=-1, SVTYPE="-", SVLEN="-", whole_coordinate_gene_1="-", whole_coordinate_gene_2="-")
somatic_data$BreakPoint_1 = as.character(somatic_data$BreakPoint_1)
somatic_data$BreakPoint_2 = as.character(somatic_data$BreakPoint_2)
somatic_data$BreakPoint_Gene_1 = as.character(somatic_data$BreakPoint_Gene_1)
somatic_data$BreakPoint_Gene_2 = as.character(somatic_data$BreakPoint_Gene_2)
somatic_data$SOMATICSCORE = as.character(somatic_data$SOMATICSCORE)
somatic_data$SVTYPE = as.character(somatic_data$SVTYPE)
somatic_data$SVLEN = as.character(somatic_data$SVLEN)
somatic_data$whole_coordinate_gene_1 = as.character(somatic_data$whole_coordinate_gene_1)
somatic_data$whole_coordinate_gene_2 = as.character(somatic_data$whole_coordinate_gene_2)

#####
# load the gene coordinates
gene_coordinates_folder = "/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/gtf_files"
gene_coordinates = read.csv(paste0(gene_coordinates_folder, "/", "EnsDb.Hsapiens.v75", ".csv"), header=TRUE, sep=",")
head(gene_coordinates); dim(gene_coordinates);
 
#
i = 1
i = 2
i = 38 # multiple
i = 430 # multiple

for (i in 1:nrow(somatic_data)) {
#for (i in 38:38) {

	print(i)
	# BreakPoint_1
	somatic_data$BreakPoint_1[i] = paste0(somatic_data$CHROM[i], ":", somatic_data$POS[i])
	
	# BreakPoint_2
	# 1). END in INFO
	tmp = as.character(somatic_data$INFO[i])
	pos = str_locate(tmp, "END=")
	pos
	
	if (!is.na(pos[1])) {
		tmp = substr(tmp, pos[2] + 1, nchar(tmp)) 
		pos = str_locate(tmp, ";")
		somatic_data$BreakPoint_2[i] = paste0(somatic_data$CHROM[i], ":", substr(tmp, 1, pos[1] - 1))
	}
	
	# 2). in "ALT"; no "chr"
	tmp = as.character(somatic_data$ALT[i])
	
	#
	pos = str_locate(tmp, ":")
	init_pos_1 = str_locate(tmp, "\\[")
	init_pos_2 = str_locate(tmp, "]")
	if (!is.na(init_pos_1[2])) {
		tmp = substr(tmp, init_pos_1[2]+1, nchar(tmp)) 
	} else if (!is.na(init_pos_2[2])) {
		tmp = substr(tmp, init_pos_2[2]+1, nchar(tmp)) 
	}
	
	#
	pos1 = str_locate(tmp, "\\[")
	pos2 = str_locate(tmp, "]")
	if (!is.na(pos1[1])) {
		somatic_data$BreakPoint_2[i] = substr(tmp, 1, pos1[1] - 1) 
	} else if (!is.na(pos2[1])) {
		somatic_data$BreakPoint_2[i] = substr(tmp, 1, pos2[1] - 1) 
	}
	
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
	
	# ;SVLEN=
	tmp = as.character(somatic_data$INFO[i])
	tmp = substr(tmp, str_locate(tmp, "SVLEN=")[2] + 1, nchar(tmp)) 
	if (!is.na(tmp)) {
		tmp
	
		# check if it is in the end
		pos = str_locate(tmp, ";")
	
		if (is.na(pos[1])) {
			somatic_data$SVLEN[i] = tmp
		} else {
			somatic_data$SVLEN[i] = substr(tmp, 1, pos[1] - 1) 
		}
	}
	
	#############
	# whole_coordinate_gene_1
	tmp = unlist(strsplit(somatic_data$BreakPoint_1[i], ":"))
	tmp
	
	check_chr = tmp[1]
	check_coordinate = as.numeric(tmp[2])
	
	tmp = subset(gene_coordinates, seqnames==check_chr); dim(tmp)
	
	# "+" strand; ssume the start is smaller
	tmp_1 = data.frame(tmp, start_difference = check_coordinate - tmp$start, end_difference = tmp$end - check_coordinate)
		
	# "-" strand; assume the end is smaller
	tmp_2 = data.frame(tmp, start_difference = tmp$start - check_coordinate, end_difference = check_coordinate - tmp$end)
	
	
	# pick up
	# "+"
	target_1 = subset(tmp_1, start_difference>=0 & end_difference>=0); dim(target_1)
	# "-"
	target_2 = subset(tmp_2, start_difference>=0 & end_difference>=0); dim(target_2)

	target = rbind(target_1, target_2)
	target
	
	if (nrow(target)>0) {
	    target = data.frame(id = paste0(target$symbol, "-", "chr", target$seqnames, ":", target$start, "-", target$end), target)
	
	    # gene name
	    somatic_data$BreakPoint_Gene_1[i] = paste(unique(target$gene_name), collapse=";")
	    somatic_data$whole_coordinate_gene_1[i] = paste(unique(target$id), collapse=";")
    }
	
	#############
	# whole_coordinate_gene_2
	tmp = unlist(strsplit(somatic_data$BreakPoint_2[i], ":"))
	tmp
	
	check_chr = tmp[1]
	check_chr
	check_coordinate = as.numeric(tmp[2])
	
	tmp = subset(gene_coordinates, seqnames==check_chr); dim(tmp)
	
	# "+" strand; ssume the start is smaller
	tmp_1 = data.frame(tmp, start_difference = check_coordinate - tmp$start, end_difference = tmp$end - check_coordinate)
		
	# "-" strand; assume the end is smaller
	tmp_2 = data.frame(tmp, start_difference = tmp$start - check_coordinate, end_difference = check_coordinate - tmp$end)
	
	
	# pick up
	# "+"
	target_1 = subset(tmp_1, start_difference>=0 & end_difference>=0)
	# "-"
	target_2 = subset(tmp_2, start_difference>=0 & end_difference>=0)

	target = rbind(target_1, target_2)
	target
	
	if (nrow(target)>0) {
	    target = data.frame(id = paste0(target$symbol, "-", "chr", target$seqnames, ":", target$start, "-", target$end), target)
	
	    # gene name
	    somatic_data$BreakPoint_Gene_2[i] = paste(unique(target$gene_name), collapse=";")
	    somatic_data$whole_coordinate_gene_2[i] = paste(unique(target$id), collapse=";")
	}
		
}

# 
head(somatic_data); dim(somatic_data)

#
write.csv(somatic_data, paste0(save_folder, "Breakpoint-gene---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)

#

####
# load the data above
somatic_data = read.csv(paste0(save_folder, "/", "Breakpoint-gene---AML---ALL_PASS_somaticSV_candidates", ".csv"), header=TRUE, sep=",")
head(somatic_data); dim(somatic_data); 


##########################################
##########################################
##########################################
# check if the SVs have a COMSMIC gene

##############
### COSMIC gene: cosmic-grch37-v84-cancer_gene_census.csv
### too long, separate them into 2 list and concatenate them
COSMIC_gene_list_1 = c("A1CF", "ABI1", "ABL1", "ABL2", "ACKR3", "ACSL3", "ACSL6", "ACVR1", "ACVR2A", "AFF1", "AFF3", "AFF4", "AKAP9", "AKT1", "AKT2", "AKT3", "ALDH2", "ALK", "AMER1", "ANK1", "APC", "APOBEC3B", "AR", "ARAF", "ARHGAP26", "ARHGAP5", "ARHGEF10", "ARHGEF10L", "ARHGEF12", "ARID1A", "ARID1B", "ARID2", "ARNT", "ASPSCR1", "ASXL1", "ASXL2", "ATF1", "ATIC", "ATM", "ATP1A1", "ATP2B3", "ATR", "ATRX", "AXIN1", "AXIN2", "B2M", "BAP1", "BARD1", "BAX", "BAZ1A", "BCL10", "BCL11A", "BCL11B", "BCL2", "BCL2L12", "BCL3", "BCL6", "BCL7A", "BCL9", "BCL9L", "BCLAF1", "BCOR", "BCORL1", "BCR", "BIRC3", "BIRC6", "BLM", "BMP5", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRD3", "BRD4", "BRIP1", "BTG1", "BTK", "BUB1B", "C15ORF65", "C2ORF44", "CACNA1D", "CALR", "CAMTA1", "CANT1", "CARD11", "CARS", "CASC5", "CASP3", "CASP8", "CASP9", "CBFA2T3", "CBFB", "CBL", "CBLB", "CBLC", "CCDC6", "CCNB1IP1", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCR4", "CCR7", "CD209", "CD274", "CD28", "CD74", "CD79A", "CD79B", "CDC73", "CDH1", "CDH10", "CDH11", "CDH17", "CDK12", "CDK4", "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2C", "CDX2", "CEBPA", "CEP89", "CHCHD7", "CHD2", "CHD4", "CHEK2", "CHIC2", "CHST11", "CIC", "CIITA", "CLIP1", "CLP1", "CLTC", "CLTCL1", "CNBD1", "CNBP", "CNOT3", "CNTNAP2", "CNTRL", "COL1A1", "COL2A1", "COL3A1", "COX6C", "CPEB3", "CREB1", "CREB3L1", "CREB3L2", "CREBBP", "CRLF2", "CRNKL1", "CRTC1", "CRTC3", "CSF1R", "CSF3R", "CSMD3", "CTCF", "CTNNA2", "CTNNB1", "CTNND1", "CTNND2", "CUL3", "CUX1", "CXCR4", "CYLD", "CYP2C8", "CYSLTR2", "DAXX", "DCAF12L2", "DCC", "DCTN1", "DDB2", "DDIT3", "DDR2", "DDX10", "DDX3X", "DDX5", "DDX6", "DEK", "DGCR8", "DICER1", "DNAJB1", "DNM2", "DNMT3A", "DROSHA", "DUX4L1", "EBF1", "ECT2L", "EED", "EGFR", "EIF1AX", "EIF3E", "EIF4A2", "ELF3", "ELF4", "ELK4", "ELL", "ELN", "EML4", "EP300", "EPAS1", "EPHA3", "EPHA7", "EPS15", "ERBB2", "ERBB3", "ERBB4", "ERC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERG", "ESR1", "ETNK1", "ETV1", "ETV4", "ETV5", "ETV6", "EWSR1", "EXT1", "EXT2", "EZH2", "EZR", "FAM131B", "FAM135B", "FAM46C", "FAM47C", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FAS", "FAT1", "FAT3", "FAT4", "FBLN2", "FBXO11", "FBXW7", "FCGR2B", "FCRL4", "FEN1", "FES", "FEV", "FGFR1", "FGFR1OP", "FGFR2", "FGFR3", "FGFR4", "FH", "FHIT", "FIP1L1", "FKBP9", "FLCN", "FLI1", "FLNA", "FLT3", "FLT4", "FNBP1", "FOXA1", "FOXL2", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FOXR1", "FSTL3", "FUBP1", "FUS", "GAS7", "GATA1", "GATA2", "GATA3", "GLI1", "GMPS", "GNA11", "GNAQ", "GNAS", "GOLGA5", "GOPC", "GPC3", "GPC5", "GPHN", "GRIN2A", "GRM3", "H3F3A", "H3F3B", "HERPUD1", "HEY1", "HIF1A", "HIP1", "HIST1H3B", "HIST1H4I", "HLA-A", "HLF", "HMGA1", "HMGA2", "HMGN2P46", "HNF1A", "HNRNPA2B1", "HOOK3", "HOXA11", "HOXA13", "HOXA9", "HOXC11", "HOXC13", "HOXD11", "HOXD13", "HRAS", "HSP90AA1", "HSP90AB1", "ID3", "IDH1", "IDH2", "IGF2BP2", "IGH", "IGK", "IGL", "IKBKB", "IKZF1", "IL2", "IL21R", "IL6ST", "IL7R", "IRF4", "IRS4", "ISX", "ITGAV", "ITK", "JAK1", "JAK2", "JAK3", "JAZF1", "JUN", "KAT6A", "KAT6B", "KAT7", "KCNJ5", "KDM5A", "KDM5C", "KDM6A", "KDR", "KDSR", "KEAP1", "KIAA1549", "KIAA1598", "KIF5B", "KIT", "KLF4", "KLF6", "KLK2", "KMT2A", "KMT2C", "KMT2D", "KNSTRN", "KRAS", "KTN1", "LARP4B", "LASP1", "LCK", "LCP1", "LEF1", "LEPROTL1", "LHFP", "LIFR", "LMNA", "LMO1", "LMO2", "LPP", "LRIG3", "LRP1B", "LSM14A", "LYL1", "LZTR1", "MAF", "MAFB", "MALAT1", "MALT1", "MAML2", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", "MAPK1", "MAX", "MB21D2", "MDM2", "MDM4", "MDS2", "MECOM", "MED12")

COSMIC_gene_list_2 = c("MEN1", "MET", "MGMT", "MITF", "MKL1", "MLF1", "MLH1", "MLLT1", "MLLT10", "MLLT11", "MLLT3", "MLLT4", "MLLT6", "MN1", "MNX1", "MPL", "MSH2", "MSH6", "MSI2", "MSN", "MTCP1", "MTOR", "MUC1", "MUC16", "MUC4", "MUTYH", "MYB", "MYC", "MYCL", "MYCN", "MYD88", "MYH11", "MYH9", "MYO5A", "MYOD1", "N4BP2", "NAB2", "NACA", "NBEA", "NBN", "NCKIPSD", "NCOA1", "NCOA2", "NCOA4", "NCOR1", "NCOR2", "NDRG1", "NF1", "NF2", "NFATC2", "NFE2L2", "NFIB", "NFKB2", "NFKBIE", "NIN", "NKX2-1", "NONO", "NOTCH1", "NOTCH2", "NPM1", "NR4A3", "NRAS", "NRG1", "NSD1", "NT5C2", "NTHL1", "NTRK1", "NTRK3", "NUMA1", "NUP214", "NUP98", "NUTM1", "NUTM2A", "NUTM2B", "OLIG2", "OMD", "P2RY8", "PABPC1", "PAFAH1B2", "PALB2", "PAX3", "PAX5", "PAX7", "PAX8", "PBRM1", "PBX1", "PCBP1", "PCM1", "PDCD1LG2", "PDE4DIP", "PDGFB", "PDGFRA", "PDGFRB", "PER1", "PHF6", "PHOX2B", "PICALM", "PIK3CA", "PIK3CB", "PIK3R1", "PIM1", "PLAG1", "PLCG1", "PML", "PMS1", "PMS2", "POLD1", "POLE", "POLG", "POLQ", "POT1", "POU2AF1", "POU5F1", "PPARG", "PPFIBP1", "PPM1D", "PPP2R1A", "PPP6C", "PRCC", "PRDM1", "PRDM16", "PRDM2", "PREX2", "PRF1", "PRKACA", "PRKAR1A", "PRKCB", "PRPF40B", "PRRX1", "PSIP1", "PTCH1", "PTEN", "PTK6", "PTPN11", "PTPN13", "PTPN6", "PTPRB", "PTPRC", "PTPRD", "PTPRK", "PTPRT", "PWWP2A", "QKI", "RABEP1", "RAC1", "RAD17", "RAD21", "RAD51B", "RAF1", "RALGDS", "RANBP2", "RAP1GDS1", "RARA", "RB1", "RBM10", "RBM15", "RECQL4", "REL", "RET", "RFWD3", "RGPD3", "RGS7", "RHOA", "RHOH", "RMI2", "RNF213", "RNF43", "ROBO2", "ROS1", "RPL10", "RPL22", "RPL5", "RPN1", "RSPO2", "RSPO3", "RUNDC2A", "RUNX1", "RUNX1T1", "S100A7", "SALL4", "SBDS", "SDC4", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SEPT5", "SEPT6", "SEPT9", "SET", "SETBP1", "SETD1B", "SETD2", "SF3B1", "SFPQ", "SFRP4", "SGK1", "SH2B3", "SH3GL1", "SIRPA", "SIX1", "SIX2", "SKI", "SLC34A2", "SLC45A3", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1", "SMARCD1", "SMARCE1", "SMC1A", "SMO", "SND1", "SOCS1", "SOX2", "SOX21", "SPECC1", "SPEN", "SPOP", "SRC", "SRGAP3", "SRSF2", "SRSF3", "SS18", "SS18L1", "SSX1", "SSX2", "SSX4", "STAG1", "STAG2", "STAT3", "STAT5B", "STAT6", "STIL", "STK11", "STRN", "SUFU", "SUZ12", "SYK", "TAF15", "TAL1", "TAL2", "TBL1XR1", "TBX3", "TCEA1", "TCF12", "TCF3", "TCF7L2", "TCL1A", "TEC", "TERT", "TET1", "TET2", "TFE3", "TFEB", "TFG", "TFPT", "TFRC", "TGFBR2", "THRAP3", "TLX1", "TLX3", "TMEM127", "TMPRSS2", "TNC", "TNFAIP3", "TNFRSF14", "TNFRSF17", "TOP1", "TP53", "TP63", "TPM3", "TPM4", "TPR", "TRA", "TRAF7", "TRB", "TRD", "TRIM24", "TRIM27", "TRIM33", "TRIP11", "TRRAP", "TSC1", "TSC2", "TSHR", "U2AF1", "UBR5", "USP44", "USP6", "USP8", "VAV1", "VHL", "VTI1A", "WAS", "WHSC1", "WHSC1L1", "WIF1", "WNK2", "WRN", "WT1", "WWTR1", "XPA", "XPC", "XPO1", "YWHAE", "ZBTB16", "ZCCHC8", "ZEB1", "ZFHX3", "ZMYM3", "ZNF198", "ZNF278", "ZNF331", "ZNF384", "ZNF429", "ZNF479", "ZNF521", "ZNRF3", "ZRSR2")

COSMIC_gene_list = c(COSMIC_gene_list_1, COSMIC_gene_list_2)
length(COSMIC_gene_list); head(COSMIC_gene_list); tail(COSMIC_gene_list)

########
# some candidates have multiple genes 
# reverse from "BreakPoint_Gene_2" to ""BreakPoint_Gene_1" 
# 1. separate "BreakPoint_Gene_2"
# COSMIC_GENE_in_Breakpoint_2:  if it is COSMIC gene, show, otherwise "-"
work_data = somatic_data
try_statistics = c()

i = 1
for (i in 1:nrow(work_data)) {

    if (i%%1000==1) print(i)
	
    tmp = as.character(work_data$BreakPoint_Gene_2[i])
	all_Involved_genes = trimws(unlist(strsplit(tmp, ";")))
	all_Involved_genes
	
	# here only consider, "all_Involved_genes" is not null, otherwise, need to check
	for (j in 1:length(all_Involved_genes)) {
		if (is.element(toupper(trimws(as.character(all_Involved_genes[j]))), toupper(trimws(as.character(COSMIC_gene_list))))) {
		    one_row = data.frame(COSMIC_GENE_in_Breakpoint_2=all_Involved_genes[j], work_data[i, ])
		} else {
		    one_row = data.frame(COSMIC_GENE_in_Breakpoint_2="-", work_data[i, ])
		}
		try_statistics = rbind(try_statistics, one_row)
	}
	
}

normalized_data = try_statistics
head(normalized_data); dim(normalized_data)


# 2. separate "BreakPoint_Gene_1"
# COSMIC_GENE_in_Breakpoint_1:  if it is COSMIC gene, show, otherwise "-"
work_data = normalized_data
try_statistics = c()

i = 1
for (i in 1:nrow(work_data)) {

    if (i%%1000==1) print(i)
	
    tmp = as.character(work_data$BreakPoint_Gene_1[i])
	all_Involved_genes = trimws(unlist(strsplit(tmp, ";")))
	all_Involved_genes
	
	# here only consider, "all_Involved_genes" is not null, otherwise, need to check
	for (j in 1:length(all_Involved_genes)) {
		if (is.element(toupper(trimws(as.character(all_Involved_genes[j]))), toupper(trimws(as.character(COSMIC_gene_list))))) {
   		    one_row = data.frame(COSMIC_GENE_in_Breakpoint_1=all_Involved_genes[j], work_data[i, ])
		} else {
		    one_row = data.frame(COSMIC_GENE_in_Breakpoint_1="-", work_data[i, ])
		}
		try_statistics = rbind(try_statistics, one_row)
	}
	
}

normalized_data = try_statistics
head(normalized_data); dim(normalized_data)

# check the COSMIC genes
tmp_1 = subset(normalized_data, is.element(toupper(trimws(as.character(COSMIC_GENE_in_Breakpoint_1))), toupper(trimws(as.character(COSMIC_gene_list)))) )
tmp_1 = data.frame(id=paste0(tmp_1$Sample, "-", tmp_1$CHROM, "-", tmp_1$POS, "-", tmp_1$BreakPoint_Gene_1, "-", tmp_1$BreakPoint_Gene_2, "-", tmp_1$COSMIC_GENE_in_Breakpoint_1, "-", tmp_1$COSMIC_GENE_in_Breakpoint_2), tmp_1)
head(tmp_1);dim(tmp_1)

#
tmp_2 = subset(normalized_data, is.element(toupper(trimws(as.character(COSMIC_GENE_in_Breakpoint_2))), toupper(trimws(as.character(COSMIC_gene_list)))) )
tmp_2 = data.frame(id=paste0(tmp_2$Sample, "-", tmp_2$CHROM, "-", tmp_2$POS, "-", tmp_2$BreakPoint_Gene_1, "-", tmp_2$BreakPoint_Gene_2, "-", tmp_2$COSMIC_GENE_in_Breakpoint_1, "-", tmp_2$COSMIC_GENE_in_Breakpoint_2), tmp_2)
head(tmp_2);dim(tmp_2)

# concatenate tmp_1 & tmp_2
tmp = rbind(tmp_1, tmp_2)
head(tmp);dim(tmp)
# unique id
length(unique(tmp$id))

# remove duplicated

res = tmp[!duplicated(tmp$id),]
dim(res)

colnames(res)[2:3] = c("COSMIC_GENE_in_Breakpoint_1", "COSMIC_GENE_in_Breakpoint_2")

#
write.csv(res, paste0(save_folder, "involved-COSMIC-gene-SVs---Breakpoint-gene---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)


##############
# check if close to CCND1
COSMIC_data = read.csv(paste0(save_folder, "/", "involved-COSMIC-gene-SVs---Breakpoint-gene---AML---ALL_PASS_somaticSV_candidates", ".csv"), header=TRUE, sep=",")
head(COSMIC_data); dim(COSMIC_data); 

# need to run some code in the beginning
CCND1_coordinate = subset(gene_coordinates, gene_name=="CCND1")
CCND1_coordinate

#     seqnames    start      end width strand         gene_id gene_name
#9611       11 69455855 69469242 13388      + ENSG00000110092     CCND1
#       gene_biotype seq_coord_system symbol
#9611 protein_coding       chromosome  CCND1

# pickup chromosome 11
target_data = data.frame(COSMIC_data, CCND1_chr = "11", CCND1_start = "69455855", CCND1_end = "69469242", chr_1="-", pos_1="-", chr_2="-", pos_2="-", same_chr1="No", same_chr2="No", distance_1="-", distance_2="-", distance_3="-", distance_4="-")
head(target_data); dim(target_data)

target_data$chr_1 = as.character(target_data$chr_1)
target_data$pos_1 = as.character(target_data$pos_1)
target_data$chr_2 = as.character(target_data$chr_2)
target_data$pos_2 = as.character(target_data$pos_2)
target_data$same_chr1 = as.character(target_data$same_chr1)
target_data$same_chr2 = as.character(target_data$same_chr2)
target_data$distance_1 = as.character(target_data$distance_1)
target_data$distance_2 = as.character(target_data$distance_2)
target_data$distance_3 = as.character(target_data$distance_3)
target_data$distance_4 = as.character(target_data$distance_4)


#
i = 1

for (i in 1:nrow(target_data)) {
	
	chr_1 = unlist(strsplit(as.character(trimws(COSMIC_data$BreakPoint_1[i])), "\\:"))[1]
	pos_1 = unlist(strsplit(as.character(trimws(COSMIC_data$BreakPoint_1[i])), "\\:"))[2]
	chr_2 = unlist(strsplit(as.character(trimws(COSMIC_data$BreakPoint_2[i])), "\\:"))[1]
	pos_2 = unlist(strsplit(as.character(trimws(COSMIC_data$BreakPoint_2[i])), "\\:"))[2]
	
	target_data$chr_1[i] = chr_1
	target_data$pos_1[i] = pos_1
	target_data$chr_2[i] = chr_2
	target_data$pos_2[i] = pos_2

    if (chr_1=="11") {
	    target_data$same_chr1[i] = "YES"
	   	target_data$distance_1[i] = as.numeric(as.character(target_data$pos_1[i])) - as.numeric(as.character(target_data$CCND1_start[i]))
		target_data$distance_2[i] = as.numeric(as.character(target_data$pos_1[i])) - as.numeric(as.character(target_data$CCND1_end[i]))
	}
    
	if (chr_2=="11") {
		target_data$same_chr2[i] = "YES"
	   	target_data$distance_3[i] = as.numeric(as.character(target_data$pos_2[i])) - as.numeric(as.character(target_data$CCND1_start[i]))
		target_data$distance_4[i] = as.numeric(as.character(target_data$pos_2[i])) - as.numeric(as.character(target_data$CCND1_end[i]))
	}

}

head(target_data); dim(target_data)


write.csv(target_data, paste0(save_folder, "compare-CCND1---involved-COSMIC-gene-SVs---Breakpoint-gene---AML---ALL_PASS_somaticSV_candidates", ".csv"), row.names = FALSE)




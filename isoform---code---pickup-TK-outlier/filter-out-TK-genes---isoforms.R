#
library(stringr)
library(plyr)

data_folder="/projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-Sample-in-ORIEN-Analysis-Workflow-August-26-2020/analysis/analysis_of_STAR-RSEM/combined_RSEM_expression_result/"

save_folder="/projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-Sample-in-ORIEN-Analysis-Workflow-August-26-2020/analysis/analysis_of_STAR-RSEM/isoform---result---pickup-TK-outlier/"

# data
#fusion_data = read.csv(paste0(data_folder, "annotate_gene--247-ORIEN---testYRIvsEU_perind_numers.counts", ".csv"), header=TRUE, sep=",")
fusion_data = read.csv(paste0(data_folder, "all_RSEM_isoforms_TPM", ".txt"), header=TRUE, sep=" ")

head(fusion_data); dim(fusion_data); 

# split gene_id
# https://stackoverflow.com/questions/40104734/split-a-whole-column-of-a-data-frame-and-keep-only-the-first-part
# not work: updated_data = data.frame(gene_name = unlist(strsplit(as.character(fusion_data$gene_id), "_"))[2], fusion_data)
updated_data = data.frame(transcipt_name = unlist(lapply(strsplit(as.character(fusion_data$transcript_id), "_"), '[[', 2)), fusion_data)

head(updated_data); dim(updated_data); 

#
updated_data = data.frame(gene_name = unlist(lapply(strsplit(as.character(updated_data$transcipt_name), "-"), '[[', 1)), updated_data)

head(updated_data); dim(updated_data); 

####### review paper gene
TK_gene_list = c("ALK", "AXL", "DDR", "DDR1", "DDR2", "EGFR", "ERB4", "ERBB4", "FGFR", "FGFR4", "FGFR3", "FGFR2", "FGFR1", "INSR", "MET", "RET", "RON", "NTRK3", "NTRK2", "NTRK1", "NTRK", "VEGFR", "MST1R", "FLT1", "VEGFR1", "KDR", "VEGFR2", "FLT4", "VEGFR3")

#######
TK_gene_list = c("LRRK1", "MLK4", "LRRK2", "SgK288", "KSR2", "ANKRD3", "TESK2", "ZAK", "HH498", "MLKL", "LIMK1", "ALK7", "LZK", "TGFbR2", "TESK1", "BMPR2", "TAK1", "MLK3", "RIPK3", "RIPK2", "RIPK1", "RAF1", "IRAK4", "MLK2", "MLK1", "MISR2", "LIMK2", "KSR1", "IRAK3", "IRAK2", "IRAK1", "ILK", "DLK", "BRAF", "ARAF", "BMPR1B", "TGFbR1", "ALK4", "BMPR1A", "ALK2", "ALK1", "ACTR2B", "ACTR2", "TYK2_b", "JAK3_b", "JAK2_b", "JAK1_b", "EphA6", "EphA10", "SuRTK106", "SRM", "PYK2", "FRK", "CTK", "BMX", "EphA7", "LMR3", "LMR2", "LMR1", "CCK4", "DDR2", "KDR", "DDR1", "ZAP70", "YES", "TYRO3", "TYK2", "TXK", "TRKC", "TRKB", "TRKA", "TNK1", "TIE1", "TIE2", "TEC", "SYK", "SRC", "RYK", "ROS", "ROR2", "ROR1", "RON", "RET", "PDGFRb", "PDGFRa", "MUSK", "MET", "MER", "LYN", "LTK", "LCK", "KIT", "JAK3", "JAK2", "JAK1", "ITK", "IRR", "INSR", "IGF1R", "ErbB4", "ErbB3", "ErbB2", "HCK", "FYN", "FLT4", "FLT1", "FLT3", "FGR", "FGFR4", "FGFR3", "FGFR2", "FGFR1", "FES", "FER", "FAK", "EphB6", "EphB4", "EphB3", "EphB2", "EphB1", "EphA8", "EphA5", "EphA4", "EphA3", "EphA2", "EphA1", "EGFR", "CSK", "FMS", "BTK", "BRK", "BLK", "AXL", "ABL2", "ALK", "ACK", "ABL1", "LRRK1", "MLK4", "LRRK2", "SgK288", "KSR2", "ANKRD3", "TESK2", "ZAK", "HH498", "MLKL", "LIMK1", "ALK7", "LZK", "TGFbR2", "TESK1", "BMPR2", "TAK1", "MLK3", "RIPK3", "RIPK2", "RIPK1", "RAF1", "IRAK4", "MLK2", "MLK1", "MISR2", "LIMK2", "KSR1", "IRAK3", "IRAK2", "IRAK1", "ILK", "DLK", "BRAF", "ARAF", "BMPR1B", "TGFbR1", "ALK4", "BMPR1A", "ALK2", "ALK1", "ACTR2B", "ACTR2", "Domain2_TYK2", "Domain2_JAK3", "Domain2_JAK2", "Domain2_JAK1", "EphA6", "EphA10", "SuRTK106", "SRM", "PYK2", "FRK", "CTK", "BMX", "EphA7", "LMR3", "LMR2", "LMR1", "CCK4", "DDR2", "KDR", "DDR1", "ZAP70", "YES", "TYRO3", "TYK2", "TXK", "TRKC", "TRKB", "TRKA", "TNK1", "TIE1", "TIE2", "TEC", "SYK", "SRC", "RYK", "ROS", "ROR2", "ROR1", "RON", "RET", "PDGFRb", "PDGFRa", "MUSK", "MET", "MER", "LYN", "LTK", "LCK", "KIT", "JAK3", "JAK2", "JAK1", "ITK", "IRR", "INSR", "IGF1R", "HER4", "HER3", "HER2", "HCK", "FYN", "FLT4", "FLT1", "FLT3", "FGR", "FGFR4", "FGFR3", "FGFR2", "FGFR1", "FES", "FER", "FAK", "EphB6", "EphB4", "EphB3", "EphB2", "EphB1", "EphA8", "EphA5", "EphA4", "EphA3", "EphA2", "EphA1", "EGFR", "CSK", "FMS", "BTK", "BRK", "BLK", "AXL", "ARG", "ALK", "ACK", "ABL", "LRRK1", "LRRK2", "ANKK1", "KSR2", "RIPK4", "TESK2", "TNNI3K", "MLKL", "LIMK1", "ACVR1C", "MAP3K13", "TGFBR2", "TESK1", "BMPR2", "MAP3K7", "MAP3K11", "RIPK3", "RIPK2", "RIPK1", "RAF1", "IRAK4", "MAP3K10", "MAP3K9", "AMHR2", "LIMK2", "KSR1", "IRAK3", "IRAK2", "IRAK1", "ILK", "MAP3K12", "BRAF", "ARAF", "BMPR1B", "TGFBR1", "ACVR1B", "BMPR1A", "ACVR1", "ACVRL1", "ACVR2B", "ACVR2A", "TYK2", "JAK3", "JAK2", "JAK1", "EPHA6", "EPHA10", "STYK1", "SRMS", "PTK2B", "FRK", "MATK", "BMX", "EPHA7", "LMTK3", "LMTK2", "AATK", "PTK7", "DDR2", "KDR", "DDR1", "ZAP70", "YES1", "TYRO3", "TYK2", "TXK", "NTRK3", "NTRK2", "NTRK1", "TNK1", "TIE1", "TEK", "TEC", "SYK", "SRC", "RYK", "ROS1", "ROR2", "ROR1", "MST1R", "RET", "PDGFRB", "PDGFRA", "MUSK", "MET", "MERTK", "LYN", "LTK", "LCK", "KIT", "JAK3", "JAK2", "JAK1", "ITK", "INSRR", "INSR", "IGF1R", "ERBB4", "ERBB3", "ERBB2", "HCK", "FYN", "FLT4", "FLT1", "FLT3", "FGR", "FGFR4", "FGFR3", "FGFR2", "FGFR1", "FES", "FER", "PTK2", "EPHB6", "EPHB4", "EPHB3", "EPHB2", "EPHB1", "EPHA8", "EPHA5", "EPHA4", "EPHA3", "EPHA2", "EPHA1", "EGFR", "CSK", "CSF1R", "BTK", "PTK6", "BLK", "AXL", "ABL2", "ALK", "TNK2", "ABL1")

#
res = subset(updated_data, is.element(trimws(toupper(gene_name)), trimws(toupper(TK_gene_list))))
head(res); dim(res); 

write.csv(res, paste0(save_folder, "/", "pickup---TK---genes---all_RSEM_isoforms_TPM", ".csv"))


#########  pickup outliers
save_name = "result--for--all_RSEM_TPM--isoforms"

# all samplles
library(plyr)

#########   top 5% maximum TPM values (95% samples remained) are remove to calculate out the mean values
#########   for example, 100 samples with TMP as 1-100, then the top 5% TMP values are removed (those are 96-100 which are removed); and the mean values are calculated out from the remained 95 samples with TMP values 1-95;
#########   if the mean value is zeros, the mean value is set as 1.
sample_number = 321
num_with_top_TMP = round(0.05*sample_number, 0)


#########    filter conditions
fold_times = 5
TMP_threshold = 10



###################
# data input
all_data = res
rownames(all_data) = all_data[,3]
all_data = all_data[,-c(1:3)]
head(all_data); dim(all_data)

#
fold_matrix = data.frame(mean="", all_data)
#head(fold_matrix); dim(fold_matrix)
fold_matrix$mean = as.character(fold_matrix$mean)
#head(fold_matrix); dim(fold_matrix)
fold_matrix = as.matrix(fold_matrix)

#
difference_matrix = data.frame(mean="", all_data)
#head(difference_matrix); dim(difference_matrix)
difference_matrix$mean = as.character(difference_matrix$mean)
#head(difference_matrix); dim(difference_matrix)
difference_matrix = as.matrix(difference_matrix)

#
i = 1

for (i in 1:nrow(fold_matrix)) {
	
	if (i%%1000==1) print(i);

	tmp = all_data[i, ]
	tmp = as.numeric(sort(tmp))
	mean(tmp);
	mean_value = mean(tmp[(1+num_with_top_TMP):(length(tmp)-num_with_top_TMP)])
	fold_matrix[i,1] = mean_value
	difference_matrix[i,1] = mean_value
	
	# remove top 5% TPM values, if others' mean is not zero
	if ( mean(tmp[(1+num_with_top_TMP):length(tmp)-num_with_top_TMP]) != 0) {
		fold_matrix[i,2:ncol(fold_matrix)] = round(as.numeric(all_data[i, ]/mean_value),2)
		difference_matrix[i,2:ncol(difference_matrix)] = round(as.numeric(all_data[i, ]-mean_value),2)
	} else {
	# remove top 5% TPM values, if others' mean IS ZERO; just keep the TPM
		fold_matrix[i,2:ncol(fold_matrix)] = as.numeric(all_data[i, ])
		difference_matrix[i,2:ncol(difference_matrix)] = round(as.numeric(all_data[i, ]-mean_value),2)
	}
	
}

#
fold_matrix_1 = fold_matrix[, -1]
difference_matrix_1 = difference_matrix[, -1]
#head(fold_matrix_1); dim(fold_matrix_1)
#head(difference_matrix_1); dim(difference_matrix_1)


### method 1: using the fold_matrix 
# extract outlier
nth_rows = nrow(fold_matrix_1); nth_cols = ncol(fold_matrix_1)
#

tmp = which(as.numeric(fold_matrix_1) > fold_times)
res = matrix("", nrow=length(tmp), ncol=7)

i = 2
for (i in 1:length(tmp)) {
	
	if (i%%2000==1) print(i);
	
	row_idx = tmp[i]%%nth_rows; 
	
	if (row_idx==0) {
		row_idx=nth_rows
		col_idx = floor(tmp[i]/nth_rows)
	} else {
		col_idx = floor(tmp[i]/nth_rows) + 1
	}
	
	fold_matrix_1[tmp[i]];
	fold_matrix_1[row_idx, col_idx]
	
	res[i, 1] = colnames(fold_matrix_1)[col_idx]
	res[i, 2] = rownames(fold_matrix_1)[row_idx]
	res[i, 3] = all_data[row_idx, col_idx] # TMP
	res[i, 4] = fold_matrix_1[row_idx, col_idx] # TMP fold times

	#
	res[i, 5] = as.numeric(fold_matrix[row_idx, 1])
	res[i, 6] = row_idx
	res[i, 7] = col_idx

	if (i%%1000==1) {
		write.csv(res, paste0(save_folder, "/", save_name, ".csv"))
		res_1 = res[as.numeric(res[,3])>TMP_threshold,]
		write.csv(res_1, paste0(save_folder, "/filtered_TMP_greater_than_", TMP_threshold, "--", save_name, ".csv"))
	}
}

colnames(res) = c("sample_ID", "transcript_ID", "TMP", "TMP_fold", "avarage_TPM", "row_idx", "col_idx")

#
head(res)

res_1 = res[as.numeric(res[,3])>TMP_threshold,]
res_1 = as.data.frame(res_1)

#
updated_res_1 = data.frame(transcipt_name = unlist(lapply(strsplit(as.character(res_1$transcript_ID), "_"), '[[', 2)), res_1)

head(updated_res_1); dim(updated_res_1); 

#
updated_res_1 = data.frame(gene_name = unlist(lapply(strsplit(as.character(updated_res_1$transcipt_name), "-"), '[[', 1)), updated_res_1)

head(updated_res_1); dim(updated_res_1); 

#
write.csv(res, paste0(save_folder, "/", save_name, ".csv"))

write.csv(updated_res_1, paste0(save_folder, "/filtered_TMP_greater_than_", TMP_threshold, "--", save_name, ".csv"))

write.csv(fold_matrix_1, paste0(save_folder, "/", "fold_change_matrix--", save_name, ".csv"))




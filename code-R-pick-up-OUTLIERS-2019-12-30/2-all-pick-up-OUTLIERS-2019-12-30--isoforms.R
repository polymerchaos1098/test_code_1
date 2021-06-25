# all samplles
library(plyr)

data_folder = "/projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-new-data-10-2019/result-RSEM-isoforms"

file_name = "all_RSEM_isoforms_TPM.txt"

save_folder = "/projectsp/foran/Team_Chan/STAR_FUSION_results/FFPE-new-data-10-2019/result-R-pick-up-OUTLIERS-2019-12-30"

save_name = "result--for--all_RSEM_TPM--isoforms"

###################
# data input
all_data = read.csv(file=paste0(data_folder,"/", file_name), header=TRUE, sep="") 
rownames(all_data) = all_data[,1]
all_data = all_data[,-1]
#head(all_data); dim(all_data)

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

	if (i%%2000==1) print(i);

	tmp = all_data[i, ]
	tmp = as.numeric(sort(tmp))
	mean(tmp);
	mean_value = mean(tmp[6:(length(tmp)-5)])
	fold_matrix[i,1] = mean_value
	difference_matrix[i,1] = mean_value
	
	# remove top 5 TPM values, if others' mean is not zero
	if ( mean(tmp[6:length(tmp)-5]) != 0) {
		fold_matrix[i,2:ncol(fold_matrix)] = round(as.numeric(all_data[i, ]/mean_value),2)
		difference_matrix[i,2:ncol(difference_matrix)] = round(as.numeric(all_data[i, ]-mean_value),2)
	} else {
	# remove top 5 TPM values, if others' mean IS ZERO; just keep the TPM
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
fold_times = 10
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
		res_1 = res[as.numeric(res[,3])>10,]
		write.csv(res_1, paste0(save_folder, "/filtered_TMP_greater_than_10--", save_name, ".csv"))
	}
}

colnames(res) = c("sample_ID", "gene_ID", "TMP", "TMP_fold", "avarage_TPM", "row_idx", "col_idx")

#
head(res)

res_1 = res[as.numeric(res[,3])>10,]

#
write.csv(res, paste0(save_folder, "/", save_name, ".csv"))

write.csv(res_1, paste0(save_folder, "/filtered_TMP_greater_than_10--", save_name, ".csv"))

write.csv(fold_matrix_1, paste0(save_folder, "/", "fold_change_matrix--", save_name, ".csv"))


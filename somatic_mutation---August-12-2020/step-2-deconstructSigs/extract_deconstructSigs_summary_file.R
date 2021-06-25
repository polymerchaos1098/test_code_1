####### data
project_folder = "/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/somatic_mutation---August-12-2020/result-deconstructSigs-pics"

filter_condition = c("filter_test_1")

file_path = paste0(project_folder, "/", filter_condition, "-summary.txt")
	
tmp = read.csv(file=file_path, header=TRUE, sep=",") 

tmp = tmp[,1:3]
colnames(tmp) = c("sample_ID", "sample_name", "mutation_signature")
head(tmp); dim(tmp); colnames(tmp)

#
sig_list = c()

i = 1
for (i in 1:nrow(tmp)) {
	
	need = as.character(tmp$mutation_signature[i])
	need = unlist(strsplit(need, " & ", fixed=TRUE))
	
	for (j in 1:length(need)) {
		a = unlist(strsplit(as.character(need[j]), " : ", fixed=TRUE))[1]
		sig_list = c(sig_list, a)
	}	
}

sig_list = sort(unique(sig_list)); sig_list

#
res = matrix("-", nrow=nrow(tmp), ncol=length(sig_list))
colnames(res) = sig_list
res = data.frame(tmp, res)
res = as.matrix(res)

head(res); dim(res)

####
i = 1
for (i in 1:nrow(res)) {
	
	need = as.character(res[i,3])
	need = unlist(strsplit(need, " & ", fixed=TRUE))
	
	j = 1
	for (j in 1:length(need)) {
		a = unlist(strsplit(as.character(need[j]), " : ", fixed=TRUE))[1]
		b = unlist(strsplit(as.character(need[j]), " : ", fixed=TRUE))[2]
		pos = which(a==colnames(res))
		res[i, pos] = b
	}	
}

#
save_name = paste0(filter_condition, "---mutation-signature-statistics")

write.csv(res, paste0(project_folder, "/", save_name, ".csv"))

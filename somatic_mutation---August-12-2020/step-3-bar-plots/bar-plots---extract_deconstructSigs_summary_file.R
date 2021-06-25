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
res = matrix("0", nrow=nrow(tmp), ncol=length(sig_list))
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
#save_name = paste0(filter_condition, "--separately")

#write.csv(res, paste0(project_folder, "/", save_name, ".csv"))



######################################
######################################
# bar plot   ----  need_data
# transpost to data with each sample as a column
library(RColorBrewer) # take a look at http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_intro.html

sig_data = res[,-c(2,3)]
rownames(sig_data) = sig_data[,1]
sig_data = sig_data[, -1]
sig_data = data.frame(sig_data)
sig_data = t(sig_data)

# update the signature names
rownames(sig_data) = gsub("nature.", " ", rownames(sig_data))
sig_data = data.frame(sig_data)
colnames(sig_data) = gsub("\\.", "-", colnames(sig_data))
colnames(sig_data) = gsub("X", "", colnames(sig_data))

# update the column names
colnames(sig_data) = gsub("IID_H", "", colnames(sig_data))
colnames(sig_data) = gsub("_T02", "", colnames(sig_data))
colnames(sig_data) = gsub("_T03", "", colnames(sig_data))
head(sig_data); dim(sig_data); 

# 
need_data = sig_data

#
head(need_data); dim(need_data); 
rownames(need_data)

# check the singture orders
need_data_1 = need_data[c(6,8:10,1:5,7),]
head(need_data_1); dim(need_data_1); 
rownames(need_data_1)

#
png(paste0(project_folder, "/", "AML-sample-component-of-signature-bar-plot.png"), width = 600, height = 480)
FD.palette <- c("#984EA3","#D1E5F0", "#92C5DE", "#FF7F00","#2166AC", "#B2182B", "#E41A1C", "#F4A582","#FDDBC7", "#053061")
options(scipen=10)
par(mar=c(8, 5, 3, 6), las=2)
barplot(as.matrix(need_data_1), ylab = "Percent",font.lab = 2.5, font.axis = 2, legend=F, beside=F,col=FD.palette, border=FD.palette, space=1)
#legend(1, 1, legend=rownames(need_data_1))
legend("topleft", legend=rownames(need_data_1), fill=FD.palette, inset=c(1,0), xpd=TRUE)
box()
dev.off()

#legend(2.8,-1,c("group A", "group B"), pch = c(1,2), lty = c(1,2), xpd=TRUE)
#legend("topright", inset=c(-0.2,0), legend=c("A","B"), pch=c(1,3), title="Group")


######
# ratio plot
ratio_thyroid = as.matrix(need_data_1); dim(ratio_thyroid)
ratio_thyroid[ratio_thyroid>0] = 1
dim(ratio_thyroid)

need = data.frame(sig=rownames(ratio_thyroid), ratio=-1)
#need = as.matrix(need)
need

for (i in 1: nrow(ratio_thyroid)) {
    need[i, 2] = sum(as.numeric(ratio_thyroid[i,]))/12*100
}
need

#
png(paste0(project_folder, "/", "AML-signature-ratio-box-plot.png"), width = 600, height = 480)

options(scipen=10)
par(mar=c(4, 5, 3, 3), las=2)
#
FD.palette <- c("#984EA3","#D1E5F0", "#92C5DE", "#FF7F00","#2166AC", "#B2182B", "#E41A1C", "#F4A582","#FDDBC7", "#053061")

barplot(need$ratio, ylab = "Percent",font.lab = 2.5, font.axis = 2, names.arg=need$sig, legend=F, beside=F,col=FD.palette, border=FD.palette)
box()

dev.off()










a = data.frame(ratio_thyroid, ratio = rowSums/12)
rowSums(ratio_thyroid)
as.numeric(need_data_1[1,1]) > 0


#
barplot(thyroid_data, main = "total revenue", names.arg = rownames(sig_data), xlab = "month", ylab = "revenue")

barplot(thyroid_data)

#
library(RColorBrewer) # take a look at http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_intro.html

#
png(paste0(project_folder, "/", "thyroid_data-mutation-signature-1---components.png"), width = 1419, height = 736)
# display.brewer.all()
FD.palette <- c("#984EA3","#377EB8","#4DAF4A","#FF7F00","#E41A1C")
options(scipen=10)
par(mar=c(6, 6, 3, 3), las=2)
data4bp <- t(thyroid_data[,c(5,4,2,3,1)])
barplot(data4bp, beside=F,col=FD.palette, border=FD.palette, space=1, legend=F, ylab="Number of People", main="Migration to the United States by Source Region (1820 - 2006)", mgp=c(4.5,1,0) )
legend( "topleft", legend=rev(rownames(data4bp)), fill=rev(FD.palette) )
box()
dev.off()



###

# bar plot
# transpost to data with each sample as a column
sig_data = res[,-c(2,3)]
rownames(sig_data) = sig_data[,1]
sig_data = sig_data[, -1]
sig_data = data.frame(sig_data)

# sig as col; sample as rows
# the format is from this websit for plotting
# http://onertipaday.blogspot.com/2009/01/statistical-visualizations-part-2.html
head(sig_data); dim(sig_data); 
colnames(sig_data) = gsub("nature.", " ", colnames(sig_data))

thyroid_data = sig_data[18:29, ]
oncocymotas_data = sig_data[1:17, ]

#
head(thyroid_data); dim(thyroid_data); 
head(oncocymotas_data); dim(oncocymotas_data); 


#
library(RColorBrewer) # take a look at http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_intro.html

#
png(paste0(project_folder, "/", "thyroid_data-mutation-signature-1---components.png"), width = 1419, height = 736)
# display.brewer.all()
FD.palette <- c("#984EA3","#377EB8","#4DAF4A","#FF7F00","#E41A1C")
options(scipen=10)
par(mar=c(6, 6, 3, 3), las=2)
data4bp <- t(thyroid_data[,c(5,4,2,3,1)])
barplot(data4bp, beside=F,col=FD.palette, border=FD.palette, space=1, legend=F, ylab="Number of People", main="Migration to the United States by Source Region (1820 - 2006)", mgp=c(4.5,1,0) )
legend( "topleft", legend=rev(rownames(data4bp)), fill=rev(FD.palette) )
box()
dev.off()


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

		j = 20
		for (j in 1:length(all_ANN_terms)) {
			one_ANN = all_ANN_terms[j]
			
			# 2nd is Annotation
			a = unlist(strsplit(one_ANN, "\\|"))[2]
			if (a!="") Annot = c(Annot, a)

			# 3rd is Annotation_Impact
			a = unlist(strsplit(one_ANN, "\\|"))[3]
			if (a!="") Annot_Impact = c(Annot_Impact, a)

			# 4th is gene
			a = unlist(strsplit(one_ANN, "\\|"))[4]
			if (a!="") gene = c(gene, a)
			
			# 10th is hgvs_c
			a = unlist(strsplit(one_ANN, "\\|"))[10]
			if (a!="") hgvs_c = c(hgvs_c, a)
			
			# 11th is hgvs_p
			a = unlist(strsplit(one_ANN, "\\|"))[11]
			if (a!="") hgvs_p = c(hgvs_p, a)
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

somatic_data = data.frame(uid=paste0(trimws(as.character(somatic_data$Sample)), "-", trimws(as.character(somatic_data$CHROM)), "-", trimws(as.character(somatic_data$POS)), "-", trimws(as.character(somatic_data$REF)), "-", trimws(as.character(somatic_data$ALT)), sep=""), somatic_data)
head(somatic_data); dim(somatic_data); 

length(unique(somatic_data$uid))

unique(somatic_data$Annot_Impact)

##############
############################################
###########    ----filter out "missense"------
############################################
missense_set = c()
i = 1
for (i in 1:nrow(somatic_data)) {
   tmp = trimws(as.character(somatic_data$Annot[i]))
   pos = str_locate(tmp, "missense_variant")
   if (!is.na(pos[1])) missense_set = rbind(missense_set, somatic_data[i,])
}
head(missense_set); dim(missense_set)

##############
############################################
###########    ----filter out "HIGH"------
############################################
HIGH_set = c()
i = 1
for (i in 1:nrow(somatic_data)) {
   tmp = trimws(as.character(somatic_data$Annot_Impact[i]))
   pos = str_locate(tmp, "HIGH")
   if (!is.na(pos[1])) HIGH_set = rbind(HIGH_set, somatic_data[i,])
}
head(HIGH_set); dim(HIGH_set)

target_set = rbind(missense_set, HIGH_set)
head(target_set); dim(target_set)

# check if duplicated
length(unique(target_set$uid))

#
write.csv(target_set, paste0(save_folder, "AML---missense_variant_and_High---ALL-Candidates", ".csv"), row.names = TRUE)



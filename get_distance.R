library("vcfR")
library("stringr")

#a = read.csv("D:/CINJ/1-16-2018/report/intersection_I-H-113385-N1-1-D1-1.vcf",sep=" ")

vcf_file = "D:/CINJ/1-16-2018/report/intersection_I-H-113385-N1-1-D1-1.vcf"
save_file = "D:/CINJ/1-16-2018/report/intersection_I-H-113385-N1-1-D1-1-get_distance.csv"

vcf <- read.vcfR( vcf_file, verbose = FALSE )


b = vcfR2loci(vcf)


# head(getFIX(vcf,getINFO = TRUE))
# https://cran.r-project.org/web/packages/vcfR/vignettes/vcf_data.html

vcf_content = getFIX(vcf,getINFO = TRUE)

add_distance = data.frame(vcf_content, End = -1, Dinstance = -1)

dim(add_distance)

add_distance = as.matrix(add_distance)

i = 1

for (i in 1:nrow(add_distance)) {

  tmp = add_distance[i,8]
  p = str_locate(tmp, "END=")
  
  if (!is.na(p[1])) {
    tmp = substr(tmp, p[2]+1, nchar(tmp))
    p = str_locate(tmp, ";")
    add_distance[i,9] = substr(tmp, 1, p[2]-1)
	add_distance[i,10] = as.numeric(substr(tmp, 1, p[2]-1))-as.numeric(add_distance[i,2])
  }
  
}

write.csv(add_distance, save_file)

dim(add_distance)

add_distance = as.data.frame(add_distance)
a = subset(add_distance, FILTER=="PASS")
dim(add_distance)






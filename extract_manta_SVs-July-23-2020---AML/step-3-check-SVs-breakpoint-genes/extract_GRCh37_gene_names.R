# gtf website
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/

# using this steps to extract genes
# https://seandavi.github.io/ITR/transcriptdb.html

#
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/

setwd("/projectsp/foran/Team_Chan/Hua_files/AML-get-summary---July-23-2020/gtf_files")

#download.file("https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz", "gencode.v24.primary_assembly.annotation.gtf.gz")

library(GenomicFeatures)
txdb = makeTxDbFromGFF('Homo_sapiens.GRCh37.75.gtf.gz')

library(AnnotationDbi)
saveDb(txdb, 'txdb.Homo_sapiens.GRCh37.75.sqlite')

library(AnnotationDbi)
txdb = loadDb(file = 'txdb.Homo_sapiens.GRCh37.75.sqlite')
a = genes(txdb)
a

#
write.table( x = data.frame(a), file = "Homo_sapiens.GRCh37.75.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )


######################
# best method
# https://biodatascience.github.io/compbio/bioc/ranges.html

library(ensembldb)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("EnsDb.Hsapiens.v75")
# this is GRCh37
library(EnsDb.Hsapiens.v75) # this pkg is about 75 Mb
edb <- EnsDb.Hsapiens.v75
g <- genes(edb)
length(g)

y = data.frame(g)
dim(y)

write.table( y[,-11], file = "EnsDb.Hsapiens.v75.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE )





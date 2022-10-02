library(DEGseq)
args = commandArgs(trailingOnly=TRUE)

intput = args[1]
out = args[2]
c = c(2)
d = c(3)
geneExpMatrix1 <- readGeneExp(file=intput, geneCol=1, valCol=c)
geneExpMatrix2 <- readGeneExp(file=intput, geneCol=1, valCol=d)

DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2), groupLabel1="Turmol",geneExpMatrix2=geneExpMatrix2, geneCol2=1, expCol2=c(2), groupLabel2="Normal",method="MARS", outputDir=out)
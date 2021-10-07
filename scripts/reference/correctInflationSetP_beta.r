#!/users/k1507306/R/bin/Rscript
##usage: file pvalueid lambda output
arg <- commandArgs(TRUE) 
lambda <- as.numeric(arg[2])
output <- arg[3]

gres <- read.table(arg[1],header=T,sep="\t")

gres$P <- pnorm(-abs((qnorm(gres$P/2)/sqrt(lambda))))*2 

gres[,"N"]<-as.integer(floor(gres[,"N"]))

gz1 <- gzfile(paste(output,".gz",sep=""),"w")
write.table(gres,gz1,quote=F,row.names=F,col.names=T,sep="\t")
close(gz1)


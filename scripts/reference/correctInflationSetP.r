#!/users/k1507306/R/bin/Rscript
##usage: file pvalueid lambda output
arg <- commandArgs(TRUE) 
pid <- arg[2]
lambda <- as.numeric(arg[3])
output <- arg[4]

gres <- read.table(arg[1],header=T,sep="\t")
minval <- min(gres$P,na.rm=T)
gres[,pid] <- pchisq((qchisq(1-as.numeric(gres[,pid]),1)/lambda), df = 1, lower = FALSE)
gres[gres$P==0]<-minval
gres[,"N"]<-floor(gres[,"N"])

gz1 <- gzfile(paste(output,".gz",sep=""),"w")
write.table(gres,gz1,quote=F,row.names=F,col.names=T,sep="\t")
close(gz1)


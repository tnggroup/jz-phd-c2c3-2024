#Create refpanel variant list

library(data.table)

nThreads<-4

filepath.bim<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.rs.CM23.qc.bim"
filepath.frq<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.rs.CM23.qc.eur.frq"
filepath.l2<-"../../ld_scores/hc1kgp3.eur.1cM/hc1kgp3.eur.l2.ldscore.gz"
filepath.newvarlist<-"hc1kgp3.eur.1cM.l2.jz2022.gz"

bim<-fread(file = filepath.bim, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim)<-c("CHR","SNP","CM","BP","A1","A2")
bim$ORDER<-1:nrow(bim)
setDT(bim)

maf <- fread(file = filepath.frq, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
maf$ORDER<-1:nrow(maf)
setDT(maf)

l2 <- fread(file = filepath.l2, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)



#data formatting
#bim$CHR<-shru::parseCHRColumn(bim$CHR) #not necessary after earlier formatting
bim$CHR<-as.integer(bim$CHR)
bim$BP<-as.integer(bim$BP)
bim[is.na(SNP),SNP:="."]
setkeyv(bim, cols = c("ORDER","CHR","BP","A2","A1","SNP"))

maf$CHR<-as.integer(maf$CHR)
maf[is.na(SNP),SNP:="."]
setkeyv(maf, cols = c("ORDER","CHR","A2","A1","SNP"))

l2$CHR<-as.integer(l2$CHR)
l2$BP<-as.integer(l2$BP)
l2$L2<-as.numeric(l2$L2)
l2[is.na(SNP),SNP:="."]
setkeyv(l2, cols = c("CHR","BP","SNP"))


#aggregate(x = testDF, by = list(testDF$by1, testDF$by2), FUN = head, 1) #get unique entries based on the chosen grouping vars and order (using parameter 1 for head)

#merge
bim[maf, on=c(ORDER="ORDER",CHR="CHR",SNP="SNP"), c("MAF","NCHROBS"):=list(i.MAF,i.NCHROBS)]

#get unique L2 sorted by absolute L2, decreasing
l2[,AL2:=abs(L2)]
l2<-l2[order(-AL2),]
#ul2<-aggregate(x = l2, by = list(l2$CHR,l2$SNP,l2$BP), FUN = head, 1)
ul2<-l2[, .(L2 = head(L2,1)), by = c("CHR","SNP","BP")]
bim[ul2, on=c(CHR="CHR",SNP="SNP",BP="BP"), c("L2.EUR"):=list(i.L2)]

#missing values
bim[is.na(L2.EUR),L2.EUR:=0]

bim[grepl(pattern = ".",x = bim$SNP,fixed = T),SNP:=paste0("somevar_",CHR,":",BP,":",ORDER)]
#sum(grepl(pattern = ".",x = bim$SNP,fixed = T)) #should be 0
fwrite(x = bim[order(ORDER),c("SNP","A1","A2","CHR","BP","CM","MAF","L2.EUR")], file = filepath.newvarlist,col.names = T, sep = "\t",nThread = nThreads, na = ".",quote = F)

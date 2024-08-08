#Create refpanel variant list - GLAD+

#revision history
#version 1 - GLAD+ version, adapted from the HC1kG variant list scripts

library(data.table)
library(shru)

nThreads<-6
setDTthreads(threads = nThreads)
getDTthreads(verbose=T)

filepath.bim<-"GLAD_EDGI_NBR.keep.CM23.bim"
filepath.frq<-"GLAD_EDGI_NBR.keep.frq"
folderpath.l2.mix<-"/scratch/prj/gwas_sumstats/ld_scores/w_ld.GLAD_EDGI_NBR.keep.b38.gcta"
filepath.newvarlist.mix<-"GLAD_EDGI_NBR.mix.l2.jz2024.gz"

allvars<-fread(file = filepath.bim, na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
setDT(allvars)
colnames(allvars)<-c("CHR","SNP","CM","BP","A1","A2")
allvars[is.na(SNP),SNP:="."] #[is.na(SNPR),SNPR:="."]
setkeyv(allvars, cols = c("CHR","BP","A2","A1","SNP"))

#rename MAF to MAF.MIX as this is the overall 1kG MAF, same for NCHROBS
allvars[,MAF.MIX:=MAF][,MAF:=NULL]
allvars[,NCHROBS.MIX:=NCHROBS][,NCHROBS:=NULL]

bim<-fread(file = filepath.bim, na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim)<-c("CHR","SNP","CM","BP","A1","A2")
#bim$ORDER<-1:nrow(bim)
setDT(bim)
bim[is.na(SNP),SNP:="."]
setkeyv(bim, cols = c("CHR","BP","A2","A1","SNP"))

l2.mix <- do.call("rbind", lapply(list.files(path = folderpath.l2.mix, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
  suppressMessages(fread(file = file.path(folderpath.l2.mix, i), na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
}))
l2.mix[is.na(SNP),SNP:="."]
setkeyv(l2.mix, cols = c("CHR","BP","SNP"))

allvars[l2.mix, on=c(CHR="CHR",BP="BP",SNP="SNP"), c("L2.MIX","MAF.MIX","SNP_NUM","MEAN_RSQ","MAX_RSQ"):=list(i.L2,i.MAF,i.SNP_NUM,i.MEAN_RSQ,i.MAX_RSQ)]
rm(l2.mix)

#missing variant ID check and fix
cat("\nN total before merge:",nrow(allvars))
cat("\nN missing SNP-values:",nrow(allvars[grepl(pattern = ".",x = allvars$SNP,fixed = T),]))
cat("\nN missing BP-values:",nrow(allvars[!is.finite(BP),]))
cat("\nN missing A1-values:",nrow(allvars[is.na(A1),]))
cat("\nN missing A2-values:",nrow(allvars[is.na(A2),]))
cat("\nN missing MAF.MIX-values:",nrow(allvars[!is.finite(MAF.MIX),]))
cat("\nN missing L2.MIX-values:",nrow(allvars[!is.finite(L2.MIX),]))
#allvars[grepl(pattern = ".",x = bim$SNP,fixed = T),SNP:=paste0("somevar_",CHR,":",BP,":",ORDER)]

#merge allele variant duplicates across SNP/rsID - NEW - a bit risky to do without strict alignment to dbSNP, but let's do it anyway
#allvars[CHR==1 & BP>(1472989-100) & BP<(1472989+100),]
setorder(allvars,-FW,-MAF.MIX,CHR,BP,A1,A2,SNP) # using MAF.MIX for determining the head SNP (highest MAF)
allvars.merged<-allvars[, .(CHR=head(CHR,1), CM = head(CM,1), BP=head(BP,1), A1=head(A1,1), A2=head(A2,1), MAF.MIX = sum(MAF.MIX), MAF.MIX.ORIG = head(MAF.MIX,1), L2.MIX=head(L2.MIX,1), SNP_NUM=head(SNP_NUM,1), MEAN_RSQ=head(MEAN_RSQ,1), MAX_RSQ=head(MAX_RSQ,1)), by = c("SNP")]

#stats showed that the above was unnecessary
cat("\n\nN total after merge:",nrow(allvars.merged))
cat("\nN missing SNP-values:",nrow(allvars.merged[grepl(pattern = ".",x = allvars.merged$SNP,fixed = T),]))
cat("\nN missing BP-values:",nrow(allvars.merged[!is.finite(BP),]))
cat("\nN missing A1-values:",nrow(allvars.merged[is.na(A1),]))
cat("\nN missing A2-values:",nrow(allvars.merged[is.na(A2),]))
cat("\nN missing MAF.MIX-values:",nrow(allvars.merged[!is.finite(MAF.MIX),]))
cat("\nN missing L2.MIX-values:",nrow(allvars.merged[!is.finite(L2.MIX),]))

#remove statistics vars
allvars.merged[,MAF.MIX.ORIG:=NULL]

setorder(allvars.merged,CHR,BP,CM,SNP,-MAF.MIX,-L2.MIX) #new for the sorted version
fwrite(x = allvars.merged[,.(SNP,CHR,BP,A1,A2,CM,MAF=MAF.MIX,L2=L2.MIX,MAF.MIX,L2.MIX)], file = filepath.newvarlist.mix,col.names = T, sep = "\t",nThread = nThreads, na = ".",quote = F)
cat("\nN final variants in variant list:",nrow(allvars.merged), "written to",filepath.newvarlist.mix)

cat("\nTHE END")

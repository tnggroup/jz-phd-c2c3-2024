#Create refpanel variant list

library(data.table)
library(shru)

nThreads<-6
setDTthreads(nThreads)

filepath.varlist<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.varlist.gz"
filepath.bim<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM23.bim"
#filepath.l2.mix<-"../../ld_scores/hc1kgp3.b38.mix.l2.jz2023/hc1kgp3.b38.mix.jz2023.1cm.250blocks.l2.ldscore.gz" #not available because of computation limits
filepath.bim.eur<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.bim"
filepath.frq.eur<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.frq"
filepath.l2.eur<-"../../ld_scores/hc1kgp3.b38.eur.l2.jz2023/hc1kgp3.b38.eur.jz2023.1cm.250blocks.l2.ldscore.gz"
filepath.newvarlist<-"hc1kgp3.b38.eur.l2.jz2023.gz"

allvars<-fread(file = filepath.varlist, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = T, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
setDT(allvars)
allvars[is.na(SNP),SNP:="."] #[is.na(SNPR),SNPR:="."]
setkeyv(allvars, cols = c("BP","A2","A1","SNP","SNPR")) #I forgot to export the CHR column in the first version of the full variant file from the edit_reference_panel.R routine. This will be added in from the BIM-file information instead.

bim<-fread(file = filepath.bim, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim)<-c("CHR","SNP","CM","BP","A1","A2")
#bim$ORDER<-1:nrow(bim)
setDT(bim)
bim[is.na(SNP),SNP:="."]
setkeyv(bim, cols = c("CHR","BP","A2","A1","SNP"))

allvars[bim, on=c(SNP="SNP",BP="BP",A1="A1",A2="A2"), c("CHR","CM"):=list(i.CHR,i.CM)]
rm(bim)
setkeyv(allvars, cols = c("CHR","BP","A2","A1","SNP","SNPR"))

# l2.mix <- fread(file = filepath.l2.mix, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)
# l2.mix[is.na(SNP),SNP:="."]
# setkeyv(l2.mix, cols = c("CHR","BP","SNP"))
# 
# allvars[l2.mix, on=c(CHR="CHR",BP="BP",SNP="SNP"), c("L2"):=list(i.L2)]
# rm(l2.mix)

bim.eur<-fread(file = filepath.bim.eur, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim.eur)<-c("CHR","SNP","CM","BP","A1","A2")
bim.eur$ORDER<-1:nrow(bim.eur)
setDT(bim.eur)
bim.eur[is.na(SNP),SNP:="."]
setkeyv(bim.eur, cols = c("ORDER","CHR","BP","A2","A1","SNP"))

maf.eur <- fread(file = filepath.frq.eur, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
maf.eur$ORDER<-1:nrow(maf.eur)
setDT(maf.eur)
maf.eur[is.na(SNP),SNP:="."]
setkeyv(maf.eur, cols = c("ORDER","CHR","A2","A1","SNP"))

bim.eur[maf.eur, on=c(ORDER="ORDER",CHR="CHR",SNP="SNP"), c("MAF","NCHROBS"):=list(i.MAF,i.NCHROBS)]
rm(maf.eur)

allvars[bim.eur, on=c(CHR="CHR",BP="BP",A1="A1",A2="A2",SNP="SNP"), c("MAF.EUR","NCHROBS.EUR"):=list(i.MAF,i.NCHROBS)]
rm(bim.eur)

l2.eur <- fread(file = filepath.l2.eur, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)
l2.eur[is.na(SNP),SNP:="."]
setkeyv(l2.eur, cols = c("CHR","BP","SNP"))

allvars[l2.eur, on=c(CHR="CHR",BP="BP",SNP="SNP"), c("L2.EUR"):=list(i.L2)]
rm(l2.eur)


#missing variant ID check and fix
cat("\nN missing SNP-values:",nrow(allvars[grepl(pattern = ".",x = allvars$SNP,fixed = T),]))
cat("\nN missing RSNP-values:",nrow(allvars[grepl(pattern = ".",x = allvars$SNPR,fixed = T),]))
cat("\nN missing BP-values:",nrow(allvars[!is.finite(BP),]))
cat("\nN missing A1-values:",nrow(allvars[is.na(A1),]))
cat("\nN missing A2-values:",nrow(allvars[is.na(A2),]))
cat("\nN missing MAF-values:",nrow(allvars[!is.finite(MAF),]))
cat("\nN missing MAF.EUR-values:",nrow(allvars[!is.finite(MAF.EUR),]))
cat("\nN missing L2.EUR-values:",nrow(allvars[!is.finite(L2.EUR),]))
#allvars[grepl(pattern = ".",x = bim$SNP,fixed = T),SNP:=paste0("somevar_",CHR,":",BP,":",ORDER)]

#duplicate variants - should not be any duplicates now after the new regference panel routine
# cat("\nN variants before removing duplicates:",nrow(allvars))
# bim<-bim[order(-MAF),]
# bim<-bim[, .(CM = head(CM,1),MAF = head(MAF,1)), by = c("CHR","SNP","BP","A1","A2")]
# cat("\nN variants after removing duplicates:",nrow(bim))

#missing L2 - what should be done with these?
#allvars[is.na(L2.EUR),L2.EUR:=0]

#sum(grepl(pattern = ".",x = bim$SNP,fixed = T)) #should be 0
setkeyv(allvars, cols = c("SNP"))
fwrite(x = allvars, file = filepath.newvarlist,col.names = T, sep = "\t",nThread = nThreads, na = ".",quote = F)

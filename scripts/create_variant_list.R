#Create refpanel variant list

#revision history
#version 3 - the mix ld-scores are available now?
#version 2 - adapted for the jz2024 version. the mix ld-scores are not available yet because they take > 2 days to run on CREATE.

library(data.table)
library(shru)

nThreads<-6
setDTthreads(threads = nThreads)
getDTthreads(verbose=T)

filepath.varlist<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.varlist.gz"
filepath.varlist.unfiltered<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.varlist.unfiltered.gz" #for creating the synonym list
filepath.bim<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM23.bim"

filepath.bim.eur<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.bim"
filepath.frq.eur<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.frq"
folderpath.l2.eur<-"/scratch/prj/gwas_sumstats/ld_scores/hc1kgp3.b38.eur.l2.jz2024.chr"
folderpath.l2.mix<-"/scratch/prj/gwas_sumstats/ld_scores/hc1kgp3.b38.mix.l2.jz2024.chr"
filepath.newvarlist.mix<-"hc1kgp3.b38.mix.l2.jz2024.gz"
filepath.newvarlist.eur<-"hc1kgp3.b38.eur.l2.jz2024.gz"
filepath.newvarlist.synonyms<-"hc1kgp3.b38.jz2024.synonyms.gz"

allvars<-fread(file = filepath.varlist, na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8", header = T, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
setDT(allvars)
allvars[is.na(SNP),SNP:="."] #[is.na(SNPR),SNPR:="."]
allvars[,CHR:=as.integer(CHR)]
setkeyv(allvars, cols = c("CHR","BP","A2","A1","SNP","SNPF","SNPR"))

#compare with existing sumstat
cSumstats <- shru::readFile(filePath = "/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/cleaned/ANXI03.gz", nThreads = nThreads)
cat("\ncSumstats (ANXI03) F+R SNPs covered: ",(sum(cSumstats$SNP %in% allvars$SNPF) + sum(cSumstats$SNP %in% allvars$SNPR))," of a total ",nrow(cSumstats),"\n")
cat("\ncSumstats (ANXI03) 'SNP' SNPs covered: ",(sum(cSumstats$SNP %in% allvars$SNP))," of a total ",nrow(cSumstats),"\n")

#rename MAF to MAF.MIX as this is the overall 1kG MAF, same for NCHROBS
allvars[,MAF.MIX:=MAF][,MAF:=NULL]
allvars[,NCHROBS.MIX:=NCHROBS][,NCHROBS:=NULL]

bim<-fread(file = filepath.bim, na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim)<-c("CHR","SNP","CM","BP","A1","A2")
#bim$ORDER<-1:nrow(bim)
setDT(bim)
bim[is.na(SNP),SNP:="."]
setkeyv(bim, cols = c("CHR","BP","A2","A1","SNP"))

allvars[bim, on=c(SNP="SNP",CHR="CHR",BP="BP",A1="A1",A2="A2"), c("CM"):=list(i.CM)] #maybe CM should have been exported from the beginning with the varlist?
rm(bim)
setkeyv(allvars, cols = c("CHR","BP","A2","A1","SNP","SNPF","SNPR"))

# l2.mix <- fread(file = filepath.l2.mix, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)
l2.mix <- do.call("rbind", lapply(list.files(path = folderpath.l2.mix, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
  suppressMessages(fread(file = file.path(folderpath.l2.mix, i), na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
}))
l2.mix[is.na(SNP),SNP:="."]
setkeyv(l2.mix, cols = c("CHR","BP","SNP"))

allvars[l2.mix, on=c(CHR="CHR",BP="BP",SNP="SNP"), c("L2.MIX"):=list(i.L2)]
rm(l2.mix)

bim.eur<-fread(file = filepath.bim.eur, na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim.eur)<-c("CHR","SNP","CM","BP","A1","A2")
bim.eur$ORDER<-1:nrow(bim.eur)
setDT(bim.eur)
bim.eur[is.na(SNP),SNP:="."]
setkeyv(bim.eur, cols = c("ORDER","CHR","BP","A2","A1","SNP"))

maf.eur <- fread(file = filepath.frq.eur, na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
maf.eur$ORDER<-1:nrow(maf.eur)
setDT(maf.eur)
maf.eur[is.na(SNP),SNP:="."]
setkeyv(maf.eur, cols = c("ORDER","CHR","A2","A1","SNP"))

bim.eur[maf.eur, on=c(ORDER="ORDER",CHR="CHR",SNP="SNP"), c("MAF","NCHROBS"):=list(i.MAF,i.NCHROBS)]
rm(maf.eur)

allvars[bim.eur, on=c(CHR="CHR",BP="BP",A1="A1",A2="A2",SNP="SNP"), c("MAF.EUR","NCHROBS.EUR"):=list(i.MAF,i.NCHROBS)]
rm(bim.eur)

l2.eur <- do.call("rbind", lapply(list.files(path = folderpath.l2.eur, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
  suppressMessages(fread(file = file.path(folderpath.l2.eur, i), na.strings =c(".",NA,"NA","<NA>",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = nThreads, data.table=T))
}))
l2.eur[is.na(SNP),SNP:="."]
setkeyv(l2.eur, cols = c("CHR","BP","SNP"))

allvars[l2.eur, on=c(CHR="CHR",BP="BP",SNP="SNP"), c("L2.EUR"):=list(i.L2)]
rm(l2.eur)


#missing variant ID check and fix
cat("\nN total before merge:",nrow(allvars))
cat("\nN missing SNP-values:",nrow(allvars[grepl(pattern = ".",x = allvars$SNP,fixed = T),]))
cat("\nN missing RSNP-values:",nrow(allvars[grepl(pattern = ".",x = allvars$SNPR,fixed = T),]))
cat("\nN missing BP-values:",nrow(allvars[!is.finite(BP),]))
cat("\nN missing A1-values:",nrow(allvars[is.na(A1),]))
cat("\nN missing A2-values:",nrow(allvars[is.na(A2),]))
cat("\nN missing MAF.MIX-values:",nrow(allvars[!is.finite(MAF.MIX),]))
cat("\nN missing MAF.EUR-values:",nrow(allvars[!is.finite(MAF.EUR),]))
cat("\nN missing L2.MIX-values:",nrow(allvars[!is.finite(L2.MIX),]))
cat("\nN missing L2.EUR-values:",nrow(allvars[!is.finite(L2.EUR),]))
cat("\nN variants with inverse mapping:",nrow(allvars[!is.na(SNPR),]))
cat("\nN inversly mapped variants (because forward mapping missing):",nrow(allvars[FW==0,]))
#allvars[grepl(pattern = ".",x = bim$SNP,fixed = T),SNP:=paste0("somevar_",CHR,":",BP,":",ORDER)]

#set reverse coded MAF - NEW
allvars[,A1.OLD:=A1][,A2.OLD:=A2]
allvars[FW==0,A1:=A2.OLD]
allvars[FW==0,A2:=A1.OLD]
allvars[,SNPF.OLD:=SNPF][,SNPR.OLD:=SNPR]
allvars[FW==0,SNPF:=SNPR.OLD]
allvars[FW==0,SNPR:=SNPF.OLD]
allvars[,A1.OLD:=NULL][,A2.OLD:=NULL][,SNPF.OLD:=NULL][,SNPR.OLD:=NULL]
allvars[FW==0,MAF.MIX:=(1-MAF.MIX)]
allvars[FW==0,MAF.EUR:=(1-MAF.EUR)]

#merge allele variant duplicates across SNP/rsID - NEW
#allvars[CHR==1 & BP>(1472989-100) & BP<(1472989+100),]
setorder(allvars,-FW,-MAF.MIX,CHR,BP,A1,A2,SNP) # using MAF.MIX for determining the head SNP (highest MAF)
allvars.merged<-allvars[, .(CHR=head(CHR,1), BP=head(BP,1), A1=head(A1,1), A2=head(A2,1), CM=head(CM,1), FW=head(FW,1), SNPF=head(SNPF,1), SNPR=head(SNPR,1), SYNF=head(SYNF,1), SYNR=head(SYNR,1), MAF.MIX = sum(MAF.MIX), MAF.MIX.ORIG = head(MAF.MIX,1), NCHROBS.MIX=head(NCHROBS.MIX,1), MAF.EUR = sum(MAF.EUR), MAF.EUR.ORIG = head(MAF.EUR,1), NCHROBS.EUR = head(NCHROBS.EUR,1), L2.MIX=head(L2.MIX,1), L2.EUR=head(L2.EUR,1)), by = c("SNP")]

cat("\n\nN total after merge:",nrow(allvars.merged))
cat("\nN missing SNP-values:",nrow(allvars.merged[grepl(pattern = ".",x = allvars.merged$SNP,fixed = T),]))
cat("\nN missing RSNP-values:",nrow(allvars.merged[grepl(pattern = ".",x = allvars.merged$SNPR,fixed = T),]))
cat("\nN missing BP-values:",nrow(allvars.merged[!is.finite(BP),]))
cat("\nN missing A1-values:",nrow(allvars.merged[is.na(A1),]))
cat("\nN missing A2-values:",nrow(allvars.merged[is.na(A2),]))
cat("\nN missing MAF.MIX-values:",nrow(allvars.merged[!is.finite(MAF.MIX),]))
cat("\nN missing MAF.EUR-values:",nrow(allvars.merged[!is.finite(MAF.EUR),]))
cat("\nN missing L2.MIX-values:",nrow(allvars.merged[!is.finite(L2.MIX),]))
cat("\nN missing L2.EUR-values:",nrow(allvars.merged[!is.finite(L2.EUR),]))
cat("\nN variants with inverse mapping:",nrow(allvars.merged[!is.na(SNPR),]))
cat("\nN inversly mapped variants (because forward mapping missing):",nrow(allvars.merged[FW==0,]))

#remove statistics vars
allvars.merged[,MAF.MIX.ORIG:=NULL][,MAF.EUR.ORIG:=NULL]
maf.max <- sqrt(max(allvars$MAF.MIX,na.rm=T))
allvars.merged[MAF.MIX>eval(maf.max),MAF.MIX:=eval(maf.max)]

cat("\ncSumstats (ANXI03) F+R SNPs covered: ",(sum(cSumstats$SNP %in% allvars.merged$SNPF) + sum(cSumstats$SNP %in% allvars.merged$SNPR))," of a total ",nrow(cSumstats),"\n")
cat("\ncSumstats (ANXI03) 'SNP' SNPs covered: ",(sum(cSumstats$SNP %in% allvars.merged$SNP))," of a total ",nrow(cSumstats),"\n")


cat("\nFiltering synonyms")
#filter synonyms
all.syn<-fread(file = filepath.varlist.unfiltered, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)
all.syn[is.na(SNP),SNP:="."] #[is.na(SNPR),SNPR:="."]
all.syn[,CHR:=as.integer(CHR)]
setkeyv(all.syn, cols = c("CHR","BP","A2","A1","SNP"))

#align synonym variants too to forward/backward direction matched
all.syn[,MAF:=NULL][,NCHROBS:=NULL][,CM:=NULL]
all.syn[,FW:=forward][,forward:=NULL]
all.syn[,A1.OLD:=A1][,A2.OLD:=A2]
all.syn[FW==0,A1:=A2.OLD]
all.syn[FW==0,A2:=A1.OLD]
all.syn[,A1.OLD:=NULL][,A2.OLD:=NULL]

#assign synonyms
all.syn[,SYN:=SNP][,SNP:=NULL]
setkeyv(all.syn, cols = c("CHR","BP","A2","A1","SYN"))
setkeyv(allvars, cols = c("CHR","BP","A2","A1","SNP"))
cat("\nFiltering synonyms - assigning to primary matches")
#merge on allvars rather than allvars.merged because we consider all allele versions
all.syn.merged<-all.syn[allvars, on=c(CHR="CHR",BP="BP",A2="A2",A1="A1"), nomatch=0] #the alleles are harmonised since previously with the reference panel - no need to do the reverse matching
#all.syn.merged.reverse<-all.syn[allvars, on=c(CHR="CHR",BP="BP",A2="A1",A1="A2"), nomatch=0]
#all.syn.merged<-rbind(all.syn.merged.forward,all.syn.merged.reverse)
# rm(all.syn.merged.forward)
# rm(all.syn.merged.reverse)

all.syn.merged<-all.syn.merged[,.(SNP,SYN,CHR,BP,A1,A2,SNPF,SYNF,SNPR,SYNR,dbSNPBuildID)]
setorder(all.syn.merged,SNP,SYN)
all.syn.merged<-all.syn.merged[SNP!=SYN,]

#exclude synonyms which are primary variants
cat("\nFiltering synonyms - non-primary only")
setkeyv(all.syn.merged, cols = c("SYN","CHR","BP","A2","A1"))
setkeyv(allvars, cols = c("SNP","CHR","BP","A2","A1"))
all.syn.merged[allvars, on=c(SYN="SNP"), c("PRIMARY"):=1]
all.syn.merged<-all.syn.merged[is.na(PRIMARY),]

#unique synonyms
all.syn.merged.unique<-all.syn.merged[, .(SNP = head(SNP,1)), by = c("SYN")]
all.syn.merged.unique<-all.syn.merged.unique[,.(SNP,SYN)]
#missing L2 - what should be done with these?
#allvars[is.na(L2.EUR),L2.EUR:=0]

cat("\nTotal number of unique synonyms: ",nrow(all.syn.merged.unique),"\n")
cat("\ncSumstats (ANXI03) SNPs in synonyms: ",sum(cSumstats$SNP %in% all.syn.merged.unique$SYN),"\n")

cat("\nWriting files")
setkeyv(all.syn.merged.unique,c("SNP","SYN"))
fwrite(x = all.syn.merged.unique, file = filepath.newvarlist.synonyms,col.names = T, sep = "\t",nThread = nThreads, na = ".",quote = F)
#sum(grepl(pattern = ".",x = bim$SNP,fixed = T)) #should be 0
#allvars.merged<-allvars.merged[,c("SNP","BID","SNPR","BIDR","CHR","BP","A1","A2","CM","MAF.MIX","NCHROBS.MIX","L2.MIX","MAF.EUR","NCHROBS.EUR","L2.EUR")]

setorder(allvars.merged,CHR,BP,CM,SNP,-NCHROBS.MIX,-MAF.MIX,-L2.MIX) #new for the sorted version
fwrite(x = allvars.merged, file = filepath.newvarlist.mix,col.names = T, sep = "\t",nThread = nThreads, na = ".",quote = F)
cat("\nN final variants in variant list:",nrow(allvars.merged), "written to",filepath.newvarlist.mix)

setorder(allvars.merged,CHR,BP,CM,SNP,-NCHROBS.EUR,-MAF.EUR,-L2.EUR)
fwrite(x = allvars.merged[,.(SNP,CHR,BP,A1,A2,CM,FW,SNPR,MAF=MAF.EUR,NCHROBS=NCHROBS.EUR,L2=L2.EUR)], file = filepath.newvarlist.eur,col.names = T, sep = "\t",nThread = nThreads, na = ".",quote = F)
cat("\nN final variants in variant list:",nrow(allvars.merged), "written to",filepath.newvarlist.eur)
cat("\nTHE END")

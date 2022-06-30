#resolve duplicate and unnamed variants in a reference panel .bim

library(data.table)

nThreads<-4

filepath.dupvar<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.dup.dupvar"
filepath.newnames<-"newnames.txt"
filepath.maf<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.dup.frq"
filepath.bim<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.bim_old"
filepath.newbim<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.bim"




#read files
dupvar<-fread(file = filepath.dupvar, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)

newnames<-fread(file = filepath.newnames, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)
colnames(newnames)<-c("oldname","newname")

maf<-fread(file = filepath.maf, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = nThreads, showProgress = T)

#important that the bim file keeps the variant order
bim<-fread(file = filepath.bim, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim)<-c("CHR","SNP","CM","BP","A1","A2")
bim$ORDER<-1:nrow(bim)
setDT(bim)


#data formatting and parsing
dupvar$ID1<-unlist(lapply(strsplit(dupvar$IDS,split=" ",fixed=T),FUN=function(x){x[[1]]}))
dupvar$ID2<-unlist(lapply(strsplit(dupvar$IDS,split=" ",fixed=T),FUN=function(x){x[[2]]}))
dupvar$CHR<-as.integer(dupvar$CHR)
dupvar$POS<-as.integer(dupvar$POS)
setkeyv(dupvar, cols = c("CHR","POS","ID1","ID2"))

setkeyv(newnames, cols = c("oldname"))

bim$CHR<-shru::parseCHRColumn(bim$CHR)
bim$CHR<-as.integer(bim$CHR)
bim$BP<-as.integer(bim$BP)
bim[is.na(SNP),SNP:="."]
setkeyv(bim, cols = c("CHR","BP","A2","A1","SNP"))

maf$CHR<-as.integer(maf$CHR)
setkeyv(maf, cols = c("CHR","SNP","A1","A2"))

bim[newnames, on=c(SNP="oldname"), c("NNAME") :=list(i.newname)]

bim[maf, on=c(CHR="CHR",SNP="SNP",A1="A1",A2="A2"), c("MAF","NCHROBS"):=list(i.MAF,i.NCHROBS)]

#check number of rs-no in NNAME
sum(grepl(pattern = "^rs.+", x = bim$NNAME))


#reserve names of duplicates
bim[dupvar,on=c(CHR="CHR",BP="POS",SNP="ID1"),NNAME2:=paste0("unknown",ORDER,"_1")]
bim[dupvar,on=c(CHR="CHR",BP="POS",SNP="ID2"),NNAME2:=paste0("unknown",ORDER,"_2")]

#reserve names of NA/. named SNPs
bim[SNP=="." & is.na(NNAME2),NNAME2:=paste0("unknown",ORDER)]

#finalise new names
bim[is.na(NNAME),NNAME:=NNAME2]


#set new names for variants with new names
bim[!is.na(NNAME),SNP:=NNAME]

#resolve name duplicates, because the rs-numbers are sometimes assigned multiple times
#VERY SLOW - DISCONTINUED
# nVar<-nrow(bim)
# #nVar<-1000000 #Test
# env <- new.env(hash = T)
# for(iVar in 1:nVar){
#   #iVar<-1
#   if(is.null(env[[bim[iVar,]$SNP]]) & !is.na(bim[iVar,]$MAF)){
#     assign(bim[iVar,]$SNP, bim[iVar,]$MAF, envir = env)
#   } else if(!is.na(bim[iVar,]$MAF)) {
#     if(bim[iVar,]$MAF > env[[bim[iVar,]$SNP]]) assign(bim[iVar,]$SNP, bim[iVar,]$MAF, envir = env)
#   }
#   #assign(bim[iVar,]$SNP, bim[iVar,]$MAF, envir = env)
# }

sum(grepl(pattern = "^rs.+", x = bim$SNP))



#output .bim format again in the right order
fwrite(x = bim[order(ORDER),c("CHR","SNP","CM","BP","A1","A2")], file = filepath.newbim,col.names = F, sep = "\t",nThread = nThreads, na = ".",quote = F)



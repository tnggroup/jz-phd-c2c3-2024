#set rsids of the reference panel .bim file from dbSNP panel. output two column filr with old and new variant names/id's to be used by plink
#run in the folder of the working reference panel

library(data.table)
library(shru)

nThreads<-5

filepath.dbSNP <- file.path("..","dbsnp.human_9606_b151_GRCh38p7","00-All.vcf") #using the unzipped version of ~131GB
#filepath.dbSNP <- file.path("..","small.vcf.sample.vcf.gz") #for testing

filepath.refpanelBim <- file.path("1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.bim")


increment<-10000000
nrows <- 56


dbSNP.dt.top<-fread(file = filepath.dbSNP, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, data.table = T,showProgress = T, nThread = nThreads, nrows = 60, skip=nrows, sep="\t")
cat("\nRead top dbSNP\n")
print(length(colnames(dbSNP.dt.top)))
print(colnames(dbSNP.dt.top))




refpanelBim.dt<-fread(file = filepath.refpanelBim, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, data.table = T,showProgress = T, nThread=nThreads, header = F)
refpanelBim.dt$V1<-parseCHRColumn(refpanelBim.dt$V1)
refpanelBim.dt$V1<-as.integer(refpanelBim.dt$V1)
setkeyv(refpanelBim.dt,cols = c("V1","V4","V5","V6"))
cat("\nRead and indexed bim\n")
cat("Overview of bim data table:\n\n")
print(refpanelBim.dt)


#loop to read the dbSNP in chunks because big
for(iIteration in 1:10000){
  #iIteration<-1
  cat("\n\n***Iteration",iIteration,"***\n")
  
  cat("\nReading dbSNP\n")
  dbSNP.dt<-fread(file = filepath.dbSNP, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, data.table = T,showProgress = T, nThread=nThreads, nrows = increment, skip = nrows, sep="\t")
  cat("\nRead dbSNP\n")
  
  if(nrow(dbSNP.dt)<1 | length(colnames(dbSNP.dt))<5){
    cat("Breaking because of empty dt!\n")
    break
  }
  
  colnames(dbSNP.dt) <- colnames(dbSNP.dt.top)
  
  dbSNP.dt$`#CHROM`<-parseCHRColumn(dbSNP.dt$`#CHROM`) #this may potentially be slow
  dbSNP.dt$`#CHROM`<-as.integer(dbSNP.dt$`#CHROM`)
  setkeyv(dbSNP.dt,cols = c("#CHROM","POS","REF","ALT"))
  cat("Indexed dbSNP\n")
  
  cat("Overview of dbSNP data table:\n\n")
  print(dbSNP.dt[,1:8])
  
  cat("\n\nUpdating rsid's\n")
  refpanelBim.dt[dbSNP.dt,on=c(V1="#CHROM",V4="POS",V5="REF",V6="ALT"),c("NNAME","CHR_REF","BP_REF") :=list(i.ID,`i.#CHROM`,i.POS)]
  cat("\nUpdating rsid's\n")
  refpanelBim.dt[dbSNP.dt,on=c(V1="#CHROM",V4="POS",V5="ALT",V6="REF"),c("NNAME","CHR_REF","BP_REF") :=list(i.ID,`i.#CHROM`,i.POS)]
  cat("\nUpdating rsid's\n")
  #refpanelBim.dt[dbSNP.dt,on=c(is.null(NNAME),V1="#CHROM",V4="POS",V6="REF"),c("NNAME","CHR_REF","BP_REF") :=list(i.ID,`i.#CHROM`,i.POS)]
  refpanelBim.dt[dbSNP.dt[refpanelBim.dt[is.na(NNAME),],on=c(`#CHROM`='V1',POS="V4",REF="V6")]
                 ,on=c(V1="#CHROM",V4="POS",V6="REF"),c("NNAME","CHR_REF","BP_REF") :=list(i.ID,`i.#CHROM`,i.POS)]
  cat("Updated rsid's\n")
  
  nrows <- nrows + increment
}

refpanelBim.dt[is.na(NNAME),NNAME:=V2]
cat("\nSet reference panel original names for variants without rs-id information\n")
refpanelBim.dt[is.na(NNAME),NNAME:=paste(V1,V4,V5,V6,sep = ":")]
cat("\nSet back-stop names for names still missing\n")
  
fwrite(x = refpanelBim.dt[NNAME!=V2,c("V2","NNAME")], file = "newnames.txt",col.names = F, sep = "\t")

cat("\nTHE END\n")

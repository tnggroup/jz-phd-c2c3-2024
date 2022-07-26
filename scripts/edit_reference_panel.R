#new edit reference panel .bim and create variant list - this outputs a new .bim rather than a name-file
# sets rsids
#resolves ambiguous naming
#deals with multiallelic variants
#replaces previous set_reference_panels_rsids, resolve_duplicate_and_unnamed_variants
#includes new concepts from supermunge and detecting systematic coordinate differences - maybe

#run in the folder of the working reference panel


library(data.table)
library(shru)

nThreads<-6

filepath.dbSNP <- file.path("..","dbsnp.human_9606_b151_GRCh38p7","00-All.vcf") #using the unzipped version of ~131GB

filepath.bim <- "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.bim_orig"
filepath.newbim <- "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.bim"
filepath.frq<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.frq"


increment<-5000000
nrows <- 56
dbSNP.nrows<-660146231 #set here manually because it is difficult for the program to know otherwise
#dbSNP.colnames<-c("CHR","BP","SNP","REF","ALT","QUAL","FILTER","INFO")
dbSNP.colnames<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

bim<-fread(file = filepath.bim, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
colnames(bim)<-c("CHR","SNP","CM","BP","A1","A2")
bim$ORDER<-1:nrow(bim)
setDT(bim)

maf <- fread(file = filepath.frq, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = F, nThread = nThreads, showProgress = T)
maf$ORDER<-1:nrow(maf)
setDT(maf)


#data formatting
#bim[,CHR:=as.integer(shru::parseCHRColumn(CHR))]
#bim[BP:=as.integer(BP)]
bim[is.na(SNP),SNP:="."]
#head(unique(bim$A1),1000)
setkeyv(bim, cols = c("ORDER","CHR","BP","A2","A1","SNP"))

#maf[,CHR:=as.integer(shru::parseCHRColumn(CHR))]
maf[is.na(SNP),SNP:="."]
setkeyv(maf, cols = c("ORDER","CHR","A2","A1","SNP"))

#The order or variants in the .bim and .frq must be the same
bim[maf, on=c(ORDER="ORDER",CHR="CHR",A2="A2",A1="A1",SNP="SNP"), c("MAF","NCHROBS"):=list(i.MAF,i.NCHROBS)]
#bim[maf, on=c(ORDER="ORDER",CHR="CHR",A2="A2",A1="A1",SNP="SNP"), c("MAF","NCHROBS"):=list(i.MAF,i.NCHROBS), by=.EACHI]
setkeyv(bim, cols = c("CHR","BP","A2","A1","SNP"))

rm(maf)


doBreak<-F
#loop to read the dbSNP in chunks because big
m<-c()
for(iIteration in 1:10000){
#for(iIteration in 1:3){ #TEST!!
  #iIteration<-1
  cat("\n\n***Iteration",iIteration,"***\n")
  
  cat("\nReading dbSNP\n")
  if(nrows+increment>dbSNP.nrows) {
    increment<-nrows+increment-dbSNP.nrows #last increment
    doBreak<-T
  }
  
  # if(increment<1) {
  #   cat("Breaking because of no increment\n")
  #   break
  # }
  
  dbSNP<-fread(file = filepath.dbSNP, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, data.table = T,showProgress = T, nThread=nThreads, nrows = increment, skip = nrows, sep="\t")
  cat("\nRead dbSNP\n")
  
  colnames(dbSNP)<-dbSNP.colnames
  
  if(doBreak | nrow(dbSNP)<1 | length(colnames(dbSNP))<5){
    cat("Breaking because of last increment or empty dt!\n")
    break
  }
  
  #colnames(dbSNP) <- colnames(dbSNP.colnames)
  
  #print(dbSNP[,1:8])
  
  #parse dbSNP INFO/metdata
  dbSNP$INFO.split<-strsplit(dbSNP$INFO, split = ";",fixed = T)
  dbSNP$dbSNPBuildID.string<-lapply(dbSNP$INFO.split,FUN = function(x){ifelse(3<=length(x),x[[3]],NA)})
  dbSNP[,dbSNPBuildID.string:=as.character(dbSNPBuildID.string)]
  dbSNP[!grepl(pattern = "^dbSNPBuildID", x = dbSNP$dbSNPBuildID.string, fixed = F),dbSNPBuildID.string:=NA]
  dbSNP$dbSNPBuildID.string.split<-strsplit(dbSNP$dbSNPBuildID.string, split = "=",fixed = T)
  dbSNP$dbSNPBuildID<-lapply(dbSNP$dbSNPBuildID.string.split,FUN = function(x){ifelse(2<=length(x),x[[2]],NA)})
  dbSNP[,dbSNPBuildID:=as.integer(dbSNPBuildID)]
  
  dbSNP[,`#CHROM`:=as.integer(shru::parseCHRColumn(`#CHROM`))]#this may potentially be slow
  #setkeyv(dbSNP,cols = c("CHR","BP","REF","ALT"))
  dbSNP<-dbSNP[,c("#CHROM","POS","ID","REF","ALT","dbSNPBuildID")]
  setkeyv(dbSNP,cols = c("#CHROM","POS","REF","ALT"))
  cat("Indexed dbSNP\n")
  
  cat("Overview of dbSNP data table:\n\n")
  print(dbSNP)
  
  
  cat("\n\nFetching rsid's\n")
  #Alleles can have alternate forms in the VCF separated by a comma
  multiREF<-dbSNP[grepl(pattern = ",", x = dbSNP$REF, fixed = T),]
  multiALT<-dbSNP[grepl(pattern = ",", x = dbSNP$ALT, fixed = T),]
  # if(nrow(multiREF)>0) {
  #   warning("Multiallele in dbSNP REF detected!")
  #   multiREF
  #   }
  # if(nrow(multiALT)>0) {
  #   warning("Multiallele in dbSNP ALT detected!")
  #   multiALT
  # }
  
  if(nrow(multiREF)>0 | nrow(multiALT)>0){
    #add alternate alleles
    
    dbSNP<-dbSNP[!(grepl(pattern = ",", x = dbSNP$REF, fixed = T) | grepl(pattern = ",", x = dbSNP$ALT, fixed = T)),][,A1:=ALT][,A2:=REF][,ALT:=NULL][,REF:=NULL]
    
    #ALT
    multiALT$asplit<-strsplit(multiALT$ALT, split = ",",fixed = T)
    asplitMaxlength<-0
    crap <- lapply(multiALT$asplit,FUN = function(x){if(length(x)>asplitMaxlength) asplitMaxlength<<- length(x)})
    rm(crap)
    for(iAllele in 1:asplitMaxlength){
      #iAllele<-3
      cMulti<-multiALT
      cMulti[,A2:=REF]
      cMulti$A1<-lapply(cMulti$asplit,FUN = function(x){ifelse(iAllele<=length(x),x[[iAllele]],NA)})
      cMulti[,ALT:=NULL][,REF:=NULL][,asplit:=NULL]
      dbSNP<-rbind(dbSNP,cMulti[!is.na(A1),])
    }
    
    dbSNP[,A1:=as.character(A1)] #not sure why A1 gets interpreted as a list sometimes
    rm(multiREF)
    rm(multiALT)
    rm(cMulti)
    
  } else {
    dbSNP[,A1:=ALT][,A2:=REF][,ALT:=NULL][,REF:=NULL]
  }
  
  
  
  #standardise dbSNP columns
  colnames(dbSNP)<-c("CHR","BP","SNP","dbSNPBuildID","A1","A2")
  setkeyv(dbSNP, cols = c("CHR","BP","A2","A1","SNP"))
  
  
  

  #capture updates update bim from dbSNP
  
  #exact matches
  m.exact<-bim[dbSNP,on=c(CHR="CHR",BP="BP",A2="A2",A1="A1"), nomatch=0][,SNP:=i.SNP]
  #m.exact.reverse<-bim[dbSNP,on=c(CHR="CHR",BP="BP",A2="A1",A1="A2"), nomatch=0] #we don't need these as dbSNP has unique rsid's for inverted allele configurations:
  #dbSNP[dbSNP,on=c(CHR="CHR",BP="BP",A2="A1",A1="A2"),nomatch=0] #check this with
  #m.tot<-rbind(m.exact,m.exact.reverse)
  
  m<-rbind(m,m.exact,fill=T)
  
  
  #nrows <- clipValues(x = (nrows + increment),min = NA, max = dbSNP.nrows)
  nrows <- nrows + increment
  
  #to save memory
  rm(dbSNP)
}

#to save memory
rm(dbSNP)

#we have possible duplicate entries over CHR,BP,A1,A2

m<-m[order(-dbSNPBuildID),]
m.unique<-m[, .(SNP = head(SNP,1),MAF = head(MAF,1),dbSNPBuildID = head(dbSNPBuildID,1)), by = c("CHR","BP","A1","A2")]
setkeyv(m.unique, cols = c("CHR","BP","A2","A1","SNP"))

#update bim
bim[m.unique,on=c(CHR="CHR",BP="BP",A2="A2",A1="A1"),c("SNP","BP.dbSNP"):=list(i.SNP,i.BP)]

cat("\nAssigned RS numbers: ",nrow(bim[grepl(pattern = "^rs.+", x = bim$SNP, fixed = F),]),"\n")

#still missing SNP
bim[is.na(SNP) | grepl(pattern = ".", x = bim$SNP, fixed = T),SNP:=paste(CHR,BP,ORDER,sep = ":")]
cat("\nSet back-stop names for names still missing\n")

#output .bim format again in the right order
fwrite(x = bim[order(ORDER),c("CHR","SNP","CM","BP","A1","A2")], file = filepath.newbim,col.names = F, sep = "\t",nThread = nThreads, na = ".",quote = F)

cat("\nTHE END\n")

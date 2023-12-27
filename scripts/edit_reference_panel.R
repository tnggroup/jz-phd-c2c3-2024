#new edit reference panel .bim and create variant list - this outputs a new .bim rather than a name-file
# sets rsids
#resolves ambiguous naming
#deals with multiallelic variants
#replaces previous set_reference_panels_rsids, resolve_duplicate_and_unnamed_variants
#removes duplicates over genomic coordinates but NOT SNP identifiers (as there may be multiple variants with the same ID)

#run in the folder of the working reference panel

#revision history
#version 2 - improved handling of variant row duplicates as seen from duplicated SNP variants later
#version 3 - improved handling of forward and backward matched rsID's and their synonyms. we need to keep both matched IDs as both may be encountered in existing GWAS sumstats. improved handling of variant row duplicates - keep SNP duplicates.


library(data.table)
library(shru)

nThreads<-6
setDTthreads(nThreads)

filepath.dbSNP <- file.path("..","dbsnp.human_9606_b151_GRCh38p7","00-All.vcf") #using the unzipped version of ~131GB

filepath.bim <- "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.bim_orig"
filepath.newbim <- "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.bim"
filepath.frq<-"1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.frq"

filepath.varlist <- "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.varlist.gz"
filepath.varlist.unfiltered <- "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.varlist.unfiltered.gz"

filepath.varlist.important <- "unmatched.sig.1kg.Rds"

increment<-10000000 #higher will lead to memory issues of CREATE
#increment<-1000000 #test!
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


#extra
importantSNPs<-c()
if(file.exists(filepath.varlist.important)){
  importantSNPs<-readRDS(filepath.varlist.important)
}

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


if(!file.exists(filepath.varlist.unfiltered)){
  
  doBreak<-F
  #loop to read the dbSNP in chunks because big
  m.forward<-c()
  m.reverse<-c()
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
    
    if(nrow(dbSNP)<1 | length(colnames(dbSNP))<5){
      cat("Breaking because of empty dt!\n")
      break
    }
    
    #colnames(dbSNP) <- colnames(dbSNP.colnames)
    
    #print(dbSNP[,1:8])
    
    #parse dbSNP INFO/metdata
    dbSNP$INFO.split<-strsplit(dbSNP$INFO, split = ";",fixed = T)
    dbSNP$dbSNPBuildID.string<-lapply(dbSNP$INFO.split,FUN = function(x){grep(pattern = "^dbSNPBuildID", x = x, fixed = F, value=T)[[1]]})
    #dbSNP$dbSNPBuildID.string<-lapply(dbSNP$INFO.split,FUN = function(x){ifelse(3<=length(x),x[[3]],NA)})
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
    
    # dbSNP.bp.max <- max(dbSNP$BP,na.rm = T)
    # dbSNP.bp.min <- min(dbSNP$BP,na.rm = T)
    
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
    #important
    dbSNP.important <- dbSNP[SNP %in% importantSNPs$SNP,]
    
    #exact matches
    m.exact<-bim[dbSNP,on=c(CHR="CHR",BP="BP",A2="A2",A1="A1"), nomatch=0][,SNP:=i.SNP][,i.SNP:=NULL]
    m.exact.reverse<-bim[dbSNP,on=c(CHR="CHR",BP="BP",A2="A1",A1="A2"), nomatch=0][,SNP:=i.SNP][,i.SNP:=NULL]
    #dbSNP has unique rsid's for inverted allele configurations:
    #dbSNP[dbSNP,on=c(CHR="CHR",BP="BP",A2="A1",A1="A2"),nomatch=0] #check this with
    
    m.forward<-rbind(m.forward,m.exact,fill=T)
    m.reverse<-rbind(m.reverse,m.exact.reverse,fill=T)
    
    #this was actually not here for the last run, so may have missed the last iteration for chr 24
    if(doBreak){
      cat("Breaking because of last increment\n")
      break
    }
    
    #nrows <- clipValues(x = (nrows + increment),min = NA, max = dbSNP.nrows)
    nrows <- nrows + increment
    
    #to save memory
    rm(dbSNP)
    
  }
  
  #to save memory
  rm(dbSNP)
  
  m.forward[,forward:=1]
  m.reverse[,forward:=0]
  
  #merge exact matches with reverse matches
  m<-rbind(m.forward,m.reverse,fill=T)
  
  rm(m.forward)
  rm(m.reverse)
  
  #output unfiltered varlist - before merging across position!
  fwrite(x = m,file = filepath.varlist.unfiltered, append = F,quote = F,sep = "\t",col.names = T,nThread=nThreads)
  
} else {
  #read in previous result
  #
  m <- shru::readFile(filePath = filepath.varlist.unfiltered, nThreads = nThreads)
}

cSumstats <- shru::readFile(filePath = "/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/cleaned/ANXI03.gz", nThreads = nThreads)
#newbim <- shru::readFile(filePath = filepath.newbim, nThreads = nThreads)

#cross version unique, including synonym lists - we need multiple SNP occurrences because the panel contains multiple instances of the same rsID SNP
#exclude possible duplicate entries over CHR,BP,A1,A2, aimed at removing duplicates due to dbSNP version updates
setkeyv(m, cols = c("CHR","BP","A2","A1","SNP"))
setorder(m,-dbSNPBuildID,-forward,-MAF,-NCHROBS,SNP)
m.xversion.unique<-m[, .(SNP = head(SNP,1), MAF = head(MAF,1), dbSNPBuildID = head(dbSNPBuildID,1), ORDER = head(ORDER,1), NCHROBS = head(NCHROBS,1), forward=head(forward,1), synonyms=paste(SNP,collapse = ",")), by = c("CHR","BP","A2","A1")] #synonyms=list(c(SNP))
# #filter to unique SNP
# setorder(m.xversion.unique,-dbSNPBuildID,-forward,-MAF,-NCHROBS,ORDER)
# m.xversion.unique<-m.xversion.unique[, .(ORDER = head(ORDER,1), CHR = head(CHR,1), BP = head(BP,1), A1 = head(A1,1), A2 = head(A2,1), MAF = head(MAF,1), dbSNPBuildID = head(dbSNPBuildID,1), NCHROBS = head(NCHROBS,1), forward=head(forward,1), synonyms=head(synonyms,1)), by = c("SNP")]

# #cross SNP unique - for bim SNP labels
# #filter to unique SNP
# setorder(m,-dbSNPBuildID,-forward,-MAF,-NCHROBS,ORDER)
# m.xsnp.unique<-m[, .(ORDER = head(ORDER,1), CHR = head(CHR,1), BP = head(BP,1), A1 = head(A1,1), A2 = head(A2,1), MAF = head(MAF,1), dbSNPBuildID = head(dbSNPBuildID,1), NCHROBS = head(NCHROBS,1), forward=head(forward,1)), by = c("SNP")]
# #exclude possible duplicate entries over CHR,BP,A1,A2, aimed at removing duplicates due to dbSNP version updates
# setorder(m.xsnp.unique,-dbSNPBuildID,-forward,-MAF,-NCHROBS,SNP)
# m.xsnp.unique<-m.xsnp.unique[, .(SNP = head(SNP,1), MAF = head(MAF,1), dbSNPBuildID = head(dbSNPBuildID,1), ORDER = head(ORDER,1), NCHROBS = head(NCHROBS,1), forward=head(forward,1)), by = c("CHR","BP","A1","A2")]


#update bim
setkeyv(bim, cols = c("ORDER","CHR","BP","A2","A1"))
bim[,SNPASSIGNED:=FALSE]
#unique forward
setkeyv(m.xversion.unique, cols = c("ORDER","CHR","BP","A2","A1"))
bim[m.xversion.unique,on=c(ORDER="ORDER",CHR="CHR",BP="BP",A2="A2",A1="A1"),c("SNP","BID","forward","SNPASSIGNED"):=list(i.SNP,i.dbSNPBuildID,i.forward,TRUE)]


#explicitly check rs-numbers in SNP
#cat("\nAssigned RS numbers: ",nrow(bim[grepl(pattern = "^rs.+", x = bim$SNP, fixed = F),]),"\n")
cat("\nAssigned SNP ID's (RS numbers): ",nrow(bim[SNPASSIGNED==TRUE,]),"\n")

#still missing SNP
bim[is.na(SNP) | grepl(pattern = ".", x = bim$SNP, fixed = T), SNP:=paste(CHR,BP,ORDER,sep = ":")]
cat("\nSet back-stop names for names still missing\n")

#output .bim format again in the right order
fwrite(x = bim[order(ORDER),c("CHR","SNP","CM","BP","A1","A2")], file = filepath.newbim,col.names = F, sep = "\t",nThread = nThreads, na = ".",quote = F)

#output general variant list + synonyms - still has the overall MAF values
#forward SNPS
bim[m.xversion.unique[forward==1,],on=c(ORDER="ORDER",CHR="CHR",BP="BP",A2="A2",A1="A1"),c("SNPF"):=list(i.SNP)]
bim[m.xversion.unique,on=c(SNPF="SNP",ORDER="ORDER",CHR="CHR",BP="BP",A2="A2",A1="A1"),c("SYNF"):=list(i.synonyms)]
#reverse SNPS
bim[m.xversion.unique[forward==0,],on=c(ORDER="ORDER",CHR="CHR",BP="BP",A2="A2",A1="A1"),c("SNPR"):=list(i.SNP)]
bim[m.xversion.unique,on=c(SNPR="SNP",ORDER="ORDER",CHR="CHR",BP="BP",A2="A2",A1="A1"),c("SYNR"):=list(i.synonyms)]


#how many important SNPs covered
cat("\nImportant SNPs covered: ",(sum(importantSNPs$SNP %in% bim$SNPF) + sum(importantSNPs$SNP %in% bim$SNPR))," of a total ",nrow(importantSNPs),"\n")

cat("\ncSumstats (ANXI03) SNPs covered: ",(sum(cSumstats$SNP %in% bim$SNPF) + sum(cSumstats$SNP %in% bim$SNPR))," of a total ",nrow(cSumstats),"\n")


fwrite(x = bim[order(ORDER),.(SNP,CHR,BP,A1,A2,MAF,NCHROBS,FW=forward,SNPF,SNPR,SYNF,SYNR)],file = filepath.varlist, append = F,quote = F,sep = "\t",col.names = T,nThread=nThreads) #added CHR also here - not included in the last run

cat("\nTHE END\n")

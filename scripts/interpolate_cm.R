#run in the genetic_recombination_mapping folder
library(optparse)
library(data.table)

clParser <- OptionParser()
clParser <- add_option(clParser, c("-c", "--chromosome"), type="character", default="22",
                       help="help")
clOptions<-parse_args(clParser)


mTot<-fread(file = file.path("genetic-map-chr-bp-rr-cm.1KGP3.grch38.gz"), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T,showProgress = F, nThread=4)

#print SHAPEIT maps per chromosome: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap
#chromosomes <- unique(mTot$CHR)

#Update the BIM!
d <- fread(file = file.path("../reference_panel/hc1kg.grch38.plink/1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.bim"), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T,showProgress = F, nThread=4, header = F)
colnames(d)<-c("CHR","SNP","CM","BP","A1","A2")
#for test!!! d<-d[1:100000,] d<-d[CHR=="chrX",]
d$NSNP <- paste0(1:nrow(d),"_interpolate_cm_UNKNOWN")
d[is.na(SNP),SNP:=NSNP]
d$CHR<-shru::parseCHRColumn(d$CHR)
d$CM<-as.numeric(d$CM)
mTot$CM<-as.numeric(mTot$CM)
mTot$RR_CM_BP<-as.numeric(mTot$RR_CM_MB)/10^6
setkeyv(d,cols = c("CHR","SNP","BP"))
chromosomes <- unique(d$CHR)
print(chromosomes)
print("Running interpolation of CM!")
#for(i in 1:length(chromosomes)){
  #i<-1
  #cChr<-chromosomes[i]
  cChr<-clOptions$chromosome
  print(paste("Processing chromosome",cChr))
  dChr<-d[CHR==eval(cChr),][order(BP,SNP),]
  #for test!!!
  dChr<-dChr[1:1000000,] #REMOVE THIS
  
  setkeyv(dChr,cols = c("SNP","BP"))
  
  minBP<-min(dChr$BP, na.rm = T)
  maxBP<-max(dChr$BP,na.rm = T)
  
  mTotChr<-mTot[CHR==eval(cChr) & BP>=eval(minBP) & BP<=eval(maxBP),][order(BP,SNP),]
  setkeyv(mTotChr,cols = c("BP"))
  
  #remove d if not using it later
  rm(d)
  
  if(nrow(mTotChr)>0){
    #dChr[mTotChr,on='BP',CM:=i.CM] #we don't have an exact match on position for all positions
    #dChr[BP %in% mTotChr$BP,]
    #mTotChr[BP %in% dChr$BP,]
    vLastBP<-0
    vLastCM<-0
    vLastRR_CM_BP<-0
    for(iMap in 1:nrow(mTotChr)){
      #iMap<-2
      #interpolation: (DCM/DBP + (RR0*(DBP-d)/DBP + RR1*(d)/DBP)*(DBP/2-abs(DBP/2-d)))d + CM0 #DOES NOT WORK!!
      
      dChr[BP>eval(vLastBP) & BP<=eval(mTotChr[iMap,]$BP),c("DBP","DCM","d","RR0","RR1","CM0"):=list((eval(mTotChr[iMap,]$BP)-eval(vLastBP)),(eval(mTotChr[iMap,]$CM)-eval(vLastCM)),(BP-eval(vLastBP)),eval(vLastRR_CM_BP),eval(mTotChr[iMap,]$RR_CM_BP),eval(vLastCM))]
      vLastBP<-mTotChr[iMap,]$BP
      vLastCM<-mTotChr[iMap,]$CM
      vLastRR_CM_BP<-mTotChr[iMap,]$RR_CM_BP
      
    }
  } else {
    #set back-stop CM values based on a genomic average
    settingGenomicAverageRecombinationRate_CM_BP<-1/10^6
    dChr[,c("DBP","DCM","d","RR0","RR1","CM0"):=list(eval(maxBP),(eval(maxBP)*eval(settingGenomicAverageRecombinationRate_CM_BP)),BP,eval(settingGenomicAverageRecombinationRate_CM_BP),eval(settingGenomicAverageRecombinationRate_CM_BP),0)]
  }
  #dChr[,CM:=(DCM/DBP + (RR0*(DBP-d)/DBP + RR1*d/DBP)*(DBP/2-abs(DBP/2-d)))*d + CM0]
  dChr[,CM:=(DCM/DBP)*d + CM0]
  
  fwrite(x = dChr[,c("CHR","SNP","CM","BP","A1","A2")], file = paste0("../reference_panel/hc1kg.grch38.plink/1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_CM_",cChr,".gz"), append = F, quote = F, sep = "\t",na = "NA", col.names = F,nThread=4)
  
  
  #d[dChr,on=c('CHR','SNP'),CM:=i.CM]
  
#}

#fwrite(x = d[,c("CHR","SNP","CM","BP","A1","A2")], file = "../reference_panel/hc1kg.grch38.plink/1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_CM.bim", append = F, quote = F, sep = "\t",na = "NA", col.names = F,nThread=4)


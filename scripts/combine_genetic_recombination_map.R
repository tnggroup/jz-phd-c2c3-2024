#read in and combine genetic recombination map
#run in the genetic_recombination_mapping folder

library(data.table)


mTot<-c()
for(iChr in 1:23){
  #iChr<-1
  m <- fread(file = file.path(paste0("genetic_map_chr",iChr,"_combined_b37.txt")), na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = T, check.names = T, fill = T, blank.lines.skip = T, key = c("position"), data.table = T,showProgress = F, nThread=4)
  m$CHR<-iChr
  mTot<-rbind(mTot,m)
}

colnames(mTot)<-c("BP","RR_CM_MB","CM","CHR")

mTot$SNP<-1:nrow(mTot)

mTot<-mTot[,c("SNP","CHR","BP","RR_CM_MB","CM")]

fwrite(x = mTot, file = "genetic-map-chr-bp-rr-cm.1KGP3.grch37.gz", append = F, quote = F, sep = "\t",na = "NA", col.names = T,nThread=4)

#if not already present
#mTot<-fread(file = file.path("genetic-map-chr-bp-rr-cm.1KGP3.grch37.gz"), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T,showProgress = F, nThread=4)

#chain file format reference: http://genome.ucsc.edu/goldenPath/help/chain.html
chain <- fread(file = "../alignment_chains/hg19ToHg38.over.chain.gz", na.strings =c(".",NA,"NA",""), encoding = "UTF-8", header = F, check.names = T, fill = T, blank.lines.skip = T, data.table = T,showProgress = F, nThread=4)

chain$row <- 1:nrow(chain)

chains.dt <- chain[V1=="chain",]
chains.dt[,V1:=NULL]
colnames(chains.dt) <- c("score","tName","tSize","tStrand","tStart","tEnd","qName","qSize","qStrand","qStart","qEnd","id","row")

#tName
chains.dt$tName<-shru::parseCHRColumn(chains.dt$tName)
chains.dt$strange_tName<-grepl(pattern = "^[^_]+_", chains.dt$tName)
indexesLengths<-regexec(pattern = "^([^_]+)_", text=chains.dt$tName)
matches<-regmatches(chains.dt$tName,indexesLengths)
chains.dt$tName[chains.dt$strange_tName] <- unlist(lapply(matches[chains.dt$strange_tName],FUN = function(x)x[2]))

#qName
chains.dt$qName<-shru::parseCHRColumn(chains.dt$qName)
chains.dt$strange_qName<-grepl(pattern = "^[^_]+_", chains.dt$qName)
indexesLengths<-regexec(pattern = "^([^_]+)_", text=chains.dt$qName)
matches<-regmatches(chains.dt$qName,indexesLengths)
chains.dt$qName[chains.dt$strange_qName] <- unlist(lapply(matches[chains.dt$strange_qName],FUN = function(x)x[2]))

#chromosome Un
#https://genome.ucsc.edu/FAQ/FAQdownloads.html#download11


setkeyv(chains.dt,cols = c("tStart","tEnd","qStart","qEnd","id","row"))
#chains.dt<-chains.dt[order("row")]

segments.dt <- chain[V1!="chain",]
segments.dt[,c("size","dt","dq"):=tstrsplit(V1,split = "\t",fixed = T)] #tstrsplit!
segments.dt[,colnames(segments.dt)[!colnames(segments.dt) %in% c("row","size","dt","dq")]:=NULL]
setkeyv(segments.dt,cols = c("row","size","dt","dq"))
chains.dt<-chains.dt[order(chains.dt$row),]

#update segments with chain id
for(i in 1:nrow(chains.dt)){
  #i<-1
  cRow<-chains.dt[i,c("row")][[1]]
  cId<-chains.dt[i,c("id")][[1]]
  segments.dt[row>cRow,chain:=cId]
}

#segments.dt[, cumsize := cumsum(size), by=list(chain)]

rm(chain) #we don't need chain anymore

#liftover mTot

#update mTot with chain id - sort on chain score!! =BLAT score?, -> overwrite lower score assignments later in the loop
chains.dt<-chains.dt[order(chains.dt$score),]
for(i in 1:nrow(chains.dt)){
  #i<-1
  cCHR<-chains.dt[i,c("tName")][[1]]
  cStartBP<-chains.dt[i,c("tStart")][[1]]
  cEndBP<-chains.dt[i,c("tEnd")][[1]]
  cNCHR<-chains.dt[i,c("qName")][[1]]
  cNStartBP<-chains.dt[i,c("qStart")][[1]]
  cId<-chains.dt[i,c("id")][[1]]
  cStrange_tName<-chains.dt[i,c("strange_tName")][[1]]
  cStrange_qName<-chains.dt[i,c("strange_qName")][[1]]
  # hits<-mTot[CHR==eval(cCHR) & BP>=eval(cStartBP) & BP<=eval(cEndBP),]
  # num<-nrow(hits)
  mTot[CHR==eval(cCHR) & BP>=eval(cStartBP) & BP<=eval(cEndBP),c('chain','strange_tName','strange_qName','CHR2','BP2') :=list(eval(cId),eval(cStrange_tName),eval(cStrange_qName),eval(cNCHR),(BP-eval(cStartBP)+eval(cNStartBP)))]
  
}

remapped <- mTot[CHR!=CHR2 | BP!= BP2,]
nrow(remapped)

unmapped <- mTot[is.na(CHR2) | is.na(BP2),]
nrow(unmapped)

mTot<-mTot[!(is.na(CHR2) | is.na(BP2)),][,c("CHR","BP"):=list(CHR2,BP2)]
nrow(mTot[CHR!=CHR2 | BP!= BP2,]) #check

mTot<-mTot[,c("SNP","CHR","BP","RR_CM_MB","CM")]
fwrite(x = mTot, file = "genetic-map-chr-bp-rr-cm.1KGP3.grch38.gz", append = F, quote = F, sep = "\t",na = "NA", col.names = T,nThread=4)


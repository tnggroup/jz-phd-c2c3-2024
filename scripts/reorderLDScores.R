#actions for reordering existing varlists
filepath.tref<-normalizePath(file.path(p$folderpath.data.variantLists,"hc1kgp3.b38.eur.l2.jz2023.gz"))
tref<-shru::readFile(filepath.tref)
setkeyv(tref, cols = c("SNP"))
setorder(tref,CHR,BP,CM,SNP,-NCHROBS,-MAF,-L2)
shru::writeFile(tref,filePath = filepath.tref)

filepath.tref<-normalizePath(file.path(p$folderpath.data.variantLists,"hc1kgp3.b38.mix.l2.jz2023.gz"))
tref<-shru::readFile(filepath.tref)
setkeyv(tref, cols = c("SNP"))
setorder(tref,CHR,BP,CM,SNP,-NCHROBS.MIX,-MAF.MIX,-L2.MIX)
shru::writeFile(tref,filePath = filepath.tref)

#reorder LD scores - LDSC-safe!
#HM3
filepath.tref<-normalizePath(file.path(p$folderpath.data.mvLDSC.ld.hm3))
tref <- do.call("rbind", lapply(list.files(path = filepath.tref, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
  suppressMessages(fread(file = file.path(filepath.tref, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = 6, data.table=T))
}))
setkeyv(tref, cols = c("SNP"))
setorder(tref,CHR,BP,CM,SNP,-MAF,-L2)
chrs<-order(unique(tref$CHR))
for(cChr in chrs){
  cFile<-file.path(paste0(filepath.tref,".ordered"),paste0(cChr,".l2.ldscore"))
  
  write.table(x = as.data.frame(tref[CHR==eval(cChr),]),file = cFile,sep="\t", quote = FALSE, row.names = F, append = F)
  nfilename.gz <- gzip(cFile,overwrite=T)
}
rm(tref)

#1kG
filepath.tref<-normalizePath(file.path(p$folderpath.data.mvLDSC.ld.1kg))
tref <- do.call("rbind", lapply(list.files(path = filepath.tref, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
  suppressMessages(fread(file = file.path(filepath.tref, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = 6, data.table=T))
}))
setkeyv(tref, cols = c("SNP"))
setorder(tref,CHR,BP,SNP,-L2)
cFile<-file.path(paste0(filepath.tref,".ordered"),paste0("1KG_Phase3.WG.CLEANED.EUR_MAF001.1cm.250blocks",".l2.ldscore"))
write.table(x = as.data.frame(tref),file = cFile,sep="\t", quote = FALSE, row.names = F, append = F)
nfilename.gz <- gzip(cFile,overwrite=T)
rm(tref)

#HC1kG EUR
filepath.tref<-normalizePath(file.path(p$folderpath.data.mvLDSC.ld.hc1kg))
tref <- do.call("rbind", lapply(list.files(path = filepath.tref, pattern = "\\.l2\\.ldscore\\.gz$"), function(i) {
  suppressMessages(fread(file = file.path(filepath.tref, i), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, showProgress = F, nThread = 6, data.table=T))
}))
setkeyv(tref, cols = c("SNP"))
setorder(tref,CHR,BP,SNP,-L2)
cFile<-file.path(paste0(filepath.tref,".ordered"),paste0("hc1kgp3.b38.eur.jz2023.1cm.250blocks",".l2.ldscore"))
write.table(x = as.data.frame(tref),file = cFile,sep="\t", quote = FALSE, row.names = F, append = F)
nfilename.gz <- gzip(cFile,overwrite=T)
rm(tref)

#HC1kG MIX - let's not do this

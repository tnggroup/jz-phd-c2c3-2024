#set rsids of the reference panel .bim file from dbSNP panel. output two column filr with old and new variant names/id's to be used by plink
#run in the folder of the working reference panel

library(data.table)
library(shru)

filepath.dbSNP <- file.path("..","dbsnp.human_9606_b151_GRCh38p7","00-All.vcf.gz")
#filepath.dbSNP <- file.path("..","small.vcf.sample.vcf.gz") #for testing

filepath.refpanelBim <- file.path("1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.bim")
cat("Reading dbSNP\n")
dbSNP.dt<-fread(file = filepath.dbSNP, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, data.table = T,showProgress = T, nThread=4)
cat("Read dbSNP\n")
#dbSNP.dt$`#CHROM`<-parseCHRColumn(dbSNP.dt$`#CHROM`) #this may potentially be slow
setkeyv(dbSNP.dt,cols = c("#CHROM","POS"))
cat("Indexed dbSNP\n")


refpanelBim.dt<-fread(file = filepath.refpanelBim, na.strings =c(".",NA,"NA",""), encoding = "UTF-8", fill = T, blank.lines.skip = T, data.table = T,showProgress = T, nThread=4, header = F)
refpanelBim.dt$V1<-parseCHRColumn(refpanelBim.dt$V1)
refpanelBim.dt$V1<-as.integer(refpanelBim.dt$V1)
setkeyv(refpanelBim.dt,cols = c("V1","V4"))
cat("Read and indexed bim\n")


refpanelBim.dt[dbSNP.dt,on=c(V1="#CHROM",V4="POS"),c("NNAME","CHR_REF","BP_REF") :=list(i.ID,`i.#CHROM`,i.POS)]
cat("Updated rsid's\n")

refpanelBim.dt[is.na(NNAME),NNAME:=V2]
cat("Set reference panel original names for variants without rs-id information\n")
refpanelBim.dt[is.na(NNAME),NNAME:=paste(V1,V4,V5,V6,sep = ":")]
cat("Set back-stop names for names still missing\n")

fwrite(x = refpanelBim.dt[NNAME!=V2,c("V2","NNAME")], file = "newnames.txt",col.names = F, sep = "\t")
cat("THE END\n")
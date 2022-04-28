library(shru)
library(data.table)
library(optparse)
library(ggplot2)
library(cowplot)
library(gt)

clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="character", default="0",
                       help="Task - used to specify trait code [default %default]")

clParser <- add_option(clParser, c("-a", "--task_argument"), type="character", default=NA,
                       help="Argument - used for specifying the level string [default %default]")

clOptions <- parse_args(clParser)

cTraitString <- clOptions$task
cLevelString <- clOptions$task_argument

#cTraitString<-"SCHI04.gz"
#cLevelString<-NA


#settings
filepathSNPReference <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt", mustWork = T)
folderpathLDscores <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/data/ld_scores/eur_w_ld_chr.1KG_Phase3", mustWork = T)
folderpathEvaluationSumstats <- normalizePath("/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/munged_1kg_eur_supermunge_noimp", mustWork = T)
folderpathEvaluationOutput <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/working_directory", mustWork = T)

#local settings
# filepathSNPReference <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt", mustWork = T)
# folderpathLDscores <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/data/ld_scores/eur_w_ld_chr.1KG_Phase3", mustWork = T)
# folderpathEvaluationSumstats <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/data/gwas_sumstats/munged_1kg_eur_supermunge_noimp", mustWork = T)
# folderpathEvaluationOutput <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/working_directory", mustWork = T)

folderpathEvaluationOutputLIMP <- file.path(folderpathEvaluationOutput,"LIMP")
#folderpathEvaluationOutputLIMP <- file.path(folderpathEvaluationOutput,"LIMP_evaluation_old_8kbwindow")
dir.create(path = folderpathEvaluationOutputLIMP)

#traitList <- c("ADHD05","ANXI02","COAD01","OBES01","SMOK04")
levelList <- c(0.05,0.1,0.15,0.2,0.25,0.3)


#munge part
# supermunge(
#   filePaths = paste0("/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/cleaned/",c("EDUC03","SCHI04","ANXI03"),".gz"),
#   refFilePath = filepathSNPReference,
#   ldDirPath=folderpathLDscores,
#   pathDirOutput = "/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/munged_1kg_eur_supermunge_noimp"
# )


filePath<-file.path(folderpathEvaluationSumstats,cTraitString)
d <- c()

d0 <- fread(file = filePath, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = 5, showProgress = F)
r0 <- fread(file = filepathSNPReference, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = 5, showProgress = F)
#limit to chr 1
d0<-d0[CHR==1,]
r0<-r0[CHR==1,]




#generate missing
if(is.na(cLevelString)){
  set.seed(2022) #to get a specific set of random values
  #which SNPs to remove
  for(iLevel in 1:length(levelList)){
    #iLevel<-1
    toremove<-sample(x = d0$SNP, size = nrow(d0)*levelList[iLevel])
    d[iLevel]<-list(d0[!SNP %in% toremove,])
    nfilepath<-file.path(folderpathEvaluationOutputLIMP,paste0(cTraitString,".",levelList[iLevel],".missing.gz"))
    fwrite(x = d[iLevel][[1]],file = nfilepath,sep="\t", quote = FALSE, row.names = F, append = F)
  }
  
  names(d)<-paste0(cTraitString,".",levelList,".imputed")
  
}

#alt read from files
if(is.na(cLevelString)){
  for(iLevel in 1:length(levelList)){
    #iLevel<-1
    nfilepath<-file.path(folderpathEvaluationOutputLIMP,paste0(cTraitString,".",levelList[iLevel],".missing.gz"))
    d[iLevel]<-list(fread(file = nfilepath, check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = 5, showProgress = F))
  }

  names(d)<-paste0(cTraitString,".",levelList,".imputed")

}
  
#rm("d0")


cOutputDir <- file.path(folderpathEvaluationOutputLIMP,cTraitString)
dir.create(path = cOutputDir)

if(is.na(cLevelString)){
  res<-supermunge(
    list_df = d,
    ref_df = r0,
    ldDirPath=folderpathLDscores,
    imputeFromLD=T,
    imputeFrameLenBp=100000,
    pathDirOutput = cOutputDir
  ) 
} else {
  cNameList <- paste0(cTraitString,".",cLevelString,".imputed")
  res<-supermunge(
    filePaths = filePath,
    ref_df = r0,
    ldDirPath=folderpathLDscores,
    traitNames = cNameList,
    imputeFromLD=T,
    imputeFrameLenBp=100000,
    pathDirOutput = cOutputDir
  ) 
}


#evaluate the result
setkeyv(d0,cols = "SNP")

dsnames<-paste0(cTraitString,'.',levelList)
ntoimpute<-list()
nimputed<-list()
ninformed<-list()
rmse<-list()
rmseInformed<-list()
cor<-list()
corInformed<-list()
corInfo<-list()
corK<-list()

scatterplots<-list()
scatterplotsInformed<-list()

for(iLevel in 1:length(levelList)){
  #iLevel<-1
  ntoimpute[iLevel]<-nrow(d0)-nrow(d[iLevel][[1]])
  
  cImp<-fread(file = file.path(folderpathEvaluationOutputLIMP,cTraitString,paste0(cTraitString,".",levelList[iLevel],".imputed.gz")), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = 5, showProgress = F)
  setkeyv(cImp,cols = "SNP")
  cImp<-cImp[!is.na(INFO.LIMP),]
  cImp[d0,on='SNP',c('EFFECT.O','SE.O'):=list(i.EFFECT,i.SE)]
  cImp<-cImp[!is.na(EFFECT.O)&!is.na(SE.O),]
  
  cImp[,MAF:=FRQ][FRQ>0.5,MAF:=1-FRQ]
  
  #evaluation parts
  
  condImputed<-cImp$MAF>0.01&cImp$MAF<0.49
  condInformed<-condImputed & cImp$INFO.LIMP>0.5
  
  nimputed[[iLevel]]<-nrow(cImp)
  cImp[,Z.O:=EFFECT.O/SE.O][,Z:=EFFECT/SE][,ABSZDIFF:=abs(Z-Z.O)][,ZDIFF2:=(Z-Z.O)^2][,IMPQ:=1/(1+abs(Z-Z.O))]
  
  rmse[[iLevel]]<-sqrt(mean(cImp[condImputed,]$ZDIFF2,na.rm=T))
  rmseInformed[[iLevel]]<-sqrt(mean(cImp[condInformed,]$ZDIFF2,na.rm=T))
  ninformed[[iLevel]]<-nrow(cImp[INFO.LIMP>0.5,])
  cor[[iLevel]]<-cor(cImp[condImputed,]$Z,cImp[condImputed,]$Z.O)
  corInformed[[iLevel]]<-cor(cImp[condInformed,]$Z,cImp[condInformed,]$Z.O)
  corInfo[[iLevel]]<-cor(-cImp[condImputed,]$ABSZDIFF,cImp[condImputed,]$INFO.LIMP)
  corK[[iLevel]]<-cor(-cImp[condImputed,]$ABSZDIFF,cImp[condImputed,]$K)
  
  scatterplots[[iLevel]]<-ggplot(cImp[condImputed,], aes(x=Z.O, Z, colour=INFO.LIMP)) +
    geom_abline(intercept =0 , slope = 1) +
    geom_point() +
    labs(x='Observed Z',y='Imputed Z', title=paste0(dsnames[iLevel],'\nN = ',nrow(cImp[condImputed,]),', cor=',round(cor[[iLevel]],digits = 2),', rmse=',round(rmse[[iLevel]],digits = 2))) +
    coord_fixed() +
    theme_half_open() +
    background_grid()
  
  
  scatterplotsInformed[[iLevel]]<-ggplot(cImp[condInformed,], aes(x=Z.O, Z, colour=INFO.LIMP)) +
    geom_abline(intercept =0 , slope = 1) +
    geom_point() +
    labs(x='Observed Z',y='Imputed Z', title=paste0(dsnames[iLevel],', informed variants','\nN = ',nrow(cImp[condInformed,]),', cor=',round(corInformed[[iLevel]],digits = 2),', rmse=',round(rmseInformed[[iLevel]],digits = 2))) +
    coord_fixed() +
    theme_half_open() +
    background_grid()
}

stats<-data.frame(
  dsnames=dsnames,
  ntoimpute=unlist(ntoimpute),
  nimputed=unlist(nimputed),
  ninformed=unlist(ninformed),
  rmse=unlist(rmse),
  rmseInformed=unlist(rmseInformed),
  cor=unlist(cor),
  corInformed=unlist(corInformed),
  corInfo=unlist(corInfo),
  corK=unlist(corK)
  )

stats$ratImputed<-stats$nimputed/stats$ntoimpute
stats$ratInformed<-stats$ninformed/stats$ntoimpute

write.csv(x = stats, file = file.path(folderpathEvaluationOutputLIMP,paste0(cTraitString,'_stats.csv')), row.names=T)
#gtsave(data = gt(data = stats), filename = file.path(folderpathEvaluationOutputLIMP,paste0(cTraitString,'_stats.rtf')))

png(file.path(folderpathEvaluationOutputLIMP,paste0(cTraitString,'_scatter.png')), units='px', height=3000, width=5000, res=300)
print(plot_grid(plotlist =scatterplots))
dev.off()

png(file.path(folderpathEvaluationOutputLIMP,paste0(cTraitString,'_scatterInformed.png')), units='px', height=3000, width=5000, res=300)
print(plot_grid(plotlist =scatterplotsInformed))
dev.off()

print("The end")

  
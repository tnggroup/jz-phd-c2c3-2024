library(shru)
library(data.table)
library(optparse)
library(ggplot2)
library(cowplot)
library(gt)

clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="character", default="generate",
                       help="Task to perform. The program can do several. generate - Generate files to apply imputation on with the specified missingness, based on the specified trait file.[default %default]")

clParser <- add_option(clParser, c("-d", "--trait"), type="character", default=NA,
                       help="Trait - specify trait file/code [default %default]")

clParser <- add_option(clParser, c("-l", "--level"), type="double", default=0.2,
                       help="Argument - used for specifying the level of missingness to impute [default %default]")

clOptions <- parse_args(clParser)

#manual settings for inline run
#clOptions$task<-"limp"
#clOptions$trait<-"ADHD05.gz"


clOptions$task
clOptions$trait
clOptions$level


#settings
# filepathSNPReference <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt", mustWork = T)
# folderpathLDscores <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/data/ld_scores/eur_w_ld_chr.1KG_Phase3", mustWork = T)
# folderpathEvaluationSumstats <- normalizePath("/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/munged_1kg_eur_supermunge", mustWork = T)
# folderpathEvaluationOutput <- normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1/working_directory", mustWork = T)

#local TEST! settings - using the old reference
filepathSNPReference <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt", mustWork = T)
folderpathLDscores <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/data/ld_scores/eur_w_ld_chr.1KG_Phase3", mustWork = T)
folderpathEvaluationSumstats <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/data/gwas_sumstats/munged_1kg_eur_supermunge", mustWork = T)
folderpathEvaluationOutput <- normalizePath("/Users/jakz/project/JZ_GED_PHD_C1/working_directory", mustWork = T)




#folder to store results in across runs
folderpathEvaluationOutput <- file.path(folderpathEvaluationOutput,"IMPTEST")
dir.create(path = folderpathEvaluationOutput)

#a string to identify this specific run, type of imputation etc
runIdentifier <-"LIMP_8kb"


filePathTrait<-ifelse(is.na(clOptions$trait),NA,file.path(folderpathEvaluationSumstats,clOptions$trait))

if(clOptions$task != "analyse"){
  d <- c()
  d0 <- fread(file = filePathTrait, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = 5, showProgress = F)
  r0 <- fread(file = filepathSNPReference, na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = 5, showProgress = F)
  
  #TEMP? limit to chr 1
  d0<-d0[CHR==1,]
  r0<-r0[CHR==1,]
  
  nrow(d0)
  
  #filter previously imputed variants if any
  if(any(colnames(d0)=="INFO.LIMP")){
    d0<-d0[is.na(INFO.LIMP),]
  }
  
  nrow(d0)

}

nfilepath<-file.path(folderpathEvaluationOutput,paste0(clOptions$trait,".",clOptions$level,".missing.gz"))
if(clOptions$task=="generate" || !file.exists(nfilepath)) {
  #generate missing
  set.seed(2022) #to get a specific set of random values
  #which SNPs to remove
  
  #iLevel<-1
  toremove<-sample(x = d0$SNP, size = nrow(d0)*clOptions$level)
  d1<-d0[!SNP %in% toremove,]
  nfilepath<-file.path(folderpathEvaluationOutput,paste0(clOptions$trait,".",clOptions$level,".missing.gz"))
  fwrite(x = d1,file = nfilepath,sep="\t", quote = FALSE, row.names = F, append = F)
  #for ssimp - needs special column formatting - ssimp can't read zipped files!!!
  nfilepath<-file.path(folderpathEvaluationOutput,paste0(clOptions$trait,".",clOptions$level,".missing.ssimp"))
  # d2 <- d1[,c("SNP","CHR","BP","A1","A2","Z","P","N")]
  # colnames(d2) <- c("SNP","chr","pos","A1","A2","Z","P","N")
  #fwrite(x = d1[,c("SNP","CHR","P","N","A1","A2","BP","Z")],file = nfilepath,sep="\t", quote = FALSE, row.names = F, append = F)
  write.table(d1[,c("SNP","CHR","P","N","A1","A2","BP","Z")], file = nfilepath, sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)
  #rm(d2)
  
  print("GENERATE task DONE!")
  if(clOptions$task=="generate") quit(save = "no")
} else {

  #alt read from files
    
  #iLevel<-1
  nfilepath<-file.path(folderpathEvaluationOutput,paste0(clOptions$trait,".",clOptions$level,".missing.gz"))
  d1<-fread(file = nfilepath, check.names = T, fill = T, blank.lines.skip = T, data.table = T, nThread = 5, showProgress = F)
}
#rm("d0")


if(clOptions$task=="limp"){
  
  #From 10.1016/j.ajhg.2008.06.005
  #via 10.1101/211821
  highld<-data.frame(matrix(data=NA,ncol=0,nrow=0))
  
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(1,48287980,52287979)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(2,86088342,101041482)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(2,134666268,138166268)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(2,183174494,190174494)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(3,47524996,50024996)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(3,83417310,96017310)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(5,44464243,50464243)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(5,97972100,100472101)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(5,128972101,131972101)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(5,135472101,138472101)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(6,25392021,33392022)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(6,56892041,63942041)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(6,139958307,142458307)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(7,55225791,66555850)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(8,7962590,11962591)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(8,42880843,49837447)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(8,111930824,114930824)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(10,36959994,43679994)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(11,46043424,57243424)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(11,87860352,90860352)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(12,33108733,41713733)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(12,111037280,113537280)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(17,31799963,33389579)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(17,40928985,42139672)
  highld[(nrow(highld)+1),c("CHR","BP1","BP2")]<-c(20,32536339,35066586)
  
  setDT(highld)
  setkeyv(highld, cols = c("CHR","BP1","BP2"))
  
  res<-supermunge(
    list_df = list(d1),
    traitNames = c(paste0(clOptions$trait,".",clOptions$level,".",runIdentifier)),
    ref_df = r0,
    ldDirPath=folderpathLDscores,
    imputeFromLD=T,
    imputeFrameLenBp=8000,
    pathDirOutput = folderpathEvaluationOutput,
    region.imputation.filter_df = highld
  ) 
  
  print("LIMP task DONE!")
  quit(save = "no")
}

if(clOptions$task=="limp"){
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
}
print("The end")

  
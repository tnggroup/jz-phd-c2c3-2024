## ----package setup, echo=FALSE, warning=F-------------------------------------------------
#install.packages("disk.frame")  
#install.packages("skimr")
#install.packages("psych")
#install.packages("Matrix")
#install.packages("tidyverse")
#install.packages("ggrepel")
#install.packages("gt")
#install.packages("kableExtra")
#remove.packages("corrplot")
#devtools::install_github("taiyun/corrplot", build_vignettes = TRUE)
#install.packages("corrplot")
#remove.packages("GenomicSEM")
#devtools::install_github("MichelNivard/GenomicSEM",ref = 'v2.1') #specify branch for better stability
#devtools::install_github("MichelNivard/GenomicSEM") #master branch as default
#remove.packages("GenomicSEM")
#devtools::install_github("johanzvrskovec/GenomicSEM",ref = 'mod-jz') #specify branch
#remove.packages("HDL")
#devtools::install_github("zhenin/HDL/HDL")
#devtools::install_github("zhenin/HDL/HDL@77cb9d0984d1302e40bfd871491e292f8f09f49d") #specify exact commit
#remove.packages("BiocManager")
#install.packages("BiocManager")
#remove.packages("MungeSumstats")
#devtools::install_github("neurogenomics/MungeSumstats")
#remove.packages("SNPlocs.Hsapiens.dbSNP144.GRCh38")
#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
#remove.packages("BSgenome.Hsapiens.NCBI.GRCh38")
#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
#install.packages("optparse")
#install.packages("stats")
#remove.packages("shru")
#devtools::install_github("johanzvrskovec/shru")
#install.packages("reticulate")
#install.packages("readr")  

#for testing supermunge
library(R.utils)
#library(disk.frame)
library(data.table)
library(optparse)
library(skimr)
library(ggrepel)
library(gt)
library(psych)
library(Matrix)
library(stats)
library(tidyverse)
library(HDL)
#library(MungeSumstats)
library(shru)
library(GenomicSEM)

#library(reticulate)



## ----command line setup-------------------------------------------------------------------
clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="character", default="0",
                help="Index of the explicit task to run separately:\n0: No task\nmvLD.mvLDSC:multivariate LDSC\nmvLD.HDL.piecewise:HDL Piecewise\nmvLD.HDL.jackknife:HDL Jackknife\nmvLD.origHDL:original HDL(jackknife)\nmvLD.origHDL.liabilityScale:original HDL with applied liability scale [default %default]")
clParser <- add_option(clParser, c("-l", "--location"), type="character", default="local",
                help="The place where the code is run [local,cluster] [default %default]")

clParser <- add_option(clParser, c("-a", "--task_argument"), type="character", default=NA,
                help="General purpose argument for tasks [default %default]")



## ----settings-----------------------------------------------------------------------------
p<-c() #create project metadata object
p$clOptions<-parse_args(clParser)
p$date.run<-Sys.Date()
p$setup.version<-7
p$setup.code<-paste0("setup",p$setup.version)
p$setup.code.date<-paste0(p$setup.code,"_",p$date.run)
p$filename.rmd<-paste0(p$setup.code,".Rmd")
p$filename.r<-paste0(p$setup.code,".R")
p$functions<-c()


p$host<-p$clOptions$location #this is the place where the code is run [local,cluster] read from command line - default local
p$seting.refreshPrepareSummaryStatistics<-FALSE
p$setting.refreshLatentFactorGWAS<-FALSE

#color theme settings
theme.color<-c()
theme.color$contrastDark1<-"#2D2D2D"
theme.color$contrastDark2<-"#CC99CC"
theme.color$contrastDark3<-"#6699CC"
theme.color$contrastDark4<-"#99CC99"
theme.color$contrastLight1<-"#66CCCC"
theme.color$contrastLight2<-"#FFCC66"
theme.color$contrastLight3<-"#F99157"
theme.color$contrastLight4<-"#F2777A"

#file path settings
##set project shared working directory **change if you have other settings**
if(p$host=="local") {
p$folderpath<-normalizePath("/Users/jakz/project/JZ_GED_PHD_C1")
} else if (p$host=="cluster") {
p$folderpath<-normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1")
}
##project working directory subfolders
p$folderpath.workingDirectory<-normalizePath(file.path(p$folderpath,"working_directory"))

p$folderpath.scripts<-normalizePath(file.path(p$folderpath,"scripts"))
#p$folderpath.includedSoftware<-normalizePath(file.path(p$folderpath,"included_software"))
p$folderpath.plots<-normalizePath(file.path(p$folderpath,"plots"))


#general data folder
if(p$host=="local") {
  p$folderpath.data<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data")
} else if (p$host=="cluster") {
  p$folderpath.data<-normalizePath(file.path(p$folderpath,"data"))
}

##cleaned sumstats folder
p$folderpath.data.sumstats.cleaned<-normalizePath(file.path(p$folderpath.data,"gwas_sumstats","cleaned"))

##munged sumstats folder
p$folderpath.data.sumstats.munged<-normalizePath(file.path(p$folderpath.data,"gwas_sumstats","munged_1kg_eur_supermunge"))
#p$folderpath.data.sumstats.munged<-normalizePath(file.path(p$folderpath.data,"gwas_sumstats","munged_1kg_eur_gSEM"))

##imputed sumstats folder
#p$folderpath.data.sumstats.imputed<-normalizePath(file.path(p$folderpath.data,"gwas_sumstat","imputed"))

##export sumstats folder - for corrected associations for secondary analyses
p$folderpath.data.sumstats.export<-normalizePath(file.path(p$folderpath.data,"gwas_sumstats","export"))


#python virtual environment folder
if(p$host=="local") {
  p$folderpath.pythonVenv<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/python-venv")
} else if (p$host=="cluster") {
  p$folderpath.pythonVenv<-normalizePath(file.path(p$folderpath,"python-venv"))
}

##Reference SNP-list (HapMap3 SNPs for example). Used for munging sumstat SNP data.
p$filepath.SNPReference.hm3<-normalizePath(file.path(p$folderpath.data,"w_hm3.noMHC.snplist")) #HapMap3 SNPs
## Used in the preparation step for performing latent factor GWAS as reference for calculating SNP variance across traits.
#p$filepath.SNPReference<-normalizePath(paste0(p$folderpath.data,"/","reference.1000G.maf.0.005.txt")) #1000 genomes phase 3
p$filepath.SNPReference.1kg<-normalizePath(file.path(p$folderpath.data,"combined.hm3_1kg.snplist.vanilla.jz2020.txt")) #custom hm3 + 1kg SNPs

p$filename.suffix.data.sumstats.munged<-".gz"
#p$filename.suffix.data.sumstats.munged<-"_noMHC.sumstats.gz"

##Reference panel folder containing individual level data reference panel. Used for GWAS sumstat imputation tasks.
#roject$folderpath.data.sumstatImp.genomeReference<-"/users/k1204688/brc_scratch/Public/1KG_Phase3/All"
p$folderpath.data.sumstatImp.genomeReference<-"/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/reference.panel.1KG_Phase3.CLEANED.EUR.cM"

##LD scores datasets folders (these strings need to have a trailing slash for the GSEM LDSC to work)
p$folderpath.data.mvLDSC.ld.1kg <- file.path(p$folderpath.data,"eur_w_ld_chr.1KG_Phase3")
p$folderpath.data.mvLDSC.ld.hm3 <- file.path(p$folderpath.data,"eur_w_ld_chr")
#Weights, if different from LD-scores
#Set weights to the same folder as ldscores
p$folderpath.data.mvLDSC.wld.1kg <- p$folderpath.data.mvLDSC.ld.1kg
p$folderpath.data.mvLDSC.wld.hm3 <- p$folderpath.data.mvLDSC.ld.hm3

##HDL LD scores reference - needs the trailing slashes!!!
if(p$host=="local") {
  #use the smallest LD reference as default for local tests
  p$folderpath.data.HDL.ld<-paste0(p$folderpath.data,"/UKB_array_SVD_eigen90_extraction/")
} else if (p$host=="cluster") {
  p$folderpath.data.HDL.ld<-paste0(p$folderpath.data,"/UKB_imputed_hm3_SVD_eigen99_extraction/")
}
  
##full script file paths
p$filepath.rmd<-normalizePath(file.path(p$folderpath.scripts,p$filename.rmd))
p$filepath.r<-normalizePath(file.path(p$folderpath.scripts,p$filename.r))

##CFA settings
p$CFA<-c()
#p$CFA$estimator=c("ML")
p$CFA$correlation<-c("COR","ORT") #ORT, OBL, and COR
p$CFA$estimator=c("ML")
p$CFA$nFactors=c(6)


##latent factor GWAS filter settings
p$lfGWAS$info.filter=.6
p$lfGWAS$maf.filter=0.01

#working directory in case of running as an R-script
setwd(dir = normalizePath(p$folderpath.workingDirectory))

#inactivated python environment until it is used
#use_virtualenv(p$folderpath.pythonVenv)



## ----plot functions-----------------------------------------------------------------------
library(corrplot) #do not run on cluster, does not work?

#c(theme.color$contrastDark3,theme.color$contrastLight3)
p$printCorr <- function(corr, SE=NULL, filename, addrect = 4, is.corr = T, number.cex = 1, number.digits=2){
  # corr <- p$mvLD$covstruct.mvLDSC.1kg$S_Stand.forPlot
  # SE <- p$mvLD$covstruct.mvLDSC.1kg$S_Stand.SE[p$sumstats.sel$code,p$sumstats.sel$code]
  # filename <- file.path(p$folderpath.plots,"rg.1kg_ld.png")
  # is.corr <- T
  
  # if(absScale){
  #   palette<-colorRampPalette(c("#FFFFFF",theme.color$contrastLight3))
  # } else {
    palette<-colorRampPalette(c(theme.color$contrastDark3,"#FFFFFF",theme.color$contrastLight3))
 # }
  
  
  
  if(is.null(SE)){
    corr.uppCI <- NULL
    corr.lowCI <- NULL
  } else {
    corr.uppCI<-clipValues(corr + 1.96 * SE, -1,1)
    corr.lowCI<-clipValues(corr - 1.96 * SE, -1,1)
  }
  
  png(filename = filename, width = 1200, height = 1000)
corrplot(
  corr = corr,
  uppCI.mat = corr.uppCI,
  lowCI.mat = corr.lowCI,
  plotCI = ifelse(is.null(SE),c("n"),"circle"),
  order = "hclust",
  hclust.method = "ward.D",
  method = "square",
  type="full",
  addCoef.col = theme.color$contrastDark1,
  addgrid.col = theme.color$contrastDark1,
  col = palette(200),
  is.corr = is.corr,
  outline = T,
  addrect = addrect,
  rect.col = theme.color$contrastLight1,
  rect.lwd = 6,
  tl.cex = 1.7,
  tl.col = theme.color$contrastDark2,
  tl.srt = 75,
  cl.cex = 2,
  cl.ratio = 0.2,
  cl.pos = "n",
  number.cex = number.cex,
  number.digits = number.digits
  )
dev.off()
  
}

p$printCorrSimplified <- function(corr, SE=NULL, filename, is.corr = T, number.cex = 1.5, number.digits=2, absScale=F){
  
  if(absScale){
    palette<-colorRampPalette(c("#FFFFFF",theme.color$contrastLight3))
  } else {
    palette<-colorRampPalette(c(theme.color$contrastDark3,"#FFFFFF",theme.color$contrastLight3))
  }
  
  if(is.null(SE)){
    corr.uppCI <- NULL
    corr.lowCI <- NULL
  } else {
    corr.uppCI<-corr + 1.96 * SE
    corr.lowCI<-corr - 1.96 * SE
  }
  
  png(filename = filename, width = 1200, height = 1000)
corrplot(
  corr = corr,
  #uppCI.mat = corr.uppCI,
  #lowCI.mat = corr.lowCI,
  #plotCI = ifelse(is.null(SE),c("n"),"circle"),
  order = "hclust",
  hclust.method = "ward.D",
  method = "color",
  type="full",
  #addCoef.col = theme.color$contrastDark1,
  addgrid.col = theme.color$contrastDark1,
  col = palette(200),
  is.corr = is.corr,
  outline = T,
  #addrect = addrect,
  #rect.col = theme.color$contrastLight1,
  #rect.lwd = 6,
  tl.cex = 1.7,
  tl.col = theme.color$contrastDark2,
  tl.srt = 35,
  cl.cex = 2,
  cl.ratio = 0.2,
  #number.cex = number.cex,
  #number.digits = number.digits
  )
dev.off()
  
}







## ----additional source setup, echo=FALSE, warning=F---------------------------------------

#source(normalizePath(file.path(p$folderpath.scripts,"sumstats.mod-jz.R")))





## ----sumstat metadata database load-------------------------------------------------------
p$filepath.sumstats<-file.path(p$folderpath.workingDirectory,paste0("sumstats.",p$setup.code,".Rds"))
if (file.exists(p$filepath.sumstats)) {
  print("Loading summary statistics metadata from previously stored file.")
  p$sumstats<-readRDS(file=p$filepath.sumstats)
} else {

#install.packages('RPostgres')
library(RPostgres)
library(DBI)

p$phenodbcon <- dbConnect(RPostgres::Postgres(),
                 dbname = 'phenodb', 
                 host = '10.200.105.5', 
                 port = 5432,
                 user = 'johan',
                 password = rstudioapi::askForPassword(prompt = "Enter database password for specified user."))

p$phenodbres <- dbSendQuery(p$phenodbcon, "SELECT \"GWAS\".*, category_id, category.name AS category_name, phenotype.name AS phenotype, phenotype.type AS phenotype_type, pmid, year FROM sumstat_old.\"GWAS\", sumstat_old.reference, sumstat_old.phenotype, sumstat_old.category 
WHERE \"GWAS\".reference_id=reference.id AND \"GWAS\".phenotype_id=phenotype.id AND phenotype.category_id = category.id
ORDER BY code,\"GWAS\".id")
p$sumstats<-dbFetch(p$phenodbres)
dbClearResult(p$phenodbres)

#fallback
#p$sumstats<-read.table(file.path(p$folderpath.data,"ukbb_sumstats_download202005.csv"), header=T, quote="\"", sep = ",", fill=T, blank.lines.skip=T,as.is = c(2), strip.white = T)

saveRDS(p$sumstats,file = p$filepath.sumstats)
write.table(p$sumstats, file = file.path(p$folderpath.workingDirectory,paste0(p$setup.code,".sumstats.tsv")), quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)
}





## ----trait setup--------------------------------------------------------------------------
#,echo=FALSE
p$trait<-data.frame(phenotype_id=c())

#ADHD
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-139
p$trait[nrow(p$trait),c("populationPrevalence")]<-.00529 #Worldwide-pooled
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1176/ajp.2007.164.6.942"

#ALCD
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-141
p$trait[nrow(p$trait),c("populationPrevalence")]<-.125 #US total measurement
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1001/archpsyc.64.7.830"

#ANOR
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-142
p$trait[nrow(p$trait),c("populationPrevalence")]<-.0245 #European average (3 secondary refs)
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1038/S41588-019-0439-2"

#ANXI
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-144
p$trait[nrow(p$trait),c("populationPrevalence")]<-.16 #Any type of anxiety disorder, Via https://doi.org/10.1038/s41380-019-0559-1 (2019), originally from https://doi.org/10.1017/S1121189X00001421 (2009)
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1017/S1121189X00001421"

#AUTI
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-145
p$trait[nrow(p$trait),c("populationPrevalence")]<-.0122
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1186/s12874-016-0280-6"

#BIPO
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-146
p$trait[nrow(p$trait),c("populationPrevalence")]<- .007 #mean of male and female global prevalence rate (2013) from https://doi.org/10.1111/bdi.12423 (2016)
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1111/bdi.12423"

#DEPR
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-149
p$trait[nrow(p$trait),c("populationPrevalence")]<-.146 #MDD, .15 from the LD-calculations in https://doi.org/10.1038/s41588-018-0090-3 (2018), but with a possible reference to https://doi.org/10.1146/annurev-publhealth-031912-114409 (2013) which states 14.6% lifetime prevalence of MDE in high-income countries.
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1146/annurev-publhealth-031912-114409"

#INSO
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-171
p$trait[nrow(p$trait),c("populationPrevalence")]<-.69 #using the 'prevalence', the higher estimate of persistent symptoms after 1y follow up.
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1093/sleep/30.3.274"

#PTSD
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-140
p$trait[nrow(p$trait),c("populationPrevalence")]<-.3 #using the moderate estimate of lifetime prevalence after trauma as used in ref
p$trait[nrow(p$trait),c("referenceDOI")]<-"10.1038/s41467-019-12576-w"

#SCHI
p$trait[nrow(p$trait)+1,c("phenotype_id")]<-151
p$trait[nrow(p$trait),c("populationPrevalence")]<-.0072 #using the lifetime morbid risk estimate, slightly more suitable for the CLOZUK sample as described in the GWAS but not using the lower point estimate of 0.4%
p$trait[nrow(p$trait),c("referenceDOI")]<-"https://doi.org/10.1093/epirev/mxn001"

p$trait



## ----setup GWAS sumstat dataset-----------------------------------------------------------
#, echo=FALSE


#View(p$sumstats)

#rename and add columns
names(p$sumstats)[names(p$sumstats)=="n_cases"]<-"n_case"
names(p$sumstats)[names(p$sumstats)=="n_controls"]<-"n_control"
p$sumstats$gwas_name.nice<-NA_character_
p$sumstats$code.trait<-NA_character_
p$sumstats$reference_doi<-NA_character_
p$sumstats$effect.logit<-as.logical(NA)
p$sumstats$dependent_variable.linprob<-as.logical(NA)
p$sumstats$se.logit<-as.logical(NA)
p$sumstats$dependent_variable.OLS<-as.logical(NA)
p$sumstats$age.min<-NA_integer_
p$sumstats$age.max<-NA_integer_
p$sumstats$age.mean<-NA_real_
p$sumstats$age.sd<-NA_real_

#add missing datasets and data
p$sumstats[nrow(p$sumstats)+1,c("code","n_case","n_control","n_total","phenotype_id","reference_id","reference_doi")]=list(
  code="DEPR05",
  n_case=16823,
  n_control=25632,
  n_total=42455,
  phenotype_id=149,
  reference_id=127,
  reference_doi=c("https://doi.org/10.1038/s41588-018-0090-3")
  )

#add missing datasets and data
p$sumstats[nrow(p$sumstats)+1,c("code","n_case","n_control","n_total","phenotype_id","reference_id","reference_doi")]=list(
  code="ANXI04",
  n_case=19012,
  n_control=58113,
  n_total=77125,
  phenotype_id=144,
  reference_id=158
  )

p$sumstats[nrow(p$sumstats)+1,c("code","name", "n_case","n_control","n_total","phenotype_id","reference_id","reference_doi")]=list(
  code="NEUR02",
  name="Neuroticism",
  n_case=449484,
  n_control=NA_integer_,
  n_total=449484,
  phenotype_id=176,
  reference_id=NA_integer_,
  reference_doi=c("https://doi.org/10.1038/s41588-018-0151-7")
  )

#add missing datasets and data
p$sumstats[nrow(p$sumstats)+1,c("code","name","n_total","phenotype_id")]=list(
  code="NOIS01",
  name="Noise 1",
  n_total=766345,
  phenotype_id=197
  )

p$sumstats[nrow(p$sumstats)+1,c("code","name","n_total","phenotype_id")]=list(
  code="NOIS02",
  name="Noise 2",
  n_total=766345,
  phenotype_id=197
  )



#reformat columns
p$sumstats$name<-as.character(p$sumstats$name)

##Add comprehensive names as in the Google sheet
p$sumstats$name[which(p$sumstats$code=="DEPR05")]="Major depressive disorder (PGC2 29) - only clinical ascertainment"

##Add trait/disorder information
p$sumstats <- p$sumstats %>%
mutate(
  code.trait=substr(x = code, start = 1, stop = 4)
       ) %>%
  left_join(p$trait[,c("phenotype_id","populationPrevalence")], by = c("phenotype_id" = "phenotype_id"))


#set code as rowname
rownames(p$sumstats)<-p$sumstats$code

##Add sumstat cleaned and munged file paths
p$sumstats$cleanedpath<-file.path(p$folderpath.data.sumstats.cleaned,paste0(p$sumstats$code,p$filename.suffix.data.sumstats.munged))

p$sumstats$mungedpath<-file.path(p$folderpath.data.sumstats.munged,paste0(p$sumstats$code,p$filename.suffix.data.sumstats.munged))


#special modification for special datasets for munging locally or at rosalind
if(p$host=="local"){
  p$sumstats["ANXI04",]$cleanedpath<-"/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/raw/GAD7_BGENIE_SexRegressed_For_LDSC.gz"
} else {
  p$sumstats["ANXI04",]$cleanedpath<-"/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/raw/GAD7_BGENIE_SexRegressed_For_LDSC.gz"
    }

if(p$host=="local"){
  p$sumstats["NEUR02",]$cleanedpath<-"/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/raw/sumstats_neuroticism_ctg_format.txt.gz"
} else {
  p$sumstats["NEUR02",]$cleanedpath<-"/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/raw/sumstats_neuroticism_ctg_format.txt.gz"
    }




##add reference year
p$sumstats["ANXI03",]$year=2019
p$sumstats["ANXI04",]$year=2019
p$sumstats["DEPR05",]$year=2018
p$sumstats["NEUR02",]$year=2018

##Add doi links for easy access to dataset publication
p$sumstats["ALCD03",]$reference_doi="https://doi.org/10.1038/s41593-018-0275-1"
p$sumstats["ANXI03",]$reference_doi="https://doi.org/10.1038/s41380-019-0559-1"
p$sumstats["ANXI04",]$reference_doi="https://doi.org/10.1038/s41380-019-0559-1"
p$sumstats["AUTI07",]$reference_doi="https://doi.org/10.1038/s41588-019-0344-8"
p$sumstats["DEPR08",]$reference_doi="https://doi.org/10.1038/s41588-018-0090-3"
p$sumstats["EDUC03",]$reference_doi="https://doi.org/10.1038/s41588-018-0147-3"
p$sumstats["HEAL01",]$reference_doi="https://doi.org/10.1093/ije/dyw219"
p$sumstats["NEUR01",]$reference_doi="https://doi.org/10.1038/ng.3552"
p$sumstats["SUBJ01",]$reference_doi="https://doi.org/10.1038/ng.3552"
p$sumstats["TIRE01",]$reference_doi="https://doi.org/10.1038/mp.2017.5"


##Add PMID
p$sumstats["ANXI03",]$pmid="31748690"
p$sumstats["ANXI04",]$pmid="31748690"
p$sumstats["DEPR05",]$pmid="29700475"

##add dependent variable type
p$sumstats$dependent_variable[which(p$sumstats$code=="DEPR05")]="binary"

##add participant numbers
p$sumstats$n_total[which(p$sumstats$code=="NEUR01")]=170911

##add ancestry details
p$sumstats$ancestry[which(p$sumstats$code=="NEUR01")]="EUR"
p$sumstats$ancestry[which(p$sumstats$code=="DEPR05")]="EUR"

##add sex details
p$sumstats$sex[which(p$sumstats$code=="ANXI04")]="both"
p$sumstats$sex[which(p$sumstats$code=="NEUR01")]="both"
p$sumstats$sex[which(p$sumstats$code=="DEPR05")]="both"

##add age range
p$sumstats$age.min[which(p$sumstats$code=="SUBJ01")]=40
p$sumstats$age.max[which(p$sumstats$code=="SUBJ01")]=73
p$sumstats$age.mean[which(p$sumstats$code=="SUBJ01")]=56.91
p$sumstats$age.sd[which(p$sumstats$code=="SUBJ01")]=7.93
p$sumstats$age.min[which(p$sumstats$code=="TIRE01")]=40
p$sumstats$age.max[which(p$sumstats$code=="TIRE01")]=73
p$sumstats$age.mean[which(p$sumstats$code=="TIRE01")]=56.91
p$sumstats$age.sd[which(p$sumstats$code=="TIRE01")]=7.93


##add sample prevalence for all datasets
p$sumstats$samplePrevalence<-p$sumstats$n_case/p$sumstats$n_total

##add number of cases or total column
p$sumstats$n_case_total<-ifelse(is.na(p$sumstats$n_case),p$sumstats$n_total,p$sumstats$n_case)

##Add nice trait names to be used in the report
p$sumstats$name.nice<-p$sumstats$name
#p$sumstats$name.nice[which(p$sumstats$code=="ADHD05")]="ADHD"
p$sumstats$name.nice[which(p$sumstats$code=="ALCD03")]="Alcohol dependence"
p$sumstats$name.nice[which(p$sumstats$code=="ANOR02")]="Anorexia nervosa"
p$sumstats$name.nice[which(p$sumstats$code=="ANXI03")]="Anxiety disorder"
p$sumstats$name.nice[which(p$sumstats$code=="ANXI04")]="Generalised anxiety symptoms"
p$sumstats$name.nice[which(p$sumstats$code=="AUTI07")]="Autism spectrum disorder"
p$sumstats$name.nice[which(p$sumstats$code=="DEPR05")]="MDD, narrow"
p$sumstats$name.nice[which(p$sumstats$code=="DEPR08")]="MDD"
p$sumstats$name.nice[which(p$sumstats$code=="EDUC03")]="Educational attainment"
p$sumstats$name.nice[which(p$sumstats$code=="EXTR01")]="Extraversion"
p$sumstats$name.nice[which(p$sumstats$code=="INCO03")]="Social deprivation"
p$sumstats$name.nice[which(p$sumstats$code=="INSO02")]="Insomnia"
p$sumstats$name.nice[which(p$sumstats$code=="INTE03")]="Cognitive ability"
p$sumstats$name.nice[which(p$sumstats$code=="LONG07")]="Longevity"
p$sumstats[c("NEUR01","NEUR02"),]$name.nice<-"Neuroticism"
p$sumstats$name.nice[which(p$sumstats$code=="RISK01")]="General risk tolerance"
p$sumstats$name.nice[which(p$sumstats$code=="RISK02")]="Risktaking, automobile"
p$sumstats$name.nice[which(p$sumstats$code=="RISK03")]="Risktaking, sex"
p$sumstats$name.nice[which(p$sumstats$code=="SCHI04")]="Schizophrenia"
p$sumstats$name.nice[which(p$sumstats$code=="SUBJ01")]="Subjective well-being"
p$sumstats$name.nice[which(p$sumstats$code=="TIRE01")]="Self-reported tiredness"

##Add a combined nice name plus code label
p$sumstats$name.nice.and_code<-paste0(p$sumstats$name.nice," (",p$sumstats$code,")")

#test
# tSumstats <- read.table(p$sumstats["NEUR01",]$cleanedpath,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# tSumstats$ptest<- 2*pnorm(abs(log(tSumstats$OR))/tSumstats$SE, lower.tail = F)
# tSumstats$ptest<- 2*pnorm(abs(tSumstats$OR)/tSumstats$SE, lower.tail = F)

##add information on wether a continous dependent variable was analysed using an OLS (linear) estimator. Used for the latent factor GWAS preparation step.
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="ADHD05")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="ALCD03")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="ANOR02")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="ANXI03")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="ANXI04")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="AUTI07")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="BIPO02")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="DEPR05")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="DEPR08")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="EDUC03")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="EXTR01")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="HEAL01")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="INCO03")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="INSO02")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="INTE03")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="LONG07")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="MIGR01")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="NEUR01")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="NEUR02")]=T #meta analysis
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="PTSD04")]=T #meta analysis
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="RISK01")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="RISK02")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="RISK03")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="SCHI04")]=F
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="SUBJ01")]=T
p$sumstats$dependent_variable.OLS[which(p$sumstats$code=="TIRE01")]=T

##add data on whether the effects were estimated with a linear regression rather than a logistic regression
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="ADHD05")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="ALCD03")]=T #from test
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="ANOR02")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="ANXI03")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="ANXI04")]=T #bgenie uses a linear estimator
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="AUTI07")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="BIPO02")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="DEPR05")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="DEPR08")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="EDUC03")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="EXTR01")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="HEAL01")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="INCO03")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="INSO02")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="INTE03")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="LONG07")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="MIGR01")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="NEUR01")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="NEUR02")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="PTSD04")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="RISK01")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="RISK02")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="RISK03")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="SCHI04")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="SUBJ01")]=F
p$sumstats$dependent_variable.linprob[which(p$sumstats$code=="TIRE01")]=F


##add data on whether the SEs are on a logistic scale or not
p$sumstats$se.logit[which(p$sumstats$code=="ADHD05")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="ALCD03")]=F #Tested
p$sumstats$se.logit[which(p$sumstats$code=="ANOR02")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="ANXI03")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="ANXI04")]=F #Tested
p$sumstats$se.logit[which(p$sumstats$code=="AUTI07")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="BIPO02")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="DEPR05")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="DEPR08")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="EDUC03")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="EXTR01")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="HEAL01")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="INCO03")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="INSO02")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="INTE03")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="LONG07")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="MIGR01")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="NEUR01")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="NEUR02")]=F #Meta analysis
p$sumstats$se.logit[which(p$sumstats$code=="PTSD04")]=F #Meta analysis
p$sumstats$se.logit[which(p$sumstats$code=="RISK01")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="RISK02")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="RISK03")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="SCHI04")]=T #Tested
p$sumstats$se.logit[which(p$sumstats$code=="SUBJ01")]=F #Continuous
p$sumstats$se.logit[which(p$sumstats$code=="TIRE01")]=F #Continuous


#p$sumstats[which(p$sumstats$code=="ANXI04"),]$populationPrevalence<-0.2 #set this to larger than for ANXI03 as this phenotype is broader
#p$sumstats[which(p$sumstats$code=="DEPR08"),]$populationPrevalence<-0.2 #set this to larger than for DEPR05 as this phenotype is broader

#set order of datasets to sorted by code
p$sumstats<-p$sumstats[order(p$sumstats$code),]

#save the project data
#saveRDS(project,file = file.path(p$folderpath.workingDirectory,paste0("project.",p$setup.code,".Rds")))

#View(p$sumstats)



## ----GWAS sumstat dataset variable selection----------------------------------------------

#selection based on specific traits
#p$sumstats.sel.code<-c("ANXI04")
#p$sumstats.sel.code<-c("ADHD05","ANXI04")
#p$sumstats.sel.code<-c("RISK02","RISK03","SCHI04","SUBJ01","TIRE01")
#p$sumstats.sel.code<-c("TIRE01")
p$sumstats.sel.code<-c("ADHD05","ALCD03","ANOR02","ANXI03","ANXI04","AUTI07","BIPO02", "DEPR05","DEPR08","EDUC03","HEAL01","INCO03","INSO02","NEUR02", "PTSD04","RISK02","RISK03","SCHI04","SUBJ01","TIRE01")
p$sumstats.sel<-p$sumstats[which(p$sumstats$code %in% p$sumstats.sel.code),]

#p$sumstats.sel$code_orig<-p$sumstats.sel$code
#p$sumstats.sel$code<-p$sumstats.sel$code.trait
p$sumstats.sel[,c("code","name","name.nice","name.nice.and_code", "year","n_case","n_control","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","dependent_variable.OLS","dependent_variable.linprob","se.logit","mungedpath")]
p$k.sel<-nrow(p$sumstats.sel)
#View(p$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")])

write.table(p$sumstats.sel[,c("code", "name.nice","year", "n_case","n_control","n_total","samplePrevalence","populationPrevalence", "reference_doi")], file = file.path(p$folderpath.workingDirectory,paste0(p$setup.code,".sumstatinfo.tsv")), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)

#View(p$sumstats.sel)



## ----GWAS sumstat munge-------------------------------------------------------------------
#test
#p$clOptions$task<-"munge"

p$filepath.lfgwas.sumstats<-file.path(p$folderpath.workingDirectory,paste0("lfGWAS.sumstats.",p$setup.code,".Rds"))
p$filepath.lfgwas.sumstats_full<-file.path(p$folderpath.workingDirectory,paste0("lfGWAS.sumstats_full.",p$setup.code,".Rds"))
p$filepath.munge.meta<-file.path(p$folderpath.workingDirectory,paste0("munge.meta.",p$setup.code,".Rds"))
p$munge<-c()

#source("../scripts/supermunge.R") #TEST OF NEW VERSION!

if(p$clOptions$task=="munge" & !all(file.exists(p$sumstats.sel$mungedpath))){
  
  #test
  #p$clOptions$task_argument<-"INSO02"
  if(!is.na(p$clOptions$task_argument)){
    #task argument sets the specific dataset to munge
    
    cat("\nSet to munge with arg",p$clOptions$task_argument,"\n")
    p$munge$filesToUse<-p$sumstats.sel$cleanedpath[which(p$sumstats.sel$code==p$clOptions$task_argument)]
    p$munge$traitNamesToUse<-p$sumstats.sel$code[which(p$sumstats.sel$code==p$clOptions$task_argument)]
    p$munge$NToUse<-p$sumstats.sel$n_total[which(p$sumstats.sel$code==p$clOptions$task_argument)]
    p$munge$OLSToUse<-p$sumstats.sel$dependent_variable.OLS[which(p$sumstats.sel$code==p$clOptions$task_argument)]
    p$munge$linprobToUse<-p$sumstats.sel$dependent_variable.linprob[which(p$sumstats.sel$code==p$clOptions$task_argument)]
    p$munge$se.logitToUse<-p$sumstats.sel$se.logit[which(p$sumstats.sel$code==p$clOptions$task_argument)]
    p$munge$propToUse<-(p$sumstats.sel$n_case/p$sumstats.sel$n_total)[which(p$sumstats.sel$code==p$clOptions$task_argument)]
    
  } else {
    #defaults
    p$munge$filesToUse<-p$sumstats.sel$cleanedpath
    p$munge$traitNamesToUse<-p$sumstats.sel$code
    p$munge$NToUse<-p$sumstats.sel$n_total
    p$munge$OLSToUse<-p$sumstats.sel$dependent_variable.OLS
    p$munge$linprobToUse<-p$sumstats.sel$dependent_variable.linprob
    p$munge$se.logitToUse<-p$sumstats.sel$se.logit
    p$munge$propToUse<-(p$sumstats.sel$n_case/p$sumstats.sel$n_total)
  }
  
  # cat("\nMunging with:\n")
  # print(p$munge$filesToUse)
  # print(p$filepath.SNPReference.1kg)
  # print(p$folderpath.data.mvLDSC.ld.1kg)
  # print(p$munge$traitNamesToUse)
  # print(p$munge$NToUse)
  # print(p$munge$OLSToUse)
  # print(p$munge$se.logitToUse)
  # print(p$munge$propToUse)
  # print(p$folderpath.data.sumstats.munged)
  
  #mask<-c(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,T,T,T,T)
  #munging with no filters applied
  p$munge$results <- supermunge(
            filePaths = p$munge$filesToUse,
            refFilePath = p$filepath.SNPReference.1kg,
            #refFilePath = p$filepath.SNPReference.hm3,
            #mask = mask,
            ldDirPath=p$folderpath.data.mvLDSC.ld.1kg,
            traitNames = p$munge$traitNamesToUse,
            setChangeEffectDirectionOnAlleleFlip = T, #T=same behaviour as genomic SEM
            imputeFromLD=T,
            #produceVariantTable = T,
            N = p$munge$NToUse,
            OLS = p$munge$OLSToUse,
            linprob=p$munge$linprobToUse,
            se.logit = p$munge$se.logitToUse,
            prop=p$munge$propToUse,
            pathDirOutput = p$folderpath.data.sumstats.munged,
            invertEffectDirectionOn = c("ANXI04") #ANXI4 has an inverted effect direction for some reason
              ) 
    
  p$munge$results$meta %>% gt()
  
  gtsave(data = (p$munge$results$meta %>% gt()), filename = file.path(p$folderpath.workingDirectory,"munge.meta.rtf"))
  
  saveRDS(object = p$munge$results$meta, file = p$filepath.munge.meta)
  
} else if(file.exists(p$filepath.munge.meta)){
  p$munge$results$meta<-readRDS(file=p$filepath.munge.meta)
  print("Read munging results metadata from previous munge result file.")
}

print(p$munge$results$meta)

if(file.exists(p$filepath.lfgwas.sumstats)) {
  # p$lfGWAS$sumstats<-readRDS(file=p$filepath.lfgwas.sumstats)
  # print("Read summary statistics for latent factor GWAS from file.")
  
  #compare with genomic SEM produced sumstats 
  # p$filepath.lfgwas.sumstats2<-file.path(p$folderpath.workingDirectory,paste0("lfGWAS.sumstats.setup4.Rds"))
  # p$lfGWAS$sumstats2<-readRDS(file=p$filepath.lfgwas.sumstats2)
  # 
  # View(p$lfGWAS$sumstats[which(p$lfGWAS$sumstats$SNP==c("rs1000000","rs10000013"))])
  
  # colBeta<-colnames(p$lfGWAS$sumstats)[grep("^BETA\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats))]
  # colSE<-colnames(p$lfGWAS$sumstats)[grep("^SE\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats))]
  # colBeta2<-colnames(p$lfGWAS$sumstats2)[grep("^BETA\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats2))]
  # colSE2<-colnames(p$lfGWAS$sumstats2)[grep("^SE\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats2))]
  # p$lfGWAS$sumstats<-as.data.frame(p$lfGWAS$sumstats)
  # for(iCol in 1:length(colBeta)){
  #   #iCol<-1
  #   print(colBeta[iCol])
  #   cBeta<-p$lfGWAS$sumstats[,colBeta[iCol]]
  #   cSE<-p$lfGWAS$sumstats[,colSE[iCol]]
  #   cBeta2<-p$lfGWAS$sumstats2[,colBeta2[iCol]]
  #   cSE2<-p$lfGWAS$sumstats2[,colSE2[iCol]]
  #   
  #   cat(mean(cBeta2,na.rm=T),"(",median(cSE2,na.rm=T),")\n")
  #   cat(mean(cBeta,na.rm=T),"(",median(cSE,na.rm=T),")\n")
  #   
  #   cat(max(abs(cBeta2),na.rm=T),"(",max(cSE2,na.rm=T),")\n")
  #   cat(max(abs(cBeta),na.rm=T),"(",max(cSE,na.rm=T),")\n")
  #   
  #   cat(min(abs(cBeta2),na.rm=T),"(",min(cSE2,na.rm=T),")\n")
  #   cat(min(abs(cBeta),na.rm=T),"(",min(cSE,na.rm=T),")\n")
  #   
  # }
  
  
  #Compare with old munged ANXI04
  #ssANXI04_old <- read.table(file.path(p$folderpath.data.sumstats.munged,"ANXI04_old.gz"),header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
  #setDT(ssANXI04_old)
  #setkeyv(ssANXI04_old, cols = c("SNP","BP","CHR"))
  #
  # ssANXI04 <- read.table(file.path(p$folderpath.data.sumstats.munged,"ANXI04.gz"),header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
  # setDT(ssANXI04)
  # setkeyv(ssANXI04, cols = c("SNP","BP","CHR"))
  #
  # ssANXI04.merged.snp<-ssANXI04[ssANXI04_old, on=c(SNP='SNP'), nomatch=0]
  # View(ssANXI04.merged.snp[1:100000,c("SNP","BP","CHR","A1","A2","i.A1","i.A2")])
  # sum(!(ssANXI04.merged.snp$A1==ssANXI04.merged.snp$i.A1 & ssANXI04.merged.snp$A2==ssANXI04.merged.snp$i.A2))
  # View(ssANXI04.merged.snp[1:100000,c("SNP","BP","CHR","EFFECT","i.EFFECT","Z","i.Z","P","i.P")])
  # View(ssANXI04.merged.snp[1:100000,c("SNP","BP","CHR","N","i.N")])

  # head(p$lfGWAS$sumstats)
  
} else if(p$clOptions$task=="munge"){
  p$munge$results<-supermunge(
    filePaths = p$sumstats.sel$mungedpath,
    refFilePath = p$filepath.SNPReference.1kg,
    traitNames = p$sumstats.sel$code,
    produceCompositeTable = T,
    process = F,
    standardiseEffectsToExposure = T,
    writeOutput = F,
    N = p$sumstats.sel$n_total,
    OLS=p$sumstats.sel$dependent_variable.OLS,
    linprob=p$sumstats.sel$dependent_variable.linprob,
    se.logit = p$sumstats.sel$se.logit,
    prop=(p$sumstats.sel$n_case/p$sumstats.sel$n_total),
    info.filter = 0.55
#    frq.filter = 0.001
    
    ) #remove rare variants only
    
  p$lfGWAS$sumstats<-p$munge$results$composite

  saveRDS(object = p$lfGWAS$sumstats, file = p$filepath.lfgwas.sumstats)
  print("Done preparing summary statistics for latent factor GWAS. The result should have been saved to a file.")
}

if(p$clOptions$task=="munge") quit(save = "no")



## ----GWAS sumstat generate noise----------------------------------------------------------
if(p$clOptions$task=="noise"){
  noiseScaffold <- read.table(p$filepath.SNPReference.1kg,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
  noiseTemplate <- read.table(p$sumstats.sel$mungedpath[which(p$sumstats.sel$code=="EDUC03")],header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
  
  noiseTemplateStats<-describe(noiseTemplate)
  
  noise<-noiseScaffold
  noise$Z<-rnorm(n=nrow(noise),mean=0,sd=1.485)
  noise$FRQ<-rnorm(n=nrow(noise),mean=0.770,sd=0.260)
  noise$FRQ<-ifelse(noise$FRQ>1, 1,noise$FRQ)
  noise$FRQ<-ifelse(noise$FRQ<0, 0,noise$FRQ)
  noise$SE<-rnorm(n=nrow(noise),mean=4.463e-03,sd=3.811e-03)
  noise$SE<-ifelse(noise$SE<1e-05, 1e-05,noise$SE)
  noise$P<-pnorm(q = noise$Z)
  #head(noise)
  
  #noiseStats<-describe(noise)
  
  nfilepath<-file.path(p$folderpath.data.sumstats.munged,"NOIS02")
  write.table(x = noise,file = nfilepath,sep="\t", quote = FALSE, row.names = F)
  nfilepath.gzip<-gzip(nfilepath)
  
  quit(save = "no")
}


## ----GWAS sumstat imputation--------------------------------------------------------------
if(p$clOptions$task=="impute"){
#Using RAISS
#https://gitlab.pasteur.fr/statistical-genetics/raiss
#https://statistical-genetics.pages.pasteur.fr/raiss/


#Using ImpG
#https://bogdan.dgsom.ucla.edu/pages/impg/
#https://github.com/huwenboshi/ImpG
  
#reticulate::import("raiss")
#py$raiss$

# system2(command = "awk",
#         args = c("BEGIN{mkdir -p ../data/gwas_sumstats/imputed/$1;} for( chr=1; chr<24; chr++){raiss --chrom $chr --gwas '$1' --ref-folder '../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM' --ld-folder '../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM' --zscore-folder '../data/gwas_sumstats/munged/$1.chr' --output-folder '../data/gwas_sumstats/imputed/$1' --l2-regularization 0.01 --eigen-threshold 0.05 --R2-threshold 0.3;}
# ",
# "project.setup2.sumstatinfo.tsv"))
  
  quit(save = "no")
}



## ----multivariate LD----------------------------------------------------------------------
print("***multivariate LD***")

p$filepath.mvLD<-file.path(p$folderpath.workingDirectory,paste0("mvLD.",p$setup.code,".Rds"))

if (file.exists(p$filepath.mvLD)) {
  print("Using existing covariance structures from previous LD computations.")
  p$mvLD<-readRDS(file=p$filepath.mvLD)
} else {
  print("Running (or reading ready intermediate results from) multivariate LD regression with different methods. This might take a while. If the procedure runs for too long you may want to abort the process.")
  
  cat("The current task is specified as:",p$clOptions$task)
  p$mvLD<-c()
  
  if(p$clOptions$task=="mvLD" || !file.exists(p$filepath.mvLD)){
    #run mvLDSC
    p$mvLD$covstruct.mvLDSC.hm3<-ldsc.mod(
      traits = p$sumstats.sel$mungedpath,
      sample.prev =  p$sumstats.sel$samplePrevalence,
      population.prev = p$sumstats.sel$populationPrevalence,
      trait.names = p$sumstats.sel$code,
      ld = p$folderpath.data.mvLDSC.ld.hm3,
      wld = p$folderpath.data.mvLDSC.ld.hm3,
      n.blocks = 200, #this was standard for the hm3 set of snps
      info.filter = 0.6,
      frq.filter = 0.01,
      mhc.filter = 37,
      N = p$sumstats.sel$n_total,
      #forceN = , p$sumstats.sel$forceN, # Consider this when some of the original N's may be untrustworthy (ANXI04!) - TODO - fix in supermunge
      ldsc.log = p$setup.code.date
      )

    p$mvLD$covstruct.mvLDSC.1kg<-ldsc.mod(
      traits = p$sumstats.sel$mungedpath,
      sample.prev =  p$sumstats.sel$samplePrevalence,
      population.prev = p$sumstats.sel$populationPrevalence,
      trait.names = p$sumstats.sel$code,
      ld = p$folderpath.data.mvLDSC.ld.1kg,
      wld = p$folderpath.data.mvLDSC.ld.1kg,
      n.blocks = 600, #a bit more here if it can support the larger reference panel used
      info.filter = 0.6,
      frq.filter = 0.01,
      mhc.filter = 37,
      N = p$sumstats.sel$n_total,
      #forceN = T, # Consider this when some of the original N's may be untrustworthy (ANXI04!) - TODO - fix in supermunge
      ldsc.log = p$setup.code.date
      )

  }
    
  # #Sys.sleep(time = 5)
  # #flip effect direction of trait
  # SModifier<-matrix(nrow = nrow(p$mvLD$covstruct.mvLDSC.1kg$S), ncol = ncol(p$mvLD$covstruct.mvLDSC.1kg$S), dimnames = list(
  #   colnames(p$mvLD$covstruct.mvLDSC.1kg$S),
  #   colnames(p$mvLD$covstruct.mvLDSC.1kg$S)
  #   )
  #   )
  # 
  # SModifier["ANXI04",]<--1
  # SModifier[,"ANXI04"]<--1
  # SModifier["ANXI04","ANXI04"]<-NA
  # # View(SModifier)
  # 
  # p$mvLD$covstruct.mvLDSC.1kg$S[!is.na(SModifier)]<-p$mvLD$covstruct.mvLDSC.1kg$S[!is.na(SModifier)]*SModifier[!is.na(SModifier)]
  
  #set the default mvLDSC object to use
  p$mvLD$covstruct.mvLDSC<-p$mvLD$covstruct.mvLDSC.1kg
  
  #Additional computations
  #saving the original S in case of smoothing experiments later stored in S
  p$mvLD$covstruct.mvLDSC$S.orig<-p$mvLD$covstruct.mvLDSC$S
  p$mvLD$covstruct.mvLDSC$S.smooth<-as.matrix((nearPD(p$mvLD$covstruct.mvLDSC$S, corr = FALSE))$mat)
  
  #save the mvLD output
  saveRDS(object = p$mvLD,file = p$filepath.mvLD)
  print("Multivariate LD correction is done now and the resulting covariance structure should have been saved to a file.")
  
}


if(p$clOptions$task=="mvLD"){
      quit(save = "no")
    }

#adding standard errors from the V matrices, unstandardised and standardised. I added these calculations to the modified ldsc.

#retrieve the standard errors of S (variances and covariances) from the diagonal of V (contains both).

# p$mvLD$covstruct.mvLDSC$S.SE<-matrix(0, p$k.sel, p$k.sel)
# rownames(p$mvLD$covstruct.mvLDSC$S.SE)<-colnames(p$mvLD$covstruct.mvLDSC$S)
# colnames(p$mvLD$covstruct.mvLDSC$S.SE)<-colnames(p$mvLD$covstruct.mvLDSC$S)
# p$mvLD$covstruct.mvLDSC$S_Stand.SE<-matrix(0, p$k.sel, p$k.sel)
# rownames(p$mvLD$covstruct.mvLDSC$S_Stand.SE)<-colnames(p$mvLD$covstruct.mvLDSC$S)
# colnames(p$mvLD$covstruct.mvLDSC$S_Stand.SE)<-colnames(p$mvLD$covstruct.mvLDSC$S)
# 
# p$mvLD$covstruct.mvLDSC$S.SE[lower.tri(p$mvLD$covstruct.mvLDSC$S.SE,diag=TRUE)] <-sqrt(diag(p$mvLD$covstruct.mvLDSC$V))
# p$mvLD$covstruct.mvLDSC$S.SE[upper.tri(p$mvLD$covstruct.mvLDSC$S.SE)]<-t(p$mvLD$covstruct.mvLDSC$S.SE)[upper.tri(p$mvLD$covstruct.mvLDSC$S.SE)]
# p$mvLD$covstruct.mvLDSC$S_Stand.SE[lower.tri(p$mvLD$covstruct.mvLDSC$S_Stand.SE,diag=TRUE)] <-sqrt(diag(p$mvLD$covstruct.mvLDSC$V_Stand))
# p$mvLD$covstruct.mvLDSC$S_Stand.SE[upper.tri(p$mvLD$covstruct.mvLDSC$S_Stand.SE)]<-t(p$mvLD$covstruct.mvLDSC$S_Stand.SE)[upper.tri(p$mvLD$covstruct.mvLDSC$S_Stand.SE)]


  
#add newly computed heritabilities to the selected summary statistics table
p$sumstats.sel$h2.liability_mvLDSC.1kg<-diag(p$mvLD$covstruct.mvLDSC.1kg$S[p$sumstats.sel$code,p$sumstats.sel$code])
p$sumstats.sel$h2.se.liability_mvLDSC.1kg<-diag(p$mvLD$covstruct.mvLDSC.1kg$S.SE[p$sumstats.sel$code,p$sumstats.sel$code])
p$sumstats.sel$h2.liability_mvLDSC.hm3<-diag(p$mvLD$covstruct.mvLDSC.hm3$S[p$sumstats.sel$code,p$sumstats.sel$code])
p$sumstats.sel$h2.se.liability_mvLDSC.hm3<-diag(p$mvLD$covstruct.mvLDSC.hm3$S.SE[p$sumstats.sel$code,p$sumstats.sel$code])


#View(p$sumstats.sel)

  






## ----EFA----------------------------------------------------------------------------------
print("***EFA***")
#saving efa results between runs as to always use the same randomised start for clustering for example
p$filepath.efa<-file.path(p$folderpath.workingDirectory,paste0("efa.",p$setup.code,".Rds"))
if (file.exists(p$filepath.efa)) {
  print("Using existing EFA results from previous run.")
  p$EFA<-readRDS(file=p$filepath.efa)
} else {

  #visualisation of up to max factors EFA models with varimax
  # for(iefa in 1:max(p$CFA$nFactors)){
  #   res <- psych::fa(r = abs(p$mvLD$covstruct.mvLDSC$S.smooth),nfactors = iefa, rotate = 'varimax', symmetric = T, warnings = T, fm='ols', max.iter = 1000)
  #   #print(res)
  #   print(res$loadings)
  #   print(res$fit)
  #   #print(kmeans(x = p$mvLD$covstruct.mvLDSC$S_Stand, centers = iefa, iter.max = 1000, nstart = 30))
  # }
  
  #visualise factors in a scree plot
  if(!p$clOptions$location=="cluster"){
      p$plots.efa.plot.scree<-fa.parallel(abs(p$mvLD$covstruct.mvLDSC$S.smooth), fa = "fa")
      png(filename = file.path(p$folderpath.plots,"efa.plot.scree.png"), width = 800, height = 500)
      fa.parallel(abs(p$mvLD$covstruct.mvLDSC$S.smooth), fa = "fa")
      dev.off()
  }
  
  
  #fit EFA models for each nFactor configuration
  print("Setting upp and computing new EFA results.")
  p$EFA<-c()
  p$EFA$PCA<-c()
  #iFactorConfiguration<-1
  for(iFactorConfiguration in 1:length(p$CFA$nFactors)) {
    cNFactors<-p$CFA$nFactors[iFactorConfiguration]
    p$EFA$PCA[[iFactorConfiguration]]<-eigen(x = abs(p$mvLD$covstruct.mvLDSC$S_Stand), symmetric = TRUE)
    rownames(p$EFA$PCA[[iFactorConfiguration]]$vectors)<-p$sumstats.sel.code
    p$EFA$PCA[[iFactorConfiguration]]$vector_values<-t(p$EFA$PCA[[iFactorConfiguration]]$vectors*p$EFA$PCA[[iFactorConfiguration]]$values)[,1:cNFactors] #has to be transposed as well!(?)
  
    p$EFA$fa.result.ORT[[iFactorConfiguration]] <- psych::fa(r = abs(p$mvLD$covstruct.mvLDSC$S.smooth),nfactors = cNFactors, rotate = 'varimax', symmetric = T, warnings = T, fm='ols', max.iter = 1000)
    
    p$EFA$fa.result.OBL[[iFactorConfiguration]] <- psych::fa(r = abs(p$mvLD$covstruct.mvLDSC$S.smooth),nfactors = cNFactors, rotate = 'oblimin', symmetric = T, warnings = T, fm='ols', max.iter = 1000)
    
    p$EFA$kmeans.result[[iFactorConfiguration]]<- kmeans(x = abs(p$mvLD$covstruct.mvLDSC$S_Stand), centers = cNFactors, iter.max = 1000, nstart = 30)
    #View(p$clustering$centers)
    #View(fitted(p$clustering))
    #View(abs(p$mvLD$covstruct.mvLDSC$S_Stand))
    #resid<-abs(p$mvLD$covstruct.mvLDSC$S_Stand)-fitted(p$clustering)
    #View(resid)
    p$EFA$kmeans.centerDistance[[iFactorConfiguration]]<-apply(X = abs(p$mvLD$covstruct.mvLDSC$S_Stand), MARGIN = 1, FUN = function(x){
      #test
      #x<-abs(p$mvLD$covstruct.mvLDSC$S_Stand)[1,]
      #cat("\nOBS:",x)
      ss<-apply(X=(abs(x)-p$EFA$kmeans.result[[iFactorConfiguration]]$centers)^2, FUN = sum, MARGIN = 1)
      #cat("\nSS:",ss)
      s<-sum(abs(x))
      #cat("\nS:",s)
      return (ss/s)
    })
    #transpose this to conform with indicator loading data frames
    p$EFA$kmeans.centerDistance[[iFactorConfiguration]]<-t(p$EFA$kmeans.centerDistance[[iFactorConfiguration]])
    rownames(p$EFA$kmeans.centerDistance[[iFactorConfiguration]])<-p$sumstats.sel.code
    
    
  cat("\nPCA vector values\n")
  print(p$EFA$PCA[[iFactorConfiguration]]$vector_values)
  cat("\nFA fa result, ORTHOGONAL rotation\n")
  print(p$EFA$fa.result.ORT[[iFactorConfiguration]])
  cat("\nFA fa result, OBLIQUE rotation\n")
  print(p$EFA$fa.result.OBL[[iFactorConfiguration]])
  #cat("\nFA factanal result\n")
  #print(p$EFA$factanal.result)
  cat("\nKmeans clustering centers\n")
  print(t(p$EFA$kmeans.result[[iFactorConfiguration]]$centers))
  cat("\nKmeans clustering residuals\n")
  print(p$EFA$kmeans.centerDistance[[iFactorConfiguration]])
    
  }
  
  
  
  saveRDS(object = p$EFA,file = p$filepath.efa)
}





## ----CFA indicator loading pattern creation-----------------------------------------------
print("***CFA indicator loading pattern creation***")

## CFA aditional settings
p$sumstats.sel$residualSizeLimitMax<-NA_real_
#p$sumstats.sel$residualSizeLimitMax[which(p$sumstats.sel$code=="ANXI03" | p$sumstats.sel$code=="DEPR05")]<-0.10
#(1/p$sumstats.sel$h2.se.liability_mvLDSC^2)/sum(1/p$sumstats.sel$h2.se.liability_mvLDSC^2)+0.01
p$CFA$nIndicators=length(p$sumstats.sel$code)

p$filepath.cfa<-file.path(p$folderpath.workingDirectory,paste0("cfa.",p$setup.code,".Rds"))
#p$filepath.cfa_converged_results<-file.path(p$folderpath.workingDirectory,paste0("cfa.",p$setup.code,".converged.txt"))
p$filepath.cfa.FactorConfiguration<-c()

if(!file.exists(p$filepath.cfa)){
  
  p$CFA$indicatorLoadingPatterns.PCA<-c()
  p$CFA$indicatorLoadingPatterns.fa.ORT<-c()
  p$CFA$indicatorLoadingPatterns.fa.OBL<-c()
  p$CFA$indicatorLoadingPatterns.kmeans<-c()
  #iFactorConfiguration<-1
  for(iFactorConfiguration in 1:length(p$CFA$nFactors)) {
  ##kmeans clustering based models
  p$CFA$indicatorLoadingPatterns.kmeans[[iFactorConfiguration]]<-semplate$generateIndicatorLoadingPatternsFromFactorLoadings(factorLoadings = t(p$EFA$kmeans.result[[iFactorConfiguration]]$centers), increment = 0.0005,forceOneIndicatorLoading = T)
  
  ##fa, orthogonal based models
  p$CFA$indicatorLoadingPatterns.fa.ORT[[iFactorConfiguration]]<-semplate$generateIndicatorLoadingPatternsFromFactorLoadings(factorLoadings = p$EFA$fa.result.ORT[[iFactorConfiguration]]$loadings, increment = 0.0005,forceOneIndicatorLoading = T)
  
  ##fa, oblique based models
  p$CFA$indicatorLoadingPatterns.fa.OBL[[iFactorConfiguration]]<-semplate$generateIndicatorLoadingPatternsFromFactorLoadings(factorLoadings = p$EFA$fa.result.OBL[[iFactorConfiguration]]$loadings, increment = 0.0005,forceOneIndicatorLoading = T)
    
  #PCA based models
  p$CFA$indicatorLoadingPatterns.PCA[[iFactorConfiguration]]<-semplate$generateIndicatorLoadingPatternsFromFactorLoadings(factorLoadings = p$EFA$PCA[[iFactorConfiguration]]$vector_values, increment = 0.0005,forceOneIndicatorLoading = T)
    
  p$CFA$sessionIndicatorLoadingPatterns[[iFactorConfiguration]]<-unique(rbind(p$CFA$indicatorLoadingPatterns.kmeans[[iFactorConfiguration]],p$CFA$indicatorLoadingPatterns.fa.ORT[[iFactorConfiguration]],p$CFA$indicatorLoadingPatterns.fa.OBL[[iFactorConfiguration]],p$CFA$indicatorLoadingPatterns.PCA[[iFactorConfiguration]]))
  
  }
  
  

}



## ----CFA model evaluation-----------------------------------------------------------------
print("***CFA model evaluation***")

#test
#p$clOptions$task<-"cfa"
#p$clOptions$task_argument<-"1"

match.row<-function(row_v,tomatch_df){
  #row_v<-lp
  #tomatch_df<-lplib
  toreturn<-vector()
  for(nTomatch in 1:nrow(tomatch_df)){
    #nTomatch<-2
    toreturn[nTomatch]<-all(row_v==tomatch_df[nTomatch,])
  }
  return(toreturn)
}

p$CFA$models<-data.frame(nModel=c(),code=c(),nFactors=c(),correlation=c(),estimator=c(), lModel=c(),lResults=c())

p$CFA$resultColumnNames<-c("chisq","df","p_chisq","AIC","CFI","SRMR")

if (file.exists(p$filepath.cfa)) {
  print("Using existing CFA results from previous run and appending to these if needed.")
  p$CFA<-readRDS(file=p$filepath.cfa)
} else {
  
  nModel<-0
  nFittingModelsFound<-0
  #compute total number of models across factor configurations
  totalNumberOfModels<-0
  nVariables<-ncol(p$mvLD$covstruct.mvLDSC$S)
  nUniqueCovariances<-nVariables*(nVariables+1)/2
  for(iFactorConfiguration in 1:length(p$CFA$nFactors)) {
    sessionPatternLength<-nrow(p$CFA$sessionIndicatorLoadingPatterns[[iFactorConfiguration]])
    totalNumberOfModels<-totalNumberOfModels+sessionPatternLength*length(p$CFA$correlation)*length(p$CFA$estimator)
  }
  
  for(iFactorConfiguration in 1:length(p$CFA$nFactors)) {
    #test
    #iFactorConfiguration<-1
    
    if(p$clOptions$task=="cfa" & length(p$clOptions$task_argument)>0 & p$clOptions$task_argument!=iFactorConfiguration) next
    
    p$filepath.cfa.FactorConfiguration[iFactorConfiguration]<-file.path(p$folderpath.workingDirectory,paste0("cfa.",p$setup.code,"_",iFactorConfiguration,".Rds"))
    
    if(file.exists(p$filepath.cfa.FactorConfiguration[iFactorConfiguration])){
      storedFactorConfigurationModelResults <- readRDS(file=p$filepath.cfa.FactorConfiguration[iFactorConfiguration])
      p$CFA$models <- rbind(p$CFA$models,storedFactorConfigurationModelResults)
      nModel<-nrow(p$CFA$models)
      next
    }
    
    
    if(p$clOptions$task!="cfa") next #do not evaluate models if not in cfa task
    
    sessionPatternLength<-nrow(p$CFA$sessionIndicatorLoadingPatterns[[iFactorConfiguration]])
    
    cat("\nAnalysing",sessionPatternLength, "patterns for this factor configuration...\n")
    for(nSessioPattern in 1:sessionPatternLength){
      #test
      #nSessioPattern<-1
      
      for(iCorrelationConfiguration in 1:length(p$CFA$correlation)){
        #iCorrelationConfiguration<-1
        cCorrelation<-p$CFA$correlation[iCorrelationConfiguration]
        
        # #avoid correlated factors if above 8 factors
        # if(p$CFA$nFactors[iFactorConfiguration]>8 & cCorrelation=="COR") next
        # 
        # #avoid oblique factors if below 9 factors
        # if(p$CFA$nFactors[iFactorConfiguration]<9 & cCorrelation=="OBL") next
        # 
        # #avoid orthogonal factors if below 7 factors
        # if(p$CFA$nFactors[iFactorConfiguration]<7 & cCorrelation=="ORT") next
        
        for(iParameterEstimator in 1:length(p$CFA$estimator)){
          #iParameterEstimator<-1
          cEstimator<-p$CFA$estimator[iParameterEstimator]
          
          #new model row
          nModel<-nModel+1
          
          if(!is.null(p$CFA$models.selected)){
            if(!nModel %in% p$CFA$models.selected$nModel) next
          }
          
          #set specific nmodel for reevaluation
          #nModel<-35
          #nSessioPattern<-35
          
          #init columns
          p$CFA$models[nModel,]<-NA
          
          #nModel
          p$CFA$models[nModel,c("nModel")]<-nModel
          
          #record code
          p$CFA$models[nModel,c("code")]<-paste0("M",p$CFA$nIndicators,"_",p$CFA$nFactors[iFactorConfiguration],
          "_",nSessioPattern,
          ".",cCorrelation,
          ".",cEstimator
          )
          
          #record nFactors
          p$CFA$models[nModel,c("nFactors")]<-p$CFA$nFactors[iFactorConfiguration]
          
          #record correlation
          p$CFA$models[nModel,c("correlation")]<-cCorrelation
          
          #record estimator
          p$CFA$models[nModel,c("estimator")]<-cEstimator
          
          #p$CFA$models$totalBitValue[nModel]<-NA #because otherwise the later assignment will crash
          p$CFA$models[nModel,c("loading_pattern","loading_pattern_pca","loading_pattern_fa","loading_pattern_kmeans")]<-NA
          p$CFA$models[nModel,c("lModel")]<-NA_character_
          p$CFA$models[nModel,c("gsemResults")]<-NA
          p$CFA$models[nModel,p$CFA$resultColumnNames]<-NA
          
          #loading pattern
          lp<-p$CFA$sessionIndicatorLoadingPatterns[[iFactorConfiguration]][nSessioPattern,]
          p$CFA$models[[nModel,c("loading_pattern")]]<-list(lp)
          cIndicatorLoadings<-matrix(data = lp, nrow = p$CFA$nIndicators, ncol = p$CFA$nFactors[iFactorConfiguration]) 
          row.names(cIndicatorLoadings)<-p$sumstats.sel$code
          
          p$CFA$models[nModel,c("loading_pattern_pca")]<-any(match.row(lp,p$CFA$indicatorLoadingPatterns.PCA[[iFactorConfiguration]]))
          p$CFA$models[nModel,c("loading_pattern_fa.ORT")]<-any(match.row(lp,p$CFA$indicatorLoadingPatterns.fa.ORT[[iFactorConfiguration]]))
          p$CFA$models[nModel,c("loading_pattern_fa.OBL")]<-any(match.row(lp,p$CFA$indicatorLoadingPatterns.fa.OBL[[iFactorConfiguration]]))
          p$CFA$models[nModel,c("loading_pattern_kmeans")]<-any(match.row(lp,p$CFA$indicatorLoadingPatterns.kmeans[[iFactorConfiguration]]))
          
          #further filter rules
          indicatorsLoadedOnFactors <- apply(cIndicatorLoadings, 1, FUN = any)
          factorsHasIndicatorsLoaded <- apply(cIndicatorLoadings, 2, FUN = any)
          ##check factor configuration is of an identified model
          if(cCorrelation!="ORT"){
            isIdentified <- (nVariables + sum(cIndicatorLoadings) + factorial(nVariables)/(factorial(2)*factorial(nVariables-2))) <= nUniqueCovariances #real condition
            p$CFA$models[nModel,c("identified")]<-isIdentified
             isIdentified <- (nVariables + sum(cIndicatorLoadings) + (factorial(nVariables)/(factorial(2)*factorial(nVariables-2)))/2) <= nUniqueCovariances #relax this condition in case my calculation is wrong
          } else { #ORT
            isIdentified<-nVariables + sum(cIndicatorLoadings) <= nUniqueCovariances
          }
          
          #allow evaluation if filter rules are met
          if(all(indicatorsLoadedOnFactors) & all(factorsHasIndicatorsLoaded)){ #& isIdentified
            
            #generate lavaan model
            p$CFA$models[nModel,c("lModel")]<- semplate$generateLavaanCFAModel(
              allow_loading.table.indicator_factor = cIndicatorLoadings,
              #indicatorArgs = p$sumstats.sel[,c("code","residualSizeLimitMax")],
              #universalResidualLimitMin = 0.0001,
              orthogonal = (cCorrelation=="ORT"),
              universalCorrelationLimitMax = ifelse((cCorrelation=="OBL"),0.3,NA)
              )
            
            #evaluate lavaan model in GenomicSEM
            cat("\n\n#Found fitting=",nFittingModelsFound,",\tevaluating new model:\t",nModel,"/",totalNumberOfModels,"\n", p$CFA$models[nModel,c("code")],"\n")
            cModelResults = tryCatch(
              usermodel.mod(covstruc = p$mvLD$covstruct.mvLDSC,
                model = p$CFA$models[nModel,c("lModel")],
                estimation = cEstimator,
                fix_resid = F,
                CFIcalc = F #ifelse(is.null(p$CFA$models.selected),F,T) #set this to true for CFI evaluation
                ), error= function(e) e
              )
            
            if(!inherits(cModelResults, "try-error") & !is.null(cModelResults$modelfit)){
              if(nrow(cModelResults$modelfit)>0 && any(p$CFA$resultColumnNames %in% colnames(cModelResults$modelfit))) {
                print(cModelResults$modelfit)
                #record results even though not fitting
                p$CFA$models[nModel,c("gsemResults")]<-NA
                p$CFA$models[[nModel,c("gsemResults")]]<-list(cModelResults)
                cRescolnames<-intersect(p$CFA$resultColumnNames,colnames(cModelResults$modelfit))
                p$CFA$models[nModel,cRescolnames]<-cModelResults$modelfit[1,cRescolnames]
                if(is.numeric(cModelResults$modelfit$chisq)){
                  #This is considered a fitting model
                  nFittingModelsFound<-nFittingModelsFound+1
                  cat("\nFITTING!:",p$CFA$models[nModel,c("code")],"\n")
                }
              } else {
                cat("\nThe model did not yield correct results.")
              }
            } else {
              cat("\nThe model did not converge.")
            }
            
          } else {warning("\nThe pattern configuration was deemed to yield an non-identified or otherwise not plausible model!")} #evaluation block
        } #for(iParameterEstimator in 1:length(p$CFA$estimator))
      } #for(iCorrelationConfiguration in 1:length(p$CFA$correlation))
    } #for(nSessioPattern in 1:sessionPatternLength)
  
    
    saveRDS(object = p$CFA$models,file = p$filepath.cfa.FactorConfiguration[iFactorConfiguration])
  print("Intermediate CFA results for this factor configuration are now done and the result should have been saved into a file.")
  
  if(p$clOptions$task=="cfa" & p$clOptions$task_argument==iFactorConfiguration){quit(save = "no")}
  
  }
  
  rownames(p$CFA$models)<-p$CFA$models$code
    
  saveRDS(object = p$CFA,file = p$filepath.cfa)
  print("CFA for this session is now done and the result should have been saved into a file.")
}



#View(p$CFA$models)
#p$CFA$models$lModel[which(p$CFA$models$nModel==57 & p$CFA$models$code=="M3-17.ML._-1107563521_389110_0_0_0_0_0_0_")]
if(p$clOptions$task=="cfa"){quit(save = "no")}



## ----CFA select---------------------------------------------------------------------------
print("***CFA select***")
#View(p$CFA$models[which(p$CFA$models$AIC<10000 & p$CFA$models$correlation=="ORT"),c("code","nModel","identified", "loading_pattern","loading_pattern_pca","loading_pattern_fa.ORT","loading_pattern_fa.OBL","loading_pattern_kmeans",p$CFA$resultColumnNames)])


p$CFA$models.selected<-p$CFA$models[c("M20_6_172.COR.ML","M20_6_30.ORT.ML"),] #this is a manual setting after choosing the best fitting models from the previous results
#View(p$CFA$models.selected)




## ----pre-analysis-------------------------------------------------------------------------

print("***Pre-analysis***")

if(p$clOptions$task=="pre"){
  p$sumstats.sel.pre.code<-c("EDUC03","NEUR02","RISK02","RISK03","SCHI04")
  p$sumstats.sel.pre<-p$sumstats[which(p$sumstats$code %in% p$sumstats.sel.pre.code),]
  
  #export modified munged files
  
  p$pre<-c()
  p$pre$filepath.mvLD<-file.path(p$folderpath.workingDirectory,paste0("mvLD.pre.",p$setup.code,".Rds"))
  p$pre$filepath.models<-file.path(p$folderpath.workingDirectory,paste0("cfa.pre.",p$setup.code,".Rds"))
  
  if (file.exists(p$pre$filepath.mvLD)) {
    print("Using existing covariance structures from previous LD computations.")
    p$pre$mvLD<-readRDS(file=p$pre$filepath.mvLD)
  } else {
      
    p$pre$mvLD<-c()
    p$pre$mvLD$covstructs<-c()
    
    p$pre$mvLD$covstructs[[1]]<-ldsc.mod(
        traits = p$sumstats.sel.pre$mungedpath,
        sample.prev =  p$sumstats.sel.pre$samplePrevalence,
        population.prev = p$sumstats.sel.pre$populationPrevalence,
        trait.names = p$sumstats.sel.pre$code,
        ld = p$folderpath.data.mvLDSC.ld.1kg,
        wld = p$folderpath.data.mvLDSC.ld.1kg,
        n.blocks = 600, #a bit more here if it can support the larger reference panel used
        info.filter = 0.6,
        frq.filter = 0.01,
        mhc.filter = 37,
        N = p$sumstats.sel.pre$n_total,
        #forceN = T, # Consider this when some of the original N's may be untrustworthy (ANXI04!) - TODO - fix in supermunge
        ldsc.log = p$setup.code.date
        )
    
    for(iChr in 1:11){
      p$pre$mvLD$covstructs[[iChr+1]]<-ldsc.mod(
        traits = p$sumstats.sel.pre$mungedpath,
        sample.prev =  p$sumstats.sel.pre$samplePrevalence,
        population.prev = p$sumstats.sel.pre$populationPrevalence,
        trait.names = p$sumstats.sel.pre$code,
        ld = p$folderpath.data.mvLDSC.ld.1kg,
        wld = p$folderpath.data.mvLDSC.ld.1kg,
        n.blocks = 600, #a bit more here if it can support the larger reference panel used
        info.filter = 0.6,
        frq.filter = 0.01,
        mhc.filter = 37,
        N = p$sumstats.sel.pre$n_total,
        #forceN = T, # Consider this when some of the original N's may be untrustworthy (ANXI04!) - TODO - fix in supermunge
        ldsc.log = p$setup.code.date,
        leave.chr=c(2*iChr-1,2*iChr)
        )
    
    }
      
    #save the mvLD output
    saveRDS(object = p$pre$mvLD,file = p$pre$filepath.mvLD)
    print("Multivariate LD correction is done now and the resulting covariance structure should have been saved to a file.")
      
  }
  
  
  #analyse the differences in genetic correlation matrices
  p$pre$mvLD$covstructs.leaveoneout<-p$pre$mvLD$covstructs[2:11]
  length(p$pre$mvLD$covstructs.leaveoneout)
  p$pre$mvLD$covstructs[[1]]$S
  p$pre$mvLD$covstructs[[1]]$S.SE
  
  p$pre$mvLD$covstructs[[2]]$S
  p$pre$mvLD$covstructs[[2]]$S.SE
  
  p$pre$mvLD$covstructs.leaveoneout.S<-matrix(NA,nrow = nrow(p$pre$mvLD$covstructs[[1]]$S), ncol = ncol(p$pre$mvLD$covstructs[[1]]$S))
  
  length(p$pre$mvLD$covstructs.leaveoneout.S)
  for(i in 1:length(p$pre$mvLD$covstructs.leaveoneout.S)){
    #i<-1
    p$pre$mvLD$covstructs.leaveoneout.S[[i]]<-mean(unlist(lapply(p$pre$mvLD$covstructs.leaveoneout,function(x){x$S[i]})))
  }
  p$pre$mvLD$covstructs.leaveoneout.S<-matrix(as.numeric(p$pre$mvLD$covstructs.leaveoneout.S),nrow = nrow(p$pre$mvLD$covstructs[[1]]$S), ncol = ncol(p$pre$mvLD$covstructs[[1]]$S))
  p$pre$mvLD$covstructs.leaveoneout.S
  p$pre$mvLD$covstructs[[1]]$S
  
  as.matrix( p$pre$mvLD$covstructs.leaveoneout.S) - as.matrix(p$pre$mvLD$covstructs[[1]]$S)
  mean(as.matrix( p$pre$mvLD$covstructs.leaveoneout.S) - as.matrix(p$pre$mvLD$covstructs[[1]]$S))
  
  #EFA
  p$pre$fa.result.OBL <- psych::fa(r = abs(p$pre$mvLD$covstructs[[1]]$S),nfactors = 2, rotate = 'oblimin', symmetric = T, warnings = T, fm='ols', max.iter = 1000)
  print(p$pre$fa.result.OBL)
  p$pre$PCA<-eigen(x = abs(p$pre$mvLD$covstructs[[1]]$S), symmetric = TRUE)
  print(p$pre$PCA)
  p$pre$PCA$vector_values<-t(p$pre$PCA$vectors*p$pre$PCA$values)[,1:2]
  
  #CFA
  p$pre$resultColumnNames<-c("chisq","df","p_chisq","AIC","CFI","SRMR")
  p$pre$indicatorLoadingPatterns<-semplate$generateIndicatorLoadingPatternsFromFactorLoadings(factorLoadings = p$pre$fa.result.OBL$loadings, increment = 0.0005,forceOneIndicatorLoading = T)
  #p$pre$indicatorLoadingPatterns<-semplate$generateIndicatorLoadingPatternsFromFactorLoadings(factorLoadings = p$pre$PCA$vector_values, increment = 0.0005,forceOneIndicatorLoading = T)
  p$pre$models<-data.frame(nModel=c(),code=c(),nFactors=c(),correlation=c(),estimator=c(), lModel=c(),lResults=c(),gsemResults=c())
  
  
  if(file.exists(p$pre$filepath.models)) {
    print("Using existing models.")
    p$pre$models<-readRDS(file=p$pre$filepath.models)
  } else {
    for(iPat in 1:nrow(p$pre$indicatorLoadingPatterns)){
      #iPat<-1
      p$pre$models[iPat,]<-NA
      p$pre$models[iPat,c("nModel")]<-iPat
      
      #record code
      p$pre$models[iPat,c("code")]<-paste0("MPRE_2",
      "_",iPat,
      ".COR",
      ".ML"
      )
      
      p$pre$models[iPat,c("nFactors")]<-2
      p$pre$models[iPat,c("correlation")]<-"COR"
      p$pre$models[iPat,c("estimator")]<-"ML"
            
      lp<-p$pre$indicatorLoadingPatterns[iPat,]
      p$pre$models[iPat,c("loading_pattern")]<-NA
      p$pre$models[[iPat,c("loading_pattern")]]<-list(lp)
      cIndicatorLoadings<-matrix(data = lp, nrow = nrow(p$pre$mvLD$covstructs[[1]]$S), ncol = 2) 
      row.names(cIndicatorLoadings)<-row.names(p$pre$mvLD$covstructs[[1]]$S)
      p$pre$models[iPat,c("lModel")]<-semplate$generateLavaanCFAModel(allow_loading.table.indicator_factor = cIndicatorLoadings)
      
      cModelResults = tryCatch(
        usermodel.mod(covstruc = p$pre$mvLD$covstructs[[1]],
          model = p$pre$models[iPat,c("lModel")],
          estimation = "ML",
          fix_resid = F,
          CFIcalc = F #ifelse(is.null(p$CFA$models.selected),F,T) #set this to true for CFI evaluation
          ), error= function(e) e
        )
      
      if(!inherits(cModelResults, "try-error") & !is.null(cModelResults$modelfit)){
                if(nrow(cModelResults$modelfit)>0 && any(p$pre$resultColumnNames %in% colnames(cModelResults$modelfit))) {
                  print(cModelResults$modelfit)
                  #record results even though not fitting
                  p$pre$models[iPat,c("gsemResults")]<-NA
                  p$pre$models[[iPat,c("gsemResults")]]<-list(cModelResults)
                  cRescolnames<-intersect(p$pre$resultColumnNames,colnames(cModelResults$modelfit))
                  p$pre$models[iPat,cRescolnames]<-cModelResults$modelfit[1,cRescolnames]
                  if(is.numeric(cModelResults$modelfit$chisq)){
                    #This is considered a fitting model
                    cat("\nFITTING!:",p$pre$models[iPat,c("code")],"\n")
                  }
                } else {
                  cat("\nThe model did not yield correct results.")
                }
              } else {
                cat("\nThe model did not converge.")
              }
      
    } #for
    
    rownames(p$pre$models)<-p$pre$models$code
    saveRDS(object = p$pre$models,file = p$pre$filepath.models)
    print("PRE CFA is done now and the resultins should have been saved to a file.")
  }
  
  #setup results to be run in the later part of the script as if they were coming out of the standard script analysis
  p$sumstats.sel<-p$sumstats.sel.pre
  p$mvLD<-c()
  p$mvLD$covstruct.mvLDSC<-p$pre$mvLD$covstructs[[1]]
  
  #Additional computations
  #saving the original S in case of smoothing experiments later stored in S
  p$mvLD$covstruct.mvLDSC$S.orig<-p$mvLD$covstruct.mvLDSC$S
  p$mvLD$covstruct.mvLDSC$S.smooth<-as.matrix((nearPD(p$mvLD$covstruct.mvLDSC$S, corr = FALSE))$mat)
  
  
  p$CFA$models.selected<-p$pre$models[c("MPRE_2_1.COR.ML"),]
  
}

if(p$clOptions$task=="pre"){
  p$clOptions$task<-"lfgwas"
  }




## ----Process CFA results------------------------------------------------------------------

p$CFA$models.selected$parsedGsemResults<-NA
p$CFA$models.selected$lModel.fixed<-NA
for(iSelected in 1:nrow(p$CFA$models.selected)){
  #test
  #iSelected<-1
  p$CFA$models.selected[[iSelected,c("parsedGsemResults")]]<-list(semplate$parseGenomicSEMResultAsMatrices(p$CFA$models.selected[iSelected,]$gsemResults[[1]][[1]]$results))
  
  # model explained variance is calculated by the parse function from the standardised residuals
  # vars<-vector()
  # for(iManifest in 1:nrow(p$mvLD$covstruct.mvLDSC$S)){
  #   vars[iManifest]<-p$mvLD$covstruct.mvLDSC$S[iManifest,iManifest]
  # }
  # 1-(p$CFA$models.selected[[iSelected,c("parsedGsemResults")]][[1]]$residualVaraiances.matrix/vars)
  
  
  #using the standardised pattern coefficients to fix the model
  #summary(p$CFA$model.bestFitting$gsemResults[[1]]$lresults, standardized=T)
  p$CFA$models.selected[[iSelected,c("parsedGsemResults")]][[1]]$patternCoefficientsSTDGenotype.matrix
  
  #prepare fixed lavaan model
      
  cIndicatorLoadings.loadingPattern<-matrix(
    data = p$CFA$models.selected[iSelected,]$loading_pattern[[1]][[1]],
    ncol = ncol(p$CFA$models.selected[iSelected,]$parsedGsemResults[[1]][[1]]$patternCoefficientsSTDGenotype.matrix),
    nrow = nrow(p$CFA$models.selected[iSelected,]$parsedGsemResults[[1]][[1]]$patternCoefficientsSTDGenotype.matrix)) 
  row.names(cIndicatorLoadings.loadingPattern)<-p$sumstats.sel$code
  cIndicatorLoadings.result<-!is.na(p$CFA$models.selected[iSelected,]$parsedGsemResults[[1]][[1]]$patternCoefficientsSTDGenotype.matrix)
  
  if(!all(cIndicatorLoadings.loadingPattern==cIndicatorLoadings.result)) stop("Difference in provided and result loading patterns detected!")
  
  cIndicatorLoadings<-cIndicatorLoadings.loadingPattern
  
  #test of free SEM parameters rather than fixed (pretend it is the fixed model)
  # p$CFA$models.selected[[iSelected,c("lModel.fixed")]]<-semplate$generateLavaanCFAModel(
  #   allow_loading.table.indicator_factor = cIndicatorLoadings,
  #   universalResidualLimitMin = NA,
  #   orthogonal = (p$CFA$models.selected[iSelected,]$correlation[[1]][[1]]=="ORT"))

  #original fixed setup
  p$CFA$models.selected[[iSelected,c("lModel.fixed")]]<-semplate$generateLavaanCFAModel(
    allow_loading.table.indicator_factor = cIndicatorLoadings,
    fix_loading.table.indicator_factor = p$CFA$models.selected[iSelected,]$parsedGsemResults[[1]][[1]]$patternCoefficientsSTDGenotype.matrix, fixResidualVariance_v = p$CFA$models.selected[iSelected,]$parsedGsemResults[[1]][[1]]$residualVaraiancesSTDGenotype.matrix,
    fix_correlation.table.factor_factor = p$CFA$models.selected[iSelected,]$parsedGsemResults[[1]][[1]]$covariancesSTDGenotype.matrix,
    orthogonal = (p$CFA$models.selected[iSelected,]$correlation[[1]][[1]]=="ORT"))

}





## ----naive variant effect meta-analysis---------------------------------------------------
#test
#p$clOptions$task<-"nmeta"
#p$clOptions$task_argument<-"M20_6_172.COR.ML:1"

cat("\n***Naive variant effect meta-analysis***\n")

if(p$clOptions$task=="nmeta"){
  
  p$nmeta<-c()
  
  p$lfGWAS$sumstats<-readRDS(file=p$filepath.lfgwas.sumstats)
  setDT(p$lfGWAS$sumstats)
  print("Read summary statistics for latent factor GWAS from file.")
  head(p$lfGWAS$sumstats)
  
  colBeta<-colnames(p$lfGWAS$sumstats)[grep("^BETA\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats))]
  colSE<-colnames(p$lfGWAS$sumstats)[grep("^SE\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats))]
  
  #sumstats QC
  ## effects
  rsum<-rowSums(abs(p$lfGWAS$sumstats[,..colBeta]),na.rm = T)
  rem<-rsum<(length(colBeta)*1e-19)
  p$lfGWAS$sumstats<-p$lfGWAS$sumstats[!rem,]
  rm("rsum")
  cat("\nRemoved ",sum(rem), " variants as part of QC before nmeta.\n")
  
  SESNP<-p$lfGWAS$sumstats[,..colSE]
  for(i in colSE)
    SESNP[is.na(get(i)), (i):=1] #assign the value 1 to NA SNP SE
  POPVARSNP<-2*p$lfGWAS$sumstats$MAF*(1-p$lfGWAS$sumstats$MAF) #2pq according to the genomic SEM publication
  
  #S<-as.matrix(p$mvLD$covstruct.mvLDSC$S)
  #V<-as.matrix(p$mvLD$covstruct.mvLDSC$V)
  I<-as.matrix(p$mvLD$covstruct.mvLDSC$I)
  #diag(I)<-ifelse(diag(I)<= 1, 1, diag(I)) #as in genomic SEM implementation. removed this because we sum all intercepts rather than only using the diagonal.
  I_diag<-diag(I)
  mI<-diag(nrow(I))
  for (x in 1:ncol(SESNP)) {
    for (y in 1:ncol(SESNP)) {
            mI[[x,y]]<-I[[x,y]]*I[[x,x]]*I[[y,y]]
        }
  }
  
  mI.sum<-colSums(mI,na.rm = T)
  
  BETASNP<-p$lfGWAS$sumstats[,..colBeta]
  #BETAVECTORS<-data.table()
  #adjust with the full ldsc intercept
  for (iCol in 1:ncol(BETASNP)){
    #iCol<-3
    col<-colnames(BETASNP)[iCol]
    set(x = BETASNP, j = col, value = BETASNP[[col]]/ifelse(mI.sum[iCol]<1,1,mI.sum[iCol]))
    #exprimental
    # if(iCol==1)
    #   set(x = BETAVECTORS, j = beta, value = BETASNP[[col]])
    # else
    #   set(x = BETAVECTORS, j = beta, value = c(BETAVECTORS[[,"beta"]],BETASNP[[col]]))
    #BETASNP[, (col):=(..col)]
    #BETASNP[, (col):=(col)(I_diag[[iCol]])]
  }
  
  #k<-apply(X = BETASNP,MARGIN = 1,FUN = function(x){sum(!is.na(x))})
  
  if(length(grep(pattern = "\\:",x = p$clOptions$task_argument))>0){
    p$nmeta$cModel<-strsplit(x = p$clOptions$task_argument, split = ':')[[1]][1]
    p$nmeta$cFactor<-strsplit(x = p$clOptions$task_argument, split = ':')[[1]][2]
    } else {
      p$nmeta$cModel<-p$clOptions$task_argument
      p$nmeta$cFactor<-NULL
    }
  
  for(iModel in 1:nrow(p$CFA$models.selected)){
    #test
    #iModel<-1
    
    cModel<-p$CFA$models.selected[iModel,]
    if(!is.null(p$nmeta$cModel) & !is.na(p$nmeta$cModel) & p$nmeta$cModel!=cModel$code) next
    cat("\nModel:  ",cModel$code)
    
    for(iFactor in 1:cModel$nFactors){
      #test
      #iFactor<-1
      if(is.null(p$nmeta$cFactor)){} else if(p$nmeta$cFactor!=iFactor) next
      factorFilepath<-file.path(p$folderpath.workingDirectory,paste0("nmeta.",cModel$code,".F",iFactor,".",p$setup.code,".Rds"))
      if(file.exists(factorFilepath)){
        modelGWAS<-readRDS(factorFilepath)
        iSNPStart<-nrow(modelGWAS[!is.na(BETA)])+1 #THIS DOES NOT DETECT THE RIGHT START POSITION!
        cat("\nContinuing nmeta from variant with index ",iSNPStart,"\n")
      } else {
        modelGWAS<-p$lfGWAS$sumstats[,c("SNP","CHR","BP","MAF","A1","A2")]
        iSNPStart<-1
      }
      cat("\nFactor:  ",iFactor,"\n")
      W<-cModel$parsedGsemResults[[1]][[1]]$relativeVarianceExplainedPerFactor[,iFactor]*sign(cModel$parsedGsemResults[[1]][[1]]$patternCoefficients.matrix[,iFactor]) #only weighted by the factor loadings here, rather than on the dataset variances
      W2<-W^2
      sqrt.W<-sqrt(abs(W))
      signed.sqrt.W<-sign(W)*sqrt.W
      k_W<-sum(!is.na(W))
      k_tot<-ncol(SESNP)
      #W0<-W
      #W0[is.na(W0)]<-0
      nSNP<-nrow(BETASNP)
      
      #mCovarTemplate<-setDT(as.data.frame(diag(k_tot)))
      #C_old<-sum(W,na.rm = T)-sum(W2,na.rm = T)/sum(W,na.rm = T)
      
      
      for(iSNP in iSNPStart:nSNP){
        #test
        #iSNP<-1L
        #iSNP<-8952400L
        #if(!iSNP %% 10 == 0) next  #for testing
        if(iSNP %% 100000 == 0){
          cat("#SNP ",iSNP, p$lfGWAS$sumstats[[iSNP,c("SNP")]],"\n")

          saveRDS(object = modelGWAS,file = factorFilepath)
          #cat("\nSaved intermediate results!\n")
          }
        
        mCovar<-diag(k_tot)
        #mI<-diag(nrow(I))
        #set snp variances adjusted for ldsc intercepts and sample overlap (through the ldsc intercepts - see genomic sem publication)
        for (x in 1:k_tot) {
          for (y in 1:k_tot) {
            #mCovar[x,y]<-(SESNP[[iSNP,y]]*SESNP[[iSNP,x]]*mI[[x,y]]*POPVARSNP[[iSNP]]^2)
            #mCovar[[x,y]]<-(SESNP[[iSNP,y]]*SESNP[[iSNP,x]]*mI[[x,y]]*POPVARSNP[[iSNP]]^2)
            mCovar[[x,y]]<-(signed.sqrt.W[x]*signed.sqrt.W[y]*SESNP[[iSNP,y]]*SESNP[[iSNP,x]]*mI[[x,y]]*POPVARSNP[[iSNP]]^2)
          }
        }
        
        abs.mCovar.sum<-abs(colSums(mCovar,na.rm = T)) #The individual variance components of each separate dataset

        #BETASNP.wmean<-weighted.mean(x = BETASNP[iSNP,], w = W0, na.rm = T)
        #BETASNP.wmean<-weighted.mean(x = BETASNP[iSNP,], w = W0/sqrt(abs.mCovar.sum), na.rm = T) #slightly corrected by the covariances
        # C<-sum(W/abs.mCovar.sum,na.rm = T)-sum(W2/(abs.mCovar.sum^2),na.rm = T)/sum(W/abs.mCovar.sum,na.rm = T)
        # Q_part<-((BETASNP[iSNP,]-BETASNP.wmean)^2)/abs.mCovar.sum
        # Q_part[is.na(Q_part)]<-0
        #Q<-sum(Q_part)
        #Q<-sum(((BETASNP[iSNP,]-BETASNP.wmean)^2)/abs.mCovar.sum, na.rm = T)
        #k<-sum(!is.na(W*BETASNP[iSNP,]))
        # T2<-max( c(
        #   (W0*Q - (k_W-1))/C,
        #   0
        #   ))
        # T2_part<-(W*Q_part - 1 + 1/k_W)/C
        # T2_part[T2_part<0]<-0

        #W_FULL<-W/(abs.mCovar.sum+T2) #as in a standard random effects meta-analysis
        #W_FULL<-W/(abs.mCovar.sum+T2_part)
        #W_FULL<-W/abs.mCovar.sum #fixed-effect meta-analysis
        #W_FULL<-W/(1+abs.mCovar.sum) #fixed-effect meta-analysis
        W_FULL<-W/(abs.mCovar.sum) #fixed-effect meta-analysis
        
        set(x = modelGWAS,i =iSNP, j = "BETA",
            #value = sum(W_FULL*BETASNP[iSNP,],na.rm = T)/sum(abs(W_FULL), na.rm = T)
            value = sum(W_FULL*BETASNP[iSNP,],na.rm = T)/sum(abs(W_FULL[!is.na(BETASNP[iSNP,])]), na.rm = T)
            )
        set(x = modelGWAS,i =iSNP, j = "SE",
            value = sqrt(sum(abs.mCovar.sum))
            #value = sqrt(sum(abs.mCovar.sum)/sqrt(k_tot)) #reduced SE correction, /k_tot sends the effects through the roof!
            )
        # set(x = modelGWAS,i =iSNP, j = "k",
        #     value = k[iSNP]
        #     )
        # set(x = modelGWAS,i =iSNP, j = "Q",
        #     value = sum(Q_part)
        #     )
        
      } #for iSNP
      
      #this is probably wrong
      #modelGWAS[, NEF:=(SE^2)/POPVARSNP] #Calculate Effective Sample Size for Factor 1 - from the Genomic SEM Wiki, citation: https://www.biorxiv.org/content/10.1101/603134v3
      modelGWAS[, P:=2*pnorm(q = abs((BETA/SE)),mean = 0, sd = 1, lower.tail = F)] #two sided!
      
      #just to conform to the use of FRQ rather than MAF when possible
      if(any(colnames(modelGWAS)=="MAF") && !any(colnames(modelGWAS)=="FRQ")) colnames(modelGWAS)[colnames(modelGWAS)=="MAF"]<-"FRQ"
      saveRDS(object = modelGWAS,file = factorFilepath) #compress = "bzip2"
      factorCode<-paste0("nmeta.",cModel$code,".F",iFactor,".",p$setup.code)
      
     
      #fwrite(x = modelGWAS, file = factorFilepath, append = F, quote = F, sep = "\t", nThread=5)
      supermunge.result<-supermunge(list_df = list(modelGWAS),traitNames = factorCode, pathDirOutput = p$folderpath.data.sumstats.munged)
      print("Naive metaanalysis of this model factor is done and the result has been exported to a file.")
    }
  } #for iModel
  
  rm("p$lfGWAS$sumstats") #to avoid oom problems 
  
}

if(p$clOptions$task=="nmeta"){quit(save = "no")}




## ----latent factor GWAS-------------------------------------------------------------------

print("***Genomic SEM latent factor GWAS***")
#library(lavaan)
#library(gdata)
#library(parallel)
#library(doParallel)

#on implicit multithreading  
#https://hpc.nih.gov/apps/R.html#threading
# p$numCores<-parallel::detectCores()
# cat("\nNumber of cores detected to be ",p$numCores)
# .Internal(setMaxNumMathThreads(p$numCores-1)) 
# .Internal(setNumMathThreads(p$numCores-1))

#test
#p$clOptions$task<-"lfgwas"
#p$clOptions$task_argument<-"M20_6_172.COR.ML:1"
#p$clOptions$task_argument<-"M20_6_198.COR.ML:1"
#p$clOptions$task_argument<-"M20_6_198.COR.ML:22:1"


 #load intermediate results
# do not remove !p$clOptions$task=="lfgwas" as this will cause any lfgwas job to generate the composite results
if(!p$clOptions$task=="lfgwas" & !file.exists(file.path(p$folderpath.workingDirectory,paste0("lfGWAS.gwas.",p$setup.code,".Rds")))){
  print("Reading in latent factor gwas intermediate results.")
  p$lfGWAS$gwas<-list()
  
  for(iModel in 1:nrow(p$CFA$models.selected)){
    #test
    #iModel<-1
    cModel<-p$CFA$models.selected[iModel,]
    
    
    p$lfGWAS$intermediateResultFiles<-list.files(path = p$folderpath.workingDirectory, pattern = paste0("^lfGWAS\\.gwas\\.",p$setup.code,"\\.",cModel$code,"\\.F.+\\..+\\.Rds"), full.names = T, ignore.case=T)
    #p$lfGWAS$intermediateResultFiles<-list.files(path = p$folderpath.workingDirectory, pattern = paste0("^lfGWAS\\.gwas\\.",p$setup.code,"\\.M.+-.+\\..+\\..+\\.F.+\\..+\\.Rds"), full.names = T, ignore.case=T)
    #lfGWAS.gwas.setup4.M25-4.74.ML.F_ALL.chr8.Rds
    nIntermediateFactors<-NULL
    if(length(p$lfGWAS$intermediateResultFiles)>0){
      p$lfGWAS$gwas[[cModel$code]]<-list()
      for(nIntermediateResultFile in 1:length(p$lfGWAS$intermediateResultFiles)){
        #nIntermediateResultFile<-1
        intermediateResult<-readRDS(file=p$lfGWAS$intermediateResultFiles[nIntermediateResultFile])
        #initialise storage
        if(is.null(nIntermediateFactors)){
          #p$lfGWAS$gwas<-list()
          nIntermediateFactors<-length(intermediateResult)
          if(nIntermediateFactors!=p$CFA$nFactors) warning(paste0("\nNumber of read lfGWAS factors (",nIntermediateFactors,") do not match settings!\n"))
          for(nFactor in 1:nIntermediateFactors){
            p$lfGWAS$gwas[[cModel$code]][[nFactor]]<-intermediateResult[[nFactor]]
          }
        } else {
          for(nFactor in 1:nIntermediateFactors){
            p$lfGWAS$gwas[[cModel$code]][[nFactor]]<-rbind(p$lfGWAS$gwas[[cModel$code]][[nFactor]],intermediateResult[[nFactor]]) 
          }
        }
      }
    }
    
  }
  
  if(length(p$lfGWAS$gwas)>0){
    saveRDS(object = p$lfGWAS$gwas,file = file.path(p$folderpath.workingDirectory,paste0("lfGWAS.gwas.",p$setup.code,".Rds")))
    print("Read latent factor gwas results and saved latent factor summary file.")
    if(p$clOptions$task=="lfgwas"){quit(save = "no")}
    
  }
  
}

if(p$clOptions$task=="lfgwas" & !file.exists(file.path(p$folderpath.workingDirectory,paste0("lfGWAS.gwas.",p$setup.code,".Rds")))) 
{
  
  p$lfGWAS$sumstats<-readRDS(file=p$filepath.lfgwas.sumstats)
  setDT(p$lfGWAS$sumstats)
  setkeyv(p$lfGWAS$sumstats, cols = c("SNP","CHR","BP"))
  print("Read summary statistics for latent factor GWAS from file.")
  
  #Filter columns of chosen traits only
  #chromosomes<- sort(as.integer(unique(p$lfGWAS$sumstats$CHR)))
  colBeta<-c()
  colSE<-c()
  colFRQ<-c()
  colINFO.LIMP<-c()
  colK<-c()
  for(iTrait in 1:nrow(p$sumstats.sel)){
    #iTrait<-1
    colBeta<-c(colBeta,colnames(p$lfGWAS$sumstats)[grep(paste0("^BETA\\.",p$sumstats.sel[iTrait,c("code")],"$"), ignore.case = TRUE,colnames(p$lfGWAS$sumstats))])
    colSE<-c(colSE,colnames(p$lfGWAS$sumstats)[grep(paste0("^SE\\.",p$sumstats.sel[iTrait,c("code")],"$"), ignore.case = TRUE,colnames(p$lfGWAS$sumstats))])
    colFRQ<-c(colFRQ,colnames(p$lfGWAS$sumstats)[grep(paste0("^FRQ\\.",p$sumstats.sel[iTrait,c("code")],"$"), ignore.case = TRUE,colnames(p$lfGWAS$sumstats))])
    colINFO.LIMP<-c(colINFO.LIMP,colnames(p$lfGWAS$sumstats)[grep(paste0("^INFO.LIMP\\.",p$sumstats.sel[iTrait,c("code")],"$"), ignore.case = TRUE,colnames(p$lfGWAS$sumstats))])
    colK<-c(colK,colnames(p$lfGWAS$sumstats)[grep(paste0("^K\\.",p$sumstats.sel[iTrait,c("code")],"$"), ignore.case = TRUE,colnames(p$lfGWAS$sumstats))])
    # colBeta<-colnames(p$lfGWAS$sumstats)[grep("^BETA\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats))]
    # colSE<-colnames(p$lfGWAS$sumstats)[grep("^SE\\.", ignore.case = TRUE,colnames(p$lfGWAS$sumstats))]
  }
  
  colAll<-c("SNP","CHR","BP","MAF","A1","A2",colBeta,colSE,colFRQ,colINFO.LIMP,colK)
  
  
  #full
  #nrow(p$lfGWAS$sumstats)
  #nrow(na.omit(p$lfGWAS$sumstats))
  
  #selected subset
  p$lfGWAS$sumstats<-p$lfGWAS$sumstats[,..colAll]
  #nrow(p$lfGWAS$sumstats)
  #nrow(na.omit(p$lfGWAS$sumstats))
  
  print("Processing summary statistics for latent factor GWAS - NA values")
  #fill NA values - should be replaced with proper imputation
  
  setDT(p$lfGWAS$sumstats)
  setkeyv(x = p$lfGWAS$sumstats, cols = "SNP")
  colIsRaw<-c()
  for(iTrait in 1:nrow(p$sumstats.sel)){
    cColBeta<-colBeta[iTrait]
    cColK<-colK[iTrait]
    cColIsRaw<-paste0("is.raw.",p$sumstats.sel[iTrait,c("code")])
    p$lfGWAS$sumstats[,(cColIsRaw):=!is.na(get(cColBeta)) & is.na(get(cColK))]
    colIsRaw<-c(colIsRaw,cColIsRaw)
  }
  p$lfGWAS$sumstats[,K:=apply(.SD,MARGIN = 1,FUN = function(x){sum(x)}), .SDcols=colIsRaw]
  p$lfGWAS$sumstats[,K.IMP:=apply(.SD,MARGIN = 1,FUN = function(x){sum(!is.na(x))}), .SDcols=colBeta]
  p$lfGWAS$sumstats[,sumK:=rowSums(.SD,na.rm = T), .SDcols=colK]
  p$lfGWAS$sumstats[,sumINFO.LIMP:=rowSums(.SD,na.rm = T), .SDcols=colINFO.LIMP]
  p$lfGWAS$sumstats[,sumBETA:=rowSums(.SD,na.rm = T), .SDcols=colBeta]
  p$lfGWAS$sumstats[,sumSE:=rowSums(.SD,na.rm = T), .SDcols=colSE]
  p$lfGWAS$sumstats[,meanBETA:=sumBETA/K.IMP][,meanSE:=sumSE/K.IMP]
  
  #process still missing values
  for(iTrait in 1:nrow(p$sumstats.sel)){
    #iTrait<-19
    cColBeta<-colBeta[iTrait]
    cColSE<-colSE[iTrait]
    #nrow(p$lfGWAS$sumstats[is.na(get(cColBeta)),])
    #p$lfGWAS$sumstats[is.na(get(cColBeta)) & K.IMP>0.5*length(colBeta) & K>0.2*length(colBeta),] 
    p$lfGWAS$sumstats[is.na(get(cColBeta)) & K.IMP>0.5*length(colBeta) & K>0.2*length(colBeta), (cColBeta):=meanBETA] 
    p$lfGWAS$sumstats[is.na(get(cColSE)) & K.IMP>0.5*length(colBeta) & K>0.2*length(colBeta), (cColSE):=meanSE] 
    #nrow(p$lfGWAS$sumstats[is.na(get(cColBeta)),])
    #nrow(p$lfGWAS$sumstats[is.na(get(cColSE)),])
    #p$lfGWAS$sumstats[is.na(get(cColBeta)), (cColBeta):=(1e-20)]
    #p$lfGWAS$sumstats[is.na(get(cColSE)), (cColSE):=1]
  }
  
  cat("\nA total of ",nrow(p$lfGWAS$sumstats), " variants before QC.\n")
  
  #nrow(na.omit(p$lfGWAS$sumstats[,..colBeta]))
  #filter columns to those required for lfGWAS
  colLfGWAS<-c("SNP","CHR","BP","MAF","A1","A2",colBeta,colSE)
  #QC - remove variants with NA effects or SE
  p$lfGWAS$sumstats<-na.omit(p$lfGWAS$sumstats[,..colLfGWAS])
  
  cat("\nA total of ",nrow(p$lfGWAS$sumstats), " variants after QC.\n")

  
  #to data frame to conform with userGWAS below
  p$lfGWAS$sumstats<-as.data.frame(p$lfGWAS$sumstats)
  
  p$lfGWAS$cModel<-NULL
  p$lfGWAS$cFn<-NULL
  p$lfGWAS$cChr<-NULL
  #test #p$lfGWAS$cChr<-22
  
  if(length(grep(pattern = "\\:",x = p$clOptions$task_argument))>0){
    p$lfGWAS$cModel<-strsplit(x = p$clOptions$task_argument, split = ':')[[1]][1]
    p$lfGWAS$cChr<-strsplit(x = p$clOptions$task_argument, split = ':')[[1]][2]
    #p$lfGWAS$cFn<-strsplit(x = p$clOptions$task_argument, split = ':')[[1]][3]
  } else if(exists(p$clOptions$task_argument)){
    p$lfGWAS$cModel<-p$clOptions$task_argument
  } else {
    p$lfGWAS$cModel<-p$CFA$models.selected[1,]$code
  }
  
  
  #if(is.null(p$lfGWAS$cFn)) p$lfGWAS$cFn<-"1"
  #if(is.null(p$lfGWAS$cChr)) p$lfGWAS$cChr<-"1"
  
  
  cat("\nPerforming latent factor GWAS")
  cat("\nSelected model:",p$lfGWAS$cModel)
  cat("\nSelected factor ",ifelse(is.null(p$lfGWAS$cFn),"ALL",paste0(p$lfGWAS$cFn))," and chromosome ",ifelse(is.null(p$lfGWAS$cChr),"ALL",paste0(p$lfGWAS$cChr)),". This will take a while!\n")
  
  #select model from code
  cModel<-p$CFA$models.selected[p$lfGWAS$cModel,]
  
  if(is.null(p$lfGWAS$cFn)){
    p$lfGWAS$lmodel<-paste0(cModel$lModel.fixed,paste0("\nF",(1:cModel$nFactors),"~SNP", collapse = ""))
  } else {
    p$lfGWAS$lmodel<-paste0(cModel$lModel.fixed,paste0("\nF",p$lfGWAS$cFn,"~SNP"))
  }
  
  cat("\nExpanded CFA model with SNP effects:")
  print(p$lfGWAS$lmodel)
  
  if(length(p$lfGWAS$cChr)>0){
    cat("\nAnalysing only the specified chromosome: ",paste0(p$lfGWAS$cChr))
    #TEST
    #p$lfGWAS$sumstats.selected<-head(p$lfGWAS$sumstats[which(as.character(p$lfGWAS$sumstats$CHR)==p$lfGWAS$cChr),])
    p$lfGWAS$sumstats<-p$lfGWAS$sumstats[which(as.character(p$lfGWAS$sumstats$CHR)==p$lfGWAS$cChr),]
  }
  
  #TEST settings - restrict to a subset of variants
  #p$lfGWAS$sumstats<-head(p$lfGWAS$sumstats, n = 100)
  
  print("Head of processed sumstats table selection")
  print(head(p$lfGWAS$sumstats))
    
  p$lfGWAS$gwas<-userGWAS.mod(
    covstruc = p$mvLD$covstruct.mvLDSC,
    SNPs = p$lfGWAS$sumstats,
    estimation = "ML",
    model = p$lfGWAS$lmodel,
    modelchi = FALSE,
    printwarn = TRUE,
    sub=paste0("F",(1:cModel$nFactors),"~SNP"),
    parallel=T,
    GC="conserv",
    smooth_check = T,
    TWAS = F,
    #parallel_outfile = file.path(p$folderpath.workingDirectory,"lfgwas.progress.txt"),
    cores = 10 #set cores as standard
    #turbo=T
    )
  
  cat("\nuser GWAS done returning results for ",length(p$lfGWAS$gwas),"factors and ",nrow(p$lfGWAS$gwas[[1]])," SNPs.")
  
  saveRDS(object = p$lfGWAS$gwas,file = file.path(p$folderpath.workingDirectory,paste0("lfGWAS.gwas.",p$setup.code,".",p$lfGWAS$cModel,
                                                                                                   ".F",ifelse(length(p$lfGWAS$cFn)>0,p$lfGWAS$cFn,"ALL"),
                                                                                                   ".",ifelse(length(p$lfGWAS$cChr)>0,paste0("chr",p$lfGWAS$cChr),"ALL"),".Rds")))
  
  print("DONE performing latent factor GWAS. The results should have been saved to a file.")

}


if(p$clOptions$task=="lfgwas"){quit(save = "no")}



## ----process latent factor GWAS results---------------------------------------------------
cat("\n***Process latent factor GWAS results***\n")
#check the first model and first factor file
if(!file.exists(file.path(p$folderpath.data.sumstats.munged,paste0(p$CFA$models.selected[1,]$code,".F",1,".gz")))){
  
  cat("\nProcessing latent factor GWAS results...\n")
  NTot<-sum(p$sumstats.sel$n_total)
  
  p$lfGWAS$gwas<-readRDS(file=file.path(p$folderpath.workingDirectory,paste0("lfGWAS.gwas.",p$setup.code,".Rds")))
    print("Read previously stored latent factor GWAS results from file.")
  
  for(iModel in 1:nrow(p$CFA$models.selected)){
    #test
    #iModel<-1
    cModel<-p$CFA$models.selected[iModel,] 
    
    inflationFactorTotal<-1/cModel$parsedGsemResults[[1]][[1]]$meanTotalRelativeVarianceExplained
    
    if(length(p$lfGWAS$gwas[[cModel$code]])>0){
      
      patCoef<-cModel$parsedGsemResults[[1]][[1]]$patternCoefficientsSTDGenotype.matrix
      patCoefSq<-patCoef^2
      disVars<-cModel$parsedGsemResults[[1]][[1]]$residualVariances.matrix
      #factorSpecificWeight<-colMeans(patCoefSq[,1:ncol(patCoefSq)]*rep((1-disVars),times=ncol(patCoefSq)),na.rm=T)
      #factorSpecificWeight<-colMeans(patCoefSq[,1:ncol(patCoefSq)],na.rm=T)
      factorSpecificWeight<-colSums(patCoefSq[,1:ncol(patCoefSq)]*rep((1-disVars),times=ncol(patCoefSq)),na.rm=T)
      
      #explainedVariance <- patCoef #TODO
      factorWeightSum<-nrow(patCoefSq)
      #factorWeightSum<-sum(factorSpecificWeight)
      effectiveN<-vector()
      
      #loop for calculating the effective N
      #inflating effect sizes - if needed - and saving factor GWAS output files
      for(iFactor in 1:cModel$nFactors){
        #test
        #iFactor<-1
        factorCode<-paste0(cModel$code,".F",iFactor)
        #factorFilepath<-file.path(p$folderpath.workingDirectory,paste0("lfGWAS.",factorCode,".",p$setup.code,".Rds"))
        
         #adjust variant effect size
        #inflationFactorFactor<-inflationFactorTotal*factorSpecificWeight[iFactor]
        #head(p$lfGWAS$gwas[[cModel$code]][[iFactor]])
        #NTot/effectiveN[iFactor]
        p$lfGWAS$gwas[[cModel$code]][[iFactor]]$Z_adj<-p$lfGWAS$gwas[[cModel$code]][[iFactor]]$Z_Estimate
        
        #adjust variant p estimate
        #p$lfGWAS$gwas[[cModel$code]][[iFactor]]$p_adj<-2*pnorm(q = p$lfGWAS$gwas[[cModel$code]][[iFactor]]$Z_Estimate, sd = (1/sqrt(nrow(p$sumstats.sel))), lower.tail = F)
        p$lfGWAS$gwas[[cModel$code]][[iFactor]]$p_adj<-2*pnorm(q = p$lfGWAS$gwas[[cModel$code]][[iFactor]]$Z_Estimate, lower.tail = F)
        #head(p$lfGWAS$gwas[[cModel$code]][[iFactor]])[,c("SNP","Z_Estimate","Z_smooth","Z_adj", "Pval_Estimate","p_adj")]
      
      
        #set which values to use for output
        #p$lfGWAS$gwas[[cModel$code]][[iFactor]]$p<-p$lfGWAS$gwas[[cModel$code]][[iFactor]]$p_adj
        p$lfGWAS$gwas[[cModel$code]][[iFactor]]$P<-p$lfGWAS$gwas[[cModel$code]][[iFactor]]$Pval_Estimate
        p$lfGWAS$gwas[[cModel$code]][[iFactor]]$BETA<-p$lfGWAS$gwas[[cModel$code]][[iFactor]]$est
        p$lfGWAS$gwas[[cModel$code]][[iFactor]]$Z<-p$lfGWAS$gwas[[cModel$code]][[iFactor]]$Z_Estimate
        p$lfGWAS$gwas[[cModel$code]][[iFactor]]$CHISQ<-p$lfGWAS$gwas[[cModel$code]][[iFactor]]$chisq
        p$lfGWAS$gwas[[cModel$code]][[iFactor]]$DF<-p$lfGWAS$gwas[[cModel$code]][[iFactor]]$chisq_df
        
        
        #save the factor GWAS to file
        #saveRDS(object = (p$lfGWAS$gwas[[cModel$code]][[iFactor]])[,c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","Z","CHISQ","DF","P","NEF")],file = factorFilepath)
        supermunge(
          list_df = list((p$lfGWAS$gwas[[cModel$code]][[iFactor]])[,c("SNP","CHR","BP","MAF","A1","A2","BETA","SE","Z","CHISQ","DF","P")]),traitNames = factorCode, pathDirOutput = p$folderpath.data.sumstats.munged
          )
        }
    }
    
  }
}


## ----setup additional traits--------------------------------------------------------------

p$sumstats$code.file<-p$sumstats$code #default as code

#GSEM pre-analysis (simple, no NA effects)
factorCode<-"PRE.F1"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"MPRE_2_1.COR.ML.F1","Risk F1",294363)
factorCode<-"PRE.F2"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"MPRE_2_1.COR.ML.F2","Risk F2",364476)


#GSEM latent factor GWAS
factorCode<-"GSEM.F1"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"M20_6_172.COR.ML.F1","Neuroticism F",364979)
factorCode<-"GSEM.F2"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"M20_6_172.COR.ML.F2","Thought F",159862)
factorCode<-"GSEM.F3"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"M20_6_172.COR.ML.F3","Depression F",448644)
factorCode<-"GSEM.F4"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"M20_6_172.COR.ML.F4","Externalising F",46293)
factorCode<-"GSEM.F5"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"M20_6_172.COR.ML.F5","SE deprivation F",271911)
factorCode<-"GSEM.F6"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"M20_6_172.COR.ML.F6","Neuropsychiatric F",91409)

#nmeta latent factor GWAS
factorCode<-"nmeta.CF1"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_172.COR.ML.F1.setup7","Anxiety CF",237392)
factorCode<-"nmeta.CF2"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_172.COR.ML.F2.setup7","Thought CF",1406685)
factorCode<-"nmeta.CF3"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_172.COR.ML.F3.setup7","Depression CF",1051408)
factorCode<-"nmeta.CF4"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_172.COR.ML.F4.setup7","Externalising CF",1131595)
factorCode<-"nmeta.CF5"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_172.COR.ML.F5.setup7","SE deprivation CF",865415)
factorCode<-"nmeta.CF6"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_172.COR.ML.F6.setup7","Neuropsychiatric CF",2091884)

factorCode<-"nmeta.OF1"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_30.ORT.ML.F1.setup7","Externalising(inv) OF",2379693)
factorCode<-"nmeta.OF2"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_30.ORT.ML.F2.setup7","Neuropsychiatric OF",1957237)
factorCode<-"nmeta.OF3"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_30.ORT.ML.F3.setup7","Fear OF",282511)
factorCode<-"nmeta.OF4"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_30.ORT.ML.F4.setup7","Depression OF",360436)
factorCode<-"nmeta.OF5"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_30.ORT.ML.F5.setup7","Distress(inv) OF",572446)
factorCode<-"nmeta.OF6"
p$sumstats[factorCode,c("code","code.file","name","n_effective")]<-c(factorCode,"nmeta.M20_6_30.ORT.ML.F6.setup7","Thought OF",1096993)


p$sumstats.sel.code<-c("PRE.F1","PRE.F2","GSEM.F1","GSEM.F2","GSEM.F3","GSEM.F4","GSEM.F5","GSEM.F6","nmeta.CF1","nmeta.CF2","nmeta.CF3","nmeta.CF4","nmeta.CF5","nmeta.CF6","nmeta.OF1","nmeta.OF2","nmeta.OF3","nmeta.OF4","nmeta.OF5","nmeta.OF6")

# set nice name for selected traits
p$sumstats[p$sumstats.sel.code,]$name.nice<-p$sumstats[p$sumstats.sel.code,]$name

##Add sumstat cleaned and munged file paths
p$sumstats[p$sumstats.sel.code,c("cleanedpath")]<-file.path(p$folderpath.data.sumstats.cleaned,paste0(p$sumstats[p$sumstats.sel.code,c("code.file")],p$filename.suffix.data.sumstats.munged))

p$sumstats[p$sumstats.sel.code,c("mungedpath")]<-file.path(p$folderpath.data.sumstats.munged,paste0(p$sumstats[p$sumstats.sel.code,c("code.file")],p$filename.suffix.data.sumstats.munged))

##Add a combined nice name plus code label
p$sumstats$name.nice.and_code<-paste0(p$sumstats$name.nice," (",p$sumstats$code,")")

#set order of datasets to sorted by code
p$sumstats<-p$sumstats[order(p$sumstats$name.nice),]

#re-calculate effective N from file information - perform manually as needed
if(F){
        factorCode<-"nmeta.OF6"
        gwasForEffectiveN<-fread(file = file.path(p$folderpath.data.sumstats.munged,paste0( p$sumstats[factorCode,]$code.file,".gz")), na.strings =c(".",NA,"NA",""), encoding = "UTF-8",check.names = T, fill = T, blank.lines.skip = T, key = c("SNP","CHR","BP"), data.table = T,showProgress = F, nThread=6)
  ##Calculate Effective Sample Size for Factor 1 - from the Genomic SEM Wiki
        #citation: https://www.biorxiv.org/content/10.1101/603134v3
        #restrict to MAF of 40% and 10%
        #change to use adjusted Z if this has been set above
        #added to supermunge also
        gwasForEffectiveN[,MAF:=ifelse(FRQ>0.5,1-FRQ,FRQ)][,VSNP:=2*FRQ*(1-FRQ)]
        gwasForEffectiveN<-gwasForEffectiveN[MAF>=.1 & MAF<=.4,]
        setDT(gwasForEffectiveN)
        #gwasForEffectiveN[,NEF:=round(((Z/EFFECT)^2)/VSNP,digits = 0)]
        effectiveN<-mean(gwasForEffectiveN$NEF,na.rm=T)
        round(effectiveN,digits = 0)
        
}


## ----GWAS sumstat dataset variable selection 2--------------------------------------------

#selection based on specific traits
#p$sumstats.sel.code<-c("ANXI04")
#p$sumstats.sel.code<-c("ADHD05","ANXI04")
#p$sumstats.sel.code<-c("RISK02","RISK03","SCHI04","SUBJ01","TIRE01")
#p$sumstats.sel.code<-c("TIRE01")
#p$sumstats.sel.uv.code<-c("EDUC03","NEUR02","RISK02","RISK03","SCHI04")
p$sumstats.sel.uv.code<-c("ADHD05","ALCD03","ANOR02","ANXI03","ANXI04","AUTI07","BIPO02", "DEPR05","DEPR08","EDUC03","HEAL01","INCO03","INSO02","NEUR02", "PTSD04","RISK02","RISK03","SCHI04","SUBJ01","TIRE01")
#p$sumstats.sel.lfgwas.code<-c("PRE.F1","PRE.F2")
p$sumstats.sel.set1_spec.code<-c("GSEM.F1","GSEM.F2","GSEM.F3","GSEM.F4","GSEM.F5","GSEM.F6")
p$sumstats.sel.set2_spec.code<-c("nmeta.CF1","nmeta.CF2","nmeta.CF3","nmeta.CF4","nmeta.CF5","nmeta.CF6")
p$sumstats.sel.set3_spec.code<-c("nmeta.OF1","nmeta.OF2","nmeta.OF3","nmeta.OF4","nmeta.OF5","nmeta.OF6")
p$sumstats.sel.set1.code<-c(p$sumstats.sel.uv.code,p$sumstats.sel.set1_spec.code)
p$sumstats.sel.set2.code<-c(p$sumstats.sel.uv.code,p$sumstats.sel.set2_spec.code)
p$sumstats.sel.set3.code<-c(p$sumstats.sel.uv.code,p$sumstats.sel.set3_spec.code)
p$sumstats.sel.code<-c(p$sumstats.sel.uv.code,p$sumstats.sel.set1_spec.code,p$sumstats.sel.set2_spec.code,p$sumstats.sel.set3_spec.code)
p$sumstats.sel<-p$sumstats[which(p$sumstats$code %in% p$sumstats.sel.code),]

p$sumstats.sel[,c("code","name","name.nice","name.nice.and_code", "year","n_case","n_control","n_total","n_effective", "pmid","reference_doi","samplePrevalence","populationPrevalence","dependent_variable.OLS","dependent_variable.linprob","se.logit","mungedpath")]
p$k.sel<-nrow(p$sumstats.sel)
#View(p$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")])
#set order of datasets to sorted by code
p$sumstats.sel[is.na(p$sumstats.sel$name.nice),c("name.nice")]<-p$sumstats.sel[is.na(p$sumstats.sel$name.nice),c("name")]
p$sumstats.sel<-p$sumstats.sel[order(p$sumstats.sel$name.nice),]

write.table(p$sumstats.sel[,c("code", "name.nice","year", "n_case","n_control","n_total","n_effective","samplePrevalence","populationPrevalence", "reference_doi")], file = file.path(p$folderpath.workingDirectory,paste0(p$setup.code,".sumstatinfo.tsv")), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)

#this should not be needed as the SNP's already have their NEF
#set n-total for new latent traits from efective n
# cond<-is.na(p$sumstats.sel$n_total)&!is.na(p$sumstats.sel$n_effective)
# p$sumstats.sel[cond,c("n_total")]<-p$sumstats.sel[cond,c("n_effective")]

#View(p$sumstats.sel)



## ----mvLD 2, including latent factors-----------------------------------------------------
print("***Multivariate LD 2***")
p$filepath.mvLD2.set1<-file.path(p$folderpath.workingDirectory,paste0("mvLD2.set1.",p$setup.code,".Rds"))
p$filepath.mvLD2.set2<-file.path(p$folderpath.workingDirectory,paste0("mvLD2.set2.",p$setup.code,".Rds"))
p$filepath.mvLD2.set3<-file.path(p$folderpath.workingDirectory,paste0("mvLD2.set3.",p$setup.code,".Rds"))

#set1 - GSEM
if (file.exists(p$filepath.mvLD2.set1)) {
  print("Using existing covariance structures from previous LD computations.")
  p$mvLD2.set1<-readRDS(file=p$filepath.mvLD2.set1)
} else {
  
  print("Running multivariate LD regression with different methods. This might take a while. If the procedure runs for too long you may want to abort the process.")
  
  cat("The current task is specified as:",p$clOptions$task)
  p$mvLD2.set1<-c()
  
  if(p$clOptions$task=="mvLD2"){
    #run mvLDSC

    p$mvLD2.set1$covstruct.mvLDSC.1kg<-ldsc.mod(
      traits = p$sumstats.sel[p$sumstats.sel.set1.code,]$mungedpath,
      sample.prev =  p$sumstats.sel[p$sumstats.sel.set1.code,]$samplePrevalence,
      population.prev = p$sumstats.sel[p$sumstats.sel.set1.code,]$populationPrevalence,
      trait.names = p$sumstats.sel[p$sumstats.sel.set1.code,]$code,
      ld = p$folderpath.data.mvLDSC.ld.1kg,
      wld = p$folderpath.data.mvLDSC.ld.1kg,
      n.blocks = 600, #a bit more here if it can support the larger reference panel used
      info.filter = 0.6,
      frq.filter = 0.01,
      mhc.filter = 37,
      #chisq.min = 1e-3,
      N = p$sumstats.sel[p$sumstats.sel.set1.code,]$n_total,
      #forceN = T, # Consider this when some of the original N's may be untrustworthy (ANXI04!) - TODO - fix in supermunge
      ldsc.log = p$setup.code.date
      )

    
    saveRDS(object = p$mvLD2.set1,file = p$filepath.mvLD2.set1)
    print("Multivariate LD correction is done now and the resulting covariance structure should have been saved to a file.")
  
  } 
}

#set2 - nmeta, best correlated factor model
if (file.exists(p$filepath.mvLD2.set2)) {
  print("Using existing covariance structures from previous LD computations.")
  p$mvLD2.set2<-readRDS(file=p$filepath.mvLD2.set2)
} else {
  
  print("Running multivariate LD regression with different methods. This might take a while. If the procedure runs for too long you may want to abort the process.")
  
  cat("The current task is specified as:",p$clOptions$task)
  p$mvLD2.set2<-c()
  
  if(p$clOptions$task=="mvLD2"){
    #run mvLDSC

    p$mvLD2.set2$covstruct.mvLDSC.1kg<-ldsc.mod(
      traits = p$sumstats.sel[p$sumstats.sel.set2.code,]$mungedpath,
      sample.prev =  p$sumstats.sel[p$sumstats.sel.set2.code,]$samplePrevalence,
      population.prev = p$sumstats.sel[p$sumstats.sel.set2.code,]$populationPrevalence,
      trait.names = p$sumstats.sel[p$sumstats.sel.set2.code,]$code,
      ld = p$folderpath.data.mvLDSC.ld.1kg,
      wld = p$folderpath.data.mvLDSC.ld.1kg,
      n.blocks = 600, #a bit more here if it can support the larger reference panel used
      info.filter = 0.6,
      frq.filter = 0.01,
      mhc.filter = 37,
      #chisq.min = 1e-3,
      N = p$sumstats.sel[p$sumstats.sel.set2.code,]$n_total,
      #forceN = T, # Consider this when some of the original N's may be untrustworthy (ANXI04!) - TODO - fix in supermunge
      ldsc.log = p$setup.code.date
      )

    
    saveRDS(object = p$mvLD2.set2,file = p$filepath.mvLD2.set2)
    print("Multivariate LD correction is done now and the resulting covariance structure should have been saved to a file.")
  
  } 
}

#set3 - nmeta, best orthogonal factor model
if (file.exists(p$filepath.mvLD2.set3)) {
  print("Using existing covariance structures from previous LD computations.")
  p$mvLD2.set3<-readRDS(file=p$filepath.mvLD2.set3)
} else {
  
  print("Running multivariate LD regression with different methods. This might take a while. If the procedure runs for too long you may want to abort the process.")
  
  cat("The current task is specified as:",p$clOptions$task)
  p$mvLD2.set3<-c()
  
  if(p$clOptions$task=="mvLD2"){
    #run mvLDSC

    p$mvLD2.set3$covstruct.mvLDSC.1kg<-ldsc.mod(
      traits = p$sumstats.sel[p$sumstats.sel.set3.code,]$mungedpath,
      sample.prev =  p$sumstats.sel[p$sumstats.sel.set3.code,]$samplePrevalence,
      population.prev = p$sumstats.sel[p$sumstats.sel.set3.code,]$populationPrevalence,
      trait.names = p$sumstats.sel[p$sumstats.sel.set3.code,]$code,
      ld = p$folderpath.data.mvLDSC.ld.1kg,
      wld = p$folderpath.data.mvLDSC.ld.1kg,
      n.blocks = 600, #a bit more here if it can support the larger reference panel used
      info.filter = 0.6,
      frq.filter = 0.01,
      mhc.filter = 37,
      #chisq.min = 1e-3,
      N = p$sumstats.sel[p$sumstats.sel.set3.code,]$n_total,
      #forceN = T, # Consider this when some of the original N's may be untrustworthy (ANXI04!) - TODO - fix in supermunge
      ldsc.log = p$setup.code.date
      )

    
    saveRDS(object = p$mvLD2.set3,file = p$filepath.mvLD2.set3)
    print("Multivariate LD correction is done now and the resulting covariance structure should have been saved to a file.")
  
  } 
}


if(p$clOptions$task=="mvLD2"){
      quit(save = "no")
}


#add newly computed heritabilities and LDSC intercepts to the selected summary statistics table
p$sumstats.sel[p$sumstats.sel.set1.code,c("h2.liability_mvLDSC")]<-diag(p$mvLD2.set1$covstruct.mvLDSC.1kg$S[p$sumstats.sel.set1.code,p$sumstats.sel.set1.code])
p$sumstats.sel[p$sumstats.sel.set1.code,c("h2.se.liability_mvLDSC")]<-diag(p$mvLD2.set1$covstruct.mvLDSC.1kg$S.SE[p$sumstats.sel.set1.code,p$sumstats.sel.set1.code])

p$sumstats.sel[p$sumstats.sel.set2_spec.code,c("h2.se.liability_mvLDSC")]<-diag(p$mvLD2.set2$covstruct.mvLDSC.1kg$S[p$sumstats.sel.set2_spec.code,p$sumstats.sel.set2_spec.code])
p$sumstats.sel[p$sumstats.sel.set2_spec.code,c("h2.se.liability_mvLDSC")]<-diag(p$mvLD2.set2$covstruct.mvLDSC.1kg$S.SE[p$sumstats.sel.set2_spec.code,p$sumstats.sel.set2_spec.code])

p$sumstats.sel[p$sumstats.sel.set3_spec.code,c("h2.se.liability_mvLDSC")]<-diag(p$mvLD2.set3$covstruct.mvLDSC.1kg$S[p$sumstats.sel.set3_spec.code,p$sumstats.sel.set3_spec.code])
p$sumstats.sel[p$sumstats.sel.set3_spec.code,c("h2.se.liability_mvLDSC")]<-diag(p$mvLD2.set3$covstruct.mvLDSC.1kg$S.SE[p$sumstats.sel.set3_spec.code,p$sumstats.sel.set3_spec.code])

colnames(p$mvLD2.set1$covstruct.mvLDSC.1kg$I)<-colnames(p$mvLD2.set1$covstruct.mvLDSC.1kg$S)
rownames(p$mvLD2.set1$covstruct.mvLDSC.1kg$I)<-rownames(p$mvLD2.set1$covstruct.mvLDSC.1kg$S)
p$sumstats.sel[p$sumstats.sel.set1.code,c("lambda_ldsc")]<-diag(p$mvLD2.set1$covstruct.mvLDSC.1kg$I[p$sumstats.sel.set1.code,p$sumstats.sel.set1.code])




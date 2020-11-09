## ----clean, include=FALSE---------------------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.


## ----settings---------------------------------------------------------------------------------------------------------
project<-c() #create project metadata object
project$date.run<-Sys.Date()
project$setup.version<-1
project$setup.code<-paste0("setup",project$setup.version)
project$setup.code.date<-paste0(project$setup.code,"_",project$date.run)
project$filename.rmd<-paste0(project$setup.code,".Rmd")
project$filename.r<-paste0(project$setup.code,".R")

project$host<-"local" #this is the place where the code is run [local,cluster]
project$seting.refreshPrepareSummaryStatistics<-FALSE
project$setting.refreshLatentFactorGWAS<-FALSE

#file path settings
##set project shared working directory **change if you have other settings**
if(project$host=="local") {
project$folderpath<-normalizePath("~/King's College London/MT-Translational Neuropsychiatric Genomics - Johan_Zvrskovec_PhD - Johan_Zvrskovec_PhD/JZ_GED_PHD_C1")
} else if (project$host=="cluster") {
project$folderpath<-normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1")
}
##project working directory subfolders
project$folderpath.workingDirectory<-normalizePath(paste0(project$folderpath,"/","workingDirectory"))
project$folderpath.data<-normalizePath(paste0(project$folderpath,"/","data"))
project$folderpath.scripts<-normalizePath(paste0(project$folderpath,"/","scripts"))
project$folderpath.include<-normalizePath(paste0(project$folderpath,"/","include"))
project$folderpath.plots<-normalizePath(paste0(project$folderpath,"/","plots"))

##cleaned datasets folder
if(project$host=="local") {
  project$folderpath.data.sumstats.cleaned<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data.sumstats.cleaned")
} else if (project$host=="cluster") {
  project$folderpath.data.sumstats.cleaned<-normalizePath("/mnt/lustre/groups/ukbiobank/sumstats/cleaned")
}

##munged datasets folder
if(project$host=="local") {
  project$folderpath.data.sumstats.munged<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data.sumstats.mungedNoMHC")
} else if (project$host=="cluster") {
  project$folderpath.data.sumstats.munged<-normalizePath("/mnt/lustre/groups/ukbiobank/sumstats/munged_noMHC")
}

##LD scores datasets folders (these strings need to have a trailing slash for the GSEM LDSC to work)
project$folderpath.data.mvLDSC.ld <- paste0(project$folderpath.data,"/eur_w_ld_chr/") #LD-scores
project$folderpath.data.mvLDSC.wld <- paste0(project$folderpath.data,"/eur_w_ld_chr/") #Weights, if different from LD-scores

##Reference file for calculating SNP variance across traits. Used in the preparation step for performing latent factor GWAS.
project$filepath.genomeReference<-normalizePath(paste0(project$folderpath.data,"/","reference.1000G.maf.0.005.txt")) #1000 genomes phase 3

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

setwd(dir = normalizePath(project$folderpath.workingDirectory))



## ----package setup, echo=FALSE, warning=F------------------------------------------------------------------------

#install.packages("skimr")
#install.packages("psych")
#install.packages("Matrix")
#install.packages("tidyverse")
#install.packages("ggrepel")
#install.packages("gt")
#install.packages("kableExtra")
#devtools::install_github("taiyun/corrplot", build_vignettes = TRUE)
#remove.packages("GenomicSEM")
#devtools::install_github("MichelNivard/GenomicSEM",ref = 'v2.1') #specify branch for better stability
#devtools::install_github("MichelNivard/GenomicSEM") #master branch as default
#remove.packages("GenomicSEM")
#devtools::install_github("johanzvrskovec/GenomicSEM",ref = 'hdl-fixes') #specify branch
#devtools::install_github("zhenin/HDL/HDL")
#devtools::install_github("zhenin/HDL/HDL@77cb9d0984d1302e40bfd871491e292f8f09f49d")
#remove.packages("HDL")
#install.packages("optparse")
#remove.packages("semPlate")
#devtools::install_github(repo = "k19049801/semPlate", host="github.kcl.ac.uk", ) #does not work
#devtools::install_local(path =  paste(project$folderpath.include,"/semPlate_0.1.tar.gz"))
#Better install semPlate from the semPlate R-studio package project

library(skimr)
library(psych)
library(Matrix)

library(tidyverse)
library(ggrepel)
library(gt)
#library(kableExtra)
library(GenomicSEM)
library(HDL)
library(optparse)
library(stats)



## ----command line setup------------------------------------------------------------------------------------------
clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="integer", default=1,
                help="Index of the explicit task to run separately:\n1: No task\n2:HDL Piecewise\n3:HDL Jackknife\n4:vLDSC\n5:original piecewise HDL\n6:original jackknife HDL [default %default]")
clParser <- add_option(clParser, c("-l", "--location"), type="character", default="local",
                help="The place where the code is run [local,cluster] [default %default]")



## ----settings----------------------------------------------------------------------------------------------------
project<-c() #create project metadata object
project$clOptions<-parse_args(clParser)
project$date.run<-Sys.Date()
project$setup.version<-1
project$setup.code<-paste0("setup",project$setup.version)
project$setup.code.date<-paste0(project$setup.code,"_",project$date.run)
project$filename.rmd<-paste0(project$setup.code,".Rmd")
project$filename.r<-paste0(project$setup.code,".R")
project$functions<-c()


project$host<-project$clOptions$location #this is the place where the code is run [local,cluster] read from command line - default local
project$seting.refreshPrepareSummaryStatistics<-FALSE
project$setting.refreshLatentFactorGWAS<-FALSE

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
if(project$host=="local") {
project$folderpath<-normalizePath("~/King's College London/MT-Translational Neuropsychiatric Genomics - Johan_Zvrskovec_PhD - Johan_Zvrskovec_PhD/JZ_GED_PHD_C1")
} else if (project$host=="cluster") {
project$folderpath<-normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1")
}
##project working directory subfolders
project$folderpath.workingDirectory<-normalizePath(paste0(project$folderpath,"/","workingDirectory"))

project$folderpath.scripts<-normalizePath(paste0(project$folderpath,"/","scripts"))
#project$folderpath.include<-normalizePath(paste0(project$folderpath,"/","include"))
project$folderpath.plots<-normalizePath(paste0(project$folderpath,"/","plots"))

#general data folder
if(project$host=="local") {
  project$folderpath.data<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data")
} else if (project$host=="cluster") {
  project$folderpath.data<-normalizePath(paste0(project$folderpath,"/","data"))
}

##cleaned sumstats folder
if(project$host=="local") {
  project$folderpath.data.sumstats.cleaned<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data.sumstats.cleaned")
} else if (project$host=="cluster") {
  project$folderpath.data.sumstats.cleaned<-normalizePath("/mnt/lustre/groups/ukbiobank/sumstats/cleaned")
}

##munged sumstats folder
if(project$host=="local") {
  project$folderpath.data.sumstats.munged<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data.sumstats.mungedNoMHC")
} else if (project$host=="cluster") {
  project$folderpath.data.sumstats.munged<-normalizePath("/mnt/lustre/groups/ukbiobank/sumstats/munged_noMHC")
}

project$filename.suffix.data.sumstats.munged<-"_noMHC.sumstats"

##LD scores datasets folders (these strings need to have a trailing slash for the GSEM LDSC to work)
project$folderpath.data.mvLDSC.ld <- paste0(project$folderpath.data,"/eur_w_ld_chr/") #LD-scores
project$folderpath.data.mvLDSC.wld <- paste0(project$folderpath.data,"/eur_w_ld_chr/") #Weights, if different from LD-scores

##HLD reference panel
if(project$host=="local") {
  #use the smallest reference panel as default for local tests
  project$folderpath.data.HLD.ld<-paste0(project$folderpath.data,"/UKB_array_SVD_eigen90_extraction/")
} else if (project$host=="cluster") {
  project$folderpath.data.HLD.ld<-paste0(project$folderpath.data,"/UKB_imputed_hm3_SVD_eigen99_extraction/")
}
  
##full script file paths
project$filepath.rmd<-normalizePath(paste0(project$folderpath.scripts,"/",project$filename.rmd))
project$filepath.r<-normalizePath(paste0(project$folderpath.scripts,"/",project$filename.r))

##Reference file for calculating SNP variance across traits. Used in the preparation step for performing latent factor GWAS.
project$filepath.genomeReference<-normalizePath(paste0(project$folderpath.data,"/","reference.1000G.maf.0.005.txt")) #1000 genomes phase 3

##CFA settings 
project$cfa.absCutoff.2F.oblq=c(0.30,0.35,0.4) #for promax
project$cfa.absCutoff.2F.orth=c(0.35,0.4,0.45) #for varimax
project$cfa.absCutoff.3F.oblq=c(0.2,0.3,0.4) #for promax
project$cfa.absCutoff.3F.orth=c(0.2,0.3,0.4,0.5) #for varimax

project$cfa.estimator=c("ML","DWLS")
project$cfa.orthogonal=c(TRUE,FALSE)
project$cfa.fixedLoadings=c(TRUE,FALSE) 

##latent factor GWAS filter settings
project$info.filter=.6
project$maf.filter=0.01

setwd(dir = normalizePath(project$folderpath.workingDirectory))







## ----additional source setup, echo=FALSE, warning=F--------------------------------------------------------------

source(normalizePath(paste0(project$folderpath.scripts,"/","shru.R")))
#source(normalizePath(paste0(project$folderpath.scripts,"/","hdl.mod.R")))



## ----trait setup-------------------------------------------------------------------------------------------------
#,echo=FALSE
project$trait<-data.frame(
  code=c("ANXI","DEPR","BIPO","ALCD"),
  populationPrevalence=c(
    .16, #Any type of anxiety disorder, Via https://doi.org/10.1038/s41380-019-0559-1 (2019), originally from https://doi.org/10.1017/S1121189X00001421 (2009)
    .15, #MDD, from the LD-calculations in https://doi.org/10.1038/s41588-018-0090-3 (2018), but with a possible reference to https://doi.org/10.1146/annurev-publhealth-031912-114409 (2013) which states 11.1% lifetime prevalence of MDE.
    .007, #mean of male and female global prevalence rate (2013) from https://doi.org/10.1111/bdi.12423 (2016)
    .159 # European measurement from https://doi.org/10.1001/archpsyc.64.7.830 via https://doi.org/10.1038/s41593-018-0275-1 (how was this calculated though?)
    ), 
  referenceDOI=c(
    "https://doi.org/10.1017/S1121189X00001421",
    "https://doi.org/10.1038/s41588-018-0090-3",
    "https://doi.org/10.1111/bdi.12423",
    "https://doi.org/10.1001/archpsyc.64.7.830"
    ))

project$trait



## ----GWAS sumstat dataset setup----------------------------------------------------------------------------------
#, echo=FALSE

project$sumstats<-read.table(paste0(project$folderpath.data,"/","ukbb_sumstats_download202005.csv"), header=T, quote="\"", sep = ",", fill=T, blank.lines.skip=T,as.is = c(2), strip.white = T)


#rename and add columns
names(project$sumstats)[names(project$sumstats)=="n_cases"]<-"n_case"
names(project$sumstats)[names(project$sumstats)=="n_controls"]<-"n_control"
project$sumstats$gwas_name.nice<-NA_character_
project$sumstats$code.trait<-NA_character_
project$sumstats$reference_doi<-NA_character_
project$sumstats$effect.logit<-as.logical(NA)
project$sumstats$se.logit<-as.logical(NA)
project$sumstats$dependent_variable.OLS<-as.logical(NA)
project$sumstats$age.min<-NA_integer_
project$sumstats$age.max<-NA_integer_
project$sumstats$age.mean<-NA_real_
project$sumstats$age.sd<-NA_real_

#add missing datasets and data
project$sumstats[nrow(project$sumstats)+1,c("code","n_case","n_control","n_total","reference_doi")]=list(
  code=c("DEPR05"),
  n_case=16823,
  n_control=25632,
  n_total=42455,
  reference_doi=c("https://doi.org/10.1038/s41588-018-0090-3")
  )

#reformat columns
project$sumstats$gwas_name<-as.character(project$sumstats$gwas_name)
project$sumstats$gwas_name.nice<-project$sumstats$gwas_name

##Add comprehensive names as in the Google sheet
project$sumstats$gwas_name[which(project$sumstats$code=="DEPR05")]="Major depressive disorder (PGC2 29)"

##Add nice trait names to be used in the report
project$sumstats$gwas_name.nice[which(project$sumstats$code=="ALCD03")]="Alcohol dependence"
project$sumstats$gwas_name.nice[which(project$sumstats$code=="ANXI03")]="Anxiety disorder"
project$sumstats$gwas_name.nice[which(project$sumstats$code=="NEUR01")]="Neuroticism"
project$sumstats$gwas_name.nice[which(project$sumstats$code=="SUBJ01")]="Subjective well-being"
project$sumstats$gwas_name.nice[which(project$sumstats$code=="TIRE01")]="Self-reported tiredness"
project$sumstats$gwas_name.nice[which(project$sumstats$code=="DEPR05")]="Major depressive disorder"


##Add trait/disorder information
project$sumstats <- project$sumstats %>%
mutate(
  code.trait=substr(x = code, start = 1, stop = 4)
       ) %>%
  left_join(project$trait[,c("code","populationPrevalence")], by = c("code.trait" = "code"))

##Add sumstat file paths
project$sumstats <- project$sumstats %>%
mutate(
  mungedpath=paste0(project$folderpath.data.sumstats.munged,"/",code,"_noMHC.sumstats.gz")
       )

project$sumstats <- project$sumstats %>%
mutate(
  cleanedpath=paste0(project$folderpath.data.sumstats.cleaned,"/",code,".gz")
       )

##add reference year
project$sumstats$reference_year[which(project$sumstats$code=="ANXI03")]=2019
project$sumstats$reference_year[which(project$sumstats$code=="NEUR01")]=2016
project$sumstats$reference_year[which(project$sumstats$code=="DEPR05")]=2018

##Add doi links for easy access to dataset publication
project$sumstats$reference_doi[which(project$sumstats$code=="ALCD03")]="https://doi.org/10.1038/s41593-018-0275-1"
project$sumstats$reference_doi[which(project$sumstats$code=="ANXI03")]="https://doi.org/10.1038/s41380-019-0559-1"
project$sumstats$reference_doi[which(project$sumstats$code=="HEAL01")]="https://doi.org/10.1093/ije/dyw219"
project$sumstats$reference_doi[which(project$sumstats$code=="NEUR01")]="https://doi.org/10.1038/ng.3552"
project$sumstats$reference_doi[which(project$sumstats$code=="SUBJ01")]="https://doi.org/10.1038/ng.3552"
project$sumstats$reference_doi[which(project$sumstats$code=="TIRE01")]="https://doi.org/10.1038/mp.2017.5"


##Add PMID
project$sumstats$pmid[which(project$sumstats$code=="ANXI03")]="31748690"
project$sumstats$pmid[which(project$sumstats$code=="NEUR01")]="27089181"
project$sumstats$pmid[which(project$sumstats$code=="SUBJ01")]="27089181"
project$sumstats$pmid[which(project$sumstats$code=="DEPR05")]="29700475"

##add dependent variable type
project$sumstats$dependent_variable[which(project$sumstats$code=="NEUR01")]="continuous"
project$sumstats$dependent_variable[which(project$sumstats$code=="DEPR05")]="binary"

##add participant numbers
project$sumstats$n_total[which(project$sumstats$code=="NEUR01")]=170911

##add ancestry details
project$sumstats$ancestry[which(project$sumstats$code=="NEUR01")]="EUR"
project$sumstats$ancestry[which(project$sumstats$code=="DEPR05")]="EUR"

##add sex details
project$sumstats$sex[which(project$sumstats$code=="NEUR01")]="both"
project$sumstats$sex[which(project$sumstats$code=="DEPR05")]="both"

##add age range
project$sumstats$age.min[which(project$sumstats$code=="SUBJ01")]=40
project$sumstats$age.max[which(project$sumstats$code=="SUBJ01")]=73
project$sumstats$age.mean[which(project$sumstats$code=="SUBJ01")]=56.91
project$sumstats$age.sd[which(project$sumstats$code=="SUBJ01")]=7.93
project$sumstats$age.min[which(project$sumstats$code=="TIRE01")]=40
project$sumstats$age.max[which(project$sumstats$code=="TIRE01")]=73
project$sumstats$age.mean[which(project$sumstats$code=="TIRE01")]=56.91
project$sumstats$age.sd[which(project$sumstats$code=="TIRE01")]=7.93

##add data on whether effects are on log odds ratio or not
project$sumstats$effect.logit[which(project$sumstats$code=="ALCD03")]=T #Did not find the original readme for this, only looked at the readme for the 2018 subsequent release, and the data. Assuming log-scale from the naming of the variables (beta) and that it follows the standard of the 2018 release - explicitly states log scale beta. Values are consistent with regression betas or log OR. THESE MIGHT HAVE BEEN CONVERTED IN THE FLASK APP LATER THOUGH!! CHECK!!               
project$sumstats$effect.logit[which(project$sumstats$code=="ANXI03")]=F # Kirstin sort of confirm that she converted both effect and s.e. to linear scales. Data and readme is available from her KCL web site. Unclear which dataset was used for ANXI03 though. Confirm this later by comparing datasets?
project$sumstats$effect.logit[which(project$sumstats$code=="HEAL01")]=F #The original data including readme downloaded to data_raw/original. States "Beta" without explanation, and I assume linear beta.
project$sumstats$effect.logit[which(project$sumstats$code=="NEUR01")]=F #The official read-me states "regression beta"
project$sumstats$effect.logit[which(project$sumstats$code=="SUBJ01")]=F #The official read-me states "regression beta"
project$sumstats$effect.logit[which(project$sumstats$code=="TIRE01")]=F #The official read-me states "beta"
project$sumstats$effect.logit[which(project$sumstats$code=="DEPR05")]=F #https://www.med.unc.edu/pgc/files/2019/04/PGC_MDD_2018_README_third_180316.pdf # The readme states daner format, but that SE is s.e. of log(OR), while the daner format does not state log scale variables here. Assuming linear scale OR, which is consistent with the data.

##add data on whether the SEs are on a logistic scale or not
project$sumstats$se.logit[which(project$sumstats$code=="ALCD03")]=T #Did not find the original readme for this, only looked at the readme for the 2018 subsequent release, and the data. SE is strangely set to 1 for top rows. Assumes this is true because of log-scale beta.
project$sumstats$se.logit[which(project$sumstats$code=="ANXI03")]=F #Assume linear becasue of Kirstin Purves assuring that both were converted to linear scale.
project$sumstats$se.logit[which(project$sumstats$code=="HEAL01")]=F #Assume linear scale because of linear beta
project$sumstats$se.logit[which(project$sumstats$code=="NEUR01")]=F #Assume linear scale because of linear beta
project$sumstats$se.logit[which(project$sumstats$code=="SUBJ01")]=F #Assume linear scale because of linear beta
project$sumstats$se.logit[which(project$sumstats$code=="TIRE01")]=F #Assume linear scale because of linear beta
project$sumstats$se.logit[which(project$sumstats$code=="DEPR05")]=F #Assume linear scale because of linear OR

##add information on wether a continous dependent variable was analysed using an OLS (linear) estimator. Used for the latent factor GWAS preparation step.
project$sumstats$dependent_variable.OLS[which(project$sumstats$code=="ALCD03")]=F
project$sumstats$dependent_variable.OLS[which(project$sumstats$code=="ANXI03")]=F
project$sumstats$dependent_variable.OLS[which(project$sumstats$code=="DEPR05")]=F
project$sumstats$dependent_variable.OLS[which(project$sumstats$code=="HEAL01")]=T
project$sumstats$dependent_variable.OLS[which(project$sumstats$code=="NEUR01")]=T
project$sumstats$dependent_variable.OLS[which(project$sumstats$code=="SUBJ01")]=T
project$sumstats$dependent_variable.OLS[which(project$sumstats$code=="TIRE01")]=T

##add sample prevalence for all datasets
project$sumstats$samplePrevalence<-project$sumstats$n_case/project$sumstats$n_total

project$sumstats<-project$sumstats[order(project$sumstats$code),]

#save the project data
saveRDS(project,file = paste0(project$folderpath.workingDirectory,"/","project.",project$setup.code,".Rds"))

#View(project$sumstats)



## ----GWAS sumstat dataset variable selection---------------------------------------------------------------------

#selection based on specific traits
project$sumstats.sel.code<-c("DEPR05","ANXI03","NEUR01","TIRE01","SUBJ01","ALCD03","HEAL01")
project$sumstats.sel<-project$sumstats[which(project$sumstats$code %in% project$sumstats.sel.code),]
project$sumstats.sel$code<-project$sumstats.sel$code.trait
project$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")]
project$k.sel<-nrow(project$sumstats.sel)
#View(project$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")])



## ----multivariate LD---------------------------------------------------------------------------------------------

project$filepath.mvLD<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".Rds")
project$filepath.mvLD.HDL.piecewise<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".HDL.piecewise.Rds")
project$filepath.mvLD.HDL.jackknife<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".HDL.jackknife.Rds")
project$filepath.mvLD.mvLDSC<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".mvLDSC.Rds")
project$filepath.mvLD.origHDL<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL.Rds")
project$filepath.mvLD.origHDL.liabilityScale<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL.liabilityScale.Rds")


if (file.exists(project$filepath.mvLD)) {
  print("Using existing covariance structure from previous HDL computation.")
  project$mvLD<-readRDS(file=project$filepath.mvLD)
} else {
  print("Running (or reading ready results from) multivariate LD regression with HDL. This might take a while. If the procedure runs for too long you may want to abort the process.")
  project$mvLD<-c()
  
  if(project$clOptions$task==1){
    
    if(!file.exists(project$filepath.mvLD.HDL.piecewise) || !file.exists(project$filepath.mvLD.HDL.jackknife) || !file.exists(project$filepath.mvLD.mvLDSC)) {
      
      shru$evercat(message = paste0(
        "\nHDL.piecewise=",file.exists(project$filepath.mvLD.HDL.piecewise),
        " HDL.jackknife=",file.exists(project$filepath.mvLD.HDL.jackknife),
        " mvLDSC=",file.exists(project$filepath.mvLD.mvLDSC),
        "\nRunning time=",round(Sys.time()-project$mvLD$startTime,3)
        ))
      Sys.sleep(time = 3)
    }
    
    project$mvLD$covstruct.HDL.piecewise<-readRDS(file=project$filepath.mvLD.HDL.piecewise)
    project$mvLD$covstruct.HDL.jackknife<-readRDS(file=project$filepath.mvLD.HDL.jackknife)
    project$mvLD$covstruct.mvLDSC<-readRDS(file=project$filepath.mvLD.mvLDSC)
    project$mvLD$covstruct.origHDL<-readRDS(file=project$filepath.mvLD.origHDL)
    project$mvLD$covstruct.origHDL.liabilityScale<-readRDS(file=project$filepath.mvLD.origHDL.liabilityScale)
  }
  
  
  if(project$clOptions$task==2 && !file.exists(project$filepath.mvLD.HDL.piecewise)){
    #run HDL piecewise
    #system(command = paste0("touch 2.txt"))
    project$mvLD$covstruct.HDL.piecewise <- hdl(traits = project$sumstats.sel$mungedpath,
                                  sample.prev = project$sumstats.sel$samplePrevalence,
                                  population.prev = project$sumstats.sel$populationPrevalence,
                                  trait.names=project$sumstats.sel$code,
                                  LD.path=project$folderpath.data.HLD.ld,
                                  method = "piecewise")
    
    saveRDS(object = project$mvLD$covstruct.HDL.piecewise,file = project$filepath.mvLD.HDL.piecewise)
    quit(save = "no")
  }
  
  if(project$clOptions$task==3 && !file.exists(project$filepath.mvLD.HDL.jackknife)){
    #run HDL jackknife
    #system(command = paste0("touch 3.txt"))
    project$mvLD$covstruct.HDL.jackknife <- hdl(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HLD.ld,
                                method = "jackknife")
  
    saveRDS(object = project$mvLD$covstruct.HDL.jackknife,file = project$filepath.mvLD.HDL.jackknife)
    quit(save = "no")
  }
  
  if(project$clOptions$task==4 && !file.exists(project$filepath.mvLD.mvLDSC)){
    #run mvLDSC
    #system(command = paste0("touch 4.txt"))
    project$mvLD$covstruct.mvLDSC<-ldsc(traits = project$sumstats.sel$mungedpath,
                                    sample.prev =  project$sumstats.sel$samplePrevalence,
                                    population.prev = project$sumstats.sel$populationPrevalence,
                                    trait.names = project$sumstats.sel$code,
                                    ld = project$folderpath.data.mvLDSC.ld,
                                    wld = project$folderpath.data.mvLDSC.ld,
                                    ldsc.log = project$setup.code.date,
                                    stand = TRUE
                                    )
    saveRDS(object = project$mvLD$covstruct.mvLDSC,file = project$filepath.mvLD.mvLDSC)
    quit(save = "no")
  }
  
  if(project$clOptions$task==5 && !file.exists(project$filepath.mvLD.origHDL)){
    #run HDL original (using the original jackknife procedure)
    project$mvLD$covstruct.origHDL<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HLD.ld, liabilityScale = FALSE)
    
    saveRDS(object = project$mvLD$covstruct.origHDL,file = project$filepath.mvLD.origHDL)
    quit(save = "no")
  }
  
  if(project$clOptions$task==6 && !file.exists(project$filepath.mvLD.origHDL.liabilityScale)){
    #run HDL original (using the original jackknife procedure) with conversion to liability scale
    project$mvLD$covstruct.origHDL.liabilityScale<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HLD.ld, liabilityScale = TRUE)
    
    saveRDS(object = project$mvLD$covstruct.origHDL.liabilityScale,file = project$filepath.mvLD.origHDL.liabilityScale)
    quit(save = "no")
  }
  
  #save the mvLD output
  saveRDS(object = project$mvLD,file = project$filepath.mvLD)
  print("Multivariate LD correction is done now and the resulting covariance structure should have been saved into a file.")
  
  #remove intermediate results
  #file.remove(project$filepath.mvLD.HDL.piecewise)
  #file.remove(project$filepath.mvLD.HDL.jackknife)
  #file.remove(project$filepath.mvLD.mvLDSC)
  
}

#Additional computations

#missing column names - should be fixed in the GSEM method by now.
#colnames(project$mvLD$covstruct.origHDL$S) <- project$sumstats.sel$code

#retrieve the standard errors of S (variances and covariances) from the diagonal of V (contains both).
  project$mvLD$covstruct.mvLDSC$S.SE<-matrix(0, project$k.sel, project$k.sel)
  project$mvLD$covstruct.mvLDSC$S_Stand.SE<-matrix(0, project$k.sel, project$k.sel)
  project$mvLD$covstruct.mvLDSC$S.SE[lower.tri(project$mvLD$covstruct.mvLDSC$S.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.mvLDSC$V))
  project$mvLD$covstruct.mvLDSC$S_Stand.SE[lower.tri(project$mvLD$covstruct.mvLDSC$S_Stand.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.mvLDSC$V_Stand))
  
  project$mvLD$covstruct.HDL.piecewise$S.SE<-matrix(0, project$k.sel, project$k.sel)
  project$mvLD$covstruct.HDL.piecewise$S_Stand.SE<-matrix(0, project$k.sel, project$k.sel)
  project$mvLD$covstruct.HDL.piecewise$S.SE[lower.tri(project$mvLD$covstruct.HDL.piecewise$S.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.HDL.piecewise$V))
  project$mvLD$covstruct.HDL.piecewise$S_Stand.SE[lower.tri(project$mvLD$covstruct.HDL.piecewise$S_Stand.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.HDL.piecewise$V_Stand))
  
  project$mvLD$covstruct.HDL.jackknife$S.SE<-matrix(0, project$k.sel, project$k.sel)
  project$mvLD$covstruct.HDL.jackknife$S_Stand.SE<-matrix(0, project$k.sel, project$k.sel)
  project$mvLD$covstruct.HDL.jackknife$S.SE[lower.tri(project$mvLD$covstruct.HDL.jackknife$S.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.HDL.jackknife$V))
  project$mvLD$covstruct.HDL.jackknife$S_Stand.SE[lower.tri(project$mvLD$covstruct.HDL.jackknife$S_Stand.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.HDL.jackknife$V_Stand))
  
#add newly computed heritabilities to the selected summary statistics table
project$sumstats.sel$h2.liability_mvLDSC<-diag(project$mvLD$covstruct.mvLDSC$S)
project$sumstats.sel$h2.se.liability_mvLDSC<-diag(project$mvLD$covstruct.mvLDSC$S.SE)
project$sumstats.sel$h2.liability_HDL.piecewise<-diag(project$mvLD$covstruct.HDL.piecewise$S)
project$sumstats.sel$h2.se.liability_HDL.piecewise<-diag(project$mvLD$covstruct.HDL.piecewise$S.SE) 
project$sumstats.sel$h2.liability_HDL.jackknife<-diag(project$mvLD$covstruct.HDL.jackknife$S)
project$sumstats.sel$h2.se.liability_HDL.jackknife<-diag(project$mvLD$covstruct.HDL.jackknife$S.SE)
project$sumstats.sel$h2.observed_origHDL<-diag(project$mvLD$covstruct.origHDL$S)
project$sumstats.sel$h2.se.observed_origHDL<-diag(project$mvLD$covstruct.origHDL$S.se)
project$sumstats.sel$h2.liability_origHDL<-diag(project$mvLD$covstruct.origHDL.liabilityScale$S)
#project$sumstats.sel$h2.se.liability_origHDL<-diag(project$mvLD$covstruct.origHDL.liabilityScale$S.se)

# #add formatted heritability strings
# project$sumstats.sel<-project$sumstats.sel %>% mutate(
#   h2.liability_mvLDSC_s = paste0(h2.liability_mvLDSC,"(",h2.se.liability_mvLDSC,")")
# )

#project$mvLD$covstruct.HDL.jackknife$complete

#View(project$sumstats.sel)

  






## ----improved annotation of chosen datasets, fig.width=10, fig.height=6, out.width="1600px", out.height="1000px"----

project$sumstats.sel.table<-project$sumstats.sel[,c("gwas_name.nice","dependent_variable","n_case","n_control","reference_year","samplePrevalence","populationPrevalence","h2.liability_mvLDSC","h2.se.liability_mvLDSC","h2.liability_HDL.piecewise","h2.se.liability_HDL.piecewise","h2.observed_origHDL","h2.se.observed_origHDL","h2.liability_origHDL")]
rownames(project$sumstats.sel.table)<-NULL #Remove the rowname column
#View(project$sumstats.sel.table)
#project$sumstats.sel.table

project$plots.sumstats.sel.table<-project$sumstats.sel.table %>% 
  gt() %>% 
  fmt_number(columns = vars(samplePrevalence, populationPrevalence), decimals = 2) %>%
  fmt_number(columns = vars(h2.liability_mvLDSC,h2.liability_HDL.piecewise,h2.observed_origHDL,h2.liability_origHDL,), decimals = 3) %>%
  fmt_number(columns = vars(h2.se.liability_mvLDSC,h2.se.liability_HDL.piecewise,h2.se.observed_origHDL), decimals = 4) %>%
  fmt_number(columns = vars(n_case,n_control), decimals = 0) %>%
  tab_header(
    title = "Selected GWAS summary statistics datasets"
  ) %>% cols_label(
    gwas_name.nice = "Trait",
    #ancestry = "Ancestry",
    #sex = "Sex",
    dependent_variable = "Dependent variable",
    n_case = "N case",
    n_control = "N control",
    reference_year = "Year",
    samplePrevalence = html("Prevalence<sub>sample</sub>"),
    populationPrevalence = html("Prevalence<sub>population</sub>"),
    h2.liability_mvLDSC = html("h<sup>2</sup><sub>mvLDSC</sub>"),
    h2.se.liability_mvLDSC = "se",
    h2.liability_HDL.piecewise = html("h<sup>2</sup><sub>HDL(pw)</sub>"),
    h2.se.liability_HDL.piecewise = "se",
    h2.observed_origHDL = html("h<sup>2</sup><sub>oHDL,OS</sub>"),
    h2.se.observed_origHDL = "se",
    h2.liability_origHDL = html("h<sup>2</sup><sub>oHDL,LS</sub>")
  ) %>%
  tab_style(
    style = cell_text(size = px(12)),
    locations = cells_column_labels(everything())       
  ) %>%
  tab_style(
    style = cell_text(size = px(12),weight = "bold"),
    locations = cells_body(everything())        
  )

project$plots.sumstats.sel.table

gtsave(data = project$plots.sumstats.sel.table, filename = paste0(project$folderpath.plots,"/sumstats.sel.table.rtf"))
















## ----prepare summary statistics----------------------------------------------------------------------------------
#eval=FALSE
#use ^this to knit without running the code in the chunk

#we need to introduce checks of the summary statistics and which scale they are on
# f.ex. Leo checks my s.e. 2 * pnorm(log(1.003) / 0.021644, mean = 0, lower.tail = FALSE)

#View(project$sumstats.sel)
#examine GWAS sumstats
# gwas_ALCD03<-read.table(file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/ALCD03_noMHC.sumstats.gz", header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# View(gwas_ALCD03)
# gwas_ANXI03<-read.table(file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/ANXI03_noMHC.sumstats.gz", header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# View(gwas_ANXI03)
# gwas_DEPR05<-read.table(file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/DEPR05_noMHC.sumstats.gz", header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# View(gwas_DEPR05)
# gwas_HEAL01<-read.table(file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/HEAL01_noMHC.sumstats.gz", header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# View(gwas_HEAL01)
# gwas_NEUR01<-read.table(file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/NEUR01_noMHC.sumstats.gz", header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# View(gwas_NEUR01)
# gwas_SUBJ01<-read.table(file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/SUBJ01_noMHC.sumstats.gz", header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# View(gwas_SUBJ01)
# gwas_TIRE01<-read.table(file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/TIRE01_noMHC.sumstats.gz", header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
# View(gwas_TIRE01)
# 
# gwas_ALCD03$SE<-NA_real_
# gwas_ANXI03$SE<-NA_real_
# gwas_DEPR05$SE<-NA_real_
# gwas_HEAL01$SE<-NA_real_
# gwas_NEUR01$SE<-NA_real_
# gwas_SUBJ01$SE<-NA_real_
# gwas_TIRE01$SE<-NA_real_
# 
# write.table(x = gwas_ALCD03, file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/ALCD03_noMHC_mod.sumstats", quote=TRUE )
# write.table(x = gwas_ANXI03, file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/ANXI03_noMHC_mod.sumstats", quote=TRUE )
# write.table(x = gwas_DEPR05, file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/DEPR05_noMHC_mod.sumstats", quote=TRUE )
# write.table(x = gwas_HEAL01, file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/HEAL01_noMHC_mod.sumstats", quote=TRUE )
# write.table(x = gwas_NEUR01, file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/NEUR01_noMHC_mod.sumstats", quote=TRUE )
# write.table(x = gwas_SUBJ01, file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/SUBJ01_noMHC_mod.sumstats", quote=TRUE )
# write.table(x = gwas_TIRE01, file = "/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data.sumstats.mungedNoMHC/TIRE01_noMHC_mod.sumstats", quote=TRUE )

project$lfGWAS<-c()

if (file.exists(paste0(project$folderpath.workingDirectory,"/","lfGWAS.sumstats.",project$setup.code,".Rds")) | project$seting.refreshPrepareSummaryStatistics) {
  project$lfGWAS<-readRDS(file=paste0(project$folderpath.workingDirectory,"/","lfGWAS.sumstats.",project$setup.code,".Rds"))
} else
{
  print("Preparing summary statistics for latent factor GWAS. This might take a while.")

  project$lfGWAS$preparedSumstats<-c()
# for(cCode in project$sumstats.sel$code)
#     project$lfGWAS$trait<-c(project$lfGWAS$trait,paste0(project$folderpath.data.sumstats.munged,"/",cCode,project$filename.suffix.data.sumstats.munged,".gz"))

project$lfGWAS$sumstats.prepared<-sumstats(
  files=project$sumstats.sel$cleanedpath,
  ref=project$filepath.genomeReference,
  trait.names=project$sumstats.sel$code,
  se.logit=project$sumstats.sel$se.logit,
  OLS=project$sumstats.sel$dependent_variable.OLS,
  linprob=NULL, #THIS SHOULD BE INVESTIGATED FURTHER, IF A LINEAR OLS ESTIMATOR ON A DICHOTOMOUS DEP. VARIABLE WAS USED FOR ANY OF THE DATSETS 
  prop=NULL,
  N=project$sumstats.sel$n_total,
  info.filter=project$info.filter,
  maf.filter=project$maf.filter,
  keep.indel=FALSE,
  parallel=FALSE,
  cores=NULL
  )

  saveRDS(object = project$lfGWAS,file = paste0(project$folderpath.workingDirectory,"/","lfGWAS.sumstats.",project$setup.code,".Rds"))
  print("Done summary statistics for latent factor GWAS. The result should have been saved to a file.")
}



## ----latent factor GWAS------------------------------------------------------------------------------------------

if(is.null(project$lfGWAS$userGWAS.correlated) | project$setting.refreshLatentFactorGWAS) {

print("Performing latent factor GWAS. This will take a while.")

#project$CFA$cfa.configurations.result[which(project$CFA$cfa.configurations.result$id==13),]
  
project$lfGWAS$modelM2fP_0.3<-"
F1 =~ NA*ALCD+ANXI+DEPR+NEUR+SUBJ
F2 =~ NA*ALCD+HEAL+TIRE
F1~~1*F1
F2~~1*F2

F1~~F2

F1~SNP

"

project$lfGWAS$modelTemplate<-"
F1 =~ NA*ALCD+ANXI+DEPR+NEUR+SUBJ
F2 =~ NA*ALCD+HEAL+TIRE
F1~~1*F1
F2~~1*F2

F1~~F2

F1~SNP
F2~SNP

"

project$lfGWAS$gwas<-userGWAS(covstruc = project$mvLD$covstruct.mvLDSC, SNPs = project$lfGWAS$sumstats.prepared, estimation = "ML", model = project$lfGWAS$modelM2fP_0.3, modelchi = FALSE, printwarn = TRUE, sub=c("F1~SNP"), GC="standard")

#project$lfGWAS$gwas<-userGWAS(covstruc = project$mvLD$covstruct.mvLDSC, SNPs = project$lfGWAS$sumstats.prepared, estimation = "ML", model = project$lfGWAS$modelM2fP_0.3, modelchi = FALSE, printwarn = TRUE, sub=c("F1~SNP","F2~SNP"), GC="standard")

saveRDS(object = project$lfGWAS$gwas,file = paste0(project$folderpath.workingDirectory,"/","lfGWAS.gwas.",project$setup.code,".Rds"))

print("DONE performing latent factor GWAS. The results should have been saved to a file.")

}



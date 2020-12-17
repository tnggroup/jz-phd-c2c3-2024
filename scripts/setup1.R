## ----package setup, echo=FALSE, warning=F-------------------------------------------------------------------------

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

library(skimr)
library(psych)
library(Matrix)

library(tidyverse)
library(ggrepel)
library(gt)
#library(kableExtra)
#library(corrplot) #do not run on cluster, does not work?
library(GenomicSEM)
library(HDL)
library(optparse)



## ----command line setup-------------------------------------------------------------------------------------------
clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="integer", default=1,
                help="Index of the explicit task to run separately:\n1: No task\n2:HDL Piecewise\n3:HDL Jackknife\n4:vLDSC\n5:original piecewise HDL\n6:original jackknife HDL [default %default]")
clParser <- add_option(clParser, c("-l", "--location"), type="character", default="local",
                help="The place where the code is run [local,cluster] [default %default]")



## ----settings-----------------------------------------------------------------------------------------------------
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







## ----additional source setup, echo=FALSE, warning=F---------------------------------------------------------------

source(normalizePath(paste0(project$folderpath.scripts,"/","shru.R")))
#source(normalizePath(paste0(project$folderpath.scripts,"/","hdl.mod.R")))



## ----trait setup--------------------------------------------------------------------------------------------------
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



## ----GWAS sumstat dataset setup-----------------------------------------------------------------------------------
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



## ----GWAS sumstat dataset variable selection----------------------------------------------------------------------

#selection based on specific traits
project$sumstats.sel.code<-c("DEPR05","ANXI03","NEUR01","TIRE01","SUBJ01","ALCD03","HEAL01")
project$sumstats.sel<-project$sumstats[which(project$sumstats$code %in% project$sumstats.sel.code),]
project$sumstats.sel$code<-project$sumstats.sel$code.trait
project$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")]

#View(project$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")])



## ----multivariate LD----------------------------------------------------------------------------------------------

project$filepath.mvLD<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".Rds")
project$filepath.mvLD.HDL.piecewise<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".HDL.piecewise.Rds")
project$filepath.mvLD.HDL.jackknife<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".HDL.jackknife.Rds")
project$filepath.mvLD.mvLDSC<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".mvLDSC.Rds")
project$filepath.mvLD.origHDL.piecewise<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL.piecewise.Rds")
project$filepath.mvLD.origHDL.jackknife<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL.jackknife.Rds")
project$filepath.mvLD.origHDL_liab.piecewise<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL_liab.piecewise.Rds")
project$filepath.mvLD.origHDL_liab.jackknife<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL_liab.jackknife.Rds")


if (file.exists(project$filepath.mvLD)) {
  print("Using existing covariance structure from previous HDL computation.")
  project$mvLD<-readRDS(file=project$filepath.mvLD)
} else {
  print("Running multivariate LD regression with HDL. This might take a while. If this runs for too long you may want to abort the process.")
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
    project$mvLD$covstruct.origHDL.jackknife<-readRDS(file=project$filepath.mvLD.origHDL.jackknife)
    project$mvLD$covstruct.origHDL.piecewise<-readRDS(file=project$filepath.mvLD.origHDL.piecewise)
    
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
  
  if(project$clOptions$task==5 && !file.exists(project$filepath.mvLD.origHDL.piecewise)){
    #run HDL original piecewise
    
    
    
    project$mvLD$covstruct.origHDL.piecewise<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HLD.ld, jackknife = FALSE)
    
    saveRDS(object = project$mvLD$covstruct.origHDL.piecewise,file = project$filepath.mvLD.origHDL.piecewise)
    quit(save = "no")
  }
  
  if(project$clOptions$task==6 && !file.exists(project$filepath.mvLD.origHDL.jackknife)){
    #run HDL original jackknife
    
    project$mvLD$covstruct.origHDL.jackknife<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HLD.ld, jackknife = TRUE)
    
    saveRDS(object = project$mvLD$covstruct.origHDL.jackknife,file = project$filepath.mvLD.origHDL.jackknife)
    quit(save = "no")
  }
  
   if(project$clOptions$task==7 && !file.exists(project$filepath.mvLD.origHDL_liab.piecewise)){
    #run HDL original piecewise using liability scale
    
    project$mvLD$covstruct.origHDL_liab.piecewise<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HLD.ld, jackknife = FALSE, liabilityScale=TRUE)
    
    saveRDS(object = project$mvLD$covstruct.origHDL_liab.piecewise,file = project$filepath.mvLD.origHDL_liab.piecewise)
    quit(save = "no")
   }
  
  if(project$clOptions$task==8 && !file.exists(project$filepath.mvLD.origHDL_liab.jackknife)){
    #run HDL original jackknife using liability scale
    
    project$mvLD$covstruct.origHDL_liab.jackknife<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HLD.ld, jackknife = TRUE, liabilityScale=TRUE)
    
    saveRDS(object = project$mvLD$covstruct.origHDL_liab.jackknife,file = project$filepath.mvLD.origHDL_liab.jackknife)
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

#some additional calculations

#prouce the standard errors of S (variances and covariances) from the diagonal of V (contains both).
  project$mvLD$covstruct.mvLDSC$S.k<-nrow(project$mvLD$covstruct.mvLDSC$S)
  project$mvLD$covstruct.mvLDSC$S.SE<-matrix(0, project$mvLD$covstruct.mvLDSC$S.k, project$mvLD$covstruct.mvLDSC$S.k)
  project$mvLD$covstruct.mvLDSC$S.SE[lower.tri(project$mvLD$covstruct.mvLDSC$S.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.mvLDSC$V))
  project$mvLD$covstruct.mvLDSC$S_Stand.SE<-matrix(0, project$mvLD$covstruct.mvLDSC$S.k, project$mvLD$covstruct.mvLDSC$S.k)
  project$mvLD$covstruct.mvLDSC$S_Stand.SE[lower.tri(project$mvLD$covstruct.mvLDSC$S_Stand.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.mvLDSC$V_Stand))
  
  # project$mvLD$covstruct.HDL.piecewise$S.k<-nrow(project$mvLD$covstruct.HDL.piecewise$S)
  # project$mvLD$covstruct.HDL.piecewise$S.SE<-matrix(0, project$mvLD$covstruct.HDL.piecewise$S.k, project$mvLD$covstruct.HDL.piecewise$S.k)
  # project$mvLD$covstruct.HDL.piecewise$S.SE[lower.tri(project$mvLD$covstruct.HDL.piecewise$S.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.HDL.piecewise$V))
  # project$mvLD$covstruct.HDL.piecewise$S_Stand.SE<-matrix(0, project$mvLD$covstruct.HDL.piecewise$S.k, project$mvLD$covstruct.HDL.piecewise$S.k)
  # project$mvLD$covstruct.HDL.piecewise$S_Stand.SE[lower.tri(project$mvLD$covstruct.HDL.piecewise$S_Stand.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.HDL.piecewise$V_Stand))

#add newly computed heritabilities to the selected summary statistics table
for(iTrait in 1:nrow(project$sumstats.sel)) {
  project$sumstats.sel$h2.liability[iTrait]<-project$mvLD$covstruct.HDL.jackknife$S[[iTrait,iTrait]]
  # for(iTrait2 in 1:iTrait) {
  #   mvLDSC$cor.test[[iTrait,iTrait2]]<-cortest(R1 = mvLDSC$output$S[[iTrait,iTrait2]], n1 = mvLDSC$output$N[iTrait,iTrait2],)
  # }
}





## ----improved annotation of chosen datasets, fig.width=9, fig.height=6, out.width="1600px", out.height="1000px"----
#skim(project$sumstats.sel)
project$sumstats.sel.table<-project$sumstats.sel[,c("gwas_name.nice","dependent_variable","n_total","reference_year","samplePrevalence","populationPrevalence","h2.liability")]
rownames(project$sumstats.sel.table)<-NULL #Remove the rowname column

#project$sumstats.sel.table

project$plots.sumstats.sel.table<-project$sumstats.sel.table %>% 
  gt() %>% 
  fmt_number(columns = vars(samplePrevalence, populationPrevalence), decimals = 2) %>%
  fmt_number(columns = vars(h2.liability), decimals = 3) %>%
  fmt_number(columns = vars(n_total), decimals = 0) %>%
  tab_header(
    title = "Selected GWAS summary statistics datasets"
  ) %>% cols_label(
    gwas_name.nice = "Trait",
    #ancestry = "Ancestry",
    #sex = "Sex",
    dependent_variable = "Dependent variable",
    n_total = "Total number of participants",
    reference_year = "Reference year",
    samplePrevalence = "Sample prevalence",
    populationPrevalence = "Population prevalence",
    h2.liability = html("h<sup>2</sup><sub>SNP,liability scale</sub>")
    
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




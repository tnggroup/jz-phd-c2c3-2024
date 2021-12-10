## ----package setup, echo=FALSE, warning=F------------------------------------------------------------------------------

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
#devtools::install_github("johanzvrskovec/GenomicSEM",ref = 'mod-jz') #specify branch
#remove.packages("HDL")
#devtools::install_github("zhenin/HDL/HDL")
#devtools::install_github("zhenin/HDL/HDL@77cb9d0984d1302e40bfd871491e292f8f09f49d") #specify exact commit
#install.packages("optparse")
#install.packages("stats")
#remove.packages("shru")
#devtools::install_github("johanzvrskovec/shru",ref = 'main') #specify branch
#install.packages("reticulate")
#install.packages("readr")

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
library(shru)

#library(reticulate)



## ----command line setup------------------------------------------------------------------------------------------------
clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="character", default="0",
                help="Index of the explicit task to run separately:\n0: No task\nmvLD.mvLDSC:multivariate LDSC\nmvLD.HDL.piecewise:HDL Piecewise\nmvLD.HDL.jackknife:HDL Jackknife\nmvLD.origHDL:original HDL(jackknife)\nmvLD.origHDL.liabilityScale:original HDL with applied liability scale [default %default]")
clParser <- add_option(clParser, c("-l", "--location"), type="character", default="local",
                help="The place where the code is run [local,cluster] [default %default]")

clParser <- add_option(clParser, c("-a", "--task_argument"), type="character", default=NA,
                help="General purpose argument for tasks [default %default]")



## ----settings----------------------------------------------------------------------------------------------------------
project<-c() #create project metadata object
project$clOptions<-parse_args(clParser)
project$date.run<-Sys.Date()
project$setup.version<-2
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
project$folderpath.workingDirectory<-normalizePath(paste0(project$folderpath,"/","working_directory"))

project$folderpath.scripts<-normalizePath(paste0(project$folderpath,"/","scripts"))
project$folderpath.includedSoftware<-normalizePath(paste0(project$folderpath,"/","included_software"))
project$folderpath.plots<-normalizePath(paste0(project$folderpath,"/","plots"))


#general data folder
if(project$host=="local") {
  project$folderpath.data<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data")
} else if (project$host=="cluster") {
  project$folderpath.data<-normalizePath(paste0(project$folderpath,"/","data"))
}

##cleaned sumstats folder
if(project$host=="local") {
  project$folderpath.data.sumstats.cleaned<-normalizePath(paste0(project$folderpath.data,"/gwas_sumstats/cleaned"))
} else if (project$host=="cluster") {
  project$folderpath.data.sumstats.cleaned<-normalizePath("/mnt/lustre/groups/ukbiobank/sumstats/cleaned")
}

##munged sumstats folder
project$folderpath.data.sumstats.munged<-normalizePath(paste0(project$folderpath.data,"/gwas_sumstats/munged"))

##imputed sumstats folder
project$folderpath.data.sumstats.imputed<-normalizePath(paste0(project$folderpath.data,"/gwas_sumstats/imputed"))

#python virtual environment folder
if(project$host=="local") {
  project$folderpath.pythonVenv<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/python-venv")
} else if (project$host=="cluster") {
  project$folderpath.pythonVenv<-normalizePath(paste0(project$folderpath,"/","python-venv"))
}

##Reference SNP-list (HapMap3 SNPs for example). Used for munging sumstat SNP data.
#project$filepath.SNPReference<-normalizePath(paste0(project$folderpath.data,"/","w_hm3.snplist.flaskapp2018")) #HapMap3 SNPs
## Used in the preparation step for performing latent factor GWAS as reference for calculating SNP variance across traits.
#project$filepath.SNPReference<-normalizePath(paste0(project$folderpath.data,"/","reference.1000G.maf.0.005.txt")) #1000 genomes phase 3
project$filepath.SNPReference<-normalizePath(paste0(project$folderpath.data,"/","combined.hm3_1kg.snplist.jz2020.txt")) #custom hm3 + 1kg SNPs

project$filename.suffix.data.sumstats.munged<-".gz"
#project$filename.suffix.data.sumstats.munged<-"_noMHC.sumstats.gz"

##Reference panel folder containing individual level data reference panel. Used for GWAS sumstat imputation tasks.
#roject$folderpath.data.sumstatImp.genomeReference<-"/users/k1204688/brc_scratch/Public/1KG_Phase3/All"
project$folderpath.data.sumstatImp.genomeReference<-"/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/reference.panel.1KG_Phase3.CLEANED.EUR.cM"

##LD scores datasets folders (these strings need to have a trailing slash for the GSEM LDSC to work)
project$folderpath.data.mvLDSC.ld <- paste0(project$folderpath.data,"/reference.panel.1KG_Phase3.CLEANED.EUR.cM/") #LD-scores
project$folderpath.data.mvLDSC.wld <- paste0(project$folderpath.data,"/reference.panel.1KG_Phase3.CLEANED.EUR.cM/") #Weights, if different from LD-scores

##HDL LD scores reference
if(project$host=="local") {
  #use the smallest LD reference as default for local tests
  project$folderpath.data.HDL.ld<-paste0(project$folderpath.data,"/UKB_array_SVD_eigen90_extraction/")
} else if (project$host=="cluster") {
  project$folderpath.data.HDL.ld<-paste0(project$folderpath.data,"/UKB_imputed_hm3_SVD_eigen99_extraction/")
}
  
##full script file paths
project$filepath.rmd<-normalizePath(paste0(project$folderpath.scripts,"/",project$filename.rmd))
project$filepath.r<-normalizePath(paste0(project$folderpath.scripts,"/",project$filename.r))

##CFA settings
project$CFA<-c()
project$CFA$estimator=c("ML")
#project$CFA$estimator=c("ML","DWLS")
project$CFA$nFactors=3

##latent factor GWAS filter settings
project$lfGWAS$info.filter=.6
project$lfGWAS$maf.filter=0.01

#working directory in case of running as an R-script
setwd(dir = normalizePath(project$folderpath.workingDirectory))

#inactivated python environment until it is used
#use_virtualenv(project$folderpath.pythonVenv)







## ----additional source setup, echo=FALSE, warning=F--------------------------------------------------------------------

source(normalizePath(file.path(project$folderpath.scripts,"sumstats.mod-jz.R")))





## ----trait setup-------------------------------------------------------------------------------------------------------
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



## ----GWAS sumstat dataset setup----------------------------------------------------------------------------------------
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
  mungedpath=paste0(project$folderpath.data.sumstats.munged,"/",code,project$filename.suffix.data.sumstats.munged)
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



## ----GWAS sumstat dataset variable selection---------------------------------------------------------------------------

#selection based on specific traits
project$sumstats.sel.code<-c("DEPR05","ANXI03","NEUR01","TIRE01","SUBJ01","ALCD03","HEAL01")
project$sumstats.sel<-project$sumstats[which(project$sumstats$code %in% project$sumstats.sel.code),]
#project$sumstats.sel$code_orig<-project$sumstats.sel$code
#project$sumstats.sel$code<-project$sumstats.sel$code.trait
project$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")]
project$k.sel<-nrow(project$sumstats.sel)
#View(project$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")])


write.table(project$sumstats.sel[,c("code", "gwas_name","reference_year", "n_case","n_control","n_total","samplePrevalence","populationPrevalence", "reference_doi")], file = file.path(project$folderpath.workingDirectory,paste0("project.",project$setup.code,".sumstatinfo.tsv")), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)



## ----GWAS sumstat munge------------------------------------------------------------------------------------------------
if(project$clOptions$task=="munge"){

  
  #construct new reference
  #library(readr)
  #project$sumstats.reference.old<-read.table(project$filepath.SNPReference, header=T, quote="\"", fill=T, blank.lines.skip=T, strip.white = T,na.strings=c(".",NA,"NA",""))
  #project$sumstats.reference.old<-data.table::fread(project$filepath.SNPReference,header=T,data.table=F)
  
  #refaddition<-as.data.frame(data.table::fread("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/munged/ref.snps.not.in.REF1KG.txt",header=F,data.table=F))
  #refaddition<-data.table::fread(normalizePath(paste0(project$folderpath.data,"/","w_hm3.snplist.flaskapp2018")),header=T,data.table=F)
  #colnames(refaddition)<-c("SNP")
  
  #diff<-setdiff(refaddition[,c("SNP","A1","A2")],project$sumstats.reference.old[,c("SNP","A1","A2")])
  #diff<-data.frame(SNP=setdiff(refaddition[,c("SNP")],project$sumstats.reference.old[,c("SNP")]),CHR=NA,BP=NA,MAF=NA,A1=NA,A2=NA)
  #nrow(diff)
  #project$sumstats.reference.new<-rbind(project$sumstats.reference.old,diff)
  #project$sumstats.reference.new<-unique(project$sumstats.reference.new2)
  #View(project$sumstats.reference.old)
  #View(project$sumstats.reference.new)
  #View(project$sumstats.reference.new2)
  #View(project$sumstats.reference.new3)
  
  

  #write.table(x = project$sumstats.reference.new,file = file.path(project$folderpath.workingDirectory,"/","combined.hm3_1kg.snplist.jz2020.txt"), quote = FALSE, row.names = F)
  #faster alt - produces a comma separated version?
  #data.table::fwrite(x = project$sumstats.reference.new, file = )
  
  #read in again as
  #project$sumstats.reference.new<-read.table(project$filepath.SNPReference, header=T, quote="\"", fill=T, blank.lines.skip=T, strip.white = T,na.strings=c(".",NA,"NA",""))
  
  #special munge of the 1KG reference SNP stats
  # munge.mod(files = c(project$filepath.genomeReference),
  #           hm3 = project$filepath.SNPReference,
  #           trait.names=c("REF1KG"),
  #           info.filter=project$info.filter,
  #           maf.filter=project$maf.filter,
  #           path.dir.output = project$folderpath.data.sumstats.munged,
  #           N = 2504
  #             ) 
  
  
  #experiment to load extended reference panel
  #Folder:
  #"/users/k19049801/project/JZ_GED_PHD_C1/data/UKB_imputed_hm3_SVD_eigen99_extraction"
  #Files:
  #UKB_snp_counter_imputed.RData
  #UKB_snp_list_imputed.vector_form.RData
  #snp.dictionary.imputed.rda
  
  #munging with no filters applied
  munge.mod(files = project$sumstats.sel$cleanedpath,
            hm3 = project$filepath.SNPReference,
            trait.names=project$sumstats.sel$code,
            info.filter=0,
            maf.filter=0,
            path.dir.output = project$folderpath.data.sumstats.munged,
            doChrSplit = FALSE
              ) 
    
    
    quit(save = "no")
  }



## ----GWAS sumstat imputation-------------------------------------------------------------------------------------------
if(project$clOptions$task=="impute"){
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



## ----multivariate LD---------------------------------------------------------------------------------------------------

project$filepath.mvLD<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".Rds")
project$filepath.mvLD.mvLDSC<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".mvLDSC.Rds")
project$filepath.mvLD.HDL.piecewise<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".HDL.piecewise.Rds")
project$filepath.mvLD.HDL.jackknife<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".HDL.jackknife.Rds")
project$filepath.mvLD.origHDL<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL.Rds")
project$filepath.mvLD.origHDL.liabilityScale<-paste0(project$folderpath.workingDirectory,"/","mvLD.",project$setup.code,".origHDL.liabilityScale.Rds")


if (file.exists(project$filepath.mvLD)) {
  print("Using existing covariance structure from previous HDL computation.")
  project$mvLD<-readRDS(file=project$filepath.mvLD)
} else {
  print("Running (or reading ready intermediate results from) multivariate LD regression with HDL. This might take a while. If the procedure runs for too long you may want to abort the process.")
  
  cat("The current task is specified as:",project$clOptions$task)
  project$mvLD<-c()
  
  
  
  if(project$clOptions$task=="mvLD.mvLDSC" && !file.exists(project$filepath.mvLD.mvLDSC)){
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
  
  if(project$clOptions$task=="mvLD.HDL.piecewise" && !file.exists(project$filepath.mvLD.HDL.piecewise)){
    #run HDL piecewise
    #system(command = paste0("touch 2.txt"))
    project$mvLD$covstruct.HDL.piecewise <- hdl.mod(traits = project$sumstats.sel$mungedpath,
                                  sample.prev = project$sumstats.sel$samplePrevalence,
                                  population.prev = project$sumstats.sel$populationPrevalence,
                                  trait.names=project$sumstats.sel$code,
                                  LD.path=project$folderpath.data.HDL.ld,
                                  method = "piecewise")
    
    saveRDS(object = project$mvLD$covstruct.HDL.piecewise,file = project$filepath.mvLD.HDL.piecewise)
    quit(save = "no")
  }
  
  if(project$clOptions$task=="mvLD.HDL.jackknife" && !file.exists(project$filepath.mvLD.HDL.jackknife)){
    #run HDL jackknife
    #system(command = paste0("touch 3.txt"))
    project$mvLD$covstruct.HDL.jackknife <- hdl.mod(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HDL.ld,
                                method = "jackknife")
  
    saveRDS(object = project$mvLD$covstruct.HDL.jackknife,file = project$filepath.mvLD.HDL.jackknife)
    quit(save = "no")
  }
  
  
  if(project$clOptions$task=="mvLD.origHDL" && !file.exists(project$filepath.mvLD.origHDL)){
    #run HDL original (using the original jackknife procedure)
    project$mvLD$covstruct.origHDL<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HDL.ld, liabilityScale = FALSE)

    saveRDS(object = project$mvLD$covstruct.origHDL,file = project$filepath.mvLD.origHDL)
    quit(save = "no")
  }
  
  if(project$clOptions$task=="mvLD.origHDL.liabilityScale" && !file.exists(project$filepath.mvLD.origHDL.liabilityScale)){
    #run HDL original (using the original jackknife procedure) with conversion to liability scale
    project$mvLD$covstruct.origHDL.liabilityScale<-hdl.original(traits = project$sumstats.sel$mungedpath,
                                sample.prev = project$sumstats.sel$samplePrevalence,
                                population.prev = project$sumstats.sel$populationPrevalence,
                                trait.names=project$sumstats.sel$code,
                                LD.path=project$folderpath.data.HDL.ld, liabilityScale = TRUE)
    
    saveRDS(object = project$mvLD$covstruct.origHDL.liabilityScale,file = project$filepath.mvLD.origHDL.liabilityScale)
    quit(save = "no")
  }
  
  if(project$clOptions$task=="0"){
    
    if(!file.exists(project$filepath.mvLD.HDL.piecewise) || !file.exists(project$filepath.mvLD.HDL.jackknife) || !file.exists(project$filepath.mvLD.mvLDSC)) {
      Sys.sleep(time = 3)
    }
    
    cat("Reading ready intermediate mvLD results and saving final mvLD result.")
    project$mvLD$covstruct.mvLDSC<-readRDS(file=project$filepath.mvLD.mvLDSC)
    project$mvLD$covstruct.HDL.piecewise<-readRDS(file=project$filepath.mvLD.HDL.piecewise)
    project$mvLD$covstruct.HDL.jackknife<-readRDS(file=project$filepath.mvLD.HDL.jackknife)
    project$mvLD$covstruct.origHDL<-readRDS(file=project$filepath.mvLD.origHDL)
    project$mvLD$covstruct.origHDL.liabilityScale<-readRDS(file=project$filepath.mvLD.origHDL.liabilityScale)
    
    #save the mvLD output
    saveRDS(object = project$mvLD,file = project$filepath.mvLD)
    print("Multivariate LD correction is done now and the resulting covariance structure should have been saved to a file.")
    
    #remove intermediate results
    #file.remove(project$filepath.mvLD.HDL.piecewise)
    #file.remove(project$filepath.mvLD.HDL.jackknife)
    #file.remove(project$filepath.mvLD.mvLDSC)
    
  }
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

  










## ----CFA---------------------------------------------------------------------------------------------------------------

#THESE USE shru::semplate and DiagrammeR !!!!!

#library(DiagrammeR)

#test
# lavaanDefinition<-semplate$generateLavaanCFAModel(code.indicator = project$sumstats.sel$code, allow_loading.table.indicator_factor = data.frame(
#   F1=c(TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE),
#   F2=c(TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
#   F3=c(FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE)),
#   fix_loading.table.indicator_factor = data.frame(
#   F1=c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE),
#   F2=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
#   F3=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)                                                                           ),
#   special.orthogonal = FALSE
#   )

# lavaanDefinitionSafe<-"F1 =~ NA*ALCD03+ANXI03+DEPR05+HEAL01+NEUR01+SUBJ01+TIRE01
# F1~~1*F1"
# 
# lavaanDefinitionSafe2 <- "
# F1 =~ NA*ALCD03+ANXI03+DEPR05+NEUR01+SUBJ01
# F2 =~ NA*ALCD03+HEAL01+TIRE01
# F1~~1*F1
# F2~~1*F2
# F1~~F2"

# lavaanDefinition<-semplate$generateLavaanCFAModel(
#   allow_loading.table.indicator_factor = data.frame(
#     F1=c(TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE),
#     F2=c(TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE),
#     F3=c(FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE)
#     ),
#   indicatorArgs = data.frame(codeIndicator=project$sumstats.sel$code,residualSizeLimit=NA),
#   orthogonal = FALSE
#   )

# lavaanDefinitionManualSafe="
# F1 =~ NA*ALCD03+ANXI03+DEPR05+NEUR01+SUBJ01
# F2 =~ NA*ALCD03+HEAL01+TIRE01
# F1~~1*F1
# F2~~1*F2
# F1~~F2"

# lavaanDefinitionManual='
# #F1 =~ NA*ANXI03+DEPR05+label(HEAL01.res)*HEAL01+start(0.5)*HEAL01+SUBJ01+TIRE01
# F1 =~ NA*ANXI03+DEPR05+label(HEAL01.res)*HEAL01+SUBJ01+TIRE01
# F2 =~ NA*ANXI03+DEPR05+HEAL01+NEUR01+SUBJ01
# F3 =~ NA*ALCD03+ANXI03+DEPR05
# 
# #HEAL01 ~~ label(HEAL01.res)*HEAL01 + start(0.5)*HEAL01
# 
# F1~~1*F1
# F2~~1*F2
# F3~~1*F3
# F1~~F2
# F1~~F3
# F2~~F3
# 
# #label(HEAL01.res)^2<1
# 
# '
# 
# lavaanDefinitionManual='
# F1 =~ NA*ALCD03+ANXI03+DEPR05+NEUR01+SUBJ01
# F2 =~ NA*ALCD03+HEAL01+TIRE01
# 
# ALCD03~~r1*ALCD03
# ANXI03~~r2*ANXI03
# DEPR05~~r3*DEPR05
# HEAL01~~r4*HEAL01
# NEUR01~~r5*NEUR01
# SUBJ01~~r6*SUBJ01
# TIRE01~~r7*TIRE01
# 
# F1~~1*F1
# F2~~1*F2
# 
# F1~~F2
# 
# abs(r2)<0.1
# r2>0.01
# abs(r3)<0.1
# r3>0.01
# r7>0.01
# '
# 
# cModelResults<-usermodel.mod(covstruc = project$mvLD$covstruct.mvLDSC,
#         model=lavaanDefinitionManual,
#                                               #model = project$CFA$models[nModel,c("lModel")],
#         estimation = project$CFA$estimator,
#         fix_resid = FALSE
#         )
# cModelResults$modelfit

# library(tidySEM)
# sem_graph_layout <- get_layout("", "F1", "","","","F2","",
#                   "V1", "V2", "V3", "V4", "V5", "V6", "V7", 
#                   "VF1", "VF2", "VF3", "VF4", "VF5", "VF6", "VF7", rows = 3)
# sem_graph_layout <- get_layout("", "F1", "","","","F2","",
#                   "V1", "V2", "V3", "V4", "V5", "V6", "V7", rows = 2)
# graph_sem(model=cModelResults$lresults, layout = sem_graph_layout)
# graph_sem(model=cModelResults$lresults)

#test
# project$clOptions$task="cfa"
# project$clOptions$task_argument="0:499"

project$sumstats.sel$residualSizeLimitMax<-NA_real_
project$sumstats.sel$residualSizeLimitMax[which(project$sumstats.sel$code=="ANXI03" | project$sumstats.sel$code=="DEPR05")]<-0.10
#(1/project$sumstats.sel$h2.se.liability_mvLDSC^2)/sum(1/project$sumstats.sel$h2.se.liability_mvLDSC^2)+0.01

project$filepath.cfa<-file.path(project$folderpath.workingDirectory,paste0("cfa.",project$setup.code,".Rds"))

if (file.exists(project$filepath.cfa)) {
print("Using existing CFA results from previous run and appending to these if needed.")
project$CFA<-readRDS(file=project$filepath.cfa)
} else {
  project$CFA$nIndicators=length(project$sumstats.sel$code)
  project$CFA$indicatorLocks<-as.data.frame(matrix(data=0, nrow = project$CFA$nIndicators, ncol = project$CFA$nFactors))
  row.names(project$CFA$indicatorLocks)<-project$sumstats.sel$code
  #project$CFA$indicatorLocks[1,]<-c(1,0,1) #ALCD03
  project$CFA$indicatorLocks[2,]<-c(1,1,1) #ANXI03
  project$CFA$indicatorLocks[3,]<-c(1,1,1) #DEPR05
  #project$CFA$indicatorLocks[4,]<-c(1,0,0) #HEAL01
  project$CFA$indicatorLocks[5,]<-c(1,1,1) #NEUR01
  #project$CFA$indicatorLocks[6,]<-c(1,0,0) #SUBJ01
  #project$CFA$indicatorLocks[7,]<-c(1,1,1) #TIRE01
  
  project$CFA$models<-data.frame(nModel=c(),code=c(),lModel=c(),lResults=c())
  
  project$CFA$resultColumnNames<-c("chisq","df","p_chisq","AIC","CFI","SRMR")
  
  project$CFA$maximumLoadingPatternNumber<-2^((project$CFA$nFactors*project$CFA$nIndicators)-sum(project$CFA$indicatorLocks))
  
  #defaults
  project$CFA$sessionLoadingPatternStartNumber<-min(c(0,project$CFA$maximumLoadingPatternNumber))
  project$CFA$sessionLoadingPatternEndNumber<-min(c(99,project$CFA$maximumLoadingPatternNumber))
  
  #load intermediate results
  if(!project$clOptions$task=="cfa"){
    project$CFA$modelsIntermediate<-list.files(path = project$folderpath.workingDirectory, pattern = paste0("^cfa\\.",project$setup.code,"\\..+\\.Rds"), full.names = T, ignore.case=T)
    for(cModelsIntermediateFilePath in project$CFA$modelsIntermediate){
      cModelsIntermediate<-readRDS(file=cModelsIntermediateFilePath)
      project$CFA$models<-rbind(project$CFA$models,cModelsIntermediate$models)
    }
    nModel=nrow(project$CFA$models)
    if(nModel>0){
      saveRDS(object = project$CFA,file = file.path(project$folderpath.workingDirectory,paste0("cfa.",project$setup.code,".Rds")))
    }
  }
}

cat("The maximum CFA loading pattern bit value for this setup is ",project$CFA$maximumLoadingPatternNumber)

nModel=nrow(project$CFA$models)
if(nModel <1 || max(project$CFA$models$searchBitValue) < project$CFA$maximumLoadingPatternNumber) {
  
if(project$clOptions$task=="cfa" & !is.na(project$clOptions$task_argument)){
  taskargs<-as.integer(unlist(noquote(strsplit(project$clOptions$task_argument,":",fixed = TRUE))))
  project$CFA$sessionLoadingPatternStartNumber<-min(c(taskargs[1],project$CFA$maximumLoadingPatternNumber))
  project$CFA$sessionLoadingPatternEndNumber<-min(c(taskargs[2],project$CFA$maximumLoadingPatternNumber))
} else if(nModel>0){
    project$CFA$sessionLoadingPatternStartNumber<-min(c(max(project$CFA$models$searchBitValue)+1,project$CFA$maximumLoadingPatternNumber))
    project$CFA$sessionLoadingPatternEndNumber<-min(c(project$CFA$sessionLoadingPatternStartNumber+499,project$CFA$maximumLoadingPatternNumber))
  }
  
  project$CFA$sessionIndicatorLoadingPatterns<-semplate$generateIndicatorLoadingPatterns(nFactors = project$CFA$nFactors, nIndicators = project$CFA$nIndicators, indicatorLocksDf = project$CFA$indicatorLocks, searchBitValues = project$CFA$sessionLoadingPatternStartNumber:project$CFA$sessionLoadingPatternEndNumber)
  #View(project$CFA$sessionIndicatorLoadingPatterns$indicatorLoadings)
  #View(project$CFA$sessionIndicatorLoadingPatterns$indicatorLoadingsMatrix)
  
  
  sessionPatternLength<-nrow(project$CFA$sessionIndicatorLoadingPatterns$indicatorLoadings)
  cat("Analysing ",sessionPatternLength, "models...")
  
  for(nSessioPattern in 1:sessionPatternLength){
    #new model row
    
    #test
    #nSessioPattern=1
    nModel=nModel+1
    project$CFA$models[nModel,]<-NA
    
    #nModel
    project$CFA$models[nModel,c("nModel")]<-nModel
    #code
    project$CFA$models[nModel,c("totalBitValue")]<-project$CFA$sessionIndicatorLoadingPatterns$totalBitValues[nSessioPattern]
    #searchBitValue
    project$CFA$models[nModel,c("searchBitValue")]<-project$CFA$sessionIndicatorLoadingPatterns$searchBitValues[nSessioPattern]
    
    project$CFA$models[nModel,c("code")]<-paste0("M",project$CFA$nFactors,"-",project$CFA$nIndicators,
    ".",project$CFA$estimator,
    ".",project$CFA$sessionIndicatorLoadingPatterns$totalBitValues[nSessioPattern])
    
    cIndicatorLoadings<-as.data.frame(project$CFA$sessionIndicatorLoadingPatterns$indicatorLoadingsMatrix[[nSessioPattern]])
    row.names(cIndicatorLoadings)<-project$sumstats.sel$code
    
    
    #further filter rules
    indicatorsLoadedOnFactors <- apply(cIndicatorLoadings, 1, FUN = any)
    factorsHasIndicatorsLoaded <- apply(cIndicatorLoadings, 2, FUN = any)
  
    #init columns
    project$CFA$models[nModel,c("lModel")]<-NA_character_
    project$CFA$models[nModel,c("lResults")]<-NA
    project$CFA$models[nModel,project$CFA$resultColumnNames]<-NA
    
    if(all(indicatorsLoadedOnFactors) & all(factorsHasIndicatorsLoaded)){
      project$CFA$models[nModel,c("lModel")]<- semplate$generateLavaanCFAModel(
        allow_loading.table.indicator_factor = cIndicatorLoadings, indicatorArgs = project$sumstats.sel[,c("code","residualSizeLimitMax")],
        universalResidualLimitMin = 0.0001,
        orthogonal = FALSE
        )
      
      #evaluate model
      cat(project$CFA$models[nModel,c("code")], "\t:",project$CFA$models[nModel,c("searchBitValue")],"/",project$CFA$maximumLoadingPatternNumber)
      cModelResults=usermodel.mod(covstruc = project$mvLD$covstruct.mvLDSC,
          model = project$CFA$models[nModel,c("lModel")],
          estimation = project$CFA$estimator,
          fix_resid = FALSE
          )
      
      project$CFA$models[[nModel,c("lResults")]]=list(cModelResults)
      if(!is.null(cModelResults$modelfit)){
        project$CFA$models[nModel,project$CFA$resultColumnNames]<-cModelResults$modelfit[1,project$CFA$resultColumnNames]
      }
      
      
    }
  }
  
  
  if(project$clOptions$task=="cfa" & !is.na(project$clOptions$task_argument)){
    saveRDS(object = project$CFA,file = file.path(project$folderpath.workingDirectory,paste0("cfa.",project$setup.code,".",project$CFA$sessionLoadingPatternStartNumber,"-",project$CFA$sessionLoadingPatternEndNumber,".Rds")))
  } else {
    saveRDS(object = project$CFA,file = file.path(project$folderpath.workingDirectory,paste0("cfa.",project$setup.code,".Rds")))
  }
  
  print("CFA for this session is now done and the result should have been saved into a file.")
  
}



#View(project$CFA$models)
#project$CFA$models$lModel[which(project$CFA$models$nModel==3)]
if(project$clOptions$task=="cfa"){quit(save = "no")}



## ----CFA select--------------------------------------------------------------------------------------------------------

project$CFA$models.selected<-project$CFA$models[which(project$CFA$models$CFI>0 & project$CFA$models$AIC<100),c("totalBitValue","code","lModel",project$CFA$resultColumnNames)]
#View(project$CFA$models.selected)
project$CFA$model.bestFitting<-project$CFA$models[which(project$CFA$models$totalBitValue==1945431),] #this is a manual setting after choosing the best fitting model from the previous results 

project$CFA$model.bestFitting





## ----prepare summary statistics----------------------------------------------------------------------------------------
#eval=FALSE
#use ^this to knit without running the code in the chunk

#we need to introduce checks of the summary statistics and which scale they are on
# f.ex. Leo checks my s.e. 2 * pnorm(log(1.003) / 0.021644, mean = 0, lower.tail = FALSE)

project$lfGWAS<-c()

if (!file.exists(file.path(project$folderpath.workingDirectory,paste0("lfGWAS.sumstats.",project$setup.code,".Rds"))) | project$seting.refreshPrepareSummaryStatistics) 
{
  print("Preparing summary statistics for latent factor GWAS. This might take a while.")


  project$lfGWAS$sumstats<-sumstats.mod(
  filenames=project$sumstats.sel$cleanedpath,
  ref=project$filepath.SNPReference,
  trait.names=project$sumstats.sel$code,
  se.logit=project$sumstats.sel$se.logit,
  OLS=project$sumstats.sel$dependent_variable.OLS,
  linprob=NULL, #THIS SHOULD BE INVESTIGATED FURTHER, IF A LINEAR OLS ESTIMATOR ON A DICHOTOMOUS DEP. VARIABLE WAS USED FOR ANY OF THE DATSETS 
  prop=NULL,
  N=project$sumstats.sel$n_total,
  info.filter=NULL,
  maf.filter=NULL,
  keep.indel=FALSE,
  parallel=FALSE, #The default = T eats lots of memory at once.
  cores=NULL
  #num = 1 #test
  )

#Error in files[[i]]$effect[[1]] : subscript out of bounds

  saveRDS(object = project$lfGWAS$sumstats,file = file.path(project$folderpath.workingDirectory,paste0("lfGWAS.sumstats.",project$setup.code,".Rds")))
  print("Done summary statistics for latent factor GWAS. The result should have been saved to a file.")
} else {
  project$lfGWAS$sumstats<-readRDS(file=file.path(project$folderpath.workingDirectory,paste0("lfGWAS.sumstats.",project$setup.code,".Rds")))
}



## ----latent factor GWAS------------------------------------------------------------------------------------------------

#inactivated
if(!is.null(project$lfGWAS$sumstats)){
  
  #specified manually here
  project$CFA$lfGWAS$lmodel<-"
  F1 =~ NA*ALCD03+ANXI03+DEPR05+NEUR01+TIRE01
  F2 =~ NA*ANXI03+DEPR05+HEAL01+NEUR01+TIRE01
  F3 =~ NA*ANXI03+DEPR05+NEUR01+SUBJ01+TIRE01
  ALCD03~~r1*ALCD03
  ANXI03~~r2*ANXI03
  DEPR05~~r3*DEPR05
  HEAL01~~r4*HEAL01
  NEUR01~~r5*NEUR01
  SUBJ01~~r6*SUBJ01
  TIRE01~~r7*TIRE01
  F1~~1*F1
  F2~~1*F2
  F3~~1*F3
  F1~~F2
  F1~~F3
  F2~~F3
  
  F1~SNP
  F2~SNP
  F3~SNP
  
  
  r1>1e-04
  abs(r2)<0.1
  r2>1e-04
  abs(r3)<0.1
  r3>1e-04
  r4>1e-04
  r5>1e-04
  r6>1e-04
  r7>1e-04
  "
  
  if (!file.exists(file.path(project$folderpath.workingDirectory,paste0("lfGWAS.gwas.",project$setup.code,".Rds"))) | project$setting.refreshLatentFactorGWAS) 
  {
  
  print("Performing latent factor GWAS. This will take a while.")
  
  #project$lfGWAS$sumstats.test<-project$lfGWAS$sumstats[which(project$lfGWAS$sumstats$CHR=='22'),]
    
  project$lfGWAS$gwas<-userGWAS(covstruc = project$mvLD$covstruct.mvLDSC, SNPs = project$lfGWAS$sumstats, estimation = "ML", model = project$CFA$lfGWAS$lmodel, modelchi = FALSE, printwarn = TRUE, sub=c("F1~SNP", "F2~SNP","F3~SNP"), GC="none")
  
  saveRDS(object = project$lfGWAS$gwas,file = file.path(project$folderpath.workingDirectory,paste0("lfGWAS.gwas.",project$setup.code,".Rds")))
  
  print("DONE performing latent factor GWAS. The results should have been saved to a file.")
  
  } else {
    project$lfGWAS$gwas<-readRDS(file=file.path(project$folderpath.workingDirectory,paste0("lfGWAS.gwas.",project$setup.code,".Rds")))
  }

}



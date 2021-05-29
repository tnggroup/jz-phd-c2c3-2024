## ----package setup, echo=FALSE, warning=F-----------------------------------------------------

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
#install.packages("optparse")
#install.packages("stats")
#remove.packages("shru")
#devtools::install_github("johanzvrskovec/shru")
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



## ----command line setup-----------------------------------------------------------------------
clParser <- OptionParser()
clParser <- add_option(clParser, c("-t", "--task"), type="character", default="0",
                help="Index of the explicit task to run separately:\n0: No task\nmvLD.mvLDSC:multivariate LDSC\nmvLD.HDL.piecewise:HDL Piecewise\nmvLD.HDL.jackknife:HDL Jackknife\nmvLD.origHDL:original HDL(jackknife)\nmvLD.origHDL.liabilityScale:original HDL with applied liability scale [default %default]")
clParser <- add_option(clParser, c("-l", "--location"), type="character", default="local",
                help="The place where the code is run [local,cluster] [default %default]")

clParser <- add_option(clParser, c("-a", "--task_argument"), type="character", default=NA,
                help="General purpose argument for tasks [default %default]")



## ----settings---------------------------------------------------------------------------------
project<-c() #create project metadata object
project$clOptions<-parse_args(clParser)
project$date.run<-Sys.Date()
project$setup.version<-3
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
project$folderpath.workingDirectory<-normalizePath(file.path(project$folderpath,"working_directory"))

project$folderpath.scripts<-normalizePath(file.path(project$folderpath,"scripts"))
#project$folderpath.includedSoftware<-normalizePath(file.path(project$folderpath,"included_software"))
project$folderpath.plots<-normalizePath(file.path(project$folderpath,"plots"))


#general data folder
if(project$host=="local") {
  project$folderpath.data<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/data")
} else if (project$host=="cluster") {
  project$folderpath.data<-normalizePath(file.path(project$folderpath,"data"))
}

##cleaned sumstats folder
project$folderpath.data.sumstats.cleaned<-normalizePath(file.path(project$folderpath.data,"gwas_sumstats","cleaned"))

##munged sumstats folder
#project$folderpath.data.sumstats.munged<-normalizePath(file.path(project$folderpath.data,"gwas_sumstats","munged_1kg_eur_supermunge"))
project$folderpath.data.sumstats.munged<-normalizePath(file.path(project$folderpath.data,"gwas_sumstats","munged_1kg_eur_gSEM"))

##imputed sumstats folder
#project$folderpath.data.sumstats.imputed<-normalizePath(file.path(project$folderpath.data,"gwas_sumstat","imputed"))

#python virtual environment folder
if(project$host=="local") {
  project$folderpath.pythonVenv<-normalizePath("~/Documents/local_db/JZ_GED_PHD_C1/python-venv")
} else if (project$host=="cluster") {
  project$folderpath.pythonVenv<-normalizePath(file.path(project$folderpath,"python-venv"))
}

##Reference SNP-list (HapMap3 SNPs for example). Used for munging sumstat SNP data.
#project$filepath.SNPReference<-normalizePath(paste0(project$folderpath.data,"/","w_hm3.snplist.flaskapp2018")) #HapMap3 SNPs
## Used in the preparation step for performing latent factor GWAS as reference for calculating SNP variance across traits.
#project$filepath.SNPReference<-normalizePath(paste0(project$folderpath.data,"/","reference.1000G.maf.0.005.txt")) #1000 genomes phase 3
project$filepath.SNPReference<-normalizePath(file.path(project$folderpath.data,"combined.hm3_1kg.snplist.vanilla.jz2020.txt")) #custom hm3 + 1kg SNPs

project$filename.suffix.data.sumstats.munged<-".gz"
#project$filename.suffix.data.sumstats.munged<-"_noMHC.sumstats.gz"

##Reference panel folder containing individual level data reference panel. Used for GWAS sumstat imputation tasks.
#roject$folderpath.data.sumstatImp.genomeReference<-"/users/k1204688/brc_scratch/Public/1KG_Phase3/All"
project$folderpath.data.sumstatImp.genomeReference<-"/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/reference.panel.1KG_Phase3.CLEANED.EUR.cM"

##LD scores datasets folders (these strings need to have a trailing slash for the GSEM LDSC to work)
project$folderpath.data.mvLDSC.ld <- file.path(project$folderpath.data,"eur_w_ld_chr.1KG_Phase3")
#project$folderpath.data.mvLDSC.ld <- file.path(project$folderpath.data,"eur_w_ld_chr")#LD-scores
#Weights, if different from LD-scores
#Set weihghts to the same folder as ldscores
project$folderpath.data.mvLDSC.wld <- project$folderpath.data.mvLDSC.ld

##HDL LD scores reference - needs the trailing slashes!!!
if(project$host=="local") {
  #use the smallest LD reference as default for local tests
  project$folderpath.data.HDL.ld<-paste0(project$folderpath.data,"/UKB_array_SVD_eigen90_extraction/")
} else if (project$host=="cluster") {
  project$folderpath.data.HDL.ld<-paste0(project$folderpath.data,"/UKB_imputed_hm3_SVD_eigen99_extraction/")
}
  
##full script file paths
project$filepath.rmd<-normalizePath(file.path(project$folderpath.scripts,project$filename.rmd))
project$filepath.r<-normalizePath(file.path(project$folderpath.scripts,project$filename.r))

##CFA settings
project$CFA<-c()
project$CFA$estimator=c("ML")
#project$CFA$estimator=c("ML","DWLS")
project$CFA$nFactors=5

##latent factor GWAS filter settings
project$lfGWAS$info.filter=.59
project$lfGWAS$maf.filter=0.009

#working directory in case of running as an R-script
setwd(dir = normalizePath(project$folderpath.workingDirectory))

#inactivated python environment until it is used
#use_virtualenv(project$folderpath.pythonVenv)







## ----additional source setup, echo=FALSE, warning=F-------------------------------------------

#source(normalizePath(file.path(project$folderpath.scripts,"sumstats.mod-jz.R")))





## ----sumstat metadata database load-----------------------------------------------------------
project$filepath.sumstats<-file.path(project$folderpath.workingDirectory,paste0("sumstats.",project$setup.code,".Rds"))
if (file.exists(project$filepath.sumstats)) {
  print("Loading summary statistics metadata from previously stored file.")
  project$sumstats<-readRDS(file=project$filepath.sumstats)
} else {

#install.packages('RPostgres')
library(RPostgres)
library(DBI)

project$phenodbcon <- dbConnect(RPostgres::Postgres(),
                 dbname = 'phenodb', 
                 host = '10.200.105.5', 
                 port = 5432,
                 user = 'johan',
                 password = rstudioapi::askForPassword(prompt = "Enter database password for specified user."))

project$phenodbres <- dbSendQuery(project$phenodbcon, "SELECT \"GWAS\".*, category_id, category.name AS category_name, phenotype.name AS phenotype, phenotype.type AS phenotype_type, pmid, year FROM sumstat_old.\"GWAS\", sumstat_old.reference, sumstat_old.phenotype, sumstat_old.category 
WHERE \"GWAS\".reference_id=reference.id AND \"GWAS\".phenotype_id=phenotype.id AND phenotype.category_id = category.id
ORDER BY code,\"GWAS\".id")
project$sumstats<-dbFetch(project$phenodbres)
dbClearResult(project$phenodbres)

#fallback
#project$sumstats<-read.table(file.path(project$folderpath.data,"ukbb_sumstats_download202005.csv"), header=T, quote="\"", sep = ",", fill=T, blank.lines.skip=T,as.is = c(2), strip.white = T)

saveRDS(project$sumstats,file = project$filepath.sumstats)
write.table(project$sumstats, file = file.path(project$folderpath.workingDirectory,paste0(project$setup.code,".sumstats.tsv")), quote = F, sep = "\t", row.names = FALSE, col.names = TRUE)
}





## ----trait setup------------------------------------------------------------------------------
#,echo=FALSE
project$trait<-data.frame(phenotype_id=c())

#ADHD
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-139
project$trait[nrow(project$trait),c("populationPrevalence")]<-.00529 #Worldwide-pooled
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1176/ajp.2007.164.6.942"

#ALCD
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-141
project$trait[nrow(project$trait),c("populationPrevalence")]<-.125 #US total measurement
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1001/archpsyc.64.7.830"

#ANXI
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-144
project$trait[nrow(project$trait),c("populationPrevalence")]<-.16 #Any type of anxiety disorder, Via https://doi.org/10.1038/s41380-019-0559-1 (2019), originally from https://doi.org/10.1017/S1121189X00001421 (2009)
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1017/S1121189X00001421"

#AUTI
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-145
project$trait[nrow(project$trait),c("populationPrevalence")]<-.0122
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1186/s12874-016-0280-6"

#BIPO
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-146
project$trait[nrow(project$trait),c("populationPrevalence")]<- .007 #mean of male and female global prevalence rate (2013) from https://doi.org/10.1111/bdi.12423 (2016)
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1111/bdi.12423"

#DEPR
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-149
project$trait[nrow(project$trait),c("populationPrevalence")]<-.146 #MDD, .15 from the LD-calculations in https://doi.org/10.1038/s41588-018-0090-3 (2018), but with a possible reference to https://doi.org/10.1146/annurev-publhealth-031912-114409 (2013) which states 14.6% lifetime prevalence of MDE in high-income countries.
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1146/annurev-publhealth-031912-114409"

#INSO
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-171
project$trait[nrow(project$trait),c("populationPrevalence")]<-.69 #using the 'prevalence', the higher estimate of persistent symptoms after 1y follow up.
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1093/sleep/30.3.274"

#SCHI
project$trait[nrow(project$trait)+1,c("phenotype_id")]<-151
project$trait[nrow(project$trait),c("populationPrevalence")]<-.0072 #using the lifetime morbid risk estimate, slightly more suitable for the CLOZUK sample as described in the GWAS but not using the lower point estimate of 0.4%
project$trait[nrow(project$trait),c("referenceDOI")]<-"https://doi.org/10.1093/epirev/mxn001"

project$trait



## ----GWAS sumstat dataset setup---------------------------------------------------------------
#, echo=FALSE


#View(project$sumstats)

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
project$sumstats[nrow(project$sumstats)+1,c("code","n_case","n_control","n_total","phenotype_id","reference_id","reference_doi")]=list(
  code=c("DEPR05"),
  n_case=16823,
  n_control=25632,
  n_total=42455,
  phenotype_id=149,
  reference_id=127,
  reference_doi=c("https://doi.org/10.1038/s41588-018-0090-3")
  )

#reformat columns
project$sumstats$name<-as.character(project$sumstats$name)

##Add comprehensive names as in the Google sheet
project$sumstats$name[which(project$sumstats$code=="DEPR05")]="Major depressive disorder (PGC2 29) - only clinical ascertainment"



##Add trait/disorder information
project$sumstats <- project$sumstats %>%
mutate(
  code.trait=substr(x = code, start = 1, stop = 4)
       ) %>%
  left_join(project$trait[,c("phenotype_id","populationPrevalence")], by = c("phenotype_id" = "phenotype_id"))

##Add sumstat munged file paths
project$sumstats <- project$sumstats %>%
mutate(
  mungedpath=file.path(project$folderpath.data.sumstats.munged,paste0(code,project$filename.suffix.data.sumstats.munged))
       )

project$sumstats <- project$sumstats %>%
mutate(
  cleanedpath=file.path(project$folderpath.data.sumstats.cleaned,paste0(code,".gz"))
       )

##add reference year
project$sumstats$year[which(project$sumstats$code=="ANXI03")]=2019
project$sumstats$year[which(project$sumstats$code=="DEPR05")]=2018

##Add doi links for easy access to dataset publication
project$sumstats$reference_doi[which(project$sumstats$code=="ALCD03")]="https://doi.org/10.1038/s41593-018-0275-1"
project$sumstats$reference_doi[which(project$sumstats$code=="ANXI03")]="https://doi.org/10.1038/s41380-019-0559-1"
project$sumstats$reference_doi[which(project$sumstats$code=="AUTI07")]="https://doi.org/10.1038/s41588-019-0344-8"
project$sumstats$reference_doi[which(project$sumstats$code=="DEPR08")]="https://doi.org/10.1038/s41588-018-0090-3"
project$sumstats$reference_doi[which(project$sumstats$code=="HEAL01")]="https://doi.org/10.1093/ije/dyw219"
project$sumstats$reference_doi[which(project$sumstats$code=="NEUR01")]="https://doi.org/10.1038/ng.3552"
project$sumstats$reference_doi[which(project$sumstats$code=="SUBJ01")]="https://doi.org/10.1038/ng.3552"
project$sumstats$reference_doi[which(project$sumstats$code=="TIRE01")]="https://doi.org/10.1038/mp.2017.5"


##Add PMID
project$sumstats$pmid[which(project$sumstats$code=="ANXI03")]="31748690"
# project$sumstats$pmid[which(project$sumstats$code=="NEUR01")]="27089181"
# project$sumstats$pmid[which(project$sumstats$code=="SUBJ01")]="27089181"
project$sumstats$pmid[which(project$sumstats$code=="DEPR05")]="29700475"
# project$sumstats$pmid[which(project$sumstats$code=="DEPR08")]="29700475"

##add dependent variable type
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

##add sample prevalence for all datasets
project$sumstats$samplePrevalence<-project$sumstats$n_case/project$sumstats$n_total

#set order of datasets to sorted by code
project$sumstats<-project$sumstats[order(project$sumstats$code),]

#save the project data
#saveRDS(project,file = file.path(project$folderpath.workingDirectory,paste0("project.",project$setup.code,".Rds")))

#View(project$sumstats)



## ----GWAS sumstat dataset variable selection--------------------------------------------------

#selection based on specific traits
project$sumstats.sel.code<-c("ADHD05","ALCD03","ANXI03","AUTI07","BIPO02", "DEPR05","DEPR08","HEAL01","INCO03","INSO02", "MIGR01","NEUR01","RISK02","RISK03","SCHI04","SUBJ01","TIRE01")
project$sumstats.sel<-project$sumstats[which(project$sumstats$code %in% project$sumstats.sel.code),]

#project$sumstats.sel$code_orig<-project$sumstats.sel$code
#project$sumstats.sel$code<-project$sumstats.sel$code.trait
project$sumstats.sel[,c("code","name","year","n_case","n_control","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")]
project$k.sel<-nrow(project$sumstats.sel)
#View(project$sumstats.sel[,c("code","n_total","pmid","reference_doi","samplePrevalence","populationPrevalence","mungedpath")])

##Add nice trait names to be used in the report
project$sumstats.sel$name.nice<-project$sumstats.sel$name
#project$sumstats.sel$name.nice[which(project$sumstats$code=="ADHD05")]="ADHD"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="ALCD03")]="Alcohol dependence"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="ANXI03")]="Anxiety disorder"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="AUTI07")]="Autism spectrum disorder"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="DEPR05")]="Major depressive disorder narrow"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="DEPR08")]="Major depressive disorder"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="INCO03")]="Social deprivation"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="INSO02")]="Insomnia"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="NEUR01")]="Neuroticism"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="RISK02")]="Risktaking, automobile"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="RISK03")]="Risktaking, sex"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="SCHI04")]="Schizophrenia"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="SUBJ01")]="Subjective well-being"
project$sumstats.sel$name.nice[which(project$sumstats.sel$code=="TIRE01")]="Self-reported tiredness"

##Add a combined nice name plus code label
project$sumstats.sel$name.nice.and_code<-paste0(project$sumstats.sel$name.nice," (",project$sumstats.sel$code,")")


##add data on whether effects are on log odds ratio or not
## Plain OR is detected by munge and sumstat procedures though

project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="ADHD05")]=T #Column OR. Munge non-OR detected. Consistent with non-OR.
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="ALCD03")]=T #Did not find the original readme for this, only looked at the readme for the 2018 subsequent release, and the data. Assuming log-scale from the naming of the variables (beta?) and that it follows the standard of the 2018 release - explicitly states "log scale beta".
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="ANXI03")]=F # Kirstin sort of confirm that she converted both effect and s.e. to linear scales. Data and readme is available from her KCL web site. Unclear which dataset was used for ANXI03 though. Confirm this later by comparing datasets?
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="AUTI07")]=T #Assuming log(OR) since not plain OR
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="BIPO02")]=T #Assuming log(OR) since not plain OR
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="DEPR05")]=F #https://www.med.unc.edu/pgc/files/2019/04/PGC_MDD_2018_README_third_180316.pdf # The readme states daner format, but that SE is s.e. of log(OR), while the daner format does not state log scale variables here. Assuming linear scale OR, which is consistent with the data, and munge.
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="DEPR08")]=F #https://www.med.unc.edu/pgc/files/2019/04/PGC_MDD_2018_README_third_180316.pdf # As DEPR05
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="HEAL01")]=F #The original data including readme downloaded to data_raw/original. States "Beta" without explanation, and I assume linear beta.
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="INCO03")]=F #As per munge
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="INSO02")]=F #As per munge
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="MIGR01")]=F #As per munge
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="NEUR01")]=F #The official read-me states "regression beta"
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="RISK02")]=F #As per munge
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="RISK03")]=F #As per munge
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="SCHI04")]=F #OR, as per munge
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="SUBJ01")]=F #The official read-me states "regression beta"
project$sumstats.sel$effect.logit[which(project$sumstats.sel$code=="TIRE01")]=F #The official read-me states "beta"


##add data on whether the SEs are on a logistic scale or not
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="ADHD05")]=T #consistent with log(OR) above

project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="ALCD03")]=T #Did not find the original readme for this, only looked at the readme for the 2018 subsequent release, and the data. SE is strangely set to 1 for top rows. Assumes this is true because of log-scale beta.
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="ANXI03")]=F #Assume linear becasue of Kirstin Purves assuring that both were converted to linear scale. Consistent with munge non-OR.
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="AUTI07")]=T #Assuming log(OR) since not plain OR
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="BIPO02")]=T #Assuming log(OR) since not plain OR
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="DEPR05")]=T #log(OR) standard errors as of documentation
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="DEPR08")]=T #Assuming log(OR) standard errors as DEPR05
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="HEAL01")]=F #Assume linear scale because of linear beta
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="INCO03")]=F #Assume linear scale because of linear beta
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="INSO02")]=T #Assuming log(OR) standard errors as some of the datasets above.
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="MIGR01")]=F #Assume linear scale because of linear beta
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="NEUR01")]=F #Assume linear scale because of linear beta
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="RISK02")]=F #Assume linear scale because of linear beta
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="RISK03")]=F #Assume linear scale because of linear beta
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="SCHI04")]=T #Assuming log(OR) standard errors as some of the datasets above.
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="SUBJ01")]=F #Assume linear scale because of linear beta
project$sumstats.sel$se.logit[which(project$sumstats.sel$code=="TIRE01")]=F #Assume linear scale because of linear beta


##add information on wether a continous dependent variable was analysed using an OLS (linear) estimator. Used for the latent factor GWAS preparation step.
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="ADHD05")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="ALCD03")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="ANXI03")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="AUTI07")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="BIPO02")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="DEPR05")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="DEPR08")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="HEAL01")]=T
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="INCO03")]=T
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="INSO02")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="MIGR01")]=T
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="NEUR01")]=T
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="RISK02")]=T
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="RISK03")]=T
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="SCHI04")]=F
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="SUBJ01")]=T
project$sumstats.sel$dependent_variable.OLS[which(project$sumstats.sel$code=="TIRE01")]=T


write.table(project$sumstats.sel[,c("code", "name","year", "n_case","n_control","n_total","samplePrevalence","populationPrevalence", "reference_doi")], file = file.path(project$folderpath.workingDirectory,paste0(project$setup.code,".sumstatinfo.tsv")), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)

#View(project$sumstats.sel)



## ----GWAS sumstat munge-----------------------------------------------------------------------
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
  
  

  #write.table(x = project$sumstats.reference.new,file = file.path(project$folderpath.workingDirectory,"combined.hm3_1kg.snplist.jz2020.txt"), quote = FALSE, row.names = F)
  #faster alt - produces a comma separated version?
  #data.table::fwrite(x = project$sumstats.reference.new, file = )
  
  # project$sumstats.reference.new.vanilla<-project$sumstats.reference.new[!is.na(project$sumstats.reference.new$A1),]
  # write.table(x = project$sumstats.reference.new.vanilla,file = file.path(project$folderpath.workingDirectory,"combined.hm3_1kg.snplist.vanilla.jz2020.txt"), quote = FALSE, row.names = F)
  
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
  
  project$munge<-c()
  #test
  #project$clOptions$task_argument<-"ADHD05"
  if(!is.na(project$clOptions$task_argument)){
    #task argument sets the specific dataset to munge
    
    cat("\nSet to munge with arg",project$clOptions$task_argument,"\n")
    project$munge$filesToUse<-project$sumstats.sel$cleanedpath[which(project$sumstats.sel$code==project$clOptions$task_argument)]
    project$munge$traitNamesToUse<-project$sumstats.sel$code[which(project$sumstats.sel$code==project$clOptions$task_argument)]
    project$munge$NToUse<-project$sumstats.sel$n_total[which(project$sumstats.sel$code==project$clOptions$task_argument)]
    
  } else {
    #defaults
    project$munge$filesToUse<-project$sumstats.sel$cleanedpath
    project$munge$traitNamesToUse<-project$sumstats.sel$code
    project$munge$NToUse<-project$sumstats.sel$n_total
  }
  
  #munging with no filters applied
  project$sumstats.meta <- supermunge(filePaths = project$munge$filesToUse,
            refFilePath = project$filepath.SNPReference,
            traitNames = project$munge$traitNamesToUse,
            N = project$munge$NToUse,
            pathDirOutput = project$folderpath.data.sumstats.munged
              ) 
    
  project$sumstats.meta <- project$sumstats.meta %>% gt()
  
  gtsave(data = project$sumstats.meta, filename = file.path(project$folderpath.workingDirectory,"sumstats.meta.rtf"))
    
    quit(save = "no")
  }



## ----GWAS sumstat imputation------------------------------------------------------------------
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



## ----multivariate LD--------------------------------------------------------------------------

project$filepath.mvLD<-file.path(project$folderpath.workingDirectory,paste0("mvLD.",project$setup.code,".Rds"))
project$filepath.mvLD.mvLDSC<-file.path(project$folderpath.workingDirectory,paste0("mvLD.",project$setup.code,".mvLDSC.Rds"))

if (file.exists(project$filepath.mvLD)) {
  print("Using existing covariance structures from previous LD computations.")
  project$mvLD<-readRDS(file=project$filepath.mvLD)
} else {
  print("Running (or reading ready intermediate results from) multivariate LD regression with different methods. This might take a while. If the procedure runs for too long you may want to abort the process.")
  
  cat("The current task is specified as:",project$clOptions$task)
  project$mvLD<-c()
  
  
  
  if(project$clOptions$task=="mvLD.mvLDSC" || !file.exists(project$filepath.mvLD.mvLDSC)){
    #run mvLDSC
    project$mvLD$covstruct.mvLDSC<-ldsc.mod(traits = project$sumstats.sel$mungedpath,
                                    sample.prev =  project$sumstats.sel$samplePrevalence,
                                    population.prev = project$sumstats.sel$populationPrevalence,
                                    trait.names = project$sumstats.sel$code,
                                    ld = project$folderpath.data.mvLDSC.ld,
                                    wld = project$folderpath.data.mvLDSC.ld,
                                    n.blocks = 600, 
                                    info.filter = 0.6,
                                    maf.filter = 0.01,
                                    N = project$sumstats.sel$n_total,
                                    ldsc.log = project$setup.code.date,
                                    stand = TRUE
                                    )
    saveRDS(object = project$mvLD$covstruct.mvLDSC,file = project$filepath.mvLD.mvLDSC)
    if(project$clOptions$task=="mvLD.mvLDSC"){
      quit(save = "no")
    }
  }
    
  if(!file.exists(project$filepath.mvLD.mvLDSC)) {
    Sys.sleep(time = 3)
  }
  
  cat("Reading ready intermediate mvLD results and saving final mvLD result.")
  project$mvLD$covstruct.mvLDSC<-readRDS(file=project$filepath.mvLD.mvLDSC)
  
  #save the mvLD output
  saveRDS(object = project$mvLD,file = project$filepath.mvLD)
  print("Multivariate LD correction is done now and the resulting covariance structure should have been saved to a file.")
  
  #remove intermediate results
  file.remove(project$filepath.mvLD.mvLDSC)
}

#Additional computations
#saving the original S in case of smoothing experiments later stored in S
project$mvLD$covstruct.mvLDSC$S.orig<-project$mvLD$covstruct.mvLDSC$S
project$mvLD$covstruct.mvLDSC$S.smooth<-as.matrix((nearPD(project$mvLD$covstruct.mvLDSC$S, corr = FALSE))$mat)

#adding standard errors from the V matrices, unstandardised and standardised. I added these calculations to the modified ldsc.

#retrieve the standard errors of S (variances and covariances) from the diagonal of V (contains both).

# project$mvLD$covstruct.mvLDSC$S.SE<-matrix(0, project$k.sel, project$k.sel)
# rownames(project$mvLD$covstruct.mvLDSC$S.SE)<-colnames(project$mvLD$covstruct.mvLDSC$S)
# colnames(project$mvLD$covstruct.mvLDSC$S.SE)<-colnames(project$mvLD$covstruct.mvLDSC$S)
# project$mvLD$covstruct.mvLDSC$S_Stand.SE<-matrix(0, project$k.sel, project$k.sel)
# rownames(project$mvLD$covstruct.mvLDSC$S_Stand.SE)<-colnames(project$mvLD$covstruct.mvLDSC$S)
# colnames(project$mvLD$covstruct.mvLDSC$S_Stand.SE)<-colnames(project$mvLD$covstruct.mvLDSC$S)
# 
# project$mvLD$covstruct.mvLDSC$S.SE[lower.tri(project$mvLD$covstruct.mvLDSC$S.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.mvLDSC$V))
# project$mvLD$covstruct.mvLDSC$S.SE[upper.tri(project$mvLD$covstruct.mvLDSC$S.SE)]<-t(project$mvLD$covstruct.mvLDSC$S.SE)[upper.tri(project$mvLD$covstruct.mvLDSC$S.SE)]
# project$mvLD$covstruct.mvLDSC$S_Stand.SE[lower.tri(project$mvLD$covstruct.mvLDSC$S_Stand.SE,diag=TRUE)] <-sqrt(diag(project$mvLD$covstruct.mvLDSC$V_Stand))
# project$mvLD$covstruct.mvLDSC$S_Stand.SE[upper.tri(project$mvLD$covstruct.mvLDSC$S_Stand.SE)]<-t(project$mvLD$covstruct.mvLDSC$S_Stand.SE)[upper.tri(project$mvLD$covstruct.mvLDSC$S_Stand.SE)]


  
#add newly computed heritabilities to the selected summary statistics table
project$sumstats.sel$h2.liability_mvLDSC<-diag(project$mvLD$covstruct.mvLDSC$S)
project$sumstats.sel$h2.se.liability_mvLDSC<-diag(project$mvLD$covstruct.mvLDSC$S.SE)


#View(project$sumstats.sel)

  








## ----clustering-------------------------------------------------------------------------------

#saving clustering results between runs as to always use the same randomised start
project$filepath.clustering<-file.path(project$folderpath.workingDirectory,paste0("clustering.",project$setup.code,".Rds"))

if (file.exists(project$filepath.clustering)) {
print("Using existing clustering results from previous run and appending to these if needed.")
project$clustering<-readRDS(file=project$filepath.clustering)
} else {
  project$clustering<- kmeans(x = abs(project$mvLD$covstruct.mvLDSC$S_Stand), centers = project$CFA$nFactors, iter.max = 1000, nstart = 30)
  #View(project$clustering$centers)
  #View(fitted(project$clustering))
  #View(abs(project$mvLD$covstruct.mvLDSC$S_Stand))
  #resid<-abs(project$mvLD$covstruct.mvLDSC$S_Stand)-fitted(project$clustering)
  #View(resid)
  project$clustering$centerDistance<-apply(X = abs(project$mvLD$covstruct.mvLDSC$S_Stand), MARGIN = 1, FUN = function(x){
    #test
    #x<-abs(project$mvLD$covstruct.mvLDSC$S_Stand)[1,]
    #cat("\nOBS:",x)
    ss<-apply(X=(abs(x)-project$clustering$centers)^2, FUN = sum, MARGIN = 1)
    #cat("\nSS:",ss)
    s<-sum(abs(x))
    #cat("\nS:",s)
    return (ss/s)
  })
  #transpose this to conform with indicator loading data frames
  project$clustering$centerDistance<-t(project$clustering$centerDistance)
  rownames(project$clustering$centerDistance)<-project$sumstats.sel.code
  saveRDS(object = project$clustering,file = file.path(project$folderpath.workingDirectory,paste0("clustering.",project$setup.code,".Rds")))
}

cat("\nClustering centers\n")
t(project$clustering$centers)
cat("\nResiduals\n")
project$clustering$centerDistance




## ----CFA--------------------------------------------------------------------------------------

# testing tidySEM
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

project$CFA$nIndicators=length(project$sumstats.sel$code)

cat("\nInitial PCA\n")
project$CFA$PCA<-eigen(x = abs(project$mvLD$covstruct.mvLDSC$S_Stand), symmetric = TRUE)
rownames(project$CFA$PCA$vectors)<-project$sumstats.sel.code
project$CFA$PCA_absolute_values<-abs(project$CFA$PCA$vectors*project$CFA$PCA$values)[,1:project$CFA$nFactors]
project$CFA$PCA_absolute_values


cat("\nInitial clustering\n")
project$CFA$clustering_absolute_values <- abs(t(project$clustering$centers))
project$CFA$clustering_absolute_values

project$CFA$loading_cutoff_lock_exclude_max<-0.25
project$CFA$loading_cutoff_lock_include_min<-0.30
cat("\nSetting loading locks based on initial absolute clustering loadings with cutoff set to lock-to-exclude:",project$CFA$loading_cutoff_lock_exclude_max,"and lock to include:",project$CFA$loading_cutoff_lock_include_min,"\n")

project$CFA$indicatorLocksExclude <- project$CFA$clustering_absolute_values < project$CFA$loading_cutoff_lock_exclude_max
row.names(project$CFA$indicatorLocksExclude)<-project$sumstats.sel$code
project$CFA$indicatorLocksInclude <- project$CFA$clustering_absolute_values > project$CFA$loading_cutoff_lock_include_min
row.names(project$CFA$indicatorLocksInclude)<-project$sumstats.sel$code
project$CFA$indicatorLocks<-project$CFA$indicatorLocksExclude | project$CFA$indicatorLocksInclude

cat("\nCFA Indicator Exclude locks set as:\n")
print(project$CFA$indicatorLocksExclude)
cat("\nCFA Indicator Include locks set as:\n")
print(project$CFA$indicatorLocksInclude)
cat("\nA total of",sum(project$CFA$indicatorLocksInclude)," locked indicator loadings and",((project$CFA$nIndicators*project$CFA$nFactors)-sum(project$CFA$indicatorLocks)),"variable indicator loadings to estimate.")

project$sumstats.sel$residualSizeLimitMax<-NA_real_
#project$sumstats.sel$residualSizeLimitMax[which(project$sumstats.sel$code=="ANXI03" | project$sumstats.sel$code=="DEPR05")]<-0.10
#(1/project$sumstats.sel$h2.se.liability_mvLDSC^2)/sum(1/project$sumstats.sel$h2.se.liability_mvLDSC^2)+0.01


# for(fi in 1:project$CFA$nFactors){
#   project$CFA$indicatorLocks[,fi]<-project$CFA$indicatorLocks[,fi]*project$CFA$eigenvalues[fi]
# }

#project$CFA$indicatorLocks<-as.data.frame(matrix(data=0, nrow = project$CFA$nIndicators, ncol = project$CFA$nFactors)) #old manual


#force fixed loadings
#project$CFA$indicatorLocks["ANXI03",]<-T
#project$CFA$indicatorLocks["DEPR05",]<-T
#project$CFA$indicatorLocks["DEPR08",]<-T

project$CFA$models<-data.frame(nModel=c(),code=c(),lModel=c(),lResults=c())

project$CFA$resultColumnNames<-c("chisq","df","p_chisq","AIC","CFI","SRMR")

project$CFA$maximumTotalLoadingPatternNumber<-2^(project$CFA$nFactors*project$CFA$nIndicators)
project$CFA$maximumLoadingPatternNumber<-2^((project$CFA$nFactors*project$CFA$nIndicators)-sum(project$CFA$indicatorLocks))

#defaults
project$CFA$sessionLoadingPatternStartNumber<-min(c(0,project$CFA$maximumLoadingPatternNumber))
project$CFA$sessionLoadingPatternEndNumber<-min(c(99,project$CFA$maximumLoadingPatternNumber))

cat("\nThe number of models to search among for this setup is ",project$CFA$maximumLoadingPatternNumber, "\nof a total of",project$CFA$maximumTotalLoadingPatternNumber," possible model configurations.")

project$filepath.cfa<-file.path(project$folderpath.workingDirectory,paste0("cfa.",project$setup.code,".Rds"))
#project$filepath.cfa_converged_results<-file.path(project$folderpath.workingDirectory,paste0("cfa.",project$setup.code,".converged.txt"))

if (file.exists(project$filepath.cfa)) {
print("Using existing CFA results from previous run and appending to these if needed.")
project$CFA<-readRDS(file=project$filepath.cfa)
} else {
  
  
  
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


nModel=nrow(project$CFA$models)
if(nModel <1 || max(project$CFA$models$searchBitValue) < project$CFA$maximumLoadingPatternNumber) {
  
if(project$clOptions$task=="cfa" & !is.na(project$clOptions$task_argument)){
  taskargs<-as.integer(unlist(noquote(strsplit(project$clOptions$task_argument,":",fixed = TRUE))))
  project$CFA$sessionLoadingPatternStartNumber<-min(c(taskargs[1],project$CFA$maximumLoadingPatternNumber))
  project$CFA$sessionLoadingPatternEndNumber<-min(c(taskargs[2],project$CFA$maximumLoadingPatternNumber))
} else if(nModel>0){
  
    project$CFA$sessionLoadingPatternStartNumber<-min(c(ifelse(is.null(project$CFA$models$searchBitValue),0,max(project$CFA$models$searchBitValue)+1),project$CFA$maximumLoadingPatternNumber))
    project$CFA$sessionLoadingPatternEndNumber<-min(c(project$CFA$sessionLoadingPatternStartNumber+49,project$CFA$maximumLoadingPatternNumber))
  }
  
  project$CFA$sessionIndicatorLoadingPatterns<-semplate$generateIndicatorLoadingPatterns(searchBitValues = (project$CFA$sessionLoadingPatternStartNumber:project$CFA$sessionLoadingPatternEndNumber), indicatorLocks.include = project$CFA$indicatorLocksInclude, indicatorLocks.exclude = project$CFA$indicatorLocksExclude)
  #View(project$CFA$sessionIndicatorLoadingPatterns$searchBitValues)
  #View(project$CFA$sessionIndicatorLoadingPatterns$indicatorLoadings[[4]][[1]]!=project$CFA$indicatorLocks)
  
  
  sessionPatternLength<-length(project$CFA$sessionIndicatorLoadingPatterns$indicatorLoadings)
  cat("\nAnalysing ",sessionPatternLength, "models...\n")
  nFittingModelsFound<-0
  for(nSessioPattern in 1:sessionPatternLength){
    #new model row
    
    #test
    #nSessioPattern=1
    nModel=nModel+1
    project$CFA$models[nModel,]<-NA
    project$CFA$models$totalBitValue<-NA #because otherwise the later assignment will crash
    
    #nModel
    project$CFA$models[nModel,c("nModel")]<-nModel
    #totalBitValue
    project$CFA$models[[nModel,c("totalBitValue")]]<-project$CFA$sessionIndicatorLoadingPatterns$totalBitValues[nSessioPattern]
    #searchBitValue
    project$CFA$models[nModel,c("searchBitValue")]<-project$CFA$sessionIndicatorLoadingPatterns$searchBitValues[nSessioPattern]
    #code
    project$CFA$models[nModel,c("code")]<-paste0("M",project$CFA$nFactors,"-",project$CFA$nIndicators,
    ".",project$CFA$estimator,
    "._",paste0(project$CFA$sessionIndicatorLoadingPatterns$totalBitValues[[nSessioPattern]],collapse="",sep="_"))
    
    cIndicatorLoadings<-project$CFA$sessionIndicatorLoadingPatterns$indicatorLoadings[[nSessioPattern]][[1]]
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
        universalResidualLimitMin = 0.0001, #The name should be changed to max!!! Error.
        orthogonal = FALSE
        )
      
      #evaluate model
      cat("\n#Found fitting=",nFittingModelsFound,"\n", project$CFA$models[nModel,c("searchBitValue")],":",project$CFA$models[nModel,c("code")], ":",nModel,"/",sessionPatternLength)
      cModelResults=usermodel.mod(covstruc = project$mvLD$covstruct.mvLDSC,
          model = project$CFA$models[nModel,c("lModel")],
          estimation = project$CFA$estimator,
          fix_resid = F,
          CFIcalc = F #set this to true for CFI evaluation
          )
      
      if(!is.null(cModelResults$modelfit)){
        if(nrow(cModelResults$modelfit)>0 && any(project$CFA$resultColumnNames %in% colnames(cModelResults$modelfit))) {
          #record results even though not fitting
          project$CFA$models[[nModel,c("lResults")]]<-list(cModelResults$lresults)
          cRescolnames<-intersect(project$CFA$resultColumnNames,colnames(cModelResults$modelfit))
          project$CFA$models[nModel,cRescolnames]<-cModelResults$modelfit[1,cRescolnames]
          if(is.numeric(cModelResults$modelfit$chisq)){
            #This is considered a fitting model
            nFittingModelsFound<-nFittingModelsFound+1
            cat("\nFITTING!:",project$CFA$models[nModel,c("searchBitValue")],"\n")
            print(cModelResults$modelfit)
          }
        }
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
#project$CFA$models$lModel[which(project$CFA$models$nModel==57 & project$CFA$models$code=="M3-17.ML._-1107563521_389110_0_0_0_0_0_0_")]
if(project$clOptions$task=="cfa"){quit(save = "no")}



## ----CFA investigate--------------------------------------------------------------------------

project$CFA$models.selected<-project$CFA$models.selected<-project$CFA$models[which(project$CFA$models$SRMR<1),c("totalBitValue","code","lModel",project$CFA$resultColumnNames)]
project$CFA$models.selected<-project$CFA$models[which(project$CFA$models$AIC<100),c("totalBitValue","code","lModel",project$CFA$resultColumnNames)]
#View(project$CFA$models.selected)




## ----CFA select-------------------------------------------------------------------------------

View(project$CFA$models.selected<-project$CFA$models[which(project$CFA$models$SRMR<1),c("totalBitValue","code","lModel",project$CFA$resultColumnNames)])
project$CFA$models.selected<-project$CFA$models[which(project$CFA$models$AIC<100),c("totalBitValue","code","lModel",project$CFA$resultColumnNames)]
#View(project$CFA$models.selected)
project$CFA$model.bestFitting<-project$CFA$models[which(project$CFA$models$totalBitValue==1945431),] #this is a manual setting after choosing the best fitting model from the previous results 

project$CFA$model.bestFitting





## ----prepare summary statistics---------------------------------------------------------------
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



## ----latent factor GWAS-----------------------------------------------------------------------

#inactivated
if(!is.null(project$lfGWAS$sumstats)){
  
  #specified manually here
  project$lfGWAS$lmodel<-"
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



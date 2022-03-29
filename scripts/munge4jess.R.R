#munge for Jess

library(shru)
library(readr)

folderpath<-normalizePath("/scratch/users/k19049801/project/JZ_GED_PHD_C1")
folderpath.data<-normalizePath(file.path(folderpath,"data"))
filepath.SNPReference.1kg<-normalizePath(file.path(folderpath.data,"combined.hm3_1kg.snplist.vanilla.jz2020.txt"))
folderpath.data.mvLDSC.ld.1kg<-file.path(folderpath.data,"eur_w_ld_chr.1KG_Phase3")
folderpath.data.sumstats.munged<-normalizePath(file.path(folderpath.data,"gwas_sumstats","munged_1kg_eur_supermunge"))

res <- supermunge(
  filePaths = c("/scratch/groups/ukbiobank/usr/jess/genetics_MDQ/results/GLADv2/MDQ_co_items/GWAS/regenie/step2/regenie_GLAD_only_co_MDQ_items_allchr_GWASv1_step2_mdq.co_sum_score.12.txt"),
  refFilePath = filepath.SNPReference.1kg,
  ldDirPath = folderpath.data.mvLDSC.ld.1kg,
  traitNames = c("MDQ"),
  pathDirOutput = folderpath.data.sumstats.munged
)

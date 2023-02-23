#Michel Nivard's GWIS of height and weight => BMI
#EDITED!!!!
#https://sites.google.com/site/mgnivard/gwis/code-examples


### Require msm, a packadge whihc contains code for delta rule aproximation, we need this to get GWIS standard errors ####

require(msm)


# BMI aproximation from height and weight GWAS # 

# From the NTR (see: Willemsen Twin Res Hum Genet. 2013 Feb; 16(1): 271-281.) we obtain mean and SD for eheight and weight, we need these because the GIANT height and Weight GWAS are based on an standardized phenotype (mean = 0 variance = 1), this is not the scale we need to compute BMI (as BMI is KG/m2):

Weight_mean <- 83.858621
Height_mean <- 1.81942192


Weight_SD <- 12.981540
Height_SD <- 0.07325379

HW_cor <- 0.370959

# Approximate the SD for BMI:
BMI_SD <- deltamethod(~x1/(x2^2) ,mean=c(Weight_mean,Height_mean),cov=diag(c(Weight_SD,Height_SD)) %*% matrix(c(1,HW_cor,HW_cor,1),2,2) %*% diag(c(Weight_SD,Height_SD )))




# read GIANT sumstats into R

# GIANT_Height_Men <- read.table(file.path(p$folderpath.data,"gwas_sumstats","raw", "GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.txt"), header=T, quote="\"")
# GIANT_Weight_Men <- read.table(file.path(p$folderpath.data,"gwas_sumstats","raw","GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.txt"), header=T, quote="\"")

GIANT_Height_Men <- as.data.frame(shru::readFile(file.path(p$folderpath.data.sumstats.munged, "HEIG04M.gz")))
GIANT_Weight_Men <- as.data.frame(shru::readFile(file.path(p$folderpath.data.sumstats.munged, "WEIG04M.gz")))



# # READ HAPMAP2 REFFREQ This file contains the Alllel frequencies in HAPMAP2, because the GIANT file seems to omit some AF's which are available and needed for GWIS.
# 
# hapmap_reffreq_PHASE2_3_CEU <- read.table("hapmap_reffreq_PHASE2_3_CEU.txt", header=T, quote="\"")



#first MERGE:
  
GIANT_Combined <- merge(GIANT_Weight_Men,GIANT_Height_Men,by=1)

# # get ref MAF From HAPMAP, and merge the files.
# 
# GIANT_Combined <- merge(GIANT_Combined,hapmap_reffreq_PHASE2_3_CEU,by=1,all.x =T)
# dim(GIANT_Combined)
# GIANT_Combined <- na.omit(GIANT_Combined)
# dim(GIANT_Combined)
# 
# 
# 
# # and repalce it for NOT flipped alleles:
# GIANT_Combined[toupper(GIANT_Combined$A1.x)==toupper(GIANT_Combined$refallele) ,]$refallele_freq<- GIANT_Combined[toupper(GIANT_Combined$A1.x)==toupper(GIANT_Combined$refallele) ,]$refallele_freq
# 
# # and replace for flipepd alleles!
# GIANT_Combined[toupper(GIANT_Combined$A1.x)!=toupper(GIANT_Combined$refallele) ,]$refallele_freq <- 1- GIANT_Combined[toupper(GIANT_Combined$A1.x)!=toupper(GIANT_Combined$refallele) ,]$refallele_freq




# compute new beta:

  # Rewscale beta and SE:
  Bw <- Weight_SD *  GIANT_Combined$BETA.x
  Bh <- Height_SD *  GIANT_Combined$BETA.y


  SEw <- Weight_SD *  GIANT_Combined$SE.x
  SEh <- Height_SD *  GIANT_Combined$SE.y

  # Frequency:
  FRQ <- as.vector(GIANT_Combined$FRQ.x)
  

  # Compute af weighted mean of 2 aproximated beta's
  AF1 <- 2*FRQ*(1-FRQ) / (2*FRQ*(1-FRQ) + FRQ^2)

  B1 <- ( ((Weight_mean  + Bw )/ (Height_mean + Bh)^2) - (Weight_mean /(Height_mean^2)))

  AF2 <- ((FRQ^2) ) / (2*FRQ*(1-FRQ) + FRQ^2)

  B2 <-(((Weight_mean  + 2*Bw) / ((Height_mean + 2*Bh)^2)) - ( Weight_mean /(Height_mean^2))) / 2


  new.beta.us <-    (AF1 * B1 + AF2 * B2 )
  new.beta <- (new.beta.us / BMI_SD) 

# compute new SE:



  new.se <- matrix(NA,nrow=length(Bh),ncol=1)
    for(i in 1:length(new.beta)){
      AF1_l <- AF1[i]
      AF2_l <- AF2[i] 
  
      new.se[i]<- deltamethod(~( AF1_l*(((Weight_mean+x1)/((Height_mean+x2)^2)) -(Weight_mean/Height_mean^2)) + (AF2_l*(((Weight_mean+2*x1)/((Height_mean+2*x2)^2)) -(Weight_mean/Height_mean^2)))/2)/BMI_SD, mean=c(Bw[i],Bh[i]),
                               
                               cov=cbind(rbind((SEw[i])^2,           (SEw[i])*HW_cor*(SEh[i])),
                                         rbind((SEw[i])*HW_cor*(SEh[i]),      (SEh[i])^2)))
      }



# build an output file:

GWIS_BMI <- as.data.frame(cbind(GIANT_Combined$SNP,GIANT_Combined$CHR.x,GIANT_Combined$BP.x,GIANT_Combined$A1.x,GIANT_Combined$A2.x,new.beta,new.se,(pchisq((new.beta/new.se)^2,1,lower.tail=F)),(GIANT_Combined$N.x+GIANT_Combined$N.y)/2))

colnames(GWIS_BMI) <- c("SNP","CHR","BP","A1","A2","BETA","SE","P","N")

GWIS_BMI$CHR<-as.integer(GWIS_BMI$CHR)
GWIS_BMI$BP<-as.integer(GWIS_BMI$BP)
GWIS_BMI$P<-as.numeric(GWIS_BMI$P)

#read previous results
#GWIS_BMI <- as.data.frame(shru::readFile(file.path(p$folderpath.data.sumstats.munged, "GWISM.gz")))
#GWISM_BMI <- as.data.frame(shru::readFile(file.path(p$folderpath.data.sumstats.munged, "GWISI.gz")))

# PLOT GWIS Mwenahttan:



#manhattan(na.omit(GWIS_BMI),chr="GIANT_Combined$chrom",bp="GIANT_Combined$pos",p="(pchisq((new.beta/new.se)^2, 1, lower.tail = F))",ylim=c(0,50))

cPlot<-plot.manhattan.custom(GWIS_BMI,var = "P")
ggsave(filename=file.path(p$folderpath.plots, paste0("manp.GWIS.GWISM.png")), plot = cPlot, width = 297, height = 210, units = "mm", dpi = 320, scale = 3)


#  Read in the real BMI results:


REFERENCE_BMI <- as.data.frame(shru::readFile(file.path(p$folderpath.data.sumstats.munged, "BODYTM.gz")))
REFERENCE_BMI2 <- as.data.frame(shru::readFile(file.path(p$folderpath.data.sumstats.munged, "BODY11.gz")))


# Merge to test:

GWIS_GWAS <- merge(GWIS_BMI,REFERENCE_BMI,by=1)
GWIS_GWAS.f <- GWIS_GWAS[GWIS_GWAS$N.x > 58000 & GWIS_GWAS$N.y > 58000,]


cor(GWIS_GWAS.f$BETA.x,GWIS_GWAS.f$BETA.y,use="pairwise.co")^2
cor(GWIS_GWAS.f$SE.x,GWIS_GWAS.f$SE.y,use="pairwise.co")^2
cor(GWIS_GWAS.f$Z.x,GWIS_GWAS.f$Z.y,use="pairwise.co")^2


tiff(file.path(p$folderpath.plots,"GWISM_REF1_corr.tiff"), width = 4, height = 8, units = 'in', res = 300)
# Make plot

par(mfrow=c(3,1))

plot(GWIS_GWAS.f$BETA.x,GWIS_GWAS.f$BETA.y, col=rgb(0,0,1,.2),pch=19,ylab="Beta's GWAS BMI", xlab="Beta's GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .935))
abline(0,1,col="red")

plot(GWIS_GWAS.f$SE.x,GWIS_GWAS.f$SE.y, col=rgb(0,0,1,.2),pch=19,ylab="SE GWAS BMI", xlab="SE GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .996))
abline(0,1,col="red")

plot(GWIS_GWAS.f$Z.x,GWIS_GWAS.f$Z.y, col=rgb(0,0,1,.2),pch=19,ylab="Z GWAS BMI", xlab="Z GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .936))
abline(0,1,col="red")

dev.off()


#`1000G_rs_ch_bp` <- read.table("C:/8. post hoc GWAS/FINAL Analysis paper/1000G_rs_ch_bp", quote="\"")


#manhattan(na.omit(GWIS_GWAS),chr="GIANT_Combined$chrom",bp="GIANT_Combined$pos",p="P.2gc",ylim=c(0,50))
cPlot<-plot.manhattan.custom(REFERENCE_BMI,var = "P")
ggsave(filename=file.path(p$folderpath.plots, paste0("manp.GWIS.BODYTM.png")), plot = cPlot, width = 297, height = 210, units = "mm", dpi = 320, scale = 3)



# Merge to test 2:

GWIS_GWAS2 <- merge(GWIS_BMI,REFERENCE_BMI2,by=1)
GWIS_GWAS2.f <- GWIS_GWAS2[GWIS_GWAS2$N.x > 58000 & GWIS_GWAS2$N.y > 58000,]


cor(GWIS_GWAS2.f$BETA.x,GWIS_GWAS2.f$BETA.y,use="pairwise.co")^2
cor(GWIS_GWAS2.f$SE.x,GWIS_GWAS2.f$SE.y,use="pairwise.co")^2
cor(GWIS_GWAS2.f$Z.x,GWIS_GWAS2.f$Z.y,use="pairwise.co")^2


tiff(file.path(p$folderpath.plots,"GWISM_REF2_corr.tiff"), width = 4, height = 8, units = 'in', res = 300)
# Make plot

par(mfrow=c(3,1))

plot(GWIS_GWAS2.f$BETA.x,GWIS_GWAS2.f$BETA.y, col=rgb(0,0,1,.2),pch=19,ylab="Beta's GWAS BMI", xlab="Beta's GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .155))
abline(0,1,col="red")

plot(GWIS_GWAS2.f$SE.x,GWIS_GWAS2.f$SE.y, col=rgb(0,0,1,.2),pch=19,ylab="SE GWAS BMI", xlab="SE GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .89))
abline(0,1,col="red")

plot(GWIS_GWAS2.f$Z.x,GWIS_GWAS2.f$Z.y, col=rgb(0,0,1,.2),pch=19,ylab="Z GWAS BMI", xlab="Z GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .185))
abline(0,1,col="red")

dev.off()


#`1000G_rs_ch_bp` <- read.table("C:/8. post hoc GWAS/FINAL Analysis paper/1000G_rs_ch_bp", quote="\"")


#manhattan(na.omit(GWIS_GWAS),chr="GIANT_Combined$chrom",bp="GIANT_Combined$pos",p="P.2gc",ylim=c(0,50))
cPlot<-plot.manhattan.custom(REFERENCE_BMI2,var = "P")
ggsave(filename=file.path(p$folderpath.plots, paste0("manp.GWIS.BODY11.png")), plot = cPlot, width = 297, height = 210, units = "mm", dpi = 320, scale = 3)




# B B Plot
#plot(GWIS_GWAS[,5],GWIS_GWAS[,12],col=rgb(0,0,1,.1),ylim=c(-.35,.35),pch=19,ylab="BMI GWAS Beta",xlab="BMI GWIS Beta")
# SE SE Plot
#plot(GWIS_GWAS[,6],GWIS_GWAS[,13],col=rgb(0,0,1,.1),ylim=c(0,.16),pch=19,ylab="BMI GWAS SE",xlab="BMI GWIS SE")
# P P Plot
plot(-log10(GWIS_GWAS$P.x),-log10(GWIS_GWAS$P.y),col=rgb(0,0,1,.1),ylim=c(0,16),xlim=c(0,16),pch=19,ylab="BMI GWAS -log10(P)",xlab="BMI GWIS -log10(P)")

# Number of Hits GWAS:
nrow(na.omit(GWIS_GWAS[GWIS_GWAS$P.x < 5e-8,]))

# Numer of hits replicated:
sum(na.omit(GWIS_GWAS[GWIS_GWAS$P.x < 5e-8,]$P.y < 5e-8))

# numer of false positives:
sum(na.omit(GWIS_GWAS[GWIS_GWAS$P.x < 5e-8,]$P.y > 5e-8))


# Table the false positives for inclusion in a supplment:

fp <-na.omit(GWIS_GWAS[GWIS_GWAS$P.x < 5e-8 & GWIS_GWAS$P.y > 5e-8 ,])

write.table(fp,"false_positive_male.txt")


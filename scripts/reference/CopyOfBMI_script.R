#Michel Nivard's GWIS of height and weight => BMI
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
BMI_SD <- deltamethod(~x1/(x2^2) ,mean=c(83.858621,1.81942192),cov=diag(c(12.981540,0.07325379 )) %*% matrix(c(1,HW_cor,HW_cor,1),2,2) %*% diag(c(12.981540,0.07325379 )))




# read GIANT sumstats into R

GIANT_Height_Men <- read.table("GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.txt", header=T, quote="\"")
GIANT_Weight_Men <- read.table("CGIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.txt", header=T, quote="\"")




# READ HAPMAP2 REFFREQ This file contains the Alllel frequencies in HAPMAP2, because the GIANT file seems to omit some AF's which are available and needed for GWIS.

hapmap_reffreq_PHASE2_3_CEU <- read.table("hapmap_reffreq_PHASE2_3_CEU.txt", header=T, quote="\"")



#first MERGE:
  
GIANT_Combined <- merge(GIANT_Weight_Men,GIANT_Height_Men,by=1)

# get ref MAF From HAPMAP, and merge the files.

GIANT_Combined <- merge(GIANT_Combined,hapmap_reffreq_PHASE2_3_CEU,by=1,all.x =T)
dim(GIANT_Combined)
GIANT_Combined <- na.omit(GIANT_Combined)
dim(GIANT_Combined)



# and repalce it for NOT flipped alleles:
GIANT_Combined[toupper(GIANT_Combined$A1.x)==toupper(GIANT_Combined$refallele) ,]$refallele_freq<- GIANT_Combined[toupper(GIANT_Combined$A1.x)==toupper(GIANT_Combined$refallele) ,]$refallele_freq

# and replace for flipepd alleles!
GIANT_Combined[toupper(GIANT_Combined$A1.x)!=toupper(GIANT_Combined$refallele) ,]$refallele_freq <- 1- GIANT_Combined[toupper(GIANT_Combined$A1.x)!=toupper(GIANT_Combined$refallele) ,]$refallele_freq




# compute new beta:

  # Rewscale beta and SE:
  Bw <- Weight_SD *  GIANT_Combined [,5]
  Bh <- Height_SD *  GIANT_Combined [,12]


  SEw <- Weight_SD *  GIANT_Combined [,6]
  SEh <- Height_SD *  GIANT_Combined [,13]

  # Frequency:
  p <- as.vector(GIANT_Combined$refallele_freq)
  

  # Compute af weighted mean of 2 aproximated beta's
  AF1 <- 2*p*(1-p) / (2*p*(1-p) + p^2)

  B1 <- ( ((Weight_mean  + Bw )/ (Height_mean + Bh)^2) - (Weight_mean /(Height_mean^2)))

  AF2 <- ((p^2) ) / (2*p*(1-p) + p^2)

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

GWIS_BMI <- cbind(GIANT_Combined[,1:4],new.beta,new.se,(pchisq((new.beta/new.se)^2,1,lower.tail=F)),apply(GIANT_Combined[,c(8,15)],1,min),GIANT_Combined$chrom,GIANT_Combined$pos)


# remove "chr" from chromose indicator and make variable numeric

GWIS_BMI[,9] <- as.numeric(substring(GWIS_BMI[,9],4))

# PLOT GWIS Mwenahttan:


manhattan(na.omit(GWIS_BMI),chr="GIANT_Combined$chrom",bp="GIANT_Combined$pos",p="(pchisq((new.beta/new.se)^2, 1, lower.tail = F))",ylim=c(0,50))


# writwe file for LD score regression
write.table(GWIS_BMI,"GWIS_BMI_MEN_LDSCORE.txt",row.names=F,quote=F )

#  Read in the real BMI results:


GIANT_BMI_Men <- read.table("C:/8. post hoc GWAS/FINAL Analysis paper/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.txt", header=T, quote="\"")



# Merge to test:

GWIS_GWAS <- merge(GWIS_BMI,GIANT_BMI_Men,by=1)


cor(GWIS_GWAS$new.beta[GWIS_GWAS[,8] > 58000 ],GWIS_GWAS$BETA[GWIS_GWAS[,8] > 58000],use="pairwise.co")^2


cor(GWIS_GWAS$new.se[GWIS_GWAS[,8] > 58000 ],GWIS_GWAS$SE.2gc[GWIS_GWAS[,8] > 58000],use="pairwise.co")^2



cor(GWIS_GWAS$new.beta[GWIS_GWAS[,8] > 58000 ]/GWIS_GWAS$new.se[GWIS_GWAS[,8] > 58000 ],GWIS_GWAS$BETA[GWIS_GWAS[,8] > 58000 ]/GWIS_GWAS$SE.2gc[GWIS_GWAS[,8] > 58000 ],use="pairwise.co")^2



tiff("Figure2 US.tiff", width = 4, height = 8, units = 'in', res = 300)
# Make plot

par(mfrow=c(3,1))


plot(GWIS_GWAS$new.beta[GWIS_GWAS[,8] > 58000 ],GWIS_GWAS$BETA[GWIS_GWAS[,8] > 58000 ], col=rgb(0,0,1,.2),pch=19,ylab="Beta's GWAS BMI", xlab="Beta's GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .937))
abline(0,1,col="red")

plot(GWIS_GWAS$new.se[GWIS_GWAS[,8] > 58000 ],GWIS_GWAS$SE.2gc[GWIS_GWAS[,8] > 58000 ], col=rgb(0,0,1,.2),pch=19,ylab="SE GWAS BMI", xlab="SE GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .995))
abline(0,1,col="red")

plot(GWIS_GWAS$new.beta[GWIS_GWAS[,8] > 58000 ]/GWIS_GWAS$new.se[GWIS_GWAS[,8] > 58000 ],GWIS_GWAS$BETA[GWIS_GWAS[,8] > 58000 ]/GWIS_GWAS$SE.2gc[GWIS_GWAS[,8] > 58000 ], col=rgb(0,0,1,.2),pch=19,ylab="Z GWAS BMI", xlab="Z GWIS BMI Male Height" )
legend("topleft",expression(R^2 ~ .939))
abline(0,1,col="red")



dev.off()


#`1000G_rs_ch_bp` <- read.table("C:/8. post hoc GWAS/FINAL Analysis paper/1000G_rs_ch_bp", quote="\"")


manhattan(na.omit(GWIS_GWAS),chr="GIANT_Combined$chrom",bp="GIANT_Combined$pos",p="P.2gc",ylim=c(0,50))





# B B Plot
#plot(GWIS_GWAS[,5],GWIS_GWAS[,12],col=rgb(0,0,1,.1),ylim=c(-.35,.35),pch=19,ylab="BMI GWAS Beta",xlab="BMI GWIS Beta")
# SE SE Plot
#plot(GWIS_GWAS[,6],GWIS_GWAS[,13],col=rgb(0,0,1,.1),ylim=c(0,.16),pch=19,ylab="BMI GWAS SE",xlab="BMI GWIS SE")
# P P Plot
plot(-log10(GWIS_GWAS[,7]),-log10(GWIS_GWAS[,14]),col=rgb(0,0,1,.1),ylim=c(0,16),xlim=c(0,16),pch=19,ylab="BMI GWAS -log10(P)",xlab="BMI GWIS -log10(P)")

# Number of Hits GWAS:
nrow(na.omit(GWIS_GWAS[GWIS_GWAS[,16] < 5e-8,]))

# Numer of hits replicated:
sum(na.omit(GWIS_GWAS[GWIS_GWAS[,16] < 5e-8,7] < 5e-8))

# numer of false positives:
sum(na.omit(GWIS_GWAS[GWIS_GWAS[,7] < 5e-8,16] > 5e-8))


# Table the false positives for inclusion in a supplment:

fp <-na.omit(GWIS_GWAS[GWIS_GWAS[,7] < 5e-8 & GWIS_GWAS[,16] > 5e-8 ,])

write.table(fp,"false_positive_male.txt")


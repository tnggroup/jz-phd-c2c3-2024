#!/bin/bash -l

#module add devtools/anaconda/2019.7-python2.7.16

PYTHONPATH=python3
#PYTHONPATH=python2.7
#PYTHONPATH=python

#MUNGEPATH=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/ldscore_2020_06/ldsc/munge_sumstats.py

LDSCPATH=/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc/ldsc.py
#LDSCPATH=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/ldscore_2020_06/ldsc/ldsc.py

RGLIST='/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/munged/DEPR05.gz,/users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/munged/DEPR08.gz'
#RGLIST='/scratch/groups/ukbiobank/sumstats/munged/DEPR05.sumstats.gz,/scratch/groups/ukbiobank/sumstats/munged/DEPR08.sumstats.gz'
#RGLIST='DEPR05.temp.sumstats.gz,DEPR08.temp.sumstats.gz'

LDREFERENCEPATH=/scratch/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/eur_w_ld_chr.1KG_Phase3/
#LDREFERENCEPATH=/scratch/users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/data/eur_w_ld_chr/
#LDREFERENCEPATH=/mnt/lustre/groups/ukbiobank/Edinburgh_Data/Software/ldscore_2020_06/eur_ref_ld_chr/

source /users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc-venv/bin/activate

#munge if needed
#$MUNGEPATH \
#--out DEPR05.temp \
#--merge-alleles /scratch/users/k19049801/project/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt \
#--sumstats /users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/cleaned/DEPR05.gz


#$MUNGEPATH \
#--out DEPR08.temp \
#--merge-alleles /scratch/users/k19049801/project/JZ_GED_PHD_C1/data/combined.hm3_1kg.snplist.vanilla.jz2020.txt \
#--sumstats /users/k19049801/project/JZ_GED_PHD_C1/data/gwas_sumstats/cleaned/DEPR08.gz


$LDSCPATH \
  --n-blocks 200 \
  --rg ${RGLIST} \
  --ref-ld-chr "${LDREFERENCEPATH}" \
  --w-ld-chr "${LDREFERENCEPATH}" \
  --out setup3.ldsc.txt \
  --samp-prev 0.3962549,0.3459495 \
  --pop-prev 0.146,0.146 \
  --no-check-alleles 



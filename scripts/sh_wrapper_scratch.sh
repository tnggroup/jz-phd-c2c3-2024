#!/bin/bash
#Load the R module
module add apps/R/3.6.0
#module add apps/R/3.6.3

srun -p brc,shared --ntasks 1 --cpus-per-task 3 --mem 8G --pty /bin/bash

sbatch --time 23:59:00 --partition brc,shared --job-name="GSEMGWAS" --ntasks 1 --cpus-per-task 6 --mem-per-cpu 8G --wrap="Rscript setup1.R -l cluster" --output "setup1_$(date +%Y%m%d).out.txt" --error "setup1_$(date +%Y%m%d).err.txt"

sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.mvLDSC" --ntasks 1 --cpus-per-task 3 --mem-per-cpu 12G --wrap="Rscript setup2.R -t mvLD.mvLDSC -l cluster" --output "setup2.mvLD.mvLDSC.out" --error "setup2.mvLD.mvLDSC.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.HDL.piecewise" --ntasks 1 --cpus-per-task 3 --mem 64G --wrap="Rscript setup2.R -t mvLD.HDL.piecewise -l cluster" --output "setup2.mvLD.HDL.piecewise.out" --error "setup2.mvLD.HDL.piecewise.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.HDL.jackknife" --ntasks 1 --cpus-per-task 3 --mem 80G --wrap="Rscript setup2.R -t mvLD.HDL.jackknife -l cluster" --output "setup2.mvLD.HDL.jackknife.out" --error "setup2.mvLD.HDL.jackknife.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.origHDL" --ntasks 1 --cpus-per-task 3 --mem 80G --wrap="Rscript setup2.R -t mvLD.origHDL -l cluster" --output "setup2.mvLD.origHDL.out" --error "setup2.mvLD.origHDL.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.origHDL.liabilityScale" --ntasks 1 --cpus-per-task 3 --mem 90G --wrap="Rscript setup2.R -t mvLD.origHDL.liabilityScale -l cluster" --output "setup2.mvLD.origHDL.liabilityScale.out" --error "setup2.mvLD.origHDL.liabilityScale.err"

sbatch --time 23:59:00 --partition brc,shared --job-name="munge" --ntasks 1 --cpus-per-task 3 --mem-per-cpu 10G --wrap="Rscript setup2.R -t munge -l cluster" --output "setup2.munge.$(date +%Y%m%d).out.txt" --error "setup2.munge.$(date +%Y%m%d).err.txt"

#raiss --gwas 'ALCD03.chr' --ref-folder '/users/k1204688/brc_scratch/Public/1KG_Phase3/All' --ld-folder '../data/eur_w_ld_chr' --zscore-folder '../data/gwas_sumstats/munged' --output-folder '../data/gwas_sumstats/imputed'
#source ../python-venv/bin/activate
raiss --chrom 1 --gwas 'ALCD03' --ref-folder '../data/reference.panel.1KG_Phase3' --ld-folder '../data/eur_w_ld_chr' --zscore-folder '../data/gwas_sumstats/munged/ALCD03.chr' --output-folder '../data/gwas_sumstats/imputed' --l2-regularization 0.01 --eigen-threshold 0.05 --R2-threshold 0.3

#python ../../software/ldsc/ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do sbatch --time 23:59:00 --partition brc,shared --job-name="LDSC" --ntasks 1 --cpus-per-task 3 --mem 15G --wrap="../../software/ldsc/ldsc.py --bfile $chr --l2 --ld-wind-cm 1 --out $chr" --output "ldsc_${chr}_$(date +%Y%m%d).out.txt" --error "ldsc_${chr}_$(date +%Y%m%d).err.txt"; done

#raiss --chrom 1 --gwas 'ALCD03' --ref-folder '../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM' --ld-folder '../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM' --zscore-folder '../data/gwas_sumstats/munged/ALCD03.chr' --output-folder '../data/gwas_sumstats/imputed' --l2-regularization 0.01 --eigen-threshold 0.05 --R2-threshold 0.3
mkdir -p ../data/gwas_sumstats/imputed/ALCD03
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23; do raiss --chrom $chr --gwas 'ALCD03' --ref-folder '../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM' --ld-folder '../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM' --zscore-folder '../data/gwas_sumstats/munged/ALCD03.chr' --output-folder '../data/gwas_sumstats/imputed/ALCD03' --l2-regularization 0.01 --eigen-threshold 0.05 --R2-threshold 0.3; done

sbatch --time 23:59:00 --partition brc,shared --job-name="impALCD03" --ntasks 1 --cpus-per-task 3 --mem 20G --wrap="source ../python-venv/bin/activate; for chr in {1..23}; do raiss --chrom \$chr --gwas ALCD03 --ref-folder ../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM --ld-folder ../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM --zscore-folder ../data/gwas_sumstats/munged/ALCD03.chr --output-folder ../data/gwas_sumstats/imputed/ALCD03 --l2-regularization 0.01 --eigen-threshold 0.05 --R2-threshold 0.3; done" --output impute.ALCD03.$(date +%Y%m%d).out.txt --error impute.ALCD03.$(date +%Y%m%d).err.txt


#awk -f ../../../../scripts/merge.imputed.dataset.awk z_${dscode}_1.txt z_${dscode}_2.txt z_${dscode}_3.txt z_${dscode}_4.txt z_${dscode}_5.txt z_${dscode}_6.txt z_${dscode}_7.txt z_${dscode}_8.txt z_${dscode}_9.txt z_${dscode}_10.txt z_${dscode}_11.txt z_${dscode}_12.txt z_${dscode}_13.txt z_${dscode}_14.txt z_${dscode}_15.txt z_${dscode}_1.txt z_${dscode}_1.txt z_${dscode}_1.txt z_${dscode}_1.txt z_${dscode}_1.txt z_${dscode}_1.txt z_${dscode}_1.txt
dscode=ALCD03; awk -v prefix="z_${dscode}_" -v suffix=".txt" -f ../../../../scripts/merge.imputed.dataset.awk > $dscode;
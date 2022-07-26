#!/bin/bash
#Load the R module
module add apps/R/3.6.0
#module add apps/R/3.6.3

#ok for smaller work
srun -p brc,shared --ntasks 1 --cpus-per-task 3 --mem 8G --pty /bin/bash
#ok for ldsc munge using 1KG
srun -p brc,shared --ntasks 1 --cpus-per-task 3 --mem 16G --pty /bin/bash
srun -p cpu --ntasks 1 --cpus-per-task 5 --mem 32G --pty /bin/bash

#misc
#merge plink datasets
plink --merge-list ~/project/JZ_GED_PHD_C1/working_directory/plink.merge.txt --make-bed --out ~/project/JZ_GED_PHD_C1/working_directory/ref4validation

rsync -avzhpt --progress /mnt/lustre/groups/ukbiobank/sumstats/ /scratch/groups/gwas_sumstats/
#sbatch --time 2-00:00:00 --partition brc,shared --job-name="rsync" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="rsync -avzhpt /mnt/lustre/groups/ukbiobank/sumstats/ /scratch/groups/gwas_sumstats/" --output "rsync.$(date +%Y%m%d).out.txt"
#This command may have added duplicate subfolders under each folder
ls /mnt/lustre/groups/ukbiobank/sumstats | xargs -n 1 -P 4 -I % rsync -avzhpt /mnt/lustre/groups/ukbiobank/sumstats/% /scratch/groups/gwas_sumstats/%
rsync -avzhpt --progress /rosalind_ro/scratch/groups/gwas_sumstats/ /scratch/prj/gwas_sumstats/
sbatch --time 2-00:00:00 --partition cpu --job-name="rsync" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="rsync -avzhpt /rosalind_ro/scratch/groups/gwas_sumstats/ /scratch/prj/gwas_sumstats/" --output "rsync.$(date +%Y%m%d).out.txt"

#Thibault's gcta command
../../JZ_GED_PHD_ADMIN_GENERAL/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta --bfile ../data/reference_panel/1KG_eur.plink/g1000_eur --mtcojo-file mtcojo_summary_data.list --ref-ld-chr ../data/ld_scores/eur_w_ld_chr.1KG_Phase3/ --w-ld-chr ../data/ld_scores/eur_w_ld_chr.1KG_Phase3/ --out test_mtcojo_result


#work

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

sbatch --time 23:59:00 --partition brc,shared --job-name="gsem1" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 0:9999 -l cluster" --output "setup2.cfa.1.out" --error "setup2.cfa.1.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem2" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 10000:19999 -l cluster" --output "setup2.cfa.2.out" --error "setup2.cfa.2.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem3" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 20000:29999 -l cluster" --output "setup2.cfa.3.out" --error "setup2.cfa.3.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem4" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 30000:39999 -l cluster" --output "setup2.cfa.4.out" --error "setup2.cfa.4.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem5" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 40000:49999 -l cluster" --output "setup2.cfa.5.out" --error "setup2.cfa.5.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem6" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 50000:59999 -l cluster" --output "setup2.cfa.6.out" --error "setup2.cfa.6.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem7" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 60000:69999 -l cluster" --output "setup2.cfa.7.out" --error "setup2.cfa.7.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem8" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 70000:79999 -l cluster" --output "setup2.cfa.8.out" --error "setup2.cfa.8.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem9" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 80000:89999 -l cluster" --output "setup2.cfa.9.out" --error "setup2.cfa.9.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem10" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 90000:99999 -l cluster" --output "setup2.cfa.10.out" --error "setup2.cfa.10.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem11" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 100000:109999 -l cluster" --output "setup2.cfa.11.out" --error "setup2.cfa.11.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem12" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 110000:119999 -l cluster" --output "setup2.cfa.12.out" --error "setup2.cfa.12.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="gsem13" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a 120000:131072 -l cluster" --output "setup2.cfa.13.out" --error "setup2.cfa.13.err"

#first run
#for gsemi in `seq 0 999 131072`; do gsemi2=$((gsemi+999)); sbatch --time 23:59:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a $gsemi:$gsemi2 -l cluster" --output "setup2.cfa.$gsemi-$gsemi2.out" --error "setup2.cfa.$gsemi-$gsemi2.err"; done
#second run
for gsemi in `seq 0 1000 32768`; do gsemi2=$((gsemi+1000)); sbatch --time 23:59:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a $gsemi:$gsemi2 -l cluster" --output "setup2.cfa.$gsemi-$gsemi2.out" --error "setup2.cfa.$gsemi-$gsemi2.err"; done
#third run
for gsemi in `seq 0 500 4096`; do gsemi2=$((gsemi+500)); sbatch --time 23:59:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup2.R -t cfa -a $gsemi:$gsemi2 -l cluster" --output "setup2.cfa.$gsemi-$gsemi2.out" --error "setup2.cfa.$gsemi-$gsemi2.err"; done

sbatch --time 23:59:00 --partition brc,shared --job-name="prepgwas" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup2.R -l cluster" --output "setup2.prepgwas.out" --error "setup2.prepgwas.err"

sbatch --time 2-00:00:00 --partition brc,shared --job-name="gsemgwas" --ntasks 1 --cpus-per-task 24 --mem-per-cpu 6G --oversubscribe --wrap="Rscript setup2.R -l cluster" --output "setup2.gsemgwas.$(date +%Y%m%d).out.txt" --error "setup2.gsemgwas.$(date +%Y%m%d).err.txt"

#setup3

sbatch --time 2-00:00:00 --partition brc,shared --job-name="mvLD.mvLDSC" --ntasks 1 --cpus-per-task 3 --mem 64G --wrap="Rscript setup3.R -t mvLD.mvLDSC -l cluster" --output "setup3.mvLD.mvLDSC.out" --error "setup3.mvLD.mvLDSC.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="LDSC" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="sh setup3.ldsc.sh" --output "setup3LDSC.out" --error "setup3LDSC.err"

sbatch --time 23:59:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 3 --mem 40G --wrap="Rscript setup3.R -t munge -l cluster" --output "setup3.munge.$(date +%Y%m%d).out.txt"
sbatch --time 23:00:00 --partition brc,shared --job-name="preplfgwas" --ntasks 1 --cpus-per-task 2 --mem 64G --wrap="Rscript setup3.R -t preplfgwas -l cluster" --output "setup3.lfGWAS.sumstats.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition brc,shared --job-name="mvLD.mvLDSC" --ntasks 1 --cpus-per-task 3 --mem 64G --wrap="Rscript setup3.R -t mvLD.mvLDSC -l cluster" --output "setup3.mvLD.mvLDSC.$(date +%Y%m%d).out.txt"
gsemi=0; gsemi2=$((gsemi+250-1)); sbatch --time 2-00:00:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t cfa -a $gsemi:$gsemi2 -l cluster" --output "setup3.cfa.$gsemi-$gsemi2.out.txt";
for gsemi in `seq 250 250 $((4096+1))`; do gsemi2=$((gsemi+250-1)); sbatch --time 2-00:00:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t cfa -a $gsemi:$gsemi2 -l cluster" --output "setup3.cfa.$gsemi-$gsemi2.out.txt"; done
sbatch --time 23:00:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t cfa -l cluster" --output "setup3.cfa.$(date +%Y%m%d).out.txt"

sbatch --time 2-00:00:00 --partition brc,shared --job-name="lfgwas" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t lfgwas -a 1 -l cluster" --output "setup3.lfgwas.chr1.$(date +%Y%m%d).out.txt"
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lfgwas" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t lfgwas -a $chr -l cluster" --output "setup3.lfgwas.chr$chr.$(date +%Y%m%d).out.txt"; done

#run lf1 chr1 first as a test
lf=1;for chr in 1; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$lf:$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t lfgwas -a $lf:$chr -l cluster" --output "setup3.lfgwas.F$lf.chr$chr.$(date +%Y%m%d).out.txt"; done
lf=1;for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$lf:$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t lfgwas -a $lf:$chr -l cluster" --output "setup3.lfgwas.F$lf.chr$chr.$(date +%Y%m%d).out.txt"; done
#lf=1;for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do Rscript setup3.R -t lfgwas -a $lf:$chr > setup3.lfgwas.chr$chr.$(date +%Y%m%d).out.txt; done

for chr in 1; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t lfgwas -a $chr -l cluster" --output "setup3.lfgwas.F_ALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup3.R -t lfgwas -a $chr -l cluster" --output "setup3.lfgwas.F_ALL.chr$chr.$(date +%Y%m%d).out.txt"; done

#setup4
sbatch --time 23:59:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup4.R -t munge -l cluster" --output "setup4.munge.$(date +%Y%m%d).out.txt"
Rscript setup4.R -t munge > "setup4.munge.$(date +%Y%m%d).out.txt" #alt
sbatch --time 2-00:00:00 --partition brc,shared --job-name="mvLD" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup4.R -t mvLD -l cluster" --output "setup4.mvLD.$(date +%Y%m%d).out.txt"
sbatch --time 4-00:00:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 3 --mem 32G --wrap="Rscript setup4.R -t cfa -l cluster" --output "setup4.cfa.$(date +%Y%m%d).out.txt"
sbatch --time 4-00:00:00 --partition brc,shared --job-name="preplfgwas" --ntasks 1 --cpus-per-task 3 --mem 64G --wrap="Rscript setup4.R -t preplfgwas -l cluster" --output "setup4.lfGWAS.sumstats.$(date +%Y%m%d).out.txt"
for chr in 1; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup4.R -t lfgwas -a M25-4.74.ML:$chr -l cluster" --output "setup4.lfgwas.M25-4.74.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup4.R -t lfgwas -a M25-4.74.ML:$chr -l cluster" --output "setup4.lfgwas.M25-4.74.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 1; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup4.R -t lfgwas -a M25-7.29.ML:$chr -l cluster" --output "setup4.lfgwas.M25-7.29.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup4.R -t lfgwas -a M25-7.29.ML:$chr -l cluster" --output "setup4.lfgwas.M25-7.29.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 1; do sbatch --time 4-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup4.R -t lfgwas -a M25-10.55.ML:$chr -l cluster" --output "setup4.lfgwas.M25-10.55.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 4-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup4.R -t lfgwas -a M25-10.55.ML:$chr -l cluster" --output "setup4.lfgwas.M25-10.55.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
sbatch --time 2-00:00:00 --partition brc,shared --job-name="mvLD2" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup4.R -t mvLD2 -l cluster" --output "setup4.mvLD2.$(date +%Y%m%d).out.txt"

#setup5
sbatch --time 2-00:00:00 --partition brc,shared --job-name="mvLD" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup5.R -t mvLD -l cluster" --output "setup5.mvLD.$(date +%Y%m%d).out.txt"
sbatch --time 4-00:00:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 3 --mem 32G --wrap="Rscript setup5.R -t cfa -l cluster" --output "setup5.cfa.$(date +%Y%m%d).out.txt"
sbatch --time 4-00:00:00 --partition brc,shared --job-name="lga1" --ntasks 1 --cpus-per-task 3 --mem 24G --wrap="Rscript setup5.R -t lfgwas -a M18-4.64.ML -l cluster" --output "setup5.lfgwas.M18-4.64.ML.FALL.chrALL.$(date +%Y%m%d).out.txt"
sbatch --time 4-00:00:00 --partition brc,shared --job-name="preplfgwas" --ntasks 1 --cpus-per-task 3 --mem 64G --wrap="Rscript setup5.R -t preplfgwas -l cluster" --output "setup5.lfGWAS.sumstats.$(date +%Y%m%d).out.txt"
for chr in 1; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup5.R -t lfgwas -a M18_4_14.C.DWLS:$chr -l cluster" --output "setup5.lfgwas.M18_4_14.C.DWLS.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup5.R -t lfgwas -a M18_4_14.C.DWLS:$chr -l cluster" --output "setup5.lfgwas.M18_4_14.C.DWLS.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
#sbatch --time 4-00:00:00 --partition brc,shared --job-name="lga1" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup5.R -t lfgwas -a M18_4_7.O.DWLS -l cluster" --output "setup5.lfgwas.M18_4_7.O.DWLS.FALL.chr$chr.$(date +%Y%m%d).out.txt" #ran out of memory?

#setup6
sbatch --time 8:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup6.R -t munge -l cluster" --output "setup6.munge.$(date +%Y%m%d).out.txt"
sbatch --time 12:00:00 --partition brc,shared --job-name="gsem1" --ntasks 1 --cpus-per-task 3 --mem 32G --wrap="Rscript setup6.R -t cfa -a 1 -l cluster" --output "setup6.cfa1.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition brc,shared --job-name="gsem2" --ntasks 1 --cpus-per-task 3 --mem 32G --wrap="Rscript setup6.R -t cfa -a 2 -l cluster" --output "setup6.cfa2.$(date +%Y%m%d).out.txt"
sbatch --time 4-00:00:00 --partition brc,shared --job-name="gsem3" --ntasks 1 --cpus-per-task 3 --mem 32G --wrap="Rscript setup6.R -t cfa -a 3 -l cluster" --output "setup6.cfa3.$(date +%Y%m%d).out.txt"
for chr in 1; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup6.R -t lfgwas -a M25_4_27.COR.ML:$chr -l cluster" --output "setup6.lfgwas.M25_4_27.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup6.R -t lfgwas -a M25_4_27.COR.ML:$chr -l cluster" --output "setup6.lfgwas.M25_4_27.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 1; do sbatch --time 4-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup6.R -t lfgwas -a M25_7_195.COR.ML:$chr -l cluster" --output "setup6.lfgwas.M25_7_195.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do sbatch --time 4-00:00:00 --partition brc,shared --job-name="lga1_$chr" --ntasks 1 --cpus-per-task 2 --mem 16G --wrap="Rscript setup6.R -t lfgwas -a M25_7_195.COR.ML:$chr -l cluster" --output "setup6.lfgwas.M25_7_195.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done

sbatch --time 2-00:00:00 --partition brc,shared --job-name="nmeta" --ntasks 1 --cpus-per-task 4 --mem 16G --wrap="Rscript setup6.R -t nmeta -a M25_7_26.ORT.ML:5 -l cluster" --output "setup6.nmeta.M25_7_26.ORT.ML:5.$(date +%Y%m%d).out.txt"
sbatch --time 3-00:00:00 --partition brc,shared --job-name="nmeta" --ntasks 1 --cpus-per-task 4 --mem 16G --wrap="Rscript setup6.R -t nmeta -a M25_7_26.ORT.ML:6 -l cluster" --output "setup6.nmeta.M25_7_26.ORT.ML:6.$(date +%Y%m%d).out.txt"
sbatch --time 3-00:00:00 --partition brc,shared --job-name="nmeta" --ntasks 1 --cpus-per-task 4 --mem 16G --wrap="Rscript setup6.R -t nmeta -a M25_10_244.ORT.ML:1 -l cluster" --output "setup6.nmeta.M25_10_244.ORT.ML:1.$(date +%Y%m%d).out.txt"

#setup7
#munging with the new basic ld-imputation procedure. slow, running datasets in parallell.
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a ADHD05 -l cluster" --output "setup7.munge.ADHD05.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a ALCD03 -l cluster" --output "setup7.munge.ALCD03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a ANOR02 -l cluster" --output "setup7.munge.ANOR02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a ANXI03 -l cluster" --output "setup7.munge.ANXI03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a ANXI04 -l cluster" --output "setup7.munge.ANXI04.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a AUTI07 -l cluster" --output "setup7.munge.AUTI07.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a BIPO02 -l cluster" --output "setup7.munge.BIPO02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a DEPR05 -l cluster" --output "setup7.munge.DEPR05.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a DEPR08 -l cluster" --output "setup7.munge.DEPR08.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a EDUC03 -l cluster" --output "setup7.munge.EDUC03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a HEAL01 -l cluster" --output "setup7.munge.HEAL01.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a INCO03 -l cluster" --output "setup7.munge.INCO03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a INSO02 -l cluster" --output "setup7.munge.INSO02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a NEUR02 -l cluster" --output "setup7.munge.NEUR02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a PTSD04 -l cluster" --output "setup7.munge.PTSD04.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a RISK02 -l cluster" --output "setup7.munge.RISK02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a RISK03 -l cluster" --output "setup7.munge.RISK03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a SCHI04 -l cluster" --output "setup7.munge.SCHI04.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a SUBJ01 -l cluster" --output "setup7.munge.SUBJ01.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 24G --wrap="Rscript setup7.R -t munge -a TIRE01 -l cluster" --output "setup7.munge.TIRE01.$(date +%Y%m%d).out.txt"
#produce summary table when all the above have been munged and imputed
sbatch --time 9:00:00 --partition brc,shared --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 48G --wrap="Rscript setup7.R -t munge -l cluster" --output "setup7.munge.$(date +%Y%m%d).out.txt"
sbatch --time 10:00:00 --partition brc,shared --job-name="pre" --ntasks 1 --cpus-per-task 4 --mem 32G --wrap="Rscript setup7.R -t pre -l cluster" --output "setup7.pre.$(date +%Y%m%d).out.txt"
for chr in 22; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t pre -a MPRE_2_1.COR.ML:$chr -l cluster" --output "setup7.lfgwas.MPRE_2_1.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 21; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t pre -a MPRE_2_1.COR.ML:$chr -l cluster" --output "setup7.lfgwas.MPRE_2_1.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 20; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t pre -a MPRE_2_1.COR.ML:$chr -l cluster" --output "setup7.lfgwas.MPRE_2_1.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 15 16 17 18 19; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t pre -a MPRE_2_1.COR.ML:$chr -l cluster" --output "setup7.lfgwas.MPRE_2_1.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t pre -a MPRE_2_1.COR.ML:$chr -l cluster" --output "setup7.lfgwas.MPRE_2_1.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done

sbatch --time 12:00:00 --partition brc,shared --job-name="mvLD" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup7.R -t mvLD -l cluster" --output "setup7.mvLD.$(date +%Y%m%d).out.txt"
sbatch --time 1-12:00:00 --partition brc,shared --job-name="gsem" --ntasks 1 --cpus-per-task 4 --mem 32G --wrap="Rscript setup7.R -t cfa -a 1 -l cluster" --output "setup7.cfa.$(date +%Y%m%d).out.txt"
#on activating implicit multithreading
#https://hpc.nih.gov/apps/R.html#threading
#for chr in 1; do sbatch --time 4-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="OMP_NUM_THREADS=3; Rscript setup7.R -t lfgwas -a M20_6_198.COR.ML:$chr:1 -l cluster" --output "setup7.lfgwas.M20_6_198.COR.ML.F1.chr$chr.$(date +%Y%m%d).out.txt"; done
#test
#for chr in 22; do sbatch --time 1:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t lfgwas -a M20_6_198.COR.ML:$chr -l cluster" --output "setup7.lfgwas.M20_6_198.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done

for chr in 20 21 22; do sbatch --time 6-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t lfgwas -a M20_6_172.COR.ML:$chr -l cluster" --output "setup7.lfgwas.M20_6_172.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; do sbatch --time 6-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t lfgwas -a M20_6_172.COR.ML:$chr -l cluster" --output "setup7.lfgwas.M20_6_172.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for chr in 1 2 3 4; do sbatch --time 6-00:00:00 --partition brc,shared --job-name="lga_$chr" --ntasks 1 --cpus-per-task 10 --mem 32G --wrap="Rscript setup7.R -t lfgwas -a M20_6_172.COR.ML:$chr -l cluster" --output "setup7.lfgwas.M20_6_172.COR.ML.FALL.chr$chr.$(date +%Y%m%d).out.txt"; done
for fac in 1 2 3 4 5 6; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="nm$fac" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup7.R -t nmeta -a M20_6_172.COR.ML:$fac -l cluster" --output "setup7.nmeta.M20_6_172.COR.ML:$fac.$(date +%Y%m%d).out.txt"; done
for fac in 1 2 3 4 5 6; do sbatch --time 2-00:00:00 --partition brc,shared --job-name="nm$fac" --ntasks 1 --cpus-per-task 3 --mem 16G --wrap="Rscript setup7.R -t nmeta -a M20_6_30.ORT.ML:$fac -l cluster" --output "setup7.nmeta.M20_6_30.ORT.ML:$fac.$(date +%Y%m%d).out.txt"; done

sbatch --time 2-00:00:00 --partition brc,shared --job-name="mvLD2" --ntasks 1 --cpus-per-task 4 --mem 64G --wrap="Rscript setup7.R -t mvLD2 -l cluster" --output "setup7.mvLD2.$(date +%Y%m%d).out.txt"

#MAGMA workshop:
#https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/instructions.txt

while read code; do
  magma --annotate window=35,10 --snp-loc ../data/gwas_sumstats/export/"$code".SNPLOC --gene-loc ../data/gene_mapping/NCBI37.3.gene.loc --out setup7.annotate."$code";
  magma --bfile ../data/reference_panel/1KG_eur.plink/g1000_eur --gene-annot setup7.annotate."$code".genes.annot --pval ../data/gwas_sumstats/export/"$code".PVAL ncol=N --out setup7.geneanalysis."$code";
done < sumstats.sel.set1.code.txt

for code in SCHI04 SUBJ01 TIRE01; do
  magma --annotate window=35,10 --snp-loc ../data/gwas_sumstats/export/"$code".SNPLOC --gene-loc ../data/gene_mapping/NCBI37.3.gene.loc --out setup7.annotate."$code";
  magma --bfile ../data/reference_panel/1KG_eur.plink/g1000_eur --gene-annot setup7.annotate."$code".genes.annot --pval ../data/gwas_sumstats/export/"$code".PVAL ncol=N --out setup7.geneanalysis."$code";
done

for code in GSEM.F4 GSEM.F5 GSEM.F6; do
  magma --annotate window=35,10 --snp-loc ../data/gwas_sumstats/export/"$code".SNPLOC --gene-loc ../data/gene_mapping/NCBI37.3.gene.loc --out setup7.annotate."$code";
  magma --bfile ../data/reference_panel/1KG_eur.plink/g1000_eur --gene-annot setup7.annotate."$code".genes.annot --pval ../data/gwas_sumstats/export/"$code".PVAL ncol=N --out setup7.geneanalysis."$code";
done

for code in HEAL01 INCO03 INSO02 NEUR02 PTSD04 RISK02 RISK03 SCHI04 SUBJ01 TIRE01 GSEM.F1 GSEM.F2 GSEM.F3 GSEM.F4 GSEM.F5 GSEM.F6; do
  magma --annotate window=35,10 --snp-loc ../data/gwas_sumstats/export/"$code".SNPLOC --gene-loc ../data/gene_mapping/NCBI37.3.gene.loc --out setup7.annotate."$code";
  magma --bfile ../data/reference_panel/1KG_eur.plink/g1000_eur --gene-annot setup7.annotate."$code".genes.annot --pval ../data/gwas_sumstats/export/"$code".PVAL ncol=N --out setup7.geneanalysis."$code";
done

while read code; do
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --gene-covar geneset.gtex8_rnaseq_median.magma --out setup7.geneset.gtex8_rnaseq_median.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --gene-covar geneset.brainspan_rnaseq_mean.magma --out setup7.geneset.brainspan_rnaseq_mean.$code";
  #sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.v7.5.1.magma --out setup7.geneset.msigdb.v7.5.1.$code";
done < sumstats.sel.set1.code.txt

while read code; do
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.h.all.v7.5.1.magma --out setup7.geneset.msigdb.h.all.v7.5.1.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c3.all.v7.5.1.magma --out setup7.geneset.msigdb.c3.all.v7.5.1.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c5.all.v7.5.1.magma --out setup7.geneset.msigdb.c5.all.v7.5.1.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c7.all.v7.5.1.magma --out setup7.geneset.msigdb.c7.all.v7.5.1.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c8.all.v7.5.1.magma --out setup7.geneset.msigdb.c8.all.v7.5.1.$code";
done < sumstats.sel.set1.code.txt

while read code; do
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c2.cp.biocarta.v7.5.1.magma --out setup7.geneset.msigdb.c2.cp.biocarta.v7.5.1.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c2.cp.kegg.v7.5.1.magma --out setup7.geneset.msigdb.c2.cp.kegg.v7.5.1.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c3.mir.mirdb.v7.5.1.magma --out setup7.geneset.msigdb.c3.mir.mirdb.v7.5.1.$code";
  sbatch --time 20:00 --partition brc,shared --job-name="m$code" --ntasks 1 --cpus-per-task 2 --mem 4G --wrap="magma --gene-results setup7.geneanalysis.$code.genes.raw --set-annot geneset.msigdb.c3.tft.gtrd.v7.5.1.magma --out setup7.geneset.msigdb.c3.tft.gtrd.v7.5.1.$code";
done < sumstats.sel.set1.code.txt

#awk 'echo $0;' < sumstats.sel.set1.code.txt


#setup8
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a ADHD05 -l cluster" --output "setup8.munge.ADHD05.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a ALCD03 -l cluster" --output "setup8.munge.ALCD03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a ANOR02 -l cluster" --output "setup8.munge.ANOR02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a ANXI03 -l cluster" --output "setup8.munge.ANXI03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a ANXI04 -l cluster" --output "setup8.munge.ANXI04.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a AUTI07 -l cluster" --output "setup8.munge.AUTI07.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a BIPO02 -l cluster" --output "setup8.munge.BIPO02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a BODY11 -l cluster" --output "setup8.munge.BODY11.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a BVOL01 -l cluster" --output "setup8.munge.BVOL01.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a DEPR05 -l cluster" --output "setup8.munge.DEPR05.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a DEPR09 -l cluster" --output "setup8.munge.DEPR09.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a EDUC03 -l cluster" --output "setup8.munge.EDUC03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a HEAL01 -l cluster" --output "setup8.munge.HEAL01.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a INCO03 -l cluster" --output "setup8.munge.INCO03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a INSO02 -l cluster" --output "setup8.munge.INSO02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a INSO03 -l cluster" --output "setup8.munge.INSO03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a NEUR02 -l cluster" --output "setup8.munge.NEUR02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a PTSD04 -l cluster" --output "setup8.munge.PTSD04.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a RISK02 -l cluster" --output "setup8.munge.RISK02.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a RISK03 -l cluster" --output "setup8.munge.RISK03.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a SCHI04 -l cluster" --output "setup8.munge.SCHI04.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a SUBJ01 -l cluster" --output "setup8.munge.SUBJ01.$(date +%Y%m%d).out.txt"
sbatch --time 24:00:00 --partition cpu --job-name="smunge" --ntasks 1 --cpus-per-task 4 --mem 10G --wrap="Rscript setup8.R -t munge -a TIRE01 -l cluster" --output "setup8.munge.TIRE01.$(date +%Y%m%d).out.txt"

sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a ADHD05 -l cluster" --output "setup8.ldimp.ADHD05.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a ALCD03 -l cluster" --output "setup8.ldimp.ALCD03.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a ANOR02 -l cluster" --output "setup8.ldimp.ANOR02.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a ANXI03 -l cluster" --output "setup8.ldimp.ANXI03.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a ANXI04 -l cluster" --output "setup8.ldimp.ANXI04.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a AUTI07 -l cluster" --output "setup8.ldimp.AUTI07.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a BIPO02 -l cluster" --output "setup8.ldimp.BIPO02.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a BODY11 -l cluster" --output "setup8.ldimp.BODY11.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a BVOL01 -l cluster" --output "setup8.ldimp.BVOL01.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a DEPR05 -l cluster" --output "setup8.ldimp.DEPR05.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a DEPR09 -l cluster" --output "setup8.ldimp.DEPR09.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a EDUC03 -l cluster" --output "setup8.ldimp.EDUC03.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a HEAL01 -l cluster" --output "setup8.ldimp.HEAL01.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a INCO03 -l cluster" --output "setup8.ldimp.INCO03.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a INSO02 -l cluster" --output "setup8.ldimp.INSO02.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a INSO03 -l cluster" --output "setup8.ldimp.INSO03.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a NEUR02 -l cluster" --output "setup8.ldimp.NEUR02.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a PTSD04 -l cluster" --output "setup8.ldimp.PTSD04.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a RISK02 -l cluster" --output "setup8.ldimp.RISK02.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a RISK03 -l cluster" --output "setup8.ldimp.RISK03.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a SCHI04 -l cluster" --output "setup8.ldimp.SCHI04.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a SUBJ01 -l cluster" --output "setup8.ldimp.SUBJ01.$(date +%Y%m%d).out.txt"
sbatch --time 2-00:00:00 --partition cpu --job-name="ld-imp" --ntasks 1 --cpus-per-task 5 --mem 20G --wrap="Rscript setup8.R -t impute -a TIRE01 -l cluster" --output "setup8.ldimp.TIRE01.$(date +%Y%m%d).out.txt"

#ssimp imputation
sbatch --time 2-00:00:00 --partition cpu --job-name="ssimpformat" --ntasks 1 --cpus-per-task 5 --mem 8G --wrap="Rscript setup8.R -t ssimpformat -l cluster" --output "setup8.ssimpformat.$(date +%Y%m%d).out.txt"


/scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --download.build.db
#test
sbatch --time 1:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 3 --mem 8G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge/ADHD05.gz --ref /scratch/prj/gwas_sumstats/ssimp_software-master/ref/small.vcf.sample.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/ADHD05.txt; date;" --output "SSIMP_evaluation.ADHD05.$(date +%Y%m%d).out.txt"


for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/SCHI04.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/SCHI04.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.SCHI04.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/BODY11.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/BODY11.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.BODY11.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/NEUR02.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/NEUR02.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.NEUR02.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/DEPR05.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/DEPR05.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.DEPR05.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/DEPR09.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/DEPR09.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.DEPR09.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/BIPO02.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/BIPO02.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.BIPO02.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/ANXI03.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/ANXI03.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.ANXI03.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/INSO02.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/INSO02.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.INSO02.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/EDUC03.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/EDUC03.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.EDUC03.$chr.$(date +%Y%m%d).out.txt"; done

for chr in 1 2 3 4 5 6 7 8 9 10 11 "12-14" "15-17" "18-23"; do sbatch --time 2-00:00:00 --partition cpu --job-name="ssimp" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="export LC_ALL=C; date; /scratch/prj/gwas_sumstats/ssimp_software-master/ssimp --gwas ../data/gwas_sumstats/munged_1kg_eur_supermunge.4ssimp/BVOL01.gz --ref ../data/reference_panel/hc1kgp3.b38.plink/hc1kgp3.b38.eur.jz2022.vcf.gz --out ../data/gwas_sumstats/imputed_1kg_eur_ssimp/BVOL01.$chr --impute.maf 0.001 --tag.maf 0.01 --impute.range $chr; date;" --output "setup8.ssimp.BVOL01.$chr.$(date +%Y%m%d).out.txt"; done

sbatch --time 2-00:00:00 --partition cpu --job-name="merge" --ntasks 1 --cpus-per-task 5 --mem 32G --wrap="Rscript setup8.R -t merge -l cluster" --output "setup8.merge.$(date +%Y%m%d).out.txt"

sbatch --time 2-00:00:00 --partition cpu --job-name="mvLD" --ntasks 1 --cpus-per-task 5 --mem 64G --wrap="Rscript setup8.R -t mvLD -l cluster" --output "setup8.mvLD.$(date +%Y%m%d).out.txt"
#for running original ldsc
sbatch --time 2-00:00:00 --partition cpu --job-name="mvLD" --ntasks 1 --cpus-per-task 5 --mem 64G --wrap="source /users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc-venv/bin/activate; Rscript setup8.R -t mvLD -l cluster;" --output "setup8.mvLD.$(date +%Y%m%d).out.txt"
Rscript setup8.R -t mvLD -l cluster > "setup8.mvLD.$(date +%Y%m%d).out.txt" #if needed to run on interactive node beacuase of python venv activation issues



#fix duplicates in PLINK refpanel for Thibault
plink2 --bfile g1000_eur --rm-dup --make-pgen --out g1000_eur.iddup



#!/bin/bash
#Load the R module
module add apps/R/3.6.0
#module add apps/R/3.6.3

#install R packages for user
#R -e 'install.packages(devtools)'
#R -e 'remove.packages("GenomicSEM")'
#R -e 'devtools::install_github(c("MichelNivard/GenomicSEM"))'
#R -e 'install.packages("skimr")'
#R -e 'install.packages("psych")'


#sbatch --time 02:59:00 --partition shared --job-name="GSEM_CFA" --mem 8G --wrap="module add apps/R/3.6.0 && Rscript analysis_plan_setup2_20200601.R" --output "analysis_plan_setup2_20200601.out" --error "analysis_plan_setup2_20200601.err" &

#sbatch --time 00:59:00 --partition shared --job-name="GSEM_GWAS_PREP" --mem 16G --wrap="module add apps/R/3.6.0 && Rscript analysis_plan_setup2_20200611.R" --output "analysis_plan_setup2_20200611.out" --error "analysis_plan_setup2_20200611.err" &

#sbatch --time 00:59:00 --partition shared --job-name="GSEM_GWAS_PREP" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 6G --wrap="module add apps/R/3.6.0 && Rscript analysis_plan_setup2_20200611.R" --output "analysis_plan_setup2_20200611.out" --error "analysis_plan_setup2_20200611.err" &

#sbatch --time 23:59:00 --partition shared --job-name="GSEMLGWAS" --ntasks 1 --cpus-per-task 6 --mem-per-cpu 6G --wrap="module add apps/R/3.6.0 && Rscript analysis_plan_setup2_20200611.R" --output "analysis_plan_setup2_20200611.out" --error "analysis_plan_setup2_20200611.err" &

#sbatch --time 23:59:00 --partition shared --job-name="GSEMLGWAS" --ntasks 1 --cpus-per-task 6 --mem-per-cpu 6G --wrap="module add apps/R/3.6.0 && Rscript analysis_plan_setup3_20200630.R" --output "analysis_plan_setup3_20200630.out" --error "analysis_plan_setup3_20200630.err" &

sbatch --time 23:59:00 --partition shared --job-name="GSEMHDL" --ntasks 1 --cpus-per-task 6 --mem-per-cpu 6G --wrap="module add apps/R/3.6.0 && Rscript setup1.R" --output "setup1.out" --error "setup1.err" &

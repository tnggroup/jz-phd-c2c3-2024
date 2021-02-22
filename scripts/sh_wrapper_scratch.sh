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

#sbatch --time 23:59:00 --partition shared --job-name="GSEMHDL" --ntasks 1 --cpus-per-task 6 --mem-per-cpu 6G --wrap="module add apps/R/3.6.0 && Rscript setup1.R" --output "setup1.out" --error "setup1.err" &

#sbatch --time 23:59:00 --partition brc --job-name="GSEMHDL" --ntasks 1 --cpus-per-task 6 --mem-per-cpu 6G --wrap="module add apps/R/3.6.0 && Rscript setup1.R" --output "setup1.out" --error "setup1.err"

# sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt6" --ntasks 1 --cpus-per-task 5 --mem-per-cpu 5G --wrap="module add apps/R/3.6.0 && Rscript setup1.R -t 6" --output "setup1_t6.out" --error "setup1_t6.err"
# sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt5" --ntasks 1 --cpus-per-task 5 --mem-per-cpu 5G --wrap="module add apps/R/3.6.0 && Rscript setup1.R -t 5" --output "setup1_t5.out" --error "setup1_t5.err"
# sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt4" --ntasks 1 --cpus-per-task 5 --mem-per-cpu 5G --wrap="module add apps/R/3.6.0 && Rscript setup1.R -t 4" --output "setup1_t4.out" --error "setup1_t4.err"
# sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt3" --ntasks 1 --cpus-per-task 5 --mem-per-cpu 5G --wrap="module add apps/R/3.6.0 && Rscript setup1.R -t 3" --output "setup1_t3.out" --error "setup1_t3.err"
# sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt2" --ntasks 1 --cpus-per-task 5 --mem-per-cpu 5G --wrap="module add apps/R/3.6.0 && Rscript setup1.R -t 2" --output "setup1_t2.out" --error "setup1_t2.err"

sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt2" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup1.R -t 2 -l cluster" --output "setup1.t2.out" --error "setup1.t2.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt3" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup1.R -t 3 -l cluster" --output "setup1.t3.out" --error "setup1.t3.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt4" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup1.R -t 4 -l cluster" --output "setup1.t4.out" --error "setup1.t4.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt5" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup1.R -t 5 -l cluster" --output "setup1.t5.out" --error "setup1.t5.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvHDLt6" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup1.R -t 6 -l cluster" --output "setup1.t6.out" --error "setup1.t6.err"

sbatch --time 23:59:00 --partition brc,shared --job-name="GSEMGWAS" --ntasks 1 --cpus-per-task 6 --mem-per-cpu 8G --wrap="Rscript setup1.R -l cluster" --output "setup1_$(date +%Y%m%d).out.txt" --error "setup1_$(date +%Y%m%d).err.txt"

sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.mvLDSC" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup2.R -t mvLD.mvLDSC -l cluster" --output "setup2.mvLD.mvLDSC.out" --error "setup2.mvLD.mvLDSC.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.HDL.piecewise" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup2.R -t mvLD.HDL.piecewise -l cluster" --output "setup2.mvLD.HDL.piecewise.out" --error "setup2.mvLD.HDL.piecewise.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.HDL.jackknife" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup2.R -t mvLD.HDL.jackknife -l cluster" --output "setup2.mvLD.HDL.jackknife.out" --error "setup2.mvLD.HDL.jackknife.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.origHDL" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup2.R -t mvLD.origHDL -l cluster" --output "setup2.mvLD.origHDL.out" --error "setup2.mvLD.origHDL.err"
sbatch --time 23:59:00 --partition brc,shared --job-name="mvLD.origHDL.liabilityScale" --ntasks 1 --cpus-per-task 4 --mem-per-cpu 5G --wrap="Rscript setup2.R -t mvLD.origHDL.liabilityScale -l cluster" --output "setup2.mvLD.origHDL.liabilityScale.out" --error "setup2.mvLD.origHDL.liabilityScale.err"





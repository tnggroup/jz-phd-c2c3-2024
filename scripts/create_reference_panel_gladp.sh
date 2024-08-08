# Commands for creating the in-sample GLAD+ reference panel and LD-scores

#run in the folder where the keep files are located (i.e. personality project /work)
sbatch --time 12:00:00 --partition cpu --job-name="keep" --ntasks 1 --cpus-per-task 4 --mem 20G --wrap="plink --bfile /users/k19049801/project/JZ_GED_PHD_C1/data/geno/GLADv3_EDGIv1_NBRv2/imputed/bfiles/GLAD_EDGI_NBR  --keep keep.all.txt --freq --list-duplicate-vars --make-bed --out GLAD_EDGI_NBR.keep" --output "GLAD_EDGI_NBR.keep.$(date +%Y%m%d).out.txt"

#create a full b38 genetic recombination map from b37 1KG (and HM2 X)
sbatch --time 1-00:00:00 --partition brc,shared --job-name="cmorgan" --ntasks 1 --cpus-per-task 4 --mem 32G --wrap="Rscript /users/k19049801/project/JZ_GED_PHD_C1/scripts/combine_genetic_recombination_map.R" --output "combine_genetic_recombination_map.$(date +%Y%m%d).out.txt"

#set new names and genomic position in cM
sbatch --time 12:00:00 --partition cpu --job-name="ref.CM" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="plink --bfile GLAD_EDGI_NBR.keep --cm-map /users/k19049801/project/JZ_GED_PHD_C1/data/genetic_recombination_mapping/genetic-map-chr-bp-rr-cm.1KGP3.b38.jz2022.SHAPEIT.chr/genetic_map_chr@_combined_b38.jz2022.txt --make-bed --out GLAD_EDGI_NBR.keep.CM" --output "GLAD_EDGI_NBR.keep.CM.plink.$(date +%Y%m%d).out.txt"
sbatch --time 12:00:00 --partition cpu --job-name="ref.CM" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="plink --bfile GLAD_EDGI_NBR.keep.CM --cm-map /users/k19049801/project/JZ_GED_PHD_C1/data/genetic_recombination_mapping/genetic-map-chr-bp-rr-cm.1KGP3.b38.jz2022.SHAPEIT.chr/genetic_map_chr23_combined_b38.jz2022.txt 23 --make-bed --out GLAD_EDGI_NBR.keep.CM23" --output "GLAD_EDGI_NBR.keep.CM23.plink.$(date +%Y%m%d).out.txt"

#LD scores
#sbatch --time 1-00:00:00 --partition cpu --job-name="ldscores" --ntasks 1 --cpus-per-task 4 --mem 70G --wrap=". /users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc-venv3.11/bin/activate; python /users/k19049801/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc/ldsc.py --bfile GLAD_EDGI_NBR.keep --l2 --ld-wind-cm 1 --out GLAD_EDGI_NBR.1cm.250blocks --n-blocks 250;" --output "GLAD_EDGI_NBR.1cm.250blocks.ldsc.$(date +%Y%m%d).out.txt" #LDSC was out of order, also we did not have cM values
#alt
sbatch --time 2-00:00:00 --partition cpu --job-name="ldscores2" --ntasks 1 --cpus-per-task 8 --mem 70G --wrap="gcta64 --bfile GLAD_EDGI_NBR.keep --ld-score --ld-wind 1000 --ld-score-adj --autosome --out GLAD_EDGI_NBR.keep.gcta.ld --thread-num 8" --output "GLAD_EDGI_NBR.1Mb.gcta.$(date +%Y%m%d).out.txt"
sbatch --time 1-00:00:00 --partition cpu --job-name="ldscores2" --ntasks 1 --cpus-per-task 6 --mem 50G --wrap="gcta64 --bfile GLAD_EDGI_NBR.keep --ld-score --ld-wind 1000 --ld-score-adj --autosome-num 23 --chr 23 --out GLAD_EDGI_NBR.keep.gcta.ld.23 --thread-num 6" --output "GLAD_EDGI_NBR.1Mb.gcta.$(date +%Y%m%d).out.txt"
#the above may be improved into a single command by using --autosome-num 23 in the first command (untested)
#In R console
# library("shru")
# library("data.table")
# ld<-readFile("GLAD_EDGI_NBR.keep.gcta.ld.score.ld")
# ld23<-readFile("GLAD_EDGI_NBR.keep.gcta.ld.23.score.ld")
# ld<-rbind(ld,ld23)
# rm(ld23)
# colnames(ld)<-c("SNP","CHR","BP","MAF","MEAN_RSQ","SNP_NUM","MAX_RSQ","L2")
# ldOut<-ld[,c("CHR","SNP","BP","L2","MAF","SNP_NUM","MEAN_RSQ","MAX_RSQ")]
# setkeyv(ldOut,cols=c("CHR","BP","SNP","MAF"))
# setorder(ldOut,CHR,BP,-MAF,SNP)
# fwrite(ldOut,file="GLAD_EDGI_NBR.keep.gcta.l2.ldscore.gz",append = F,quote = F,sep = "\t",col.names = T)

#this was not actually run, but rather run interactively
sbatch --time 1:00:00 --partition cpu --job-name="varlist" --ntasks 1 --cpus-per-task 6 --mem 80G --wrap="Rscript /users/k19049801/project/JZ_GED_PHD_C1/scripts/create_variant_list_gladp.R" --output "create_variant_list_gladp.$(date +%Y%m%d).out.txt"

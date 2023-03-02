# Commands for creating the hc1kgp3 reference panel and LD-scores

#new high coverage 1KG reference panel
sbatch --time 12:00:00 --partition brc,shared --job-name="wget" --ntasks 1 --cpus-per-task 4 --mem 8G --wrap="wget -r –level=0 -E –ignore-length -x -k -p -erobots=off -np -N http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV" --output "wget.hc1kg.$(date +%Y%m%d).out.txt"
vcf-concat 1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | gzip > 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.vcf.gz #change this to use bgzip -c instead, for getting the correct zip format to work with tabix.
#alternatively with bcftools(not tested) to get the correct BGZF zip format:
#bcftools concat -o 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
bcftools view -Oz -o 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_2.vcf.gz 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.vcf.gz #used to get the correct zip-format. rename second version to original file-name and continue. discard the first file. this step is supposedly not needed if you used bcftools for the concatenation above.
tabix -p vcf 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
sbatch --time 12:00:00 --partition brc,shared --job-name="refpan" --ntasks 1 --cpus-per-task 4 --mem 60G --wrap="plink --vcf ../hc1kgp3.b38.vcf/1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --out 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.plink.$(date +%Y%m%d).out.txt"
sbatch --time 12:00:00 --partition brc,shared --job-name="refpan" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="plink --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel --freq --list-duplicate-vars --make-bed --out 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.plink.$(date +%Y%m%d).out.txt" #this extra step makes the chromosome encoding numeric and align the .bim and .frq variants.


#run combine_genetic_recombination_map.R - does the liftover!!!!
#otherwise
#Liftover needs the UCSC tools from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
#UCSC LiftOver
#https://genome.sph.umich.edu/wiki/LiftOver
#liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed

#create a full b38 genetic recombination map from b37 1KG (and HM2 X)
sbatch --time 1-00:00:00 --partition brc,shared --job-name="cmorgan" --ntasks 1 --cpus-per-task 4 --mem 32G --wrap="Rscript ../../../JZ_GED_PHD_C1/scripts/combine_genetic_recombination_map.R" --output "combine_genetic_recombination_map.$(date +%Y%m%d).out.txt"

#dbSNP
sbatch --time 12:00:00 --partition brc,shared --job-name="wget" --ntasks 1 --cpus-per-task 4 --mem 8G --wrap="wget -r –level=0 -E –ignore-length -x -k -p -erobots=off -np -N https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz" --output "wget.dbsnp.human_9606_b151_GRCh38p7.$(date +%Y%m%d).out.txt"
wget -r –level=0 -E –ignore-length -x -k -p -erobots=off -np -N https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz.tbi
#there was potentially some indexing done here also
#gunzip the dbSNP file before running this
sbatch --time 2-00:00:00 --partition cpu --job-name="rsids" --ntasks 1 --cpus-per-task 7 --mem 160G --wrap="Rscript /users/k19049801/project/JZ_GED_PHD_C1/scripts/edit_reference_panel.R" --output "edit_reference_panel.$(date +%Y%m%d).out.txt"


#set new names and genomic position in cM
sbatch --time 12:00:00 --partition brc,shared --job-name="ref.CM" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="plink --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq --cm-map ../../genetic_recombination_mapping/genetic-map-chr-bp-rr-cm.1KGP3.b38.jz2022.SHAPEIT.chr/genetic_map_chr@_combined_b38.jz2022.txt --make-bed --out 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.TCM" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.TCM.plink.$(date +%Y%m%d).out.txt"
sbatch --time 12:00:00 --partition brc,shared --job-name="ref.CM" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="plink --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.TCM --cm-map ../../genetic_recombination_mapping/genetic-map-chr-bp-rr-cm.1KGP3.b38.jz2022.SHAPEIT.chr/genetic_map_chr23_combined_b38.jz2022.txt 23 --make-bed --out 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.plink.$(date +%Y%m%d).out.txt"
#sbatch --time 12:00:00 --partition brc,shared --job-name="ref.qc" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="plink --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.rs.CM23 --mac 3 --make-bed --out 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.rs.CM23.qc" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.rs.CM23.qc.plink.$(date +%Y%m%d).out.txt"
plink2 --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM --rm-dup force-first --make-pgen --out tmp.iddup #check id-duplicates with plink2

#Extract European subset - non-QC'd!!
sbatch --time 12:00:00 --partition brc,shared --job-name="ref.eur" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="plink --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM --keep ../1KG_eur.plink/g1000_eur.fam --freq --make-bed --out 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.plink.$(date +%Y%m%d).out.txt"
#produce vcf for subset
sbatch --time 12:00:00 --partition brc,shared --job-name="ref.eur" --ntasks 1 --cpus-per-task 4 --mem 50G --wrap="plink --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur --recode vcf-iid --out 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.vcf" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.plink.vcf.$(date +%Y%m%d).out.txt"
bcftools view -Oz9 -o hc1kgp3.b38.eur.jz2022.vcf.gz 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.vcf.vcf
tabix -p vcf hc1kgp3.b38.eur.jz2022.vcf.gz


#create LD score libraries
#source ~/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc-venv/bin/activate
#run in the reference panel folder
#python ~/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc/ldsc.py --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.rs.CM23.eur --l2 --ld-wind-cm 1 --out hc1kgp3.eur --n-blocks 2000 > ldsc.hc1kgp3.eur.txt
sbatch --time 1-00:00:00 --partition brc,shared --job-name="ldsc" --ntasks 1 --cpus-per-task 4 --mem 80G --wrap="source ~/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc-venv/bin/activate; python ~/project/JZ_GED_PHD_ADMIN_GENERAL/software/ldsc/ldsc.py --bfile 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur --l2 --ld-wind-cm 1 --out hc1kgp3.b38.eur.1cm.400blocks --n-blocks 400;" --output "1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel.frq.CM.eur.1cm.400blocks.ldsc.$(date +%Y%m%d).out.txt"
sbatch --time 2:00:00 --partition brc,shared --job-name="varlist" --ntasks 1 --cpus-per-task 4 --mem 120G --wrap="Rscript ../../../../JZ_GED_PHD_C1/scripts/create_variant_list.R" --output "create_variant_list.$(date +%Y%m%d).out.txt"






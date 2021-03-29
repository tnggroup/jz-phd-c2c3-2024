#!/usr/bin/awk -f

BEGIN{
folderpath_reference="../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM";
folderpath_ld="../data/reference.panel.1KG_Phase3.CLEANED.EUR.cM";
folderpath_parent_zscore="../data/gwas_sumstats/munged";
folderpath_parent_output="../data/gwas_sumstats/imputed";
filepath_python_venv="../python-venv/bin/activate"

}

#Adapt this wether there are headers or not
FNR>1 && NF>3{
	#unquote
	code=$1;
	gsub(/"/,"",code);
  	cmd = "mkdir -p " folderpath_parent_output "/" code; 
	#print cmd;
	system(cmd);
	#split("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23",chromosomes);
	cmd="source " filepath_python_venv "; for chr in {1..23}; do raiss --chrom \\$chr --gwas " code " --ref-folder " folderpath_reference " --ld-folder " folderpath_ld " --zscore-folder "folderpath_parent_zscore"/"code".chr --output-folder "folderpath_parent_output"/"code" --l2-regularization 0.01 --eigen-threshold 0.05 --R2-threshold 0.3; done";
  	#print cmd;
	cmd2="sbatch --time 23:59:00 --partition brc,shared --job-name=\"imp" code "\" --ntasks 1 --cpus-per-task 3 --mem 20G --wrap=\"" cmd "\" --output impute."code".$(date +%Y%m%d).out.txt --error impute."code".$(date +%Y%m%d).err.txt"
	#print cmd2;
	system(cmd2);

  	
}

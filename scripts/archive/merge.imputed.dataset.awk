#set the code, prefix and suffix variables for awk with the -v flag when running awk

BEGIN{
	#print "Start";
	#nfile=0;
	folderpath_parent="../data/gwas_sumstats/imputed/";
	chromosomes_len=split("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22",chromosomes);
	#print "Chromosomes length =" chromosomes_len;
	for(i=1;i<=chromosomes_len;i++)
	{
		ARGV[i]=folderpath_parent code "/" prefix chromosomes[i] suffix;
		#print ARGV[i];
		ARGC++;
	}
}

NR==1 && FNR==1{
	print $0;
}

#FNR==1{
	#nfile++;
	#print ARGV[nfile];
#}

FNR>1{
	print $0;
}


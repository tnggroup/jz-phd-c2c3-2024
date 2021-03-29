#use awk -v chr=[YOUR VALUE] to set which chromosome to extract
BEGIN{columnNumbers[0]=0;}
NF>3{
	if(CHR){chr=CHR;}

	if(length(columnNumbers)<3)
	{
		for(n=1;n<=NF;n++)
		{
			columnNumbers[$(n)]=n;
		}
		print $0;
	}
	else 
	{
		if(chr==$(columnNumbers["CHR"])){print $0;}
	}
}

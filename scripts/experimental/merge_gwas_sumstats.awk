#!/usr/bin/awk -f

BEGIN{
  tagSNPrsID="SNP";
  columnNames1[1]=0;
  columnNames2[1]=0;
  columnIndices1[1]=-1;
  columnIndices2[1]=-1;
}


NR==FNR && NF>3{
  if(length(columnNames1)<2)
  {
    for(i=0; i<NF; i++)
    {
      columnNames1[i]=toupper($(i+1));
      columnIndices1[toupper($(i+1))]=i;
      if(!iSNPrsID1&&columnNames1[i]==tagSNPrsID){iSNPrsID1=i;}
    }
    #print "Found header 1:", columnNames1;
    #for(i=0; i<length(columnNames1); i++)
    #{
     # 	printf "%s",columnIndices1[columnNames1[i]],",";
    #}
    #print "";
  
  }
  else
  {
    #print "Stored ";
    for(i=0; i<NF; i++)
    {
      d1[$(iSNPrsID1+1)][i]=$(i+1);
      #print "",$(i+1);
    }
  }
  
} 

NR!=FNR && NF>3{
  
  if(length(columnNames2)<2)
  {
    for(i=0; i<NF; i++)
    {
      columnNames2[i]=toupper($(i+1));
      columnIndices2[toupper($(i+1))]=i;  
      if(!iSNPrsID2&&columnNames2[i]==tagSNPrsID){iSNPrsID2=i;}
    }

    for(i=0; i<length(columnNames1)+length(columnNames2);i++)
    {
	#Renaming columns
	if(i<length(columnNames1)){c=columnNames1[i];}
	else{c=columnNames2[i-length(columnNames1)];}
        
        cn=0;
	if(c=="SNP"){cn="rsID";}
	if(c=="ORIGBP"){cn="pos";}
	
	if(cn)
	{
		#print cn;	
		if(i<length(columnNames1)){columnNames1[i]=cn;}
		else{columnNames2[i-length(columnNames1)]=cn;}
	}
    }    

    for(i=0; i<length(columnNames1); i++)
    {
	printed[columnNames1[i]]=(i+1);
      	printf "%s%s",columnNames1[i],OFS;
    }  

    for(i=0;i<length(columnNames2);i++)
    {
	if(!printed[columnNames2[i]])
	{
      		printf "%s%s",columnNames2[i],OFS;
	}
    }
    print "";

  }
  else
  {
    if(length(d1[$(iSNPrsID2+1)])>2)
    {
      #print "Found ", $(iSNPrsID2+1);
      #printf "%s:",$(iSNPrsID2+1);
      for(i=0; i<length(columnNames1); i++)
      {
	#print columnNames1[i];
      	printf "%s%s",d1[$(iSNPrsID2+1)][i],OFS;
      }  

      for(i=0;i<length(columnNames2);i++)
      {
	if(!printed[columnNames2[i]])
	{
      		printf "%s%s",$(i+1),OFS;
	}
      }
      print "";
    }
  }
}

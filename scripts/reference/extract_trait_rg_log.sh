#!/bin/bash -l

###Extract rG for trait vs other traits in correlation logs
SCRIPTPATH=${BASH_SOURCE}

SCRIPTNAME=extract_trait_rg
TRAITCODE=ANXI03

#either normal or parallel
RUNTYPE=parallel


SCRIPTDIRPATH=/scratch/users/k19049801/project/JZ_GED_RESEARCHPROJECT/scripts
MUNGEDNOMHCDIRPATH=/mnt/lustre/groups/ukbiobank/sumstats/munged_noMHC
#CORRELATIONSPATH=/mnt/lustre/groups/ukbiobank/sumstats/correlations

LDREFERENCEPATH=/scratch/users/k19049801/project/JZ_GED_RESEARCHPROJECT/data_raw/eur_w_ld_chr
WORKDIRPATH="/scratch/users/k19049801/project/JZ_GED_RESEARCHPROJECT/working_directory/$SCRIPTNAME"

PYTHONPATH=/scratch/groups/ukbiobank/Edinburgh_Data/Software/anaconda2/bin/python2.7
LDSCPATH=/scratch/groups/ukbiobank/usr/abi/ldsc/ldsc.py

#make work dir in work path
mkdir -p "$WORKDIRPATH" && cd $WORKDIRPATH

TASK=$1
RGLIST=$2
IDSTRING=$3

echo "RUNNING: $SCRIPTPATH"
echo "TASK IS $TASK"

if [ "$TASK" == "START" ];
then
#read the file names in mungednomhcpath to a file
ls "$MUNGEDNOMHCDIRPATH" > "$WORKDIRPATH"/munged_noMHC_files.txt
#remove '_noMHC.sumstats.gz' from the list
#sed 's/_noMHC.sumstats.gz//g' $WORKPATH/munged_noMHC.txt > $WORKPATH/munged_noMHC.txt
gawk 'match($1, /^(.*)_noMHC\.sumstats\.gz$/, matches) {print matches[1]}' "$WORKDIRPATH"/munged_noMHC_files.txt > $WORKDIRPATH/munged_noMHC_list.txt


while read p; do
  echo "$p"
  
  if [ ! -z "$p" ];
  then
    #set -x
    TRAITCODE_REL=$p
    RGLIST="${MUNGEDNOMHCDIRPATH}/${TRAITCODE}_noMHC.sumstats.gz,${MUNGEDNOMHCDIRPATH}/${TRAITCODE_REL}_noMHC.sumstats.gz"
    #RGLIST="${MUNGEDNOMHCDIRPATH}/${TRAITCODE}_noMHC.sumstats.gz,${MUNGEDNOMHCDIRPATH}/ALCO01_noMHC.sumstats.gz,${MUNGEDNOMHCDIRPATH}/INSO01_noMHC.sumstats.gz"
    
    
    if [ "$RUNTYPE" == "parallel" ];
    then
      sbatch \
      --time 00:10:00 \
      --partition shared \
      --job-name=JZ_LDSC \
      --wrap="sh ${SCRIPTDIRPATH}/extract_trait_rg_log.sh PROC ${RGLIST} ${TRAITCODE}_${TRAITCODE_REL}" \
      --output "$WORKDIRPATH/$SCRIPTNAME.${TRAITCODE}_${TRAITCODE_REL}.out" \
      --error "$WORKDIRPATH/$SCRIPTNAME.${TRAITCODE}_${TRAITCODE_REL}.err"
      sleep 1s
    else
      #does not work?
      sh ${SCRIPTDIRPATH}/extract_trait_rg_log.sh PROC ${RGLIST} ${TRAITCODE}_${TRAITCODE_REL} > "$WORKDIRPATH/$SCRIPTNAME.${TRAITCODE}_${TRAITCODE_REL}.out"
      
    fi
    #set +x
    #--wrap="sh ${SCRIPTDIRPATH}/extract_trait_rg_log.sh PROC ${RGLIST} ${TRAITCODE}_MULTI" \
      #--output "$WORKDIRPATH/$SCRIPTNAME.${TRAITCODE}_MULTI.out" \
      #--error "$WORKDIRPATH/$SCRIPTNAME.${TRAITCODE}_MULTI.err"
    
    #--mail-user=k19049801@kcl.ac.uk
    #--test-only \
  fi
  
done < "$WORKDIRPATH"/munged_noMHC_list.txt


elif [ "$TASK" == "PROC" ];
then
  
  echo "PROC[${IDSTRING}]:$RGLIST"
  
  set -x
  
  $PYTHONPATH $LDSCPATH \
  --n-blocks 200 \
  --rg ${RGLIST} \
  --ref-ld-chr "${LDREFERENCEPATH}/" \
  --w-ld-chr "${LDREFERENCEPATH}/" \
  --out "${WORKDIRPATH}/${IDSTRING}.ldsc"
  
  #--print-delete-vals \
  set +x
  
elif [ "$TASK" == "EXTRACT" ]
then

  #touch "${WORKDIRPATH}/${SCRIPTNAME}.${TRAITCODE}_gc.txt"
  #echo -e 'code\tgc\n'>"${WORKDIRPATH}/${TRAITCODE}_gc.txt" #use -e for interpreting escape characters
  echo "Performing gawk"
  gawk 'BEGIN {FS="\t"; OFS="\t"; print "p1","p2","rg","se","z","p","h2_obs","h2_obs_se","h2_int","h2_int_se","gcov_int","gcov_int_se";} ' > "${WORKDIRPATH}/${TRAITCODE}_gc.txt"
  echo "Done gawking"
  
  while read p; do
    echo "$p"
    
    if [ ! -z "$p" ];
    then
      #set -x
      TRAITCODE_REL=$p
      
      echo "Matched GC[${TRAITCODE_REL}]:${matchedGC}"
      echo "In file ${WORKDIRPATH}/${TRAITCODE}_${TRAITCODE_REL}.ldsc.log"
      #This gawk step can probably be improved further since it captures a lot of unnecessary whitespace.
      gawk 'BEGIN {FS="\t"; OFS="\t"} { if (match($0,/gcov_int_se/) ) {found=1;} if (found==1 && match($0,/(_noMHC\.sumstats\.gz)()/)) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12} }' "${WORKDIRPATH}/${TRAITCODE}_${TRAITCODE_REL}.ldsc.log" >> "${WORKDIRPATH}/${TRAITCODE}_gc.txt"
      echo -e '\n' >> "${WORKDIRPATH}/${TRAITCODE}_gc.txt"
    fi
    
  done < "$WORKDIRPATH"/munged_noMHC_list.txt
  
  
else
  echo "Use task = START to start"
fi
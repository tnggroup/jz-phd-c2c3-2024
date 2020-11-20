#!/usr/bin/awk -f
BEGIN{
  FS="|"
}

{
  id[$1]= 1;
  data[$1,FILENAME]=$2
}

END {
  for (i in id) print trim(i)"|"trim(data[i,"scientific"])"|"trim(data[i,"genbank"])
}

function trim (x) {
  sub(/^[ \t]*/,"",x);
  sub(/[ \t]*$/,"",x);
  return x
} 

#!/usr/bin/awk -f

NR == FNR {
  printf("%s  ",$0)
  while(getline line < ARGV[2]){
    split(line,data)
    if(data[1]==$1) printf("%s  ",data[2])
  }
  close(ARGV[2])
  print ""
}

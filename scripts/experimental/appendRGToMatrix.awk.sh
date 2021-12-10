#!/usr/bin/awk -f

BEGIN 
{
  #settings
  fileFolderPath=ARGV[1]
  outputFolderPath=ARGV[2]
  #system("ls outputFolderPath > _files.txt")
  "ls outputFolderPath" | getline filesList
  print filesList
}

#Merge phenotypic data template, Johan Zvrskovec 2023

#packages
library(data.table) #this merge routine relies heavily on the data.table package (which is great)

#settings - configure these for your environment!!
folderpath.data<-normalizePath("/Users/jakz/Library/CloudStorage/OneDrive-SharedLibraries-King'sCollegeLondon/MT-TNG BioResource EDIT - ilovedata - ilovedata",mustWork = T) #path to the folder with cleaned data


#read data to be merged
folderpath.data.current<-file.path(folderpath.data,"data","latest_freeze_2021","glad","demographics")
d.dem.age <- readRDS(
  file = file.path(folderpath.data.current,"age_glad_clean.rds")
)
d.dem.sex_gender_sexuality<- readRDS(
  file = file.path(folderpath.data.current,"sex_gender_sexuality_glad_clean.rds")
)

folderpath.data.current<-file.path(folderpath.data,"data","latest_freeze_2021","glad","environmental_factors")
d.dem.ats <- readRDS(
  file = file.path(folderpath.data.current,"ats_glad_clean.rds")
)

folderpath.data.current<-file.path(folderpath.data,"data","latest_freeze_2021","coping_glad","personality")
d.dem.HEXACO <- readRDS(
  file = file.path(folderpath.data.current,"HEXACO_coping_glad_clean.rds")
)

#set data.table nature
setDT(d.dem.age)
setkeyv(d.dem.age,cols = "ID")

setDT(d.dem.sex_gender_sexuality)
setkeyv(d.dem.sex_gender_sexuality,cols = "ID")

setDT(d.dem.ats)
setkeyv(d.dem.ats,cols = "ID")

setDT(d.dem.HEXACO)
setkeyv(d.dem.HEXACO,cols = "ID")

# Merge all data
merged <- data.table(
  ID=unique(c(
    d.dem.age[,ID],
    d.dem.sex_gender_sexuality[,ID],
    d.dem.ats[,ID],
    d.dem.HEXACO[,ID]
  ))
)
setkeyv(merged,cols = "ID")
dim(merged)

merged<-merge(x = merged, y = d.dem.age, by = "ID", all = T)
merged<-merge(x = merged, y = d.dem.sex_gender_sexuality, by = "ID", all = T)
merged<-merge(x = merged, y = d.dem.ats, by = "ID", all = T)
merged<-merge(x = merged, y = d.dem.HEXACO, by = "ID", all = T)
dim(merged)

merged<-as.data.frame(merged) #remove data.table nature



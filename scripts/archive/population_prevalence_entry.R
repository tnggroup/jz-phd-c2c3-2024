#DO NOT FORGET TO SET WORKING DIRECTORY
#setwd("~/King's College London/MT-Translational Neuropsychiatric Genomics - Johan_Zvrskovec_PhD - Johan_Zvrskovec_PhD/JZ_GED_PHD_C1/workingDirectory")
library(DBI)
#install.packages('RPostgres')
#library(RPostgres)

# Connect to the database
con <- dbConnect(RPostgres::Postgres(),
                 dbname = 'phenodb', 
                 host = '10.200.105.5', 
                 port = 5432,
                 user = 'johan',
                 password = 'hej123')


popprev<-read.table("population-prevalence_20210311.csv", header=T, quote="\"", sep = ",", fill=T, blank.lines.skip=T,as.is = c(2), strip.white = T)

View(popprev)

dbWriteTable(con, "popprev20210311", popprev)


# dbListTables(con)
# dbWriteTable(con, "mtcars", mtcars)
# dbListTables(con)
# 
# dbListFields(con, "mtcars")
# dbReadTable(con, "mtcars")

# Disconnect from the database
dbDisconnect(con)

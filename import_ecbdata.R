####
library(dplyr)    # to use %>% -- loads also "last"
library(tidyverse)
library(data.table)
library(xts)

# library(magrittr) # to use %>%
# library(purrr)


###### IMPORT ######
rm(list=ls())
getwd()
setwd("/Users/gianlorenzo/Desktop/TESI")
df <- list.files(path = "./ECB_data/", pattern = "*.csv", full.names = TRUE) %>% # Identify all CSV files
  .[ !grepl("GDP_MNA.Q.Y.I9.W2.S1.S1.B.B1GQ._Z._Z._Z.EUR.V.N.csv", .) ] %>%   # Exclude csv files
  lapply(read_csv) %>%                            # Store all files in list
  reduce(full_join, by = "DATE")                  # Full-join data sets into one data set
df <- na.omit(df)

#Drop col using subset method and rename col
df <- subset(select(df, -c("TIME PERIOD.x", "TIME PERIOD.y", "TIME PERIOD")))
setnames(df, colnames(df), 
         new = c("DATE","HICP", "M1", "SX5E")) # rename col

# transform data.frame into xts
Date_df <-strptime(df$DATE, "%Y-%m-%d",tz="GMT")
df.xts <-as.xts(df[,2:ncol(df)],Date_df)

# import GDP.csv
gdp<- list.files(path = "./ECB_data/", pattern = "*.csv", full.names = TRUE) %>% # Identify all CSV files
  .[ grepl("GDP_MNA.Q.Y.I9.W2.S1.S1.B.B1GQ._Z._Z._Z.EUR.V.N.csv", .) ] %>%   # select csv files
  read.csv()
gdp <- gdp[, -1]
colnames(gdp) <- c("Quarter", "GDP")
z <- read.zoo(gdp, FUN = function(x) as.yearqtr(gdp$Quarter))
GDP <- zooreg(na.approx(c(t(cbind(z, NA, NA)))), 
                start = as.yearmon(start(z)), freq = 12)
gdp <- fortify.zoo(GDP)[-1:-24,]; rm(GDP, z)

# transform into xts
Date_gdp <- Date_df[1:328,]
gdp.xts <-as.xts(gdp[,2:ncol(gdp)],Date_gdp)

# resizing
range(time(gdp.xts))
df.xts <- df.xts['1997-01-31/2024-04-30']

# merge
df.xts <- merge.xts(df.xts, gdp.xts);
names(df.xts)[names(df.xts) == 'gdp.xts'] <- 'GDP'   

# reorder
df.xts <- df.xts[,c(3,4,1,2)]; df.xts = na.omit(df.xts)

# first differences
# diff_xts = diff(df.xts)
# diff_df = as.data.frame(diff_xts)

df.xts <- log(df.xts) # log levels
lev <- as.data.frame(df.xts)

# cleaning
my_list <- c("df.xts", "lev") 
rm(list = setdiff(ls(), my_list))


# save enviroment
save(list = ls(), file = "obj.rData")

rm(list = ls())
# load("obj.rData")



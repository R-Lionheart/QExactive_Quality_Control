library(dplyr)
library(GGally)
library(readr)
library("reshape2")
library(stringr)
library(tidyr)
library(tidyverse)

areas.raw.test <- read.csv("./datafiles/ExampleSkylineOutput.csv", header = TRUE, na.strings = "#N/A", stringsAsFactors = F) %>%
  filter(Precursor.Ion.Name != "GTP",
         Precursor.Ion.Name != "Argininosuccinic Acid",
         Precursor.Ion.Name != "Cys-Gly")
areas.raw.noIS.test <- areas.raw.test %>%
  filter(!Protein.Name == "Internal Stds_neg")

run.type <- as.factor(tolower(str_extract(areas.raw.noIS.test$Replicate.Name, "(?<=_)[^_]+(?=_)")))
areas.raw.noIS.test[ , "Sample.Type"] <- run.type


before <- sapply(areas.raw.noIS.test[-1], class)
cols.to.change <- c(7:9, 12)
areas.raw.noIS.test[cols.to.change] <- sapply(areas.raw.noIS.test[cols.to.change], as.numeric)
after <- sapply(areas.raw.noIS.test[-1], class)

areas.split.test <- split(areas.raw.noIS.test, areas.raw.noIS.test$Sample.Type)
run.type.options.test <- names(areas.split.test)

RT.range.test <- sapply(split(areas.split.test[["std"]]$Retention.Time,
                         areas.split.test[["std"]]$Precursor.Ion.Name),
                   range, na.rm = T)

blk.range.test <- sapply(split(areas.split.test[["blk"]]$Area,
                          areas.split.test[["blk"]]$Precursor.Ion.Name),
                    range, na.rm = T)

###

blank.data.test <- areas.split.test[["blk"]]
samp.data.test <- areas.split.test[["smp"]]
if (any(run.type.options.test == "poo")) {
  poo.data.test <- areas.split.test[["poo"]]
  samp.data.test <- rbind(samp.data.test, poo.data.test)
}



cmpds.test <- unique(samp.data.test$Precursor.Ion.Name)
cmpd.samp.dfs.test <- split(samp.data.test, samp.data.test$Precursor.Ion.Name)
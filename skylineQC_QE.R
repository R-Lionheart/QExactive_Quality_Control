## Quality Control for Skyline output

## Note that the csv can have blanks but not characters in 
## the RT, MZ, and Area columns


# Description of code -----------------------------------------------------

# This code takes one file input (the output of a Skyline processing pipeline) and applies a custom quality control.

# Reads in the csv, accounting for NA values.
# Remove particular compounds.
# Create one dataframe of just internal standards by extracting Protein.Names that include Internal standards, and rename original dataframe to reflect that data has been extracted.
# Set numerical values for the parameters of the quality control, including suggestions.
# Identify run types by string splitting between standards, samples, blanks and pools, and renaming them to "std", "blk", "smp", "poo"
# Change variable types to factors or numeric and create variables for standard and blank runs. Split areas.raw.noIS by sample type.
# Create range of retention times and blank ranges. Include pooled samples in "sample" list. Isolate compounds and replicates. Create dataframes of samples grouped by compound.
# Set up acceptable Retention Time range matrix using predetermined RT.flex, and obtain sample RTs.
# Make new table for matrix format.
# Calculate signal to noise ratio and add areas, heights, ion ratios, names, etc by 
# Add warning if peak height is overloaded.
# Check for peak height.
# Check for signal to noise.
# Check for blank.
# Check for retention time range.
# Read data from internal standards.
# Attach blank data to output.
# Add comment to output and re-save as a new csv.

# TODO (rlionheart): Are the same compounds always removed? 
# TODO (rlionheart): Is ppm parts per million?
# TODO (rlionheart): When creating int.stds.data, should the filter parameter be a regex to account for various inputs indicating internal standards?
# TODO (rlionheart): No Mass.Error.PPM column to add to areas.raw.noIS? Also presents problems in lines 167 - 178.

# TODO (lionhearts): What specifically is happening in the ID run type section, as well as the RT.matrix section? 
# TODO (lionhearts): Where is the range coming from in RT.range under Range of Retention Times?
# TODO (lionhearts): The if statement on line 137 seems strange. Why pull "poo" from run.type.options?

library(dplyr)
library(GGally)
library(readr)
library("reshape2")
library(stringr)
library(tidyr)
library(tidyverse)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  filename <- args[1]
  # actions <- args[-1]
  # stopifnot(action %in% c("--filter", "--extract"))
  
  if (length(filename) == 0) {
    process(file("stdin"))
  } else {
    for (f in filename) {
      QC_parameters()
      process(f)
    }
  }
}

QC_parameters <- function() {
  cat("Pick an overload value (QE suggestion: 5.0e8): ");
  max.height <- readLines("stdin", n = 1);
  cat("Pick the minimum height to be counted as a 'real' peak (QE suggestion: HILIC - 1000, Cyano - 5000): " );
  min.height <- readLines("stdin", n = 1);
  cat("Pick retention time (RT) flexibility (QE suggestion: +/- 0.4 min for HILIC, +/- 0.2 min for Cyano): ")
  RT.flex <- readLines("stdin", n = 1);
  cat("Pick signal size comparison between sample and blank to merit inclusion (QE suggestion: +/- 20%): ")
  blk.thresh <- readLines("stdin", n = 1);
  cat("Pick acceptable signal to noise ratio value. Note: broader peaks create more background noise(QE suggestion: 5 for Cyano, 4 for HILIC. Negotiable.): ")
  SN.thresh <- readLines("stdin", n = 1);
  cat("Pick an absolute value for a cutoff for parts per million (ppm) (QE suggestion: 7): ")
  ppm.thresh <- readLines("stdin", n = 1);
  
  cat("Your max height value is:", max.height, "\n")
  cat("Your minimum height value is:", min.height, "\n")
  cat("Your retention time flexibility value is:", RT.flex, "\n")
  cat("Your blank threshold value is:", blk.thresh, "\n")
  cat("Your signal to noise ratio value is:", SN.thresh, "\n")
  cat("Your parts per million value is:", ppm.thresh, "\n")
  
  cat("Ready to move on? (Y/N):" );
  answer <- readLines("stdin", n = 1);
  if (regexpr(answer, 'y', ignore.case = TRUE) == 1) {
    continue = TRUE
  } else if (regexpr(answer, 'n', ignore.case = TRUE) == 1) {
    cat("Ok, let's try again!")
    break
  }
}

process <- function(filename) {
  areas.raw <- read.csv(file = filename, header = TRUE) %>%
    filter(Precursor.Ion.Name != "GTP",
           Precursor.Ion.Name != "Argininosuccinic Acid",
           Precursor.Ion.Name != "Cys-Gly")
  areas.raw.noIS <- areas.raw %>%
    filter(!Protein.Name == "Internal Stds_neg")
  print(head(areas.raw))
  print(head(areas.raw.noIS))
}

main()



# ## Set the parameters for the QC ----------------------------------------
# # Pick overload value. 
# # QE suggestion: 5e8
# max.height <- 5.0e8
# 
# # Pick the minimum height to be counted as a 'real' peak.
# # QE suggestion: HILIC - 1000, Cyano - 5000
# min.height <- 1000
# 
# # Pick retention time (RT) flexibility.
# # QE suggestion: +/- 0.4 min for HILIC, +/- 0.2 min for Cyano. 
# RT.flex <- 0.4
# 
# # Pick signal size comparison between sample and blank to merit inclusion.
# # QE suggestion: +/- 20%
# blk.thresh <- 0.5
# 
# # Pick acceptable signal to noise ratio value. Note: broader peaks create more background noise.
# # QE suggestion: 5 for Cyano, 4 for HILIC. Negotiable.
# SN.thresh <- 4
# 
# # Pick an absolute value for a cutoff for parts per million (ppm)
# # QE suggestion: 7
# ppm.thresh <- 7

# 
# ## ID run types ---------------------------
# # Standards (std), Samples (smp), Blanks (blk), Pooled (poo)
# t1 <- strsplit(areas.raw.noIS$Replicate.Name, "[_]")
# run.type <- vector()
# for (i in 1:length(t1)) {
#      this.one = t1[[i]][2]
#      run.type <- c(run.type, this.one)  
# }
# run.type <- tolower(run.type)
# 
# 
# ## Change variable types ---------------------------------------------------
# areas.raw.noIS$sample.type        <- as.factor(run.type)
# areas.raw.noIS$Precursor.Ion.Name <- as.factor(areas.raw.noIS$Precursor.Ion.Name)
# areas.raw.noIS$Area               <- as.numeric(areas.raw.noIS$Area)
# areas.raw.noIS$Retention.Time     <- as.numeric(areas.raw.noIS$Retention.Time)
# areas.raw.noIS$Background         <- as.numeric(areas.raw.noIS$Background)
# areas.raw.noIS$Height             <- as.numeric(areas.raw.noIS$Height)
# areas.raw.noIS$Mass.Error.PPM     <- as.numeric(areas.raw.noIS$Mass.Error.PPM)
# 
# stdRows <- areas.raw.noIS$sample.type == "std"
# blkRows <- areas.raw.noIS$sample.type == "blk"
# 
# areas.split <- split(areas.raw.noIS, areas.raw.noIS$sample.type)
# 
# run.type.options <- names(areas.split)
# 
# cmpd.blk.list <- split(areas.split[["blk"]],
#                        areas.split[["blk"]]$Precursor.Ion.Name)
# 
# ## Check the range of Retention Times and ion ratio in Standards ---------------------
# # Range of Retention Times (RTs) and pooled sample inclusion
# RT.range <- sapply(split(areas.split[["std"]]$Retention.Time, 
#                          areas.split[["std"]]$Precursor.Ion.Name),
#                          range, na.rm = T)
# 
# blk.range <- sapply(split(areas.split[["blk"]]$Area, 
#                           areas.split[["blk"]]$Precursor.Ion.Name),
#                           range, na.rm = T) 
# 
# samp.data <- areas.split[["smp"]]
# blank.data <- areas.split[["blk"]]
# 
# # If there are pooled samples, include them in the 'sample' list so only blanks and standards are excluded.
# if (any(run.type.options == "poo")) {
#      poo.data <- areas.split[["poo"]]
#      samp.data <- rbind(samp.data, poo.data)
# }
# 
# cmpds <- unique(samp.data$Precursor.Ion.Name)
# samples <- unique(samp.data$Replicate.Name)
# cmpd.samp.dfs <- split(samp.data, samp.data$Precursor.Ion.Name)
# 
# ## RT range matrix ------------------------------
# RT.ok <- RT.range
# RT.ok[1, ] <- RT.range[1, ] - RT.flex
# RT.ok[2, ] <- RT.range[2, ] + RT.flex
# 
# # Get the sample RTs
# RT.matrix <- c()
# for (i in 1:length(cmpds)) {
#      RT <- sapply(split(cmpd.samp.dfs[[cmpds[i]]]$Retention.Time,
#                         cmpd.samp.dfs[[cmpds[i]]]$Replicate.Name),
#                         mean, na.rm = T)
#      RT.matrix <- cbind(RT.matrix, RT)
#      colnames(RT.matrix)[ncol(RT.matrix)] <- as.character(cmpds[i])
# }
# 
# ## Convert short format matrix to long form ----------------------
# output <- melt(RT.matrix, value.name = "Retention Time")
# colnames(output) <- c("Replicate.Name", "Compound.Name", "Retention.Time")
# 
# ## Add areas, heights, and ion ratios to output ----------------------
# samp.data$S.N <- ((samp.data$Area + samp.data$Background) / samp.data$Background)
# 
# output <- full_join(output, samp.data[, c("Replicate.Name", "Precursor.Ion.Name",
#                                          "Area", "Height", "Background",
#                                          "Mass.Error.PPM", "S.N")], 
#                     by = c("Replicate.Name", "Compound.Name" = "Precursor.Ion.Name"))
# output <- output %>%
#      rename(ppm = Mass.Error.PPM) %>%
#      mutate(Notes = "",
#             rawArea = Area,
#             AreaBlkSub = Area,
#             BlkRatio = NA)
# 
# ## Is it overloaded? ---------------------------------------
# # Check peak height, if it is > max.height then keep the area data but make a note that it may be overloaded.
# for (i in 1:nrow(output)) {
#      if (!is.na(output$Height[i]) & output$Height[i] > max.height) {
#           output$Notes[i] <- paste(output$Notes[i], "overloaded?", sep = "")
#      }
# }
# 
# ## Absolute height check  ------
# for (i in 1:nrow(output)) {
#      if (!is.na(output$Height[i]) & output$Height[i] < min.height) {
#           output$Area[i] <- NA
#           output$Notes[i] <- paste(output$Notes[i], "too small", sep = "")
#      }
# }
# 
# ## Signal to Noise check -----------------------------------------------
# # If S/N is below threshold, throw it out!
# for (i in 1:nrow(output)) {
#   if (!is.na(output$S.N[i]) & output$S.N[i] < SN.thresh) {
#     output$Area[i] <- NA
#     output$Notes[i] <- paste(output$Notes[i], "bad S/N", sep = "")
#     
#   }
# }
# 
# ## Parts per million check -----------------------------------------------
# # If ppm is greater than threshold, throw it out!
# for (i in 1:nrow(output)) {
#   if (!is.na(output$ppm[i]) & abs(output$ppm[i]) > ppm.thresh) {
#     output$Notes[i] <- paste(output$Notes[i], "bad ppm", sep = "")
#     output$Area[i] <- NA
#   }
# }
# 
# ## Blank check -----------------------------------------------
# # If area is not much greater than blank, throw it out!
# for (i in 1:nrow(output)) {
#      key <- as.character(output$Compound.Name[i])
#      if (!is.na(output$Area[i]) & (output$Area[i] * blk.thresh) < mean(blk.range[, key])) {
#           output$Notes[i] <- paste(output$Notes[i], "comparable to blank", sep = "")
#           output$AreaBlkSub[i] <- output$Area[i] - mean(blk.range[, key])
#           output$BlkRatio[i] <-  mean(blk.range[, key]) / output$Area[i]
#           output$Area[i] <- NA
#      }
# }
# 
# ## RT range check --------------
# for (i in 1:nrow(output)) { 
#      key <- as.character(output$Compound.Name[i])
#      if (!is.na(output$Retention.Time[i]) & output$Retention.Time[i] > max(RT.ok[ , key])) {
#          output$Notes[i] <- paste(output$Notes[i], "Bad RT", sep = "")
#          output$Area[i] <- NA
#      }
# }
# 
# for (i in 1:nrow(output)) { 
#   key <- as.character(output$Compound.Name[i])
#   if (!is.na(output$Retention.Time[i]) & 
#       output$Retention.Time[i] < min(RT.ok[, key])) {
#       output$Notes[i] <- paste(output$Notes[i], "Bad RT", sep = "")
#       output$Area[i] <- NA
#   }
# }
# 
# ## Reinput all data from the internal standards-------------
# output <- full_join(output, samp.data[, c("Replicate.Name", "Precursor.Ion.Name", "Protein.Name")],
#                     by = c("Replicate.Name", "Compound.Name" = "Precursor.Ion.Name")) %>%
#           rename(Compound.Type = Protein.Name) %>%
#           mutate(Compound.Type = ifelse(grepl("Internal", Compound.Type), "Internal Stds", Compound.Type))
# 
# for (i in 1:nrow(output)) {
#   if (!is.na(output$Compound.Type[i]) & output$Compound.Type[i] == "Internal Stds") {
#     output$Notes[i] <- paste(output$Notes[i], "internal standard", sep = "")
#     if (!grepl(pattern = "loaded", output$Notes[i])) {
#       output$Area[i] <- output$rawArea[i]
#     }
#   }
# }
# 
# ## Attach blank data to output ---------------------------
# blank.data$Compound.Name <- blank.data$Precursor.Ion.Name
# blank.data$Notes <- rep("Blank used for comparison", nrow(blank.data))
# blank.data$S.N <- (blank.data$Area+blank.data$Background) / blank.data$Background
# blank.data$rawArea <- blank.data$Area
# blank.data$ppm <- blank.data$Mass.Error.PPM
# blank.data$Compound.Type <- blank.data$Protein.Name
# blank.data$AreaBlkSub <- blank.data$Area
# blank.data$BlkRatio <- NA
# final.output <- rbind(output, blank.data[, colnames(output)])
# 
#      
# ## Output with comment-------------------------
# # Ion name, area, was a peak removed?
# comment.text <- paste("# Hello! welcome to your data! ", "Overload height: ", 
#                       max.height, ". ", "RT flexibility: ", RT.flex, ". ",
#                       "Blank can be this fraction of a sample: ",blk.thresh, ". ", 
#                       "S/N threshold: " , SN.thresh, ". ",
#                       "Minimum peak height: ", min.height, ". ",
#                       "Processed on: ", Sys.time(), sep = "")
# new.filename <- paste("QEQC_output", filename, sep = "")
# con <- file(new.filename, open = "wt")
# writeLines(paste(comment.text), con)
# write.csv(final.output, con)
# close(con)
# 

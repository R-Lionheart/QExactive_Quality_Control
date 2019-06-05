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
# TODO (rlionheart): Figure out how to interactively enter columns that need variable changes (locations may change!)

# TODO (lionhearts): What specifically is happening in the ID run type section, as well as the RT.matrix section? 
# TODO (lionhearts): Where is the range coming from in RT.range under Range of Retention Times?
# TODO (lionhearts): Why does print(paste Variables) line in testing.R return two instances of the $Replicate.Name column?

library(dplyr)
library(GGally)
library(readr)
library("reshape2")
library(stringr)
library(tidyr)
library(tidyverse)

## Filter out unnecessary compounds---------------------------
filter_compounds <- function(area.data) {
  area.data %>%
    filter(Precursor.Ion.Name != "GTP",
           Precursor.Ion.Name != "Argininosuccinic Acid",
           Precursor.Ion.Name != "Cys-Gly",
           Protein.Name       != "Internal Stds_neg")
  return(area.data)
}

## Identify run types ---------------------------
# Standards (std), Samples (smp), Blanks (blk), Pooled (poo)
identify_runtypes <- function(area.data) {
  run.type <- as.factor(tolower(str_extract(area.data$Replicate.Name, "(?<=_)[^_]+(?=_)")))
  area.data[ , "Sample.Type"] <- run.type
  print("New column of run types:", quote = FALSE)
  print(head(run.type, n = 20))
  return(area.data)
}

## Change variable types ---------------------------------------------------
transform_variables <- function(area.data) {
  before <- lapply(area.data[-1], class)
  print("Original class variables ", quote = FALSE)
  print(paste(colnames(area.data)[-1], ":", before))
  
  area.data$Precursor.Ion.Name <- as.factor(area.data$Precursor.Ion.Name)
  cols.to.change <- c(4,5,6,8)
  area.data[cols.to.change] <- lapply(area.data[cols.to.change], as.numeric)
  
  after <- lapply(area.data[-1], class)
  print("New class variables ", quote = FALSE)
  print(paste(colnames(area.data)[-1], ":", after))
  
  return(area.data)
}


args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

if (length(filename) == 0) {
  stop("Please enter a csv of Skyline output.")
}

cat("Pick an overload value (QE suggestion: 5.0e8): ");
max.height <- as.double(readLines("stdin", n = 1));
cat("Pick the minimum height to be counted as a 'real' peak (QE suggestion: HILIC - 1000, Cyano - 5000): " );
min.height <- as.double(readLines("stdin", n = 1));
cat("Pick retention time (RT) flexibility (QE suggestion: +/- 0.4 min for HILIC, +/- 0.2 min for Cyano): ")
RT.flex <- as.double(readLines("stdin", n = 1));
cat("Pick signal size comparison between sample and blank to merit inclusion (QE suggestion: +/- 0.2): ")
blk.thresh <- as.double(readLines("stdin", n = 1));
cat("Pick acceptable signal to noise ratio value. Note: broader peaks create more background noise(QE suggestion: 5 for Cyano, 4 for HILIC. Negotiable.): ")
SN.thresh <- as.double(readLines("stdin", n = 1));
cat("Pick an absolute value for a cutoff for parts per million (ppm) (QE suggestion: 7): ")
ppm.thresh <- as.double(readLines("stdin", n = 1));

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

areas.raw <- read.csv(file = filename, header = TRUE, na.strings = "#N/A", stringsAsFactors = F)
areas.filtered <- filter_compounds(areas.raw)
areas.runtypes <- identify_runtypes(areas.filtered)
areas.transformed <- transform_variables(areas.runtypes)
areas.split <- split(areas.transformed, areas.transformed$Sample.Type)


## Check the range of Retention Times and ion ratio in Standards ---------------------
# Range of Retention Times (RTs) and pooled sample inclusion

RT.range <- lapply(split(areas.split[["std"]]$Retention.Time,
                         areas.split[["std"]]$Precursor.Ion.Name),
                         range, na.rm = T)

blank.range <- lapply(split(areas.split[["blk"]]$Area,
                          areas.split[["blk"]]$Precursor.Ion.Name),
                          range, na.rm = T)

samp.data <- rbind(areas.split[["smp"]], areas.split[["poo"]])
blank.data <- areas.split[["blk"]]


cmpds <- unique(samp.data$Precursor.Ion.Name)
cmpd.samp.dfs <- split(samp.data, samp.data$Precursor.Ion.Name)


## RT range matrix ------------------------------
RT.ok <- RT.range
RT.ok[1, ] <- RT.range[1, ] - RT.flex
RT.ok[2, ] <- RT.range[2, ] + RT.flex

# Get the sample RTs, new
RT.matrix <- c()
for (i in 1:length(cmpds)) {
  RT <- sapply(split(cmpd.samp.dfs[[cmpds[i]]]$Retention.Time,
                     cmpd.samp.dfs[[cmpds[i]]]$Replicate.Name),
               mean, na.rm = T)
  RT.matrix <- cbind(RT.matrix, RT)
  colnames(RT.matrix)[ncol(RT.matrix)] <- as.character(cmpds[i])
}

## Convert short format matrix to long form ----------------------
output <- melt(RT.matrix, value.name = "Retention Time")
colnames(output) <- c("Replicate.Name", "Compound.Name", "Retention.Time")

## Add areas, heights, and ion ratios to output ----------------------
samp.data$S.N <- ((samp.data$Area + samp.data$Background) / samp.data$Background)

output <- full_join(output, samp.data[, c("Replicate.Name", "Precursor.Ion.Name",
                                          "Area", "Height", "Background",
                                          "Mass.Error.PPM", "S.N")],
                    by = c("Replicate.Name", "Compound.Name" = "Precursor.Ion.Name"))
output <- output %>%
  rename(ppm = Mass.Error.PPM) %>%
  mutate(Notes = "",
         rawArea = Area,
         AreaBlkSub = Area,
         BlkRatio = NA)

## Is it overloaded? ---------------------------------------
# Check peak height, if it is > max.height then keep the area data but make a note that it may be overloaded.
for (i in 1:nrow(output)) {
  if (!is.na(output$Height[i]) & output$Height[i] > max.height) {
    output$Notes[i] <- paste(output$Notes[i], "overloaded?", sep = "")
  }
}

## Absolute height check  ------
for (i in 1:nrow(output)) {
  if (!is.na(output$Height[i]) & output$Height[i] < min.height) {
    output$Area[i] <- NA
    output$Notes[i] <- paste(output$Notes[i], "too small", sep = "")
  }
}

## Signal to Noise check -----------------------------------------------
# If S/N is below threshold, throw it out!
for (i in 1:nrow(output)) {
  if (!is.na(output$S.N[i]) & output$S.N[i] < SN.thresh) {
    output$Area[i] <- NA
    output$Notes[i] <- paste(output$Notes[i], "bad S/N", sep = "")
    
  }
}

## Parts per million check -----------------------------------------------
# If ppm is greater than threshold, throw it out!
for (i in 1:nrow(output)) {
  if (!is.na(output$ppm[i]) & abs(output$ppm[i]) > ppm.thresh) {
    output$Notes[i] <- paste(output$Notes[i], "bad ppm", sep = "")
    output$Area[i] <- NA
  }
}

## Blank check -----------------------------------------------
# If area is not much greater than blank, throw it out!
for (i in 1:nrow(output)) {
  key <- as.character(output$Compound.Name[i])
  if (!is.na(output$Area[i]) & (output$Area[i] * blk.thresh) < mean(blank.range[, key])) {
    output$Notes[i] <- paste(output$Notes[i], "comparable to blank", sep = "")
    output$AreaBlkSub[i] <- output$Area[i] - mean(blank.range[, key])
    output$BlkRatio[i] <-  mean(blank.range[, key]) / output$Area[i]
    output$Area[i] <- NA
  }
}

## RT range check --------------
for (i in 1:nrow(output)) {
  key <- as.character(output$Compound.Name[i])
  if (!is.na(output$Retention.Time[i]) & output$Retention.Time[i] > max(RT.ok[ , key])) {
    output$Notes[i] <- paste(output$Notes[i], "Bad RT", sep = "")
    output$Area[i] <- NA
  }
}

for (i in 1:nrow(output)) {
  key <- as.character(output$Compound.Name[i])
  if (!is.na(output$Retention.Time[i]) &
      output$Retention.Time[i] < min(RT.ok[, key])) {
    output$Notes[i] <- paste(output$Notes[i], "Bad RT", sep = "")
    output$Area[i] <- NA
  }
}

## Reinput all data from the internal standards-------------
output <- full_join(output, samp.data[, c("Replicate.Name", "Precursor.Ion.Name", "Protein.Name")],
                    by = c("Replicate.Name", "Compound.Name" = "Precursor.Ion.Name")) %>%
  rename(Compound.Type = Protein.Name) %>%
  mutate(Compound.Type = ifelse(grepl("Internal", Compound.Type), "Internal Stds", Compound.Type))

for (i in 1:nrow(output)) {
  if (!is.na(output$Compound.Type[i]) & output$Compound.Type[i] == "Internal Stds") {
    output$Notes[i] <- paste(output$Notes[i], "internal standard", sep = "")
    if (!grepl(pattern = "loaded", output$Notes[i])) {
      output$Area[i] <- output$rawArea[i]
    }
  }
}

## Attach blank data to output ---------------------------
blank.data$Compound.Name <- blank.data$Precursor.Ion.Name
blank.data$Notes         <- rep("Blank used for comparison", nrow(blank.data))
blank.data$S.N           <- (blank.data$Area+blank.data$Background) / blank.data$Background
blank.data$rawArea       <- blank.data$Area
blank.data$ppm           <- blank.data$Mass.Error.PPM
blank.data$Compound.Type <- blank.data$Protein.Name
blank.data$AreaBlkSub    <- blank.data$Area
blank.data$BlkRatio      <- NA
final.output             <- rbind(output, blank.data[, colnames(output)])


## Output with comment-------------------------
# Ion name, area, was a peak removed?
comment.text <- paste("# Hello! welcome to your data! ", "Overload height: ",
                      max.height, ". ", "RT flexibility: ", RT.flex, ". ",
                      "Blank can be this fraction of a sample: ",blk.thresh, ". ",
                      "S/N threshold: " , SN.thresh, ". ",
                      "Minimum peak height: ", min.height, ". ",
                      "Processed on: ", Sys.time(), sep = "")

new.filename <- paste("QEQC_output", filename, sep = "")
con <- file(new.filename, open = "wt")
writeLines(paste(comment.text), con)
write.csv(final.output, con)
close(con)

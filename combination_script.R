library(rlist)
library(tidyverse)

# This script is intended to be used for .csv files output from Skyline (https://skyline.ms/project/home/begin.view?).
# Run using the command line. When in the proper directory, enter Rscript <this script> <skyline output csv> <blank sample matching csv>
# The script will then guide you through the steps to obtaining your results, filtered for quality control!

# TODO (kheal): Anything to add in the general description?

PromptToContinue <- function() {
  repeat {
    cat("Ready to move on? (Y/N): " )
    answer <- tolower(readLines("stdin", n = 1))
    if (answer == "y")  {
      break
    } else if (answer == "n") {
      stop("Not ready to move on.")
    }
  }
}

## Set your flags to find standards to base RT off of
PromptForStdFlags <- function() {
  std.flags <- list()
  repeat {
    cat("Sample flag (enter blank if finished): ")
    flag <- readLines("stdin", n = 1);
    if (flag == "") {
      return(std.flags)
    }
    std.flags <- list.append(std.flags, flag)
  }
}

## Begin data processing here ----------------------------------------------

## Transform variable classes
TransformVariables <- function(skyline.output) {
  before <- lapply(skyline.output, class)
  print("Original class variables ", quote = FALSE)
  print(paste(colnames(skyline.output), ":", before))

  skyline.output <- skyline.output %>%
    mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
    mutate(Area           = as.numeric(as.character(Area))) %>%
    mutate(Background     = as.numeric(as.character(Background))) %>%
    mutate(Mass.Error.PPM = as.numeric(as.character(Mass.Error.PPM)))

  after <- lapply(skyline.output, class)
  print("New class variables ", quote = FALSE)
  print(paste(colnames(skyline.output), ":", after))

  return(skyline.output)
}

## Create Signal to noise (SN), Parts per million (ppm), and Area minimum (area.min) flags.
CreateFirstFlags <- function(skyline.output, area.min, SN.min, ppm.flex) {
  first.flags <- skyline.output %>%
    filter(Replicate.Name %in% blank.matcher$Replicate.Name) %>%
    mutate(SNFlag      = ifelse(((Area / Background) < SN.min), "SNFlag", NA)) %>%
    mutate(ppmFlag     = ifelse((Mass.Error.PPM > ppm.flex), "ppmFlag", NA)) %>%
    mutate(areaminFlag = ifelse((Area < area.min), "areaminFlag", NA))

  return(first.flags)
}

## Create dataset with Retention Time (RT) flags, when compared to a standard or template.
CreateRTFlags <- function(skyline.output, std.flags) {
  retention.time.flags <- skyline.output %>%
    filter(Replicate.Name %in% std.flags) %>%
    select(Precursor.Ion.Name, Retention.Time) %>%
    group_by(Precursor.Ion.Name) %>%
    summarise(RT_ref = mean((Retention.Time), na.rm = TRUE))
  
  return(retention.time.flags)
}

## Create dataset with blank flags
CreateBlankFlags <- function(skyline.output, blank.matcher) {
  blank.flags <- skyline.output %>%
    filter(Replicate.Name %in% blank.matcher$Blank.Name) %>%
    rename(Blank.Name = Replicate.Name,
           Blank.Area = Area) %>%
    select(Blank.Name, Precursor.Ion.Name, Blank.Area) %>%
    left_join(blank.matcher) %>% select(-Blank.Name)
  
  return(blank.flags)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please enter a csv of Skyline output and sample blanks to compare.")
}

skyline.output <- read.csv(file = args[1], header = TRUE)
blank.matcher <- read.csv(file = args[2], header = TRUE)

## Set parameters for quality control.
cat("Pick the minimum height to be counted as a 'real' peak (QE suggestion: HILIC - 1000, Cyano - 5000): " )
area.min        <- as.double(readLines("stdin", n = 1))
cat("Pick retention time (RT) flexibility (QE suggestion: +/- 0.4 min for HILIC, +/- 0.2 min for Cyano): ")
RT.flex         <- as.double(readLines("stdin", n = 1))
cat("Pick signal size comparison between sample and blank to merit inclusion (QE suggestion: +/- 0.2): ")
blank.ratio.max <- as.double(readLines("stdin", n = 1))
cat("Pick acceptable signal to noise ratio value. Note: broader peaks create more background noise(QE suggestion: 5 for Cyano, 4 for HILIC.): ")
SN.min          <- as.double(readLines("stdin", n = 1))
cat("Pick an absolute value for a cutoff for parts per million (ppm) (QE suggestion: 7): ")
ppm.flex        <- as.double(readLines("stdin", n = 1))

cat("Your minimum height value is:", area.min, "\n")
cat("Your retention time flexibility value is:", RT.flex, "\n")
cat("Your blank threshold value is:", blank.ratio.max, "\n")
cat("Your signal to noise ratio value is:", SN.min, "\n")
cat("Your parts per million value is:", ppm.flex, "\n")

PromptToContinue()
std.flags <- PromptForStdFlags()

skyline.columns.dropped <- skyline.output %>%
  select(-Protein.Name, -Protein)
skyline.classes.transformed <- TransformVariables(skyline.columns.dropped)

SNPPMAM.flags <- CreateFirstFlags(skyline.classes.transformed, area.min, SN.min, ppm.flex)
RT.flags      <- CreateRTFlags(skyline.classes.transformed, std.flags)
blank.flags   <- CreateBlankFlags(skyline.classes.transformed, blank.matcher)

## Join those datasets!

## Add RT flags in to dataset.
first.join <- SNPPMAM.flags %>%
  left_join(RT.flags) %>%
  mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
  mutate(RTFlag = ifelse((abs(Retention.Time - RT_ref) > RT.flex), "RTFlag", NA))

## Add Blank flags into dataset
second.join <- first.join %>%
  left_join(blank.flags) %>%
  mutate(BlankFlag = ifelse((Area / Blank.Area) < blank.ratio.max, "BlankFlag", NA))

## Finally, combine all the flags and throw out any peak with a flag
last.join <- second.join %>%
  mutate(Flags = paste(SNFlag, ppmFlag, areaminFlag, RTFlag, BlankFlag, sep = ", ")) %>%
  mutate(Flags = as.character(Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
  mutate(QC_area = ifelse(str_detect(Flags, "Flag"), NA, Area))

write.csv(last.join, "RML_QCoutput.csv")
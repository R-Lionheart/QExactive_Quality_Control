library(rlist)
library(tidyverse)

# This script is intended to be used for .csv files output from Skyline (https://skyline.ms/project/home/begin.view?).
# Run using the command line. When in the proper directory, enter Rscript <this script> <skyline output csv> <blank sample matching csv>
# The script will then guide you through the steps to obtaining your results, filtered for quality control!

# TODO (kheal): Anything to add in the general description?

args <- commandArgs(trailingOnly = TRUE)
filename1 <- args[1]
filename2 <- args[2]

if (filename1 == 0 | filename2 == 0) {
  stop("Please enter a csv of Skyline output and sample blanks to compare.")
}

skyline.output <- read.csv(file = filename1)
blank.matcher <- read.csv(file = filename2)

## Set parameters for quality control.
cat("Pick the minimum height to be counted as a 'real' peak (QE suggestion: HILIC - 1000, Cyano - 5000): " )
Areamin <- as.double(readLines("stdin", n = 1))
cat("Pick retention time (RT) flexibility (QE suggestion: +/- 0.4 min for HILIC, +/- 0.2 min for Cyano): ")
RTflex <- as.double(readLines("stdin", n = 1))
cat("Pick signal size comparison between sample and blank to merit inclusion (QE suggestion: +/- 0.2): ")
BlankRatiomax <- as.double(readLines("stdin", n = 1))
cat("Pick acceptable signal to noise ratio value. Note: broader peaks create more background noise(QE suggestion: 5 for Cyano, 4 for HILIC. Negotiable.): ")
SNmin <- as.double(readLines("stdin", n = 1))
cat("Pick an absolute value for a cutoff for parts per million (ppm) (QE suggestion: 7): ")
ppmflex <- as.double(readLines("stdin", n = 1))

cat("Your minimum height value is:", Areamin, "\n")
cat("Your retention time flexibility value is:", RTflex, "\n")
cat("Your blank threshold value is:", BlankRatiomax, "\n")
cat("Your signal to noise ratio value is:", SNmin, "\n")
cat("Your parts per million value is:", ppmflex, "\n")

cat("Ready to move on? (Y/N):" )
answer <- readLines("stdin", n = 1)

if (regexpr(answer, 'y', ignore.case = TRUE) == 1) {
  continue = TRUE
} else if (regexpr(answer, 'n', ignore.case = TRUE) == 1) {
  stop("Ok, let's try again!")
  break
}

# TODO (lionhearts)
random.sample <- function(x) {
  repeat {
    # do something
    i <- sample(nrow(df), 1)
    x <- df[sample(nrow(df), 1), ]
    # exit if the condition is met
    if (x$SCORE > 0) break
  }
  return(x)
}

std_flag_input <- function() {
  StdFlags <- list()
  repeat {
    cat("Sample flag (leave blank if finished): ")
    flag <- readLines("stdin", n = 1);
    if (flag == "") {
      #print(StdFlags)
      return(StdFlags)
    }
    StdFlags <- list.append(StdFlags, flag)
  }
}

std_flag_input()


# ## Begin data processing here ----------------------------------------------
# 
# ## Drop columns
# drop_columns <- function(skyline.output) {
#   # Removes Protein Name and Protein columns from skyline output.
#   skyline.output %>%
#     select(-Protein.Name, -Protein)
# 
#   return(skyline.output)
# }
# 
# ## Transform variable classes
# transform_variables <- function(skyline.output) {
#   before <- lapply(skyline.output, class)
#   print("Original class variables ", quote = FALSE)
#   print(paste(colnames(skyline.output), ":", before))
# 
#   mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
#   mutate(Area = as.numeric(as.character(Area))) %>%
#   mutate(Background = as.numeric(as.character(Background))) %>%
#   mutate(Mass.Error.PPM = as.numeric(as.character(Mass.Error.PPM)))
# 
#   after <- lapply(skyline.output, class)
#   print("New class variables ", quote = FALSE)
#   print(paste(colnames(skyline.output), ":", after))
# 
#   return(skyline.output)
# }
# 
# ## Create Signal to noise (SN), Parts per million (ppm), and Area minimum (areamin) flags. Add to dataset.
# first_flags <- function(skyline.output) {
#   first_flags_added <- skyline.output %>%
#     filter(Replicate.Name %in% BlankMatcher$Replicate.Name) %>%
#     mutate(SNFlag = ifelse((Area/Background < SNmin), "SNFlag", NA)) %>%
#     mutate(ppmFlag = ifelse((Mass.Error.PPM > ppmflex), "ppmFlag", NA)) %>%
#     mutate(areaminFlag = ifelse((Area < Areamin), "areaminFlag", NA))
# 
#   return(SNPPMAM.flags.made)
# }
# 
# ## Create dataset with Retention Time (RT) flags, when compared to a standard or template.
# RT_flags <- function(skyline.output) {
#   Retention <- skyline.output %>%
#     filter(Replicate.Name %in% StdFlag) %>%
#     select(Precursor.Ion.Name, Retention.Time) %>%
#     group_by(Precursor.Ion.Name) %>%
#     summarise(RT_ref = mean((Retention.Time), na.rm = TRUE))
# 
#   return(RT_flags)
# }
# 
# ## Create dataset with blank flags
# Blank_flags <- function(skyline.output) {
#   Area_of_Blanks <- skyline.output %>%
#     filter(Replicate.Name %in% BlankMatcher$Blank.Name) %>%
#     rename(Blank.Name = Replicate.Name,
#            Blank.Area = Area) %>%
#     select(Blank.Name, Precursor.Ion.Name, Blank.Area) %>%
#     left_join(BlankMatcher) %>% select(-Blank.Name)
# }
# 
# 
# columns.dropped <- drop_columns(skyline.output)
# classes.transformed <- transform_variables(columns.dropped)
# SNPPMAM.flags.made <- first_flags(skyline.output)
# RT.flags.made <- RT_flags(skyline.output)
# Blank.flags.made <- Blank_flags(skyline.output)
# 
# ## Join those datasets!
# 
# ## Add RT flags in to dataset.
# first_join <- SNPPMAM.flags.made %>%
#   left_join(RT_flags) %>%
#   mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
#   mutate(RTFlag = ifelse((abs(Retention.Time - RT_ref) > RTflex), "RTFlag", NA))
# 
# 
# ## Add Blank flags into dataset
# second_join <- first_join %>%
#   left_join(Area_of_Blanks) %>%
#   mutate(BlankFlag = ifelse(Area/Blank.Area < BlankRatiomax, "BlankFlag", NA))
# 
# 
# ## Finally, combine all the flags and throw out any peak with a flag
# last_join <- second_join %>%
#   mutate(Flags = paste(SNFlag, ppmFlag, areaminFlag, RTFlag, BlankFlag, sep = ", ")) %>%
#   mutate(Flags = as.character(Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
#   mutate(QC_area = ifelse(str_detect(Flags, "Flag"), NA, Area))
# 
# 
# write.csv(Dat4, "QCoutput.csv")
# 
# 
# 
# 

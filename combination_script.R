library(rlist)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
filename1 <- args[1]
filename2 <- args[2]

if (filename1 == 0 | filename2 == 0) {
  stop("Please enter a csv of Skyline output and sample blanks to compare.")
}

machine.output <- read.csv(file = filename1)
blank.matcher <- read.csv(file = filename2)

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

# TODO (lionhearts) 
if (regexpr(answer, 'y', ignore.case = TRUE) == 1) {
  continue = TRUE
} else if (regexpr(answer, 'n', ignore.case = TRUE) == 1) {
  stop("Ok, let's try again!")
  break
}

# TODO (lionhearts)


std_flag_input <- function() {
  last_one <- readlines("stdin", n = 1)
  StdFlags = list()
  while(last_one != 'y') {
    n <- cat(prompt = "Sample flag: ")
    StdFlags = c(StdFlags, n)
  }
  if (regexpr(last_one, 'y', ignore.case = TRUE) == 1) {
  } return (StdFlags)
}

std_flag_input()

## Begin data processing here ----------------------------------------------

# ## Drop columns
# drop_columns <- function(machine.output) {
#   machine.output %>%
#     select(-Protein.Name, -Protein)
# 
#   return(machine.output)
# }
# 
# ## Transform variable classes
# transform_variables <- function(machine.output) {
#   before <- lapply(machine.output, class)
#   print("Original class variables ", quote = FALSE)
#   print(paste(colnames(machine.output), ":", before))
# 
#   mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
#   mutate(Area = as.numeric(as.character(Area))) %>%
#   mutate(Background = as.numeric(as.character(Background))) %>%
#   mutate(Mass.Error.PPM = as.numeric(as.character(Mass.Error.PPM)))
# 
#   after <- lapply(machine.output, class)
#   print("New class variables ", quote = FALSE)
#   print(paste(colnames(machine.output), ":", after))
# 
#   return(machine.output)
# }
# 
# ## Create Signal to noise (SN), Parts per million (ppm), and Area minimum (areamin) flags. Add to dataset.
# first_flags <- function(machine.output) {
#   first_flags_added <- machine.output %>%
#     filter(Replicate.Name %in% BlankMatcher$Replicate.Name) %>%
#     mutate(SNFlag = ifelse((Area/Background < SNmin), "SNFlag", NA)) %>%
#     mutate(ppmFlag = ifelse((Mass.Error.PPM > ppmflex), "ppmFlag", NA)) %>%
#     mutate(areaminFlag = ifelse((Area < Areamin), "areaminFlag", NA))
# 
#   return(SNPPMAM.flags.made)
# }
# 
# ## Create dataset with Retention Time (RT) flags, when compared to a standard or template.
# RT_flags <- function(machine.output) {
#   Retention <- machine.output %>%
#     filter(Replicate.Name %in% StdFlag) %>%
#     select(Precursor.Ion.Name, Retention.Time) %>%
#     group_by(Precursor.Ion.Name) %>%
#     summarise(RT_ref = mean((Retention.Time), na.rm = TRUE))
# 
#   return(RT_flags)
# }
# 
# ## Create dataset with blank flags
# Blank_flags <- function(machine.output) {
#   Area_of_Blanks <- machine.output %>%
#     filter(Replicate.Name %in% BlankMatcher$Blank.Name) %>%
#     rename(Blank.Name = Replicate.Name,
#            Blank.Area = Area) %>%
#     select(Blank.Name, Precursor.Ion.Name, Blank.Area) %>%
#     left_join(BlankMatcher) %>% select(-Blank.Name)
# }
# 
# 
# columns.dropped <- drop_columns(machine.output)
# classes.transformed <- transform_variables(columns.dropped)
# SNPPMAM.flags.made <- first_flags(machine.output)
# RT.flags.made <- RT_flags(machine.output)
# Blank.flags.made <- Blank_flags(machine.output)
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
# ##  Finally, combine all the flags and throw out any peak with a flag
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

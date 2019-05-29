## Quality Control for Skyline output

## Note that the csv can have blanks but not characters in 
## the RT, MZ, and Area columns
## Load needed libraries ------

library(dplyr)
library(GGally)
library(githubinstall)
library("plyr")
library(readr)
library("reshape2")
library(stringr)
library(tidyr)

## Load your data here---------
# filename <- "HILICNeg_QE_SkylineResults.csv"
filename <- "HILICPos_QE_SkylineResults.csv"

areas.raw <- read.csv(filename, na.strings = "#N/A", stringsAsFactors = F)

## ditch compounds that are missing data but are also not present and are confusing the code:
areas.raw <- areas.raw %>% 
     filter(Precursor.Ion.Name!="GTP",
            Precursor.Ion.Name!="Argininosuccinic Acid",
            Precursor.Ion.Name!="Cys-Gly")

## extract internal standards data ------------------
int.stds.data <- areas.raw %>% filter(Protein.Name=="Internal Stds_neg")
areas.raw.noIS <- areas.raw #%>% filter(Protein.Name!="Internal Stds_neg")

## Set the parameters for the QC ----------------------------------------
## pick overload value (suggestion: 5e8 for QE)
max.height <- 5.0e8
## pick the minimum height to be counted as a 'real' peak
## (suggestions on QE: HILIC - 1000, cyano - 5000)
min.height <- 1000
## pick RT flexibility 
## (suggeestions: +/- 0.4 min for HILIC, +/- 0.2 min for cyano)
RT.flex <- 0.4
## pick how big of a signal the samples need to have in comparison
## to a blank in order to 'count' (suggestion: +/- 20%)
blk.thresh <- 0.5
## pick the value that will be the acceptable signal to noise ratio
## (suggestion maybe 5 for cyano and 4 for HILIC?  (broader peaks make for more background))
SN.thresh <- 4
## pick a value for a cutoff for ppm (absolute value)
ppm.thresh <- 7

## ID run types ---------------------------
## Standards (std), Samples (smp), Blanks (blk), Pooled (poo)
t1<-strsplit(areas.raw.noIS$Replicate.Name, "[_]")
type <- vector()
for(i in 1:length(t1)){
     this.one = t1[[i]][2]
     type <- c(type, this.one)  
}
type<- tolower(type)

areas.raw.noIS$sample.type <- as.factor(type)
#areas.raw.noIS$Precursor.Ion.Name <- as.factor(areas.raw.noIS$Precursor.Ion.Name)

areas.raw.noIS$Area <- as.numeric(areas.raw.noIS$Area)
areas.raw.noIS$Retention.Time <- as.numeric(areas.raw.noIS$Retention.Time)
areas.raw.noIS$Background <- as.numeric(areas.raw.noIS$Background)
areas.raw.noIS$Height <- as.numeric(areas.raw.noIS$Height)
areas.raw.noIS$Mass.Error.PPM <- as.numeric(areas.raw.noIS$Mass.Error.PPM)

stdRows <- areas.raw.noIS$sample.type == "std"
blkRows <- areas.raw.noIS$sample.type == "blk"

areas.split<-split(areas.raw.noIS,areas.raw.noIS$sample.type)

type.options <- names(areas.split)

cmpd.blk.list<- split(areas.split[["blk"]],
                      areas.split[["blk"]]$Precursor.Ion.Name)

## Check the range of RT and ion ratio in Standards ---------------------
## Range of RTs
RT.range<-sapply(split(areas.split[["std"]]$Retention.Time,areas.split[["std"]]$Precursor.Ion.Name),
                 range,na.rm = T)

blk.range<-sapply(split(areas.split[["blk"]]$Area,areas.split[["blk"]]$Precursor.Ion.Name),
                 range,na.rm = T)
samp.data <- areas.split[["smp"]]
blank.data <- areas.split[["blk"]]

## if there are pooled samples, include them in the 'sample' list
if(any(type.options=="poo")){
     poo.data <- areas.split[["poo"]]
     samp.data <- rbind(samp.data,poo.data)
}
cmpds <- unique(samp.data$Precursor.Ion.Name)
samples <- unique(samp.data$Replicate.Name)

cmpd.samp.dfs <- split(samp.data,samp.data$Precursor.Ion.Name)

## RT range matrix ------------------------------
RT.ok <- RT.range
RT.ok[1,] <- RT.range[1,] - RT.flex
RT.ok[2,] <- RT.range[2,] + RT.flex

# Get the sample RTs
RT.matrix <- c()
for (i in 1:length(cmpds)){
     RT <- sapply(split(cmpd.samp.dfs[[cmpds[i]]]$Retention.Time,
                        cmpd.samp.dfs[[cmpds[i]]]$Replicate.Name),
                  mean, na.rm = T)
     RT.matrix <- cbind(RT.matrix, RT)
     colnames(RT.matrix)[ncol(RT.matrix)]<-as.character(cmpds[i])
}

#convert short format matrix to long form ----------------------
output <- melt(RT.matrix, value.name = "Retention Time")

colnames(output) <- c("Replicate.Name","Compound.Name","Retention.Time")

## Add areas, heights, and Ion Ratios to output ----------------------
# S.N <- c()
samp.data$S.N <- ((samp.data$Area +samp.data$Background) / samp.data$Background)

output <- full_join(output,samp.data[,c("Replicate.Name","Precursor.Ion.Name",
                                         "Area","Height","Background",
                                         "Mass.Error.PPM","S.N")], 
                    by = c("Replicate.Name",
                           "Compound.Name"="Precursor.Ion.Name"))
output <- output %>%
     rename(ppm  = Mass.Error.PPM) %>%
     mutate(Notes = "",
            rawArea = Area,
            AreaBlkSub = Area,
            BlkRatio = NA)

## Is it overloaded? ---------------------------------------
## Check peak height
## If it is > specified value (max.height) then keep the area data
## but make a note that it may be overloaded
for (i in 1:nrow(output)){
     if (!is.na(output$Height[i]) & output$Height[i]>max.height){
          output$Notes[i]<-paste(output$Notes[i],"overloaded?",sep="")
     }
}

## Absolute height check  ------
for (i in 1:nrow(output)){
     if (!is.na(output$Height[i]) & output$Height[i]<min.height){
          output$Area[i] <- NA
          output$Notes[i] <- paste(output$Notes[i],"too small",sep="")
     }
}

## Signal to Noise check -----------------------------------------------
## If S/N is below threshold, throw it out!
for (i in 1:nrow(output)){
  if (!is.na(output$S.N[i]) & output$S.N[i]<SN.thresh){
    output$Area[i] <- NA
    output$Notes[i] <- paste(output$Notes[i],"bad S/N",sep="")
    
  }
}

## ppm  check -----------------------------------------------
## If ppm is greater than threshold, throw it out!
for (i in 1:nrow(output)){
  if (!is.na(output$ppm[i]) & abs(output$ppm[i])>ppm.thresh){
    output$Notes[i] <- paste(output$Notes[i],"bad ppm",sep="")
    output$Area[i] <- NA
  }
}

## Blank check -----------------------------------------------
## If area is not much greater than blank, throw it out!
for (i in 1:nrow(output)){
     key <- as.character(output$Compound.Name[i])
     if (!is.na(output$Area[i]) & 
         (output$Area[i]*blk.thresh)<mean(blk.range[ , key])){
          output$Notes[i] <- paste(output$Notes[i],"comparable to blank",sep="")
          output$AreaBlkSub[i] <- output$Area[i] - mean(blk.range[ , key])
          output$BlkRatio[i] <-  mean(blk.range[ , key])/output$Area[i]
          output$Area[i] <- NA
     }
}

##RT range check --------------
for (i in 1:nrow(output)) { 
     key <- as.character(output$Compound.Name[i])
     if (!is.na(output$Retention.Time[i]) & 
         output$Retention.Time[i] > max(RT.ok[ , key])){
          output$Notes[i] <- paste(output$Notes[i],"Bad RT",sep="")
          output$Area[i] <- NA
     }
}
for (i in 1:nrow(output))
{ key <- as.character(output$Compound.Name[i])
if (!is.na(output$Retention.Time[i]) & 
    output$Retention.Time[i] < min(RT.ok[ , key])){
  output$Notes[i] <- paste(output$Notes[i],"Bad RT",sep="")
  output$Area[i] <- NA
}
}

## Reinput all data from the internal standards-------------
output <- full_join(output, samp.data[,c("Replicate.Name","Precursor.Ion.Name",
                                         "Protein.Name")],
                    by = c("Replicate.Name", "Compound.Name" = "Precursor.Ion.Name")) %>%
     rename(Compound.Type = Protein.Name) %>%
     mutate(Compound.Type = ifelse(grepl("Internal",Compound.Type),"Internal Stds", Compound.Type))

for (i in 1:nrow(output)){
  if (!is.na(output$Compound.Type[i]) & output$Compound.Type[i] == "Internal Stds"){
    output$Notes[i] <- paste(output$Notes[i],"internal standard",sep="")
    if (!grepl(pattern = "loaded",output$Notes[i])){
      output$Area[i] <- output$rawArea[i]
    }
  }
}

## attach blank data to output ---------------------------
blank.data$Compound.Name <- blank.data$Precursor.Ion.Name
blank.data$Notes <- rep("Blank used for comparison",nrow(blank.data))
blank.data$S.N <- (blank.data$Area+blank.data$Background)/blank.data$Background
blank.data$rawArea <- blank.data$Area
blank.data$ppm <- blank.data$Mass.Error.PPM
blank.data$Compound.Type <- blank.data$Protein.Name
blank.data$AreaBlkSub <- blank.data$Area
blank.data$BlkRatio <- NA
final.output <- rbind(output, blank.data[,colnames(output)])

     
## Output with: comment -------------------------
## Ion name, area, was a peak removed?
comment.text <- paste("# Hello! welcome to your data! ","Overload height: ", 
                      max.height, ". ", "RT flexibility: ", RT.flex, ". ",
                      "Blank can be this fraction of a sample: ",blk.thresh, ". ", 
                      "S/N threshold: " , SN.thresh, ". ",
                      "Minimum peak height: ", min.height, ". ",
                      "Processed on: ", Sys.time(), sep="")
new.filename <- paste("QEQC_output",filename,sep="")
con <- file(new.filename, open="wt")
writeLines(paste(comment.text), con)
write.csv(final.output, con)
close(con)


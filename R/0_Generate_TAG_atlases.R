#### 0. Install and load R libraries ####
# install.packages("rgl")
library(rgl)
library(geomorph)
library(devtools)
install_github("marta-vidalgarcia/morpho.tools.GM")
library(morpho.tools.GM)
library(Morpho)
library(Rvcg)
library(magick)
library(ggplot2)
library(vegan)
library(factoextra)
library(Evomorph)


#### 1. Convert files FCSV to TAG ####

# Temporary function to automatically convert fcsv to tag, will add it to morpho.tools.GM

fcsv2tag <- function(input_file, output_file){
  LMs <- as.matrix(read.csv(file = input_file, 
                            sep = ",", skip = 3, header = FALSE)[, 2:4])
  labels_raw <- as.matrix(read.csv(file = input_file, 
                                   sep = ",", skip = 3, header = FALSE)[, 12])
  labels <- substr(labels_raw, 1, nchar(labels_raw)-2)
  
  mat_tag <- c()
  for (i in 1:nrow(LMs)) {
    if (i == nrow(LMs)) {
      point <- paste0(" ", LMs[i, 1], " ", LMs[i, 2], " ", 
                      LMs[i, 3], " ", 1, " ", 1, " ", 1, " ", paste0("\"", labels[i], "\";"))
      mat_tag <- append(mat_tag, point)
    }
    else {
      point <- paste0(" ", LMs[i, 1], " ", LMs[i, 2], " ", 
                      LMs[i, 3], " ", 1, " ", 1, " ", 1, " ", paste0("\"", labels[i], "\""))
      mat_tag <- append(mat_tag, point)
    }
  }
  file_tag <- c("MNI Tag Point File", "Volumes = 1;", paste0("% Volume 1 points converted with morpho.tools.GM::fcsv2tag(). Input file = ", 
                                                             paste0(getwd(), gsub("./", "/", input_file))), "", "Points =", mat_tag)
  writeLines(file_tag, output_file)
}


# Generate the TAG atlases
setwd("./data/atlas/")

dir()
input_string <- dir(pattern = "*.fcsv")
output_string <- gsub(pattern = "Merged.fcsv", replacement = "ATLAS.tag", x = input_string)

for (i in 1:length(input_string)){
  fcsv2tag(input_file = input_string[i], output_file = output_string[i])
}


# New TAG2LM function to accommodate the tag files generated in 3DSlicer, and also give a list of LM labels
# Will add it to morpho.tools.GM soon too

tag2lm_m <- function(file, processed = TRUE) {
  if (processed == TRUE) {
    number_skip <- 5
  }
  else {
    number_skip <- 4
  }
  lm_matrix <- suppressWarnings(read.table(file = file, skip = number_skip, 
                                           sep = " ", header = F))[, 2:4]
  labels <- as.factor(suppressWarnings(read.table(file = file, skip = number_skip, 
                                           sep = " ", header = F))[, 8])
  colnames(lm_matrix) <- c("x", "y", "z")
  return(list(LMs = lm_matrix, labels = labels))
}

# Example on how to import one of the tag files
dir(pattern = "*.tag")

DMSO_LM_ATLAS <- tag2lm_m("./DMSO_LM_ATLAS.tag")

DMSO_LM_ATLAS$labels
DMSO_LM_ATLAS$LMs

setwd("../../")

getwd()


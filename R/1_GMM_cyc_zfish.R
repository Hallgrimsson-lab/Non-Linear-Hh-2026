#### 0. Intro to objects and common functions ####

# check working directory
getwd() # different to terminal $ pwd

# Set working directory
# setwd("/Users/craig.jacobs/Documents/GitRepos/Zfish_GMM/Zfish_Shh_GMM")
# in the terminal $ cd /Users/craig.jacobs/Documents/GitRepos/Zfish_GMM/Zfish_Shh_GMM

# check files in working directory
dir() # in the terminal $ ls

# Check files with specific pattern
dir(pattern = "*.md") # terminal-equivalent of $ ls *.md
dir(pattern = "R")

# Check objects in R Global Environment
ls()


#### 1. Install and load R libraries ####
install.packages("rgl")
library(rgl)
library(geomorph)
library(devtools)
install_github("marta-vidalgarcia/morpho.tools.GM")
library(morpho.tools.GM)
library(Morpho)
library(Rvcg)
library(magick)
library(ggplot2)
install.packages("ggExtra")
library(ggExtra)
library(grid)
library(glue)
library(vegan)
library(dplyr)
library(factoextra)
library(Evomorph)
install.packages("gtools")
library(gtools)
library(patchwork)

# Create directories
getwd() # where are we? We should always be in the project main directory, where the .git file is (check with $ ls -la)

folder_structure <- c("./figs", "./figs/pc_morphs",
                      "./data", "./data/atlas", "./data/Prop_LMs", "./output")

# I am not creating R folder

folder_structure 

# To check structure of object
str(folder_structure)

# Check class of object
class(folder_structure)

dir.exists(folder_structure) # it should all be false unless you created folders manually. I recommend doing it this way

# Little loop to create directories locally
for (i in 1:length(folder_structure)){
  dir.create(folder_structure[i], mode = "0777") # create these folders, equivalent to chmod 777 in bash
}

# Then save this script inside the R folder
# Copy landmark data to the Prop_LMs folder
# Copy atlases inside the data/atlas folder. We will need:
# PLY of each atlas (ascii)
# CURVESLIDE FILES (just 1 that can apply to all atlases)
# ATLAS TAG file for each atlas
# Copy classifiers file inside data and make sure it is a csv and there are no issues

# I want to push these empty directories to github, how do I do that?
# Using .gitkeep files in each empty subdirectory
# I use the function touch in the terminal to create a file (any extension I want), and also hidden files
# $ touch .gitignore # to create gitignore file if I dont have one
# $ touch .gitkeep # needs to be inside each empty folder

# $ find . -type d -empty -not -path "./.git/*" -exec touch {}/.gitkeep \;

#### 2. Import data ####
# 2.1 Classifiers ####

classifiers_no_order <- read.csv(file = "./data/Zfish_Classifiers.csv", header = TRUE, sep = ",")

head(classifiers_no_order)

?setdiff # vectors you are looking at here will be - match the row names with the dimnames (specimen name)
# (1) row.names(classifiers_no_order)
# (2) dimnames(example_array_LMs)[[3]] # [3] is the specimen name in the array (row, column, specimen)

setdiff(row.names(classifiers_no_order), dimnames(head_array)[[3]])
setdiff(dimnames(head_array)[[3]], row.names(classifiers_no_order))

# do both and remove carefully whatever does not match

row.names(classifiers_no_order) <- classifiers_no_order$File.Name
classifiers_no_order$Treatment <- as.factor(classifiers_no_order$Treatment)
str(classifiers_no_order)

row.names(classifiers_no_order)

dimnames(classifiers_no_order)

dimnames(head_array)

head(head_array)

classifiers_no_order$Treatment


# 2.2. Propagated Landmarks ####

# make sure all specimens have the same number of LMs, otherwise the tag2array function will not work

?tag2array

setwd("./data/Prop_LMs/")
dir()

head_array <- tag2array(string_del = "_Landmarks", propagated = TRUE) # removes _Landmarks from the end of the file to match with the classifiers

classifiers_no_order$File.Name

# MAKE SURE THE ORDER OF THE FILES IN THE CSV AS THE SAME AS THE ARRAY - Otherwise the treatments etc will be grouped incorrectly


?tag2lm

dim(head_array)
str(head_array)
dimnames(head_array)[[3]]

head_array[,,116]

setwd("../atlas")
dir()


setwd("../../")

# 3. LM scheme & semiLMs ####
# 3.1 Import data ####
# make sure to run the tag2lm_m function found in 0_Generate_TAG_atlases.R before this step

DMSO_LM_ATLAS <- tag2lm_m("./data/atlas/DMSO_LM_ATLAS.tag")
Cyc20_LM_ATLAS <- tag2lm_m("./data/atlas/Cyc20_LM_ATLAS.tag")
Cyc40_LM_ATLAS <- tag2lm_m("./data/atlas/Cyc40_LM_ATLAS.tag")
Cyc60_LM_ATLAS <- tag2lm_m("./data/atlas/Cyc60_LM_ATLAS.tag")
Cyc80_LM_ATLAS <- tag2lm_m("./data/atlas/Cyc80_LM_ATLAS.tag")
Cyc100_LM_ATLAS <- tag2lm_m("./data/atlas/Cyc100_LM_ATLAS.tag")

length(levels(DMSO_LM_ATLAS$labels))
length(levels(Cyc20_LM_ATLAS$labels))
length(levels(Cyc40_LM_ATLAS$labels))
length(levels(Cyc60_LM_ATLAS$labels))
length(levels(Cyc80_LM_ATLAS$labels))
length(levels(Cyc100_LM_ATLAS$labels))


# Make sure to fix the last label as it comes with ; from TAG format
levels(DMSO_LM_ATLAS$labels)

str(DMSO_LM_ATLAS$labels)
DMSO_LM_ATLAS$labels[which(DMSO_LM_ATLAS$labels == "ventral_R_6;")] <- "ventral_R_6"

DMSO_LM_ATLAS$labels <- as.factor(as.character(DMSO_LM_ATLAS$labels))
levels(DMSO_LM_ATLAS$labels)

# Cyc20 LMs

Cyc20_LM_ATLAS$labels

Cyc20_LM_ATLAS$labels[which(Cyc20_LM_ATLAS$labels == "ventral_R_6;")] <- "ventral_R_6"

Cyc20_LM_ATLAS$labels <- as.factor(as.character(Cyc20_LM_ATLAS$labels))
levels(Cyc20_LM_ATLAS$labels)

# Cyc40 LMs

Cyc40_LM_ATLAS$labels

Cyc40_LM_ATLAS$labels[which(Cyc40_LM_ATLAS$labels == "ventral_R_6;")] <- "ventral_R_6"

Cyc40_LM_ATLAS$labels <- as.factor(as.character(Cyc40_LM_ATLAS$labels))
levels(Cyc40_LM_ATLAS$labels)

# Cyc60 LMs

Cyc60_LM_ATLAS$labels

Cyc60_LM_ATLAS$labels[which(Cyc60_LM_ATLAS$labels == "ventral_R_6;")] <- "ventral_R_6"

Cyc60_LM_ATLAS$labels <- as.factor(as.character(Cyc60_LM_ATLAS$labels))
levels(Cyc60_LM_ATLAS$labels)

# Cyc80 LMs

Cyc80_LM_ATLAS$labels

Cyc80_LM_ATLAS$labels[which(Cyc80_LM_ATLAS$labels == "ventral_R_6;")] <- "ventral_R_6"

Cyc80_LM_ATLAS$labels <- as.factor(as.character(Cyc80_LM_ATLAS$labels))
levels(Cyc80_LM_ATLAS$labels)

# Cyc100 LMs

Cyc100_LM_ATLAS$labels

Cyc100_LM_ATLAS$labels[which(Cyc100_LM_ATLAS$labels == "ventral_R_6;")] <- "ventral_R_6"

Cyc100_LM_ATLAS$labels <- as.factor(as.character(Cyc100_LM_ATLAS$labels))
levels(Cyc100_LM_ATLAS$labels)


# Need to make sure Meshes are in ASCII format and not binary
# $ meshlabserver # or open directly on your computer
# Jaron working on new meshes

Zfish_Atlas_DMSO_atlas_mesh <- file2mesh("./data/atlas/Zfish_Atlas_DMSO_ascii.ply")
Zfish_Atlas_Cyc20_atlas_mesh <- file2mesh("./data/atlas/Zfish_Atlas_Cyc20_ascii.ply")
Zfish_Atlas_Cyc40_atlas_mesh <- file2mesh("./data/atlas/Zfish_Atlas_Cyc40_ascii.ply")
Zfish_Atlas_Cyc60_atlas_mesh <- file2mesh("./data/atlas/Zfish_Atlas_Cyc60_ascii.ply")
Zfish_Atlas_Cyc80_atlas_mesh <- file2mesh("./data/atlas/Zfish_Atlas_Cyc80_ascii.ply")
Zfish_Atlas_Cyc100_atlas_mesh <- file2mesh("./data/atlas/Zfish_Atlas_Cyc100_ascii.ply")

# Make sure that there is no internal mesh as that will affect morphs and heatmaps.
# Better to play with density thresholding so that you do not have internal 'sub-meshes' that will affect morphing later
# work on DMSO, worry about severe morph meshes later
# Cut neck too.
# Craig & Jaron will do that

# 3.2. Plotting LM schemes ####

# create positions for easy viewing of each atlas
# new meshes need to be generated and uploaded (Jaron working on it)

# DMSO

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700)) 
shade3d(Zfish_Atlas_DMSO_atlas_mesh, color = "white", alpha = 1, specular = 1) # alpha is transparency %

lateral <- par3d()$userMatrix # navigate mesh display to lateral view then run this
dorsal <- par3d()$userMatrix # dorsal view
ventral <- par3d()$userMatrix # ventral view
transverse <- par3d()$userMatrix # transverse view
save(lateral, dorsal, ventral, transverse, file = "./data/atlas/RGL_positions_DMSO_atlas.Rdata")

# once saved these can be loaded without needing to be recreated

load("./data/atlas/RGL_positions_DMSO_atlas.Rdata")

# save in all views

morph_views <- list(dorsal = dorsal, lateral = lateral, ventral = ventral, transverse = transverse)

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(Zfish_Atlas_DMSO_atlas_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("DMSO_atlas_", view, ".png")
  rgl.snapshot(file.path("./figs", output_file), fmt = "png", top = TRUE)
  close3d()
}

rgl.close() # closes RGL window

# Cyc20

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = ventral) 
shade3d(Zfish_Atlas_Cyc20_atlas_mesh, color = "white", alpha = 1, specular = 1) # alpha is transparency %

lateral <- par3d()$userMatrix # lateral view
dorsal <- par3d()$userMatrix # dorsal view
ventral <- par3d()$userMatrix # ventral view
transverse <- par3d()$userMatrix # transverse view
save(lateral, dorsal, ventral, transverse, file = "./data/atlas/RGL_positions_Cyc20_atlas.Rdata")

# once saved these can be loaded without needing to be recreated

load("./data/atlas/RGL_positions_Cyc20_atlas.Rdata")

rgl.close() # closes RGL window

# Cyc40

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700)) 
shade3d(Zfish_Atlas_Cyc40_atlas_mesh, color = "gray", alpha = 0.75) # alpha is transparency %

lateral <- par3d()$userMatrix # lateral view
dorsal <- par3d()$userMatrix # dorsal view
ventral <- par3d()$userMatrix # ventral view
transverse <- par3d()$userMatrix # transverse view
save(lateral, dorsal, ventral, transverse, file = "./data/atlas/RGL_positions_Cyc40_atlas.Rdata")

# once saved these can be loaded without needing to be recreated

load("./data/atlas/RGL_positions_Cyc40_atlas.Rdata")

rgl.close() # closes RGL window

# Cyc60

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700)) 
shade3d(Zfish_Atlas_Cyc60_atlas_mesh, color = "gray", alpha = 0.75) # alpha is transparency %

lateral <- par3d()$userMatrix # lateral view
dorsal <- par3d()$userMatrix # dorsal view
ventral <- par3d()$userMatrix # ventral view
transverse <- par3d()$userMatrix # transverse view
save(lateral, dorsal, ventral, transverse, file = "./data/atlas/RGL_positions_Cyc60_atlas.Rdata")

# once saved these can be loaded without needing to be recreated

load("./data/atlas/RGL_positions_Cyc60_atlas.Rdata")

rgl.close() # closes RGL window

# Cyc80

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal) 
shade3d(Zfish_Atlas_Cyc80_atlas_mesh, color = "gray", alpha = 0.75) # alpha is transparency %

lateral <- par3d()$userMatrix # lateral view
dorsal <- par3d()$userMatrix # dorsal view
ventral <- par3d()$userMatrix # ventral view
transverse <- par3d()$userMatrix # transverse view
save(lateral, dorsal, ventral, transverse, file = "./data/atlas/RGL_positions_Cyc80_atlas.Rdata")

# once saved these can be loaded without needing to be recreated

load("./data/atlas/RGL_positions_Cyc80_atlas.Rdata")

rgl.close() # closes RGL window

# Cyc100

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal) 
shade3d(Zfish_Atlas_Cyc100_atlas_mesh, color = "gray", alpha = 0.75) # alpha is transparency %

lateral <- par3d()$userMatrix # lateral view
dorsal <- par3d()$userMatrix # dorsal view
ventral <- par3d()$userMatrix # ventral view
transverse <- par3d()$userMatrix # transverse view
save(lateral, dorsal, ventral, transverse, file = "./data/atlas/RGL_positions_Cyc100_atlas.Rdata")

# once saved these can be loaded without needing to be recreated

load("./data/atlas/RGL_positions_Cyc100_atlas.Rdata")

rgl.close() # closes RGL window

# Generate vectors of positions of LM, curve semiLM and surface semiLM. Can use any atlas as LMs are the same

# Vector of positions of LMs
vec_pos_LMs <- which(Cyc20_LM_ATLAS$labels == "Homologous_LM")

# Vector of positions of curve semiLMs
# vec_pos_curves <- which(Cyc20_LM_ATLAS$labels == "mouth_curve_1" | Cyc20_LM_ATLAS$labels == "mouth_curve_2" | 
                        #  Cyc20_LM_ATLAS$labels == "ventral_curve" | Cyc20_LM_ATLAS$labels == "forebrain_midline_curve_1")

vec_pos_curves <- c(which(Cyc20_LM_ATLAS$labels == "mouth_curve_1"), which(Cyc20_LM_ATLAS$labels == "mouth_curve_2"),
                    which(Cyc20_LM_ATLAS$labels == "ventral_curve"), which(Cyc20_LM_ATLAS$labels == "forebrain_midline_curve_1"))

# Vector of positions of surface semiLMs

vec_pos_surface_semis <- c(1:length(Cyc20_LM_ATLAS$labels))[-c(vec_pos_LMs, vec_pos_curves)]


# Generating LM scheme figures

# load appropriate userMatrix so positions are correct
load("./data/atlas/RGL_positions_DMSO_atlas.Rdata")

# open mesh
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal) 
shade3d(Zfish_Atlas_DMSO_atlas_mesh, color = "white", alpha = 1, specular = 1)

# plot homologous LMs
plot3d(DMSO_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
# plot curve semi LMs
plot3d(DMSO_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
# plot surface semi LMs
plot3d(DMSO_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
# save image
rgl.snapshot("./figs/ATLAS_LM_scheme_DMSO_dorsal.png", top = TRUE)
close3d()

# other views
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral) 
shade3d(Zfish_Atlas_DMSO_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(DMSO_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(DMSO_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(DMSO_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_DMSO_lateral.png", top = TRUE)
close3d()
# now do it for other example atlases

load("./data/atlas/RGL_positions_Cyc20_atlas.Rdata")

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
shade3d(Zfish_Atlas_Cyc20_atlas_mesh, color = "white", alpha = 1, specular = 1) # alpha is transparency %

# Plot homologous landmarks
plot3d(Cyc20_LM_ATLAS$LMs[which(Cyc20_LM_ATLAS$labels == "Homologous_LM"),], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)

# Plot curve semiLMs
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)

# Plot surface semiLMs
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)

# Save image of LM scheme
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc20_lateral.png", top = TRUE)


rgl.close()

# Now plot it on dorsal view
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal)
shade3d(Zfish_Atlas_Cyc20_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc20_dorsal.png", top = TRUE)

rgl.close()

# Now plot it on ventral view
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = ventral)
shade3d(Zfish_Atlas_Cyc20_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc20_ventral.png", top = TRUE)

rgl.close()

# Now plot it on transverse view 
open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = transverse)
shade3d(Zfish_Atlas_Cyc20_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc20_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc20_transverse.png", top = TRUE)

rgl.close()

# Generate figure for other atlases (only doing Cyc60 and Cyc100 to demonstrate how the scheme was adapted to severe phenotypes)

# Cyc60
load("./data/atlas/RGL_positions_Cyc60_atlas.Rdata")

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
shade3d(Zfish_Atlas_Cyc60_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc60_lateral.png", top = TRUE)

rgl.close()

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal)
shade3d(Zfish_Atlas_Cyc60_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc60_dorsal.png", top = TRUE)

rgl.close()

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = ventral)
shade3d(Zfish_Atlas_Cyc60_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc60_ventral.png", top = TRUE)

rgl.close()

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = transverse)
shade3d(Zfish_Atlas_Cyc60_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc60_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc60_transverse.png", top = TRUE)

rgl.close()

# Cyc100
load("./data/atlas/RGL_positions_Cyc100_atlas.Rdata")

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = lateral)
shade3d(Zfish_Atlas_Cyc100_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc100_lateral.png", top = TRUE)

rgl.close()

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = dorsal)
shade3d(Zfish_Atlas_Cyc100_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc100_dorsal.png", top = TRUE)

rgl.close()

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = ventral)
shade3d(Zfish_Atlas_Cyc100_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc100_ventral.png", top = TRUE)

rgl.close()

open3d(zoom = 0.75, windowRect = c(0, 0, 1000, 700), userMatrix = transverse)
shade3d(Zfish_Atlas_Cyc100_atlas_mesh, color = "white", alpha = 1, specular = 1)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_LMs,], aspect = "iso", type = "s", size = 1.75,
       col ="darkblue", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_curves,], aspect = "iso", type = "s", size = 0.75,
       col ="red", add = TRUE)
plot3d(Cyc100_LM_ATLAS$LMs[vec_pos_surface_semis,], aspect = "iso", type = "s", size = 0.5,
       col ="magenta", add = TRUE)
rgl.snapshot("./figs/ATLAS_LM_scheme_Cyc100_transverse.png", top = TRUE)

rgl.close()


# 3.3. curve semiLMs - curveslide files ####
vec_pos_LMs
mouth_curve_1_pos <- which(Cyc20_LM_ATLAS$labels == "mouth_curve_1") # sliding between LM83 and LM81
mouth_curve_2_pos <- which(Cyc20_LM_ATLAS$labels == "mouth_curve_2") # sliding between LM83 and LM82
ventral_curve_pos <- which(Cyc20_LM_ATLAS$labels == "ventral_curve") # sliding between LM84 and LM85
forebrain_midline_curve_1_pos <- which(Cyc20_LM_ATLAS$labels == "forebrain_midline_curve_1") # sliding between LM80 and LM77


# We need to create a curveslide matrix to know where from to where to the semis slide
curveslide_mouth_curve_1_pos <- c(83, mouth_curve_1_pos, 81)
curveslide_mouth_curve_2_pos <- c(83, mouth_curve_2_pos, 82)

ventral_curve_pos_left <- c(84, ventral_curve_pos[1:c(length(ventral_curve_pos)-1)])
ventral_curve_pos_right <- c(ventral_curve_pos[2:c(length(ventral_curve_pos))], 85)
curveslide_ventral_curve_pos <- cbind(ventral_curve_pos_left, ventral_curve_pos, ventral_curve_pos_right)

forebrain_midline_curve_1_pos_left <- c(80, forebrain_midline_curve_1_pos[1:c(length(forebrain_midline_curve_1_pos)-1)])
forebrain_midline_curve_1_pos_right <- c(forebrain_midline_curve_1_pos[2:c(length(forebrain_midline_curve_1_pos))], 77)
curveslide_forebrain_midline_curve_1_pos <- cbind(forebrain_midline_curve_1_pos_left, forebrain_midline_curve_1_pos, forebrain_midline_curve_1_pos_right)


ls(pattern = "curveslide*") # Search within Global Environment

curveslide_list <- lapply(ls(pattern = "curveslide*"), get)
str(curveslide_list)
curveslide_all <- do.call(rbind, curveslide_list)
colnames(curveslide_all) <- c("before", "sliding", "after")

curveslide_all

write.csv(curveslide_all, "./data/atlas/curveslide_semis_zfish.csv", row.names = FALSE)

load("./data/atlas/curveslide_semis_zfish.csv")




# 4. Generalised Procrustes Analysis (GPA) ####

# The reason why we made the curveslide file is to have the info on how curve semis slide during GPA
getwd()

setwd("./data/Prop_LMs/")
?tag2array

dir()

# we previously made an array of all the propLMs
head_array

dim(head_array)

dim(classifiers_no_order)
dimnames(classifiers_no_order)
classifiers_no_order$Treatment
str(classifiers_no_order$Treatment)

head_array <- tag2array(string_del = "_Landmarks", propagated = TRUE)

head_array

head(head_array)

# trying to create arrays for each treatment - NEED TO ANALYSE WITHIN THE SPECIFIC GROUPS OR DISCARD CERTAIN TREATMENTS FROM GPA ANALYSIS

dimnames(head_array)

classifiers_no_order$Treatment


#subset_indices <- which(classifiers_no_order$Treatment == "DMSO")
#head_array[subset_indices]

#tag_files <- list.files(pattern = "\\.tag$")
#dmso_tag <- tag_files[grep("_DMSO_", tag_files)]

#DMSO_array <- tag2array(ID = dmso_tag, string_del = "_Landmarks", propagated = TRUE)


?tag2array

?subset.data.frame()

# Make sure to run/load curveslide and vector code above before trying to run the GPA

?gpagen

gpagen

GPA_All <- gpagen(A = head_array, curves = curveslide_all, surfaces = vec_pos_surface_semis) # most simple way to run a GPA 

GPA_All$coords
GPA_All$Csize

setwd("../../")
getwd()

# If we wanted to include symmetry (for later)
?Morpho::procSym


# Plot outliers (good to show distribution of data)


?morpho.tools.GM::plotOutliers_percentile # Use shape coordinates
?plotOutliers


plotOutliers_percentile(GPA_All$coords, groups = NULL, inspect.outliers = FALSE, percentile = 0.95, save.plot = FALSE)

plotOutliers_percentile(GPA_All$coords, groups = classifiers_no_order$Treatment, inspect.outliers = FALSE, percentile = 0.95, save.plot = FALSE)

levels(classifiers_no_order$Treatment)

outliers_95 <- plotOutliers_percentile(GPA_All$coords, groups = NULL, inspect.outliers = FALSE, percentile = 0.95, save.plot = FALSE)

outliers_95$Proc_d_percentile

outliers_99 <- plotOutliers_percentile(GPA_All$coords, groups = NULL, inspect.outliers = FALSE, percentile = 0.99, save.plot = FALSE)

outliers_90 <- plotOutliers_percentile(GPA_All$coords, groups = NULL, inspect.outliers = FALSE, percentile = 0.90, save.plot = FALSE)

outliers_90$Proc_d_percentile

outliers_85 <- plotOutliers_percentile(GPA_All$coords, groups = NULL, inspect.outliers = FALSE, percentile = 0.85, save.plot = FALSE)
outliers_85$Proc_d_percentile

outliers_80 <- plotOutliers_percentile(GPA_All$coords, groups = NULL, inspect.outliers = FALSE, percentile = 0.80, save.plot = FALSE)
outliers_80$Proc_d_percentile

# Identified outliers 


# Then get rid of that specimen (s) in both ARRAY & classifiers files
# Only remove specimen if you are SURE it is an outlier (I did not remove any here, as the "outliers" were all the in the high concentration group - a phenotype)


# 5. Principle Component Analysis ####

PCA_head <- gm.prcomp(A = GPA_All$coords)

summary(PCA_head)

class(PCA_head)

PCA_var_comp <- PCA_head
class(PCA_var_comp) <- "princomp"

pdf("./figs/PCA_skull_shape_scree_plot.pdf", height = 5, width = 5)
PCA_scree <- fviz_eig(PCA_var_comp, addlabels = TRUE, barfill = "darkgrey",
                      barcolor = "black", linecolor = "blue") + theme_classic()

PCA_scree


#PCA_Cyc60 <- gm.prcomp(A = Cyc60_coords)

#summary(PCA_Cyc60)
#class(PCA_Cyc60)

#PCA_var_comp_Cyc60 <- PCA_Cyc60
#class(PCA_var_comp_Cyc60) <- "princomp"

#PCA_scree_Cyc60 <- fviz_eig(PCA_var_comp_Cyc60, addlabels = TRUE, barfill = "darkgrey",
 #                     barcolor = "black", linecolor = "blue") + theme_classic()

#PCA_scree_Cyc60

# 5.1 PCA by treatment ####

install.packages("ggordiplots")
library(ggordiplots)

str(classifiers_no_order$Treatment)

palette()
levels (classifiers_no_order$Treatment)
palette(c("black", "blue", "purple", "pink", "orange", "red"))

classifiers_no_order$Treatment <- factor(classifiers_no_order$Treatment, levels = c("DMSO", "Cyc20", "Cyc40", "Cyc60", "Cyc80", "Cyc100"))

# PCA plot of PC1 ~ PC2

plot(PCA_head, axis1 = 1, axis2 = 2, col = classifiers_no_order$Treatment, pch = 19, cex = 1.25)
ordiellipse(PCA_head, classifiers_no_order$Treatment, choices = c(1,2), kind = "sd", conf = 0.95, col = palette(), 
            draw = "polygon", alpha = 0.15, lty = 0)
legend("topright", pch = 19, cex = 1.5, col = palette(), legend = levels(classifiers_no_order$Treatment))
title("Head Shape - PC1 ~ PC2")

# PCA plots of different PCs

plot(PCA_head, axis1 = 2, axis2 = 3, col = classifiers_no_order$Treatment, pch = 19, cex = 1.25)
ordiellipse(PCA_head, classifiers_no_order$Treatment, choices = c(2,3), kind = "sd", conf = 0.95, col = palette(), 
            draw = "polygon", alpha = 0.15, lty = 0)
legend("topright", pch = 19, cex = 1.5, col = palette(), legend = levels(classifiers_no_order$Treatment))
title("Head Shape - PC2 ~ PC3")

plot(PCA_head, axis1 = 1, axis2 = 3, col = classifiers_no_order$Treatment, pch = 19, cex = 1.25)
ordiellipse(PCA_head, classifiers_no_order$Treatment, choices = c(1,3), kind = "sd", conf = 0.95, col = palette(), 
            draw = "polygon", alpha = 0.15, lty = 0)
legend("topright", pch = 19, cex = 1.5, col = palette(), legend = levels(classifiers_no_order$Treatment))
title("Head Shape - PC1 ~ PC3")

plot(PCA_head, axis1 = 1, axis2 = 4, col = classifiers_no_order$Treatment, pch = 19, cex = 1.25)
ordiellipse(PCA_head, classifiers_no_order$Treatment, choices = c(1,4), kind = "sd", conf = 0.95, col = palette(), 
            draw = "polygon", alpha = 0.15, lty = 0)
legend("topright", pch = 19, cex = 1.5, col = palette(), legend = levels(classifiers_no_order$Treatment))
title("Head Shape - PC1 ~ PC4")

plot(PCA_head, axis1 = 3, axis2 = 4, col = classifiers_no_order$Treatment, pch = 19, cex = 1.25)
ordiellipse(PCA_head, classifiers_no_order$Treatment,  choices = c(3,4), kind = "sd", conf = 0.95, col = palette(), 
            draw = "polygon", alpha = 0.15, lty = 0)
legend("topright", pch = 19, cex = 1.5, col = palette(), legend = levels(classifiers_no_order$Treatment))
title("Head Shape - PC3 ~ PC4")

# 5.1.2 PCA plots with marginal densities ####

# translated from code from Aponte

library(ggplot2)
library(grid)
install.packages("ggExtra")
library(ggExtra) #lets you easily plot in the margins
library(glue) #lets you create strings in a pretty handy way

# put PCA into dataframe, make sure columns are PC1, PC2 etc and the data are grouped correctly into treatments

PCA_df <- as.data.frame(PCA_head$x)
  colnames(PCA_df) <- paste0("PC", 1:ncol(PCA_df))
  groups <- classifiers_no_order$Treatment
  PCA_df$Group <- groups

head(PCA_df)

# Calculate the variance explained by each principal component to use for axis labels
pca_variance <- (PCA_head$sdev^2) / sum(PCA_head$sdev^2) * 100

# Create the axis labels with variance explained
x_label <- paste0("PC1 (", round(pca_variance[1], 2), "%)")
y_label <- paste0("PC2 (", round(pca_variance[2], 2), "%)")

# Calculate maximum limit based on data range for symmetry (with padding). I don't like the PC1 symmetrical, but PC2 looks good with the increased padding
max_limitPC1 <- max(abs(c(PCA_df$PC1))) * 1.2
max_limitPC2 <- max(abs(c(PCA_df$PC2))) * 2

# create vector for PCA plot

PC1_PC2_plot <- ggplot(PCA_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = palette()) +
  scale_fill_manual(values = palette()) +
  stat_ellipse(aes(fill = Group), type = "norm", level = 0.95, geom = "polygon", alpha = 0.15, color = NA) +
  scale_y_continuous(limits = c(-max_limitPC2, max_limitPC2)) +
  theme_minimal() +
  labs(x = x_label, y = y_label, color = "Treatment") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", size = 0.2),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.ticks.length = unit(5, "pt")
        ) +
  guides(fill = "none")

# Add marginal density plots using newly created vector

ggMarginal(PC1_PC2_plot, type = "density", groupColour = TRUE, groupFill = TRUE)

# large plot for poster

PC1_PC2_plotL <- ggplot(PCA_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 5) +
  scale_color_manual(values = palette()) +
  scale_fill_manual(values = palette()) +
  stat_ellipse(aes(fill = Group), type = "norm", level = 0.95, geom = "polygon", alpha = 0.15, color = NA) +
  scale_y_continuous(limits = c(-max_limitPC2, max_limitPC2)) +
  theme_minimal() +
  labs(x = x_label, y = y_label, color = "Treatment") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 28),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.key.height = unit(32, "pt"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 1),
        axis.ticks.length = unit(5, "pt"),
        axis.title = element_text(size = 26),
        axis.text = element_text(size = 20)
  ) +
  guides(fill = "none")

ggMarginal(PC1_PC2_plotL, type = "density", groupColour = TRUE, groupFill = TRUE)


# Now what if you want to generate/save out PC1/2, PC3/4, PC5/6 without always manually changing the code?
# get is some janky ass R function that will let you construct a string to call a variable. Check it out

pcPairs <- matrix(1:6, ncol = 2, byrow = T)

for(i in 1:nrow(pcPairs)){
  
  p <- ggplot(data, aes(x = get(glue("PC{pcPairs[i,1]}")), y = get(glue("PC{pcPairs[i,2]}")), color = group, alpha = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Group 1" = "#3F51B5", "Group 2" = "#F44336", "Group 3" = "darkgreen")) +
    scale_alpha_manual(values = c("Group 1" = 0.5, "Group 2" = .5, "Group 3" = .5)) +
    theme_minimal(base_line_size = 0) +
    labs(x = glue("PC{pcPairs[i,1]}"), y = glue("PC{pcPairs[i,2]}")) +
    theme(legend.position = "none")
  
  # Add marginal density plots and save to pdf
  pdf(file = glue("PC{pcPairs[i,1]}_{pcPairs[i,2]}.pdf"), height = 5, width = 5)
  print(ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE))
  dev.off()
  
}

# 5.1.3 PC Scores plots and curve ####

# PC scores box plots

ggplot(PCA_df, aes(x= Group, y = PC1, fill = Group)) +
  scale_fill_manual(values = palette()) +
  geom_boxplot(alpha = 0.75) + 
  theme_classic() + geom_jitter(width = 0.2, size = 1.25) +
  labs(x = "Treatment", y = "PC1 Score", title = "PC1 Scores by Hh Inhibition") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
        axis.title = element_text(size = 26), axis.text = element_text(size = 20),
        axis.line = element_line(size = 1), axis.ticks.length = unit(5, "pt")
  )

PC1_anova_result <- aov(PC1 ~ Group, data = PCA_df)
summary(PC1_anova_result)

TukeyHSD(PC1_anova_result)

# testing concentration ~ PC1 score linearity


# make concentration numeric
dose_map <- c("DMSO" = 0, "Cyc20" = 20, "Cyc40" = 40, "Cyc60" = 60, "Cyc80" = 80, "Cyc100" = 100)
PCA_df$conc <- dose_map[PCA_df$Group]
PCA_summary$conc <- dose_map[PCA_summary$Group]

# Fit linear model to group means
lin_mod <- lm(mean_PC1 ~ conc, data = PCA_summary)

# Fit quadratic model to check for nonlinearity
quad_mod <- lm(mean_PC1 ~ poly(conc, 2), data = PCA_summary)

# Extract goodness of fit
lin_r2 <- summary(lin_mod)$r.squared
quad_r2 <- summary(quad_mod)$r.squared
delta_AIC <- AIC(lin_mod) - AIC(quad_mod)
LinearityDeparture <- anova(lin_mod, quad_mod)

# RMSE-based linearity score (0 = not linear, 1 = perfectly linear)
rmse <- sqrt(mean(residuals(lin_mod)^2))
lin_score <- 1 - rmse / sd(PCA_summary$mean_PC1)

list(
  Linear_R2 = lin_r2,
  Quadratic_R2 = quad_r2,
  Delta_AIC = delta_AIC,
  Linearity_Score_0to1 = lin_score,
  Depature_from_Linearity = LinearityDeparture)

# using raw values, not just means - USE THIS IT IS MORE REPRESENTATIVE
lin_mod_raw  <- lm(PC1 ~ conc, data = PCA_df)
quad_mod_raw <- lm(PC1 ~ conc + poly(conc, 2), data = PCA_df)

# Extract goodness of fit
lin_r2_raw <- summary(lin_mod_raw)$r.squared
quad_r2_raw <- summary(quad_mod_raw)$r.squared
delta_AIC_raw <- AIC(lin_mod_raw) - AIC(quad_mod_raw)
LinearityDeparture_raw <- anova(lin_mod_raw, quad_mod_raw)

rmse <- sqrt(mean(residuals(lin_mod_raw)^2))
lin_score_raw <- 1 - rmse / sd(PCA_summary$mean_PC1)

linearity_tests <- list(
                    Linear_R2 = lin_r2_raw,
                    Quadratic_R2 = quad_r2_raw,
                    Delta_AIC = delta_AIC_raw,
                    Linearity_Score_0to1 = lin_score_raw,
                    Depature_from_Linearity = LinearityDeparture_raw)
capture.output(linearity_tests, file = "linearity_tests.txt")


# plot linear vs quadratic fit

# assign palette

my_palette <- c(
  "DMSO"   = "black",
  "Cyc20"  = "blue",
  "Cyc40"  = "purple",
  "Cyc60"  = "pink",
  "Cyc80"  = "orange",
  "Cyc100" = "red"
)

# Prediction grid on numeric dose column
xgrid <- data.frame(conc = seq(min(PCA_summary$conc), max(PCA_summary$conc), length.out = 200))
xgrid$linear_fit <- predict(lin_mod_raw, newdata = xgrid)
xgrid$quad_fit   <- predict(quad_mod_raw, newdata = xgrid)

# Labels
lin_label  <- paste0("Linear R² = ", round(lin_r2_raw, 3))
quad_label <- paste0("Quadratic R² = ", round(quad_r2_raw, 3))
depart_label <- paste0("p < 0.0001")

# Breaks = actual doses (0,20,40,60,80,100)
dose_ticks <- sort(unique(PCA_summary$conc))

ggplot() +
  # Raw scatter
  geom_jitter(data = PCA_df,
              aes(x = conc, y = PC1, colour = Group),
              width = 1, alpha = 0.6, size = 2, shape = 16) +
  # means 
  geom_segment(
    data = PCA_summary,
    aes(x = conc - 2, xend = conc + 2, y = mean_PC1, yend = mean_PC1, colour = Group),
    size = 2) +
  
  # Linear fit (dashed)
  geom_line(data = xgrid, aes(x = conc, y = linear_fit),
            linetype = "dashed", size = 1.5, colour = "darkgrey") +
  
  # Quadratic fit (solid)
  geom_line(data = xgrid, aes(x = conc, y = quad_fit),
            size = 2, colour = "darkblue") +
  labs(
    title = "Linear vs Quadratic Fit to PC1 Scores",
    x = "Hh Inhibitor Concentration (µM)",
    y = "PC1 Score"
  ) +
  scale_x_continuous(breaks = dose_ticks, labels = dose_ticks, expand = c(0,0)) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_blank(), axis.text = element_text(size = 20),
        axis.line = element_line(color = "black", size = 2), axis.ticks.length = unit(5, "pt"),
        legend.position = "none") +
  scale_colour_manual(values = my_palette) +
  scale_fill_manual(values = my_palette)

# PC1 drawn as a curve

# Compute summary statistics (mean and standard deviation)
PCA_summary <- PCA_df %>%
  group_by(Group) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1 = sd(PC1, na.rm = TRUE),
    n = n(),
    se_PC1 = sd_PC1/sqrt(n),
    ci_lower = mean_PC1-1.96 * se_PC1,
    ci_upper = mean_PC1+1.96 * se_PC1
  )

PCA_summary$Group <- factor(PCA_summary$Group, levels = c("DMSO", "Cyc20", "Cyc40", "Cyc60", "Cyc80", "Cyc100"))

# Create the plot with SD shaded
ggplot(PCA_summary, aes(x = Group, y = mean_PC1, group = 1)) +
  geom_ribbon(aes(ymin = mean_PC1 - 2*sd_PC1, ymax = mean_PC1 + 2*sd_PC1), 
              fill = "blue", alpha = 0.2) +  # Shaded area for ±1 SD
  geom_line(color = "blue", size = 1) +     # Line through means
  geom_point(color = "red", size = 3) +     # Points for mean values
  labs(x = "Treatment", y = "Mean PC1 Score", 
       title = "PC1 Scores by Treatment with ±2 SD") +
  scale_y_reverse() +
  theme_minimal()

# SE shaded

ggplot(PCA_summary, aes(x = Group, y = mean_PC1, group = 1)) +
  geom_ribbon(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), 
              fill = "blue", alpha = 0.2) +  # Shaded area for ±1 SD
  geom_line(color = "blue", size = 1) +     # Line through means
  geom_point(color = "red", size = 3) +     # Points for mean values
  labs(x = "Treatment", y = "Mean PC1 Score", 
       title = "PC1 Scores by Treatment with ±1 SE") +
  theme_minimal()

# 95% confidence intervals shaded 

ggplot(PCA_summary, aes(x = Group, y = mean_PC1, group = 1)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), 
              fill = "blue", alpha = 0.2) +  # Shaded area for ±1 SD
  geom_line(color = "blue", size = 1) +     # Line through means
  geom_point(color = "red", size = 3) +     # Points for mean values
  labs(x = "Treatment", y = "Mean PC1 Score", 
       title = "PC1 Scores by Treatment with 95% CI") +
  theme_minimal()

# trying out smoothing of curve

PCA_summary <- PCA_df %>%
  group_by(Group) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1 = sd(PC1, na.rm = TRUE),
  ) %>%
  mutate(Dose = as.numeric(Group))  # Ensure Dose is numeric

# Plot with smoothing
ggplot(PCA_summary, aes(x = Dose, y = mean_PC1)) +
  geom_ribbon(aes(ymin = mean_PC1 - sd_PC1, ymax = mean_PC1 + sd_PC1), 
              fill = "blue", alpha = 0.2) +  # Shaded area for ±1 SD
  geom_smooth(method = "loess", color = "blue", se = FALSE) +  # Smoothed trend line
  geom_point(color = "red", size = 3) +  # Points for mean values
  labs(x = "Level of Inhibition", y = "Mean PC1 Score", 
       title = "PC1 Scores by Hh Inhibition with Smoothed Curve") +
  scale_y_reverse() +
  theme_minimal()

# smoothing of curve and shaded region

PCA_summary <- PCA_df %>%
  group_by(Group) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1 = sd(PC1, na.rm = TRUE)
  ) %>%
  mutate(Dose = as.numeric(Group))  # Ensure Dose is numeric

# Create fine-grid dose values for smooth interpolation
dose_grid <- data.frame(Dose = seq(min(PCA_summary$Dose), max(PCA_summary$Dose), length.out = 100))

# Fit LOESS models for smooth interpolation
loess_mean <- loess(mean_PC1 ~ Dose, data = PCA_summary, span = 0.7)
loess_upper <- loess((mean_PC1 + 3*sd_PC1) ~ Dose, data = PCA_summary, span = 0.7)
loess_lower <- loess((mean_PC1 - 3*sd_PC1) ~ Dose, data = PCA_summary, span = 0.7)

# Predict smoothed values
dose_grid$mean_PC1 <- predict(loess_mean, newdata = dose_grid)
dose_grid$upper <- predict(loess_upper, newdata = dose_grid)
dose_grid$lower <- predict(loess_lower, newdata = dose_grid)

# Plot with smoothed shaded region
ggplot(dose_grid, aes(x = Dose, y = mean_PC1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#80CBC4", alpha = 0.4) +  # Smoothed shaded region
  geom_line(color = "#00695C", size = 2) +  # Smoothed mean curve
  #geom_point(data = PCA_summary, aes(x = Dose, y = mean_PC1), color = "black", size = 3) +  # Original points
  labs(x = "Hh inhibition", y = "Mean PC1 Score", 
       title = "PC1 Scores by Hh Inhibition ±3 SD") +
  scale_y_reverse() +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
      legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "black", size = 2),
      axis.ticks = element_blank(),
      #axis.text = element_blank()
      )

# fit linear model to the same LOESS curve

ggplot(dose_grid, aes(x = Dose, y = mean_PC1)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#80CBC4", alpha = 0.4) +
  geom_line(color = "#00695C", size = 2) +
  
  # Add linear model line using the mean values
  geom_smooth(
    data = PCA_summary, aes(x = Dose, y = mean_PC1),
    method = "lm", se = FALSE, linetype = "dashed", size = 1.5, colour = "darkgrey"
  ) +
  
  labs(x = "Hh inhibition", y = "Mean PC1 Score",
       title = "PC1 Scores by Hh Inhibition ±3 SD") +
  scale_y_reverse() +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(), legend.text = element_text(size = 12),
    legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(color = "black", size = 2),
    axis.ticks = element_blank()
  )

# 5.2 Box plots - Csize, Pdist, Trace ####

class(GPA_All)

# define Pdist

Pdist <- ShapeDist(shapes = GPA_All$coords, reference = GPA_All$consensus)

# create data frame so data can be subset into treatments etc

?geomorph.data.frame

gdf_head <- geomorph.data.frame(GPA_All, Treatment = classifiers_no_order$Treatment, Pdist = Pdist)

class(gdf_head)

# Subset data frame into separate treatments 

gdf_head

gdf_head$Treatment
gdf_head$coords

class(gdf_head$coords)

DMSO_fish <- which(gdf_head$Treatment == "DMSO")
DMSO_fish # make sure this lines up with correct specimens in the classifier

DMSO_coords <- gdf_head$coords[,,DMSO_fish]

DMSO_coords # should be all coords of DMSO treated embryos

head(DMSO_coords)
tail(DMSO_coords)

class(DMSO_coords)

DMSO_coords

# calculate mean shape of DMSO fish to use as a reference later
DMSO_mean_coords <- mshape(DMSO_coords)

# calculate mean shape of other treatment group. This can probably be turned into a loop, but I am not proficient enough

Cyc20_fish <- which(gdf_head$Treatment == "Cyc20")
Cyc20_fish # make sure this lines up with correct specimens in the classifier

Cyc20_coords <- gdf_head$coords[,,Cyc20_fish]
Cyc20_mean_coords <- mshape(Cyc20_coords)

Cyc40_fish <- which(gdf_head$Treatment == "Cyc40")
Cyc40_fish # make sure this lines up with correct specimens in the classifier

Cyc40_coords <- gdf_head$coords[,,Cyc40_fish]
Cyc40_mean_coords <- mshape(Cyc40_coords)

Cyc60_fish <- which(gdf_head$Treatment == "Cyc60")
Cyc60_fish # make sure this lines up with correct specimens in the classifier

Cyc60_coords <- gdf_head$coords[,,Cyc60_fish]
Cyc60_mean_coords <- mshape(Cyc60_coords)

Cyc80_fish <- which(gdf_head$Treatment == "Cyc80")
Cyc80_fish # make sure this lines up with correct specimens in the classifier

Cyc80_coords <- gdf_head$coords[,,Cyc80_fish]
Cyc80_mean_coords <- mshape(Cyc80_coords)

Cyc100_fish <- which(gdf_head$Treatment == "Cyc100")
Cyc100_fish # make sure this lines up with correct specimens in the classifier

Cyc100_coords <- gdf_head$coords[,,Cyc100_fish]
Cyc100_mean_coords <- mshape(Cyc100_coords)

# create a matrix for GGPlot with everything except shape coords

gdf_head$Treatment <- factor(gdf_head$Treatment, levels = c("DMSO", "Cyc20", "Cyc40", "Cyc60", "Cyc80", "Cyc100"))
gdf_head$Csize <- as.numeric(gdf_head$Csize)

ggplot_df <- as.data.frame(cbind(as.character(gdf_head$Csize),
                                 as.character(gdf_head$Treatment),
                                 as.character(gdf_head$Pdist)))


row.names(ggplot_df) <- dimnames(GPA_All$coords)[[3]]
colnames(ggplot_df) <- c("Csize", "Treatment", "Pdist")

ggplot_df$Csize <- as.numeric(ggplot_df$Csize)
ggplot_df$Treatment <- as.factor(ggplot_df$Treatment)
ggplot_df$Pdist <- as.numeric(ggplot_df$Pdist)

ggplot_df$Treatment <- factor(ggplot_df$Treatment, levels = c("DMSO", "Cyc20", "Cyc40", "Cyc60", "Cyc80", "Cyc100"))

# Csize differences

ggplot(ggplot_df, aes(x=Treatment, y = Csize, fill = Treatment)) +
  scale_fill_manual(values = palette()) +
  geom_boxplot(alpha = 0.75) + 
  theme_classic() + geom_jitter(width = 0.2, size = 1.25) +
  ggtitle("Centroid Size") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
        axis.title = element_text(size = 26), axis.text = element_text(size = 20),
        axis.line = element_line(size = 1), axis.ticks.length = unit(5, "pt")
  )

# Procrustes distance dispersion

Pdist <- ShapeDist(shapes = GPA_All$coords, reference = GPA_All$consensus)

DMSO_mean_coords

Pdist_from_DMSO <- ShapeDist(shapes = GPA_All$coords, reference = DMSO_mean_coords)

# Pdist dispersion by treatment from general consensus
ggplot(ggplot_df, aes(Pdist, fill = Treatment)) +
  scale_fill_manual(values = palette()) + geom_density(alpha =0.75) + 
  theme_classic() + 
  ggtitle("Procrustes Distance Dispersion") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Pdist dispersion by treatment from DMSO
ggplot(ggplot_df, aes(Pdist_from_DMSO, fill = Treatment)) +
  scale_fill_manual(values = palette()) + geom_density(alpha =0.75) + 
  theme_classic() + 
  ggtitle("Procrustes Distance Dispersion from DMSO") + 
  labs(x = "Procrustes distance from DMSO", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
        axis.title = element_text(size = 26), axis.text = element_text(size = 20),
        axis.line = element_line(size = 1), axis.ticks.length = unit(5, "pt")
  )

# var() plot variance

ggplot_df

?var()

var(ggplot_df$Pdist)
var(ggplot_df$Csize)

# Calculate the variances for Pdist and Csize within each Treatment group
variance_by_Treatment <- aggregate(. ~ Treatment, data = ggplot_df, FUN = var)

print(variance_by_Treatment)

# show variances grouped by treatment on a plot

?geom_bar

palette(c("black", "blue", "purple", "pink", "orange", "red"))

# Pdist
variances_Pdist <- ggplot_df %>%
  group_by(Treatment) %>%
  summarise(variance = var(Pdist))

ggplot(variances_Pdist, aes(x = Treatment, y = variance)) + 
  scale_fill_manual(values = palette()) +
  geom_bar(stat = "identity", fill = palette(), colour = "black", size = 0.5) + 
  theme_classic() +
  labs(title = "Variance of Procrustes Distance", x = "Treatment", y = "Pdist Variance") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Csize
variances_Csize <- ggplot_df %>%
  group_by(Treatment) %>%
  summarise(variance = var(Csize))

ggplot(variances_Csize, aes(x = Treatment, y = variance)) + 
  scale_fill_manual(values = palette()) +
  geom_bar(stat = "identity", fill = palette(), colour = "black", size = 0.5) + 
  theme_classic() +
  labs(title = "Variance of Centroid Size", x = "Treatment", y = "Csize Variance") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



# variance-covariance trace to display shape variance (gives true indication of shape variance based on the GPA)

# calculate trace

gdf_head

gdf_head$Treatment

treatment <- as.character(unique(gdf_head$Treatment))

treatment_combinations <- combinations(treatment, n = length(treatment), r = 2)

treatment_levels <- treatment

# create empty data frame to receive loop data
trace_results <- data.frame(array(NA, dim = c(length(treatment_levels),2)))
colnames(trace_results) <- c("Group", "Trace")


# for loop to generate trace of each treatment group

for (i in 1:length(treatment_levels)) {
  # Subset the data for the current treatment level
  group_i <- treatment_levels[i]
  treatment_vec <- gdf_head$Treatment
  
  # create 2d array that can be used to calculate trace (cannot interpret a 3d array)
  coord_data <- gdf_head$coords
  two_d_coords <- two.d.array(coord_data)
  
  data_i <- two_d_coords[which(treatment_vec == group_i),]
  trace_i <- sum(diag(cov(data_i)))
  
  trace_results[i,1] <- paste0(group_i)
  trace_results[i,2] <- trace_i
  }

trace_results$Group <- factor(trace_results$Group, levels = c("DMSO", "Cyc20", "Cyc40", "Cyc60", "Cyc80", "Cyc100")) # not working?

print(trace_results)

# graph for trace, colours kept being wrong unless I add in specific line

palette(c("black", "blue", "purple", "pink", "orange", "red"))



ggplot(trace_results, aes(x = factor(Group, levels = c("DMSO", "Cyc20", "Cyc40", "Cyc60", "Cyc80", "Cyc100")), y = Trace)) + 
  scale_fill_manual(values = palette()) +
  geom_bar(stat = "identity", fill = palette(), colour = "black", size = 0.5) + 
  theme_classic() +
  labs(title = "Shape Variance (Trace)", x = "Treatment", y = "Trace") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

trace_plot <- ggplot(trace_results, aes(x = Group, y = Trace, fill = Group)) +
  geom_bar(stat = "identity", colour = "black", size = 0.5) +
  scale_fill_manual(values = c(DMSO = "black", Cyc20 = "blue", Cyc40 = "purple", Cyc60 = "pink", 
                               Cyc80 = "orange", Cyc100 = "red")) +
  theme_classic() +
  labs(title = "Shape Variance (Trace)", x = "Treatment", y = "Trace") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
        axis.title = element_text(size = 26), axis.text = element_text(size = 20),
        axis.line = element_line(size = 1), axis.ticks.length = unit(5, "pt"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
        )

trace_plot

ggsave("Trace_plot_transparentbg.png", trace_plot, bg = "transparent", width = 9, height = 5)



# 5.3 ANOVAs ####

summary(gdf_head)

# packages to display results as table?
install.packages("broom")
install.packages("knitr")
install.packages("kableExtra")
install.packages("webshot")
library(broom)
library(knitr)
library(kableExtra)
library(webshot)

# Are specimens in different treatment groups significantly different in size?
Csize_ANOVA <- aov(gdf_head$Csize ~ gdf_head$Treatment, data = gdf_head)
summary(Csize_ANOVA)

# Are specimens in different treatment groups significantly different in shape?
shape_ANOVA <- procD.lm(GPA_All$coords ~ Treatment, data = gdf_head, RRPP = TRUE)
summary(shape_ANOVA)

# Are these shape differences correlated to size?
allometry <- procD.lm(GPA_All$coords ~ Csize, data = gdf_head, RRPP = TRUE)
summary(allometry)

# Are these shape differences correlated to size and the interaction with treatment?
allometry_treatment <- procD.lm(GPA_All$coords ~ Csize*Treatment, data = gdf_head, RRPP = TRUE)
summary(allometry_treatment)

# extract DMSO allometry - regression of shape on size in DMSO - set of coefficients allows shape prediction at specific sizes
# multiply shape coords of each fish of score from 

# regress out the DMSO allometry vector - regression score per fish - take residuals


# 5.4 Post-hoc tests ####

class(Csize_ANOVA)

?TukeyHSD

tukey_Csize <- TukeyHSD(Csize_ANOVA)

tukey_Csize

tukey_Csize_table <- as.data.frame(tukey_Csize$`gdf_head$Treatment`)
tukey_Csize_table <- tukey_Csize_table[, c("diff", "lwr", "upr", "p adj")]

setwd("./figs")
dir()

# under this currently doesn't work well, I made tables manually

kable(tukey_Csize_table, format = "markdown") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  save_kable(file = "tukey_Csize_results_table.html")

html_Csize_table <- kable(tukey_Csize_table, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

webshot::install_phantomjs()
webshot(html_Csize_table, file = "tukey_Csize_results_table.png", zoom = 2)

tukey_Shape <- TukeyHSD(shape_ANOVA)


plot(tukey_Csize, las = 1)

# for the procD.lm analyses - I can't get any of this to work. I want to be able to find out which groups are significantly different from one another, like I did with Csize
# loop to compare all combinations to run separate MANOVAs and combine and display into one table - example code from Benedikt

gdf_head

gdf_head$Treatment

treatment <- as.character(unique(gdf_head$Treatment))

treatment_combinations <- combinations(treatment, n = length(treatment), r = 2)

class(gdf_head$Treatment)
class(treatment)
class(treatment_combinations)

dim(treatment_combinations)

# now loop with treatment_combinations to do multiple pairwise comparisons

summary(GPA_All$coords)

?procD.lm

# Loop through combinations and perform pairwise shape ANOVA

for (i in 1:nrow(treatment_combinations)) {
  group1 <- treatment_combinations[i, 1]
  group2 <- treatment_combinations[i, 2]
  
  # Subset data for the two groups
  subset_data <- gdf_head$coords[,,gdf_head$Treatment %in% c(group1, group2)]
  subset_Treatment <- gdf_head$Treatment[gdf_head$Treatment %in% c(group1, group2)]
  
  # Perform pairwise ANOVA
  pairwise_ANOVA <- procD.lm(subset_data ~ subset_Treatment, RRPP = TRUE)
  
  # Print or store the results as needed
  print(paste("Pairwise ANOVA for groups", group1, "and", group2))
  print(summary(pairwise_ANOVA))
}


# asymmetric variance in geomorph - use procsym, specify symmetric landmarks then find the asymmetric component
# haven't got anywhere with this yet

?pairwise

posthoc_shape <- pairwise(shape_ANOVA, ~ Treatment)

class(shape_ANOVA)


install.packages("lsmeans")
library(lsmeans)

lsmeans_shape <- lsmeans(shape_ANOVA)

shape_residuals <- residuals(shape_ANOVA)
shape_model_matrix <- model.matrix(shape_ANOVA)

class(shape_residuals)

tukey_shape <- TukeyHSD(shape_ANOVA, data = gdf_head())



# 6. Morphs and Heatmaps ####
# find code in the R_intro_morphometrics git

# Make sure meshes have been imported
DMSO_atlas_mesh <- file2mesh("./data/atlas/Zfish_Atlas_DMSO_ascii.ply")

open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(DMSO_atlas_mesh, color = "gray", alpha =0.9)
close3d()

ref_LMs <- as.matrix(DMSO_LM_ATLAS$LMs)

# 6.1. Generating PC and mean morphs ####

# generate PC1 min/max landmarks

PC1_min_shape <- as.matrix(PCA_head$shapes[[1]]$min)
PC1_max_shape <- as.matrix(PCA_head$shapes[[1]]$max)

# generate PC2 min/max landmarks

PC2_min_shape <- as.matrix(PCA_head$shapes[[2]]$min)
PC2_max_shape <- as.matrix(PCA_head$shapes[[2]]$max)

?tps3d

# create and check min/max morphs

# PC1
PC1_min_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = PC1_min_shape)
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(PC1_min_mesh, color = "white", alpha = 1, specular = 1)
close3d()

PC1_max_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = PC1_max_shape)
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(PC1_max_mesh, color = "white", alpha = 1, specular = 1)
close3d()

# PC2
PC2_min_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = PC2_min_shape)
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(PC2_min_mesh, color = "white", alpha = 1, specular = 1)
close3d()

PC2_max_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = PC2_max_shape)
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(PC2_max_mesh, color = "white", alpha = 1, specular = 1)
close3d()

# these created meshes did not have the same orientation as the atlas mesh, create new userMatrix positions and save.

lateral <- par3d()$userMatrix # navigate mesh display to lateral view then run this
dorsal <- par3d()$userMatrix # dorsal view
ventral <- par3d()$userMatrix # ventral view
transverse <- par3d()$userMatrix # transverse view
save(lateral, dorsal, ventral, transverse, file = "./data/atlas/RGL_positions_morphs.Rdata")

# LOAD NEW VIEWS FOR MORPHS
load("./data/atlas/RGL_positions_morphs.Rdata")

# make these views a list to use in a for loop later
morph_views <- list(dorsal = dorsal, lateral = lateral, ventral = ventral, transverse = transverse)

# loops to save PC1 morphs in all views

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(PC1_min_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("PC1_min_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps/pc_morphs", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(PC1_max_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("PC1_max_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps/pc_morphs", output_file), fmt = "png", top = TRUE)
  close3d()
}

# loops to save PC2 morphs in all views

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(PC2_min_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("PC2_min_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps/pc_morphs", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(PC2_max_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("PC2_max_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps/pc_morphs", output_file), fmt = "png", top = TRUE)
  close3d()
}


# create morphs using mean coords - mean coords were subset in 5.2 

DMSO_mean_coords
DMSO_mean_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = as.matrix(DMSO_mean_coords))
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(DMSO_mean_mesh, color = "white", alpha = 1, specular = 1)
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(DMSO_mean_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("DMSO_mean_", view, ".png")
  rgl.snapshot(file.path("./figs", output_file), fmt = "png", top = TRUE)
  close3d()
}

Cyc20_mean_coords
Cyc20_mean_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = as.matrix(Cyc20_mean_coords))
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(Cyc20_mean_mesh, color = "white", alpha = 1, specular = 1)
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(Cyc20_mean_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("Cyc20_mean_", view, ".png")
  rgl.snapshot(file.path("./figs", output_file), fmt = "png", top = TRUE)
  close3d()
}

Cyc40_mean_coords
Cyc40_mean_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = as.matrix(Cyc40_mean_coords))
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(Cyc40_mean_mesh, color = "white", alpha = 1, specular = 1)
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(Cyc40_mean_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("Cyc40_mean_", view, ".png")
  rgl.snapshot(file.path("./figs", output_file), fmt = "png", top = TRUE)
  close3d()
}

Cyc60_mean_coords
Cyc60_mean_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = as.matrix(Cyc60_mean_coords))
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(Cyc60_mean_mesh, color = "white", alpha = 1, specular = 1)
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(Cyc60_mean_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("Cyc60_mean_", view, ".png")
  rgl.snapshot(file.path("./figs", output_file), fmt = "png", top = TRUE)
  close3d()
}

Cyc80_mean_coords
Cyc80_mean_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = as.matrix(Cyc80_mean_coords))
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(Cyc80_mean_mesh, color = "white", alpha = 1, specular = 1)
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(Cyc80_mean_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("Cyc80_mean_", view, ".png")
  rgl.snapshot(file.path("./figs", output_file), fmt = "png", top = TRUE)
  close3d()
}

Cyc100_mean_coords
Cyc100_mean_mesh <- tps3d(x = DMSO_atlas_mesh, refmat = ref_LMs, tarmat = as.matrix(Cyc100_mean_coords))
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal)
shade3d(Cyc100_mean_mesh, color = "white", alpha = 1, specular = 1)
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  shade3d(Cyc100_mean_mesh, color = "white", alpha = 1, specular = 1)
  output_file <- paste0("Cyc100_mean_", view, ".png")
  rgl.snapshot(file.path("./figs", output_file), fmt = "png", top = TRUE)
  close3d()
}

# 6.2. Making heatmaps ####

# heat maps between treatment mean shapes

# compute distance ranges so we can create a standardised, symmetrical colour ramp across the full distance range
dist_DMSO_Cyc20 <- meshDist(DMSO_mean_mesh, Cyc20_mean_mesh, sign = TRUE)
mindist_DMSO_Cyc20 <- min(dist_DMSO_Cyc20$dists, na.rm = TRUE)
maxdist_DMSO_Cyc20 <- max(dist_DMSO_Cyc20$dists, na.rm = TRUE)
print(range(c(mindist_DMSO_Cyc20, maxdist_DMSO_Cyc20)))

dist_DMSO_Cyc40 <- meshDist(DMSO_mean_mesh, Cyc40_mean_mesh, sign = TRUE)
mindist_DMSO_Cyc40 <- min(dist_DMSO_Cyc40$dists, na.rm = TRUE)
maxdist_DMSO_Cyc40 <- max(dist_DMSO_Cyc40$dists, na.rm = TRUE)
print(range(c(mindist_DMSO_Cyc20, maxdist_DMSO_Cyc20)))

dist_DMSO_Cyc60 <- meshDist(DMSO_mean_mesh, Cyc60_mean_mesh, sign = TRUE)
mindist_DMSO_Cyc60 <- min(dist_DMSO_Cyc60$dists, na.rm = TRUE)
maxdist_DMSO_Cyc60 <- max(dist_DMSO_Cyc60$dists, na.rm = TRUE)
print(range(c(mindist_DMSO_Cyc60, maxdist_DMSO_Cyc60)))

dist_DMSO_Cyc80 <- meshDist(DMSO_mean_mesh, Cyc80_mean_mesh, sign = TRUE)
mindist_DMSO_Cyc80 <- min(dist_DMSO_Cyc80$dists, na.rm = TRUE)
maxdist_DMSO_Cyc80 <- max(dist_DMSO_Cyc80$dists, na.rm = TRUE)
print(range(c(mindist_DMSO_Cyc80, maxdist_DMSO_Cyc80)))

dist_DMSO_Cyc100 <- meshDist(DMSO_mean_mesh, Cyc100_mean_mesh, sign = TRUE)
mindist_DMSO_Cyc100 <- min(dist_DMSO_Cyc100$dists, na.rm = TRUE)
maxdist_DMSO_Cyc100 <- max(dist_DMSO_Cyc100$dists, na.rm = TRUE)
print(range(c(mindist_DMSO_Cyc100, maxdist_DMSO_Cyc100)))

dist_Cyc20_Cyc40 <- meshDist(Cyc20_mean_mesh, Cyc40_mean_mesh, sign = TRUE)
mindist_Cyc20_Cyc40 <- min(dist_Cyc20_Cyc40$dists, na.rm = TRUE)
maxdist_Cyc20_Cyc40 <- max(dist_Cyc20_Cyc40$dists, na.rm = TRUE)
print(range(c(mindist_Cyc20_Cyc40, maxdist_Cyc20_Cyc40)))

dist_Cyc20_Cyc60 <- meshDist(Cyc20_mean_mesh, Cyc60_mean_mesh, sign = TRUE)
mindist_Cyc20_Cyc60 <- min(dist_Cyc20_Cyc60$dists, na.rm = TRUE)
maxdist_Cyc20_Cyc60 <- max(dist_Cyc20_Cyc60$dists, na.rm = TRUE)
print(range(c(mindist_Cyc20_Cyc60, maxdist_Cyc20_Cyc60)))

dist_Cyc20_Cyc80 <- meshDist(Cyc20_mean_mesh, Cyc80_mean_mesh, sign = TRUE)
mindist_Cyc20_Cyc80 <- min(dist_Cyc20_Cyc80$dists, na.rm = TRUE)
maxdist_Cyc20_Cyc80 <- max(dist_Cyc20_Cyc80$dists, na.rm = TRUE)
print(range(c(mindist_Cyc20_Cyc80, maxdist_Cyc20_Cyc80)))

dist_Cyc20_Cyc100 <- meshDist(Cyc20_mean_mesh, Cyc100_mean_mesh, sign = TRUE)
mindist_Cyc20_Cyc100 <- min(dist_Cyc20_Cyc100$dists, na.rm = TRUE)
maxdist_Cyc20_Cyc100 <- max(dist_Cyc20_Cyc100$dists, na.rm = TRUE)
print(range(c(mindist_Cyc20_Cyc100, maxdist_Cyc20_Cyc100)))

dist_Cyc40_Cyc60 <- meshDist(Cyc40_mean_mesh, Cyc60_mean_mesh, sign = TRUE)
mindist_Cyc40_Cyc60 <- min(dist_Cyc40_Cyc60$dists, na.rm = TRUE)
maxdist_Cyc40_Cyc60 <- max(dist_Cyc40_Cyc60$dists, na.rm = TRUE)
print(range(c(mindist_Cyc40_Cyc60, maxdist_Cyc40_Cyc60)))

dist_Cyc40_Cyc80 <- meshDist(Cyc40_mean_mesh, Cyc80_mean_mesh, sign = TRUE)
mindist_Cyc40_Cyc80 <- min(dist_Cyc40_Cyc80$dists, na.rm = TRUE)
maxdist_Cyc40_Cyc80 <- max(dist_Cyc40_Cyc80$dists, na.rm = TRUE)
print(range(c(mindist_Cyc40_Cyc80, maxdist_Cyc40_Cyc80)))

dist_Cyc40_Cyc100 <- meshDist(Cyc40_mean_mesh, Cyc100_mean_mesh, sign = TRUE)
mindist_Cyc40_Cyc100 <- min(dist_Cyc40_Cyc100$dists, na.rm = TRUE)
maxdist_Cyc40_Cyc100 <- max(dist_Cyc40_Cyc100$dists, na.rm = TRUE)
print(range(c(mindist_Cyc40_Cyc100, maxdist_Cyc40_Cyc100)))

dist_Cyc60_Cyc80 <- meshDist(Cyc60_mean_mesh, Cyc80_mean_mesh, sign = TRUE)
mindist_Cyc60_Cyc80 <- min(dist_Cyc60_Cyc80$dists, na.rm = TRUE)
maxdist_Cyc60_Cyc80 <- max(dist_Cyc60_Cyc80$dists, na.rm = TRUE)
print(range(c(mindist_Cyc60_Cyc80, maxdist_Cyc60_Cyc80)))

dist_Cyc60_Cyc100 <- meshDist(Cyc60_mean_mesh, Cyc100_mean_mesh, sign = TRUE)
mindist_Cyc60_Cyc100 <- min(dist_Cyc60_Cyc100$dists, na.rm = TRUE)
maxdist_Cyc60_Cyc100 <- max(dist_Cyc60_Cyc100$dists, na.rm = TRUE)
print(range(c(mindist_Cyc60_Cyc100, maxdist_Cyc60_Cyc100)))

dist_Cyc80_Cyc100 <- meshDist(Cyc80_mean_mesh, Cyc100_mean_mesh, sign = TRUE)
mindist_Cyc80_Cyc100 <- min(dist_Cyc80_Cyc100$dists, na.rm = TRUE)
maxdist_Cyc80_Cyc100 <- max(dist_Cyc80_Cyc100$dists, na.rm = TRUE)
print(range(c(mindist_Cyc80_Cyc100, maxdist_Cyc80_Cyc100)))

# combine all min and max distances
all_minmax <- c(mindist_DMSO_Cyc20, maxdist_DMSO_Cyc20, mindist_DMSO_Cyc40, maxdist_DMSO_Cyc40, mindist_DMSO_Cyc60, maxdist_DMSO_Cyc60, mindist_DMSO_Cyc80, maxdist_DMSO_Cyc80, mindist_DMSO_Cyc100, maxdist_DMSO_Cyc100,
                mindist_Cyc20_Cyc40, maxdist_Cyc20_Cyc40, mindist_Cyc20_Cyc60, maxdist_Cyc20_Cyc60, mindist_Cyc20_Cyc80, maxdist_Cyc20_Cyc80, mindist_Cyc20_Cyc100, maxdist_Cyc20_Cyc100,
                mindist_Cyc40_Cyc60, maxdist_Cyc40_Cyc60, mindist_Cyc40_Cyc80, maxdist_Cyc40_Cyc80, mindist_Cyc40_Cyc100, maxdist_Cyc40_Cyc100,
                mindist_Cyc60_Cyc80, maxdist_Cyc60_Cyc80, mindist_Cyc60_Cyc100, maxdist_Cyc60_Cyc100, mindist_Cyc80_Cyc100, maxdist_Cyc80_Cyc100)

# Make all values positive
all_minmax_positive <- abs(all_minmax)

# Find and print the largest positive value - this will be the "to" argument in meshDist
max_value <- max(all_minmax_positive)
print(max_value)

# make the max value negative - this will be the "from" argument in meshDist
max_value_negative <- -max_value
print(max_value_negative)

# define colour ramp
custom_ramp <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100)

# create and check heat maps using standardised colour ramp based on global symmetrical scale

# DMSO-Cyc20 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc20_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# DMSO-Cyc40 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc40_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# DMSO-Cyc60 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc60_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# DMSO-Cyc80 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc80_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# DMSO-Cyc100 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc100_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# Cyc20-Cyc40 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(Cyc20_mean_mesh, Cyc40_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# Cyc20-Cyc60 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(Cyc20_mean_mesh, Cyc50_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# Cyc40-Cyc60 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(Cyc40_mean_mesh, Cyc60_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# for loops to save heat maps in each view (Can I make this one loop?)

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc20_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc20_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc40_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc40_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc60_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc60_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc80_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc80_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc100_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc100_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(Cyc20_mean_mesh, Cyc40_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("Cyc20-Cyc40_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(Cyc20_mean_mesh, Cyc60_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("Cyc20-Cyc60_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(Cyc40_mean_mesh, Cyc60_mean_mesh, from = max_value_negative, to = max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("Cyc40-Cyc60_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

# non-global heat map scale to exaggerate differences

# DMSO-Cyc20 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc20_mean_mesh, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc20_mean_mesh, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc20_heatmap_defaultscale_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

# DMSO-Cyc40 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc40_mean_mesh, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc40_mean_mesh, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc40_heatmap_defaultscale_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

# Cyc20-Cyc40 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(Cyc20_mean_mesh, Cyc40_mean_mesh, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(Cyc20_mean_mesh, Cyc40_mean_mesh, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("Cyc20-Cyc40_heatmap_defaultscale_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

# same thing, but make symmetrical scale using just the minmax distances associated with Cyc20 and Cyc40

# combine min and max distances for relevant groups
condensed_minmax <- c(mindist_DMSO_Cyc20, maxdist_DMSO_Cyc20, mindist_DMSO_Cyc40, maxdist_DMSO_Cyc40, mindist_Cyc20_Cyc40, maxdist_Cyc20_Cyc40)

# Make all values positive
condensed_minmax_positive <- abs(condensed_minmax)

# Find and print the largest positive value - this will be the "to" argument in meshDist
condensed_max_value <- max(condensed_minmax_positive)
print(condensed_max_value)

# make the max value negative - this will be the "from" argument in meshDist
condensed_max_value_negative <- -condensed_max_value
print(condensed_max_value_negative)

# DMSO-Cyc20 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc20_mean_mesh, from = condensed_max_value_negative, to = condensed_max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc20_mean_mesh, from = condensed_max_value_negative, to = condensed_max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc20_heatmap_condensedscale_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

# DMSO-Cyc40 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(DMSO_mean_mesh, Cyc40_mean_mesh, from = condensed_max_value_negative, to = condensed_max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(DMSO_mean_mesh, Cyc40_mean_mesh, from = condensed_max_value_negative, to = condensed_max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc40_heatmap_condensedscale_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

# Cyc20-Cyc40 means
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(Cyc20_mean_mesh, Cyc40_mean_mesh, from = condensed_max_value_negative, to = condensed_max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(Cyc20_mean_mesh, Cyc40_mean_mesh, from = condensed_max_value_negative, to = condensed_max_value, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("DMSO-Cyc40_heatmap_condensedscale_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps", output_file), fmt = "png", top = TRUE)
  close3d()
}

# 6.2.1 PC morph heatmaps ####

# PC1 heat maps

# calculate distance min/max for colour ramp scale
dist_PC1min_PC1max <- meshDist(PC1_min_mesh, PC1_max_mesh, sign = TRUE)
mindist_PC1min_PC1max <- min(dist_PC1min_PC1max$dists, na.rm = TRUE)
maxdist_PC1min_PC1max <- max(dist_PC1min_PC1max$dists, na.rm = TRUE)
print(range(c(mindist_PC1min_PC1max, maxdist_PC1min_PC1max)))

dist_PC1max_PC1min <- meshDist(PC1_max_mesh, PC1_min_mesh, sign = TRUE)
mindist_PC1max_PC1min <- min(dist_PC1max_PC1min$dists, na.rm = TRUE)
maxdist_PC1max_PC1min <- max(dist_PC1max_PC1min$dists, na.rm = TRUE)
print(range(c(mindist_PC1max_PC1min, maxdist_PC1max_PC1min)))

# combine all min and max distances
PC1_minmax <- c(mindist_PC1min_PC1max, maxdist_PC1min_PC1max, mindist_PC1max_PC1min, maxdist_PC1max_PC1min)

# Make all values positive
PC1_minmax_positive <- abs(PC1_minmax)

# Find and print the largest positive value - this will be the "to" argument in meshDist
max_value_PC1only <- max(PC1_minmax_positive)
print(max_value_PC1only)

# make the max value negative - this will be the "from" argument in meshDist
max_value_PC1only_negative <- -max_value_PC1only
print(max_value_PC1only_negative)

# create and check heat maps of PC1_min and max using common, symmetrical scale

# PC1_min-PC1_max
open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(PC1_min_mesh, PC1_max_mesh, from = max_value_PC1only_negative, to = max_value_PC1only, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

open3d(windowRect = c(0,0,1000,700), userMatrix = dorsal, zoom =0.75)
meshDist(PC1_max_mesh, PC1_min_mesh, from = max_value_PC1only_negative, to = max_value_PC1only, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
close3d()

# for loop to save them in all views
for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(PC1_min_mesh, PC1_max_mesh, from = max_value_PC1only_negative, to = max_value_PC1only, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("PC1min-PC1max_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps/pc_morphs", output_file), fmt = "png", top = TRUE)
  close3d()
}

for (view in names(morph_views)) 
{
  open3d(windowRect = c(0, 0, 1000, 700), userMatrix = morph_views[[view]], zoom = 0.75)
  meshDist(PC1_max_mesh, PC1_min_mesh, from = max_value_PC1only_negative, to = max_value_PC1only, steps = 100, rampcolors = custom_ramp, sign = TRUE, titleplot = "Common Point Distance")
  output_file <- paste0("PC1max-PC1min_heatmap_", view, ".png")
  rgl.snapshot(file.path("./figs/heatmaps/pc_morphs", output_file), fmt = "png", top = TRUE)
  close3d()
}



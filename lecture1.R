# Lecture 1. Introduction and Setup

# Installing the current release of Bioconductor
BiocManager::install(version = "3.16")

# Installing packages from Bioconductor follows the command:
BiocManager::install(c("GWASTools", "GENESIS", "gdsfmt"), 
                     force = TRUE)

# Install 'generalize' and 'MetaCor' packages from Tamar Sofer GitHub:
library("devtools")
install_github("tamartsi/generalize@Package_update",
               subdir = "generalize")
install_github("tamartsi/MetaCor")

# Install the remaining packages
install.packages("mvtnorm")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("RColorBrewer")
####

# Ok! So with all packages installed, lets finish our setup:
## Saving our working directory in a object
## That's the way they use it to read files and stuff..
dir <- "/home/chpassos/association_mapping/sisg2017"

# Downloading datasets:
## The datasets can be found on the SISG course module website: 
## https://si.biostat.washington.edu/suminst/archives/SISG2017/SM1712
# wget https://si.biostat.washington.edu/sites/default/files/modules/datasets_1.zip
# unzip datasets_1.zip

# Installing SUGEN software
## We'll do this in the command line, using the following lines of code:
# wget https://github.com/dragontaoran/SUGEN/archive/master.zip
# unzip master.zip
# cd SUGEN-master
# make


# Instaling EasyStrata
install.packages("EasyStrata")
## Apparently I'm not able to install EasyStrata. 
## Fow now, lets proceed without it.

# Lecture 3. Simulating a dataset similar to HCHS/SOL

# Quick review on HCHS/SOL dataset
## Yes, I'm writing again
## But just because we're going to need in order to simulate some data

# HCHS/SOL review:
## 1. Longitudinal Study
## 2. Four locations: Chicago, Bronx, Miami, San Diego
## 3. Admixture: 3-way ancestry - African, European and Native American
## 4. Self-Identification: 
#### 4.a Mexican, Dominican, Puerto Rican, South Americna, Cuban and Central American
## 5. Local Ancestry was estimated in HCHS/SOL

# OK! That was our quick review.
# Why simulate data?
## Due to privacy restrictions, I don't have access to HCHS/SOL data
## (neither do the students who'd taken this course :) )
## 'Cause of this, we'll generate data that is somewhat similar to HCHS/SOL

# Simulating Data
## The data is ready-to-go for us, but how they generated it?
### 1. Hapgen software 
### 2. Simulate data from two populations: CEU and MEX
##### 2.a Not the same references to HCHS/SOL (3-way admixture, remember?)
### 3. Assumed that each individual had two parents
#### 3.a Of the two parental chr, one was entirely CEU and the other MEX
### 4. They randomly assigned intervals inherited from each parent 
##### 4.a to be those from the first or the second chromosome
### 5. The probability of CEU ancestry was:
##### 5.a Either 0.8 for the "UW genetic analysis group";
##### 5.b Or 0.5 for the "UNC genetic analysis group".
####### Not sure exactly what is UW and UNC groups
####### It is referring to the unis, but what is its meaning in the data?
#### 6. Individuals have both ancestries in each chromosome (mosaic). 

# Files that we're going to use
dir <- "/home/chpassos/association_mapping/sisg2017"
list.files(dir)

# Load GWASTools package
library("GWASTools", quietly=TRUE)
## Manual for GWASTools package
### https://www.bioconductor.org/packages/devel/bioc/manuals/GWASTools/man/GWASTools.pdf

# GWASTools works with GDS files
## GDS files have an "attached" Variant Annotation File

# When working with genotype data:
## 1. We first define a Genotype Reader Object [GdsGenotypeReader]
## 2. Then a genotype data object [Genotype Data]
## 3. The latter could be associated with the SNP annotation


# Using GWASTools
## GdsGenotypeReader is used to read genotype data stored in GDS files
gds <- GdsGenotypeReader(file.path(dir,
                                   "SISG_snp_dosages.gds"))
gds
## genotype -> Matrix with dimensions (snp, sample)
## sample.id -> a unique integer vector of scan ids
## snp.id -> a unique integer vector of snp ids
## snp.chromosome -> integer chromosome codes
## snp.position -> integer position values

# Looking at our Data 
## getChromosome returns a vector of chromosomes present in GDS
head(getChromosome(gds))
## Sample IDs of the first 5 individuals
### getScanID is basically returning a vector of our sample ids
getScanID(gds)[1:5]
## getSnpID returns a unique integer vector of SNP ids 
head(getSnpID(gds))
## getPosition returns a vector of base pair positions
head(getPosition(gds))

# Connecting GDS file with SNP annotation file 
## Via genotypeData object (just like we said before)

## getobj returns an R object stored in in .RData file
snpAnnot <- getobj(file.path(dir,
                             "SISG_snp_dosages_snpAnnot.RData"))
## pData returns ALL annotation variables as a dataframe
pData(snpAnnot)
dim(pData(snpAnnot)) # Dimensions of that dataframe
head(pData(snpAnnot)[,c(1:5)]) # First five columns 
head(pData(snpAnnot)[,c(6:9)]) # Remaining columns
## varMetadata returns metadata describing the annotation variables
varMetadata(snpAnnot) # So we can understand what each columns is

# Cool! Lets connect Genotype Reader Object to snpAnnot object
## Creating a genotypeData object
### Basically, the GenotypeData object is a CLASS of object!
### It is used for storing genotype data from a genome-wide association study
### together with the metadata associated with the subjects and SNPs
### used in the study
#### Create an GenotypeData object with GenotypeData()
#### providing the GDS file and the SNP annotation file
genoData <- GenotypeData(gds, snpAnnot=snpAnnot)

getAlleleA(genoData)[1:5] # Vector of A alleles

#  Returns the snp annotation variable varname
## I can replace the "rsID" below with any column on genoData
rsIDs <- getSnpVariable(genoData, "rsID") 
rsIDs[1:5]

# getGenotypeSelection Extracts genotype Values -> NUMBER OF A ALLELES
## We can extract for a single SNP (just like below)
## Or omit snp argument, then return for all SNPs
## scan argument is returning the first 10 individuals (?)
## We can change the scan
getGenotypeSelection(genoData, 
                     snp = (rsIDs == "rs2977656"), 
                     scan = 1:10)
getGenotypeSelection(genoData)

# Connecting Genotype Data with SAMPLE ANNOTATIONS
## I'm assuming for now that those are the phenotypes!
### We already know what getobj does.
scanAnnot <- getobj(file.path(dir,
                              "SISG_phenotypes.RData"))
scanAnnot
# varLabels returns a vector with the name of all columns 
varLabels(scanAnnot)[1:4]
varLabels(scanAnnot)

# Apparently, the sex column is important
## We must have the sex column in the data
## It needs to be named "sex"
## Males are "M"; Females are "F"


# Now connecting it to our Genotype Data
genoData <- GenotypeData(gds,
                         snpAnnot=snpAnnot, scanAnnot = scanAnnot)
## See now how we get the original GDS file, SNP annotation AND SCAN annotation?
genoData 

# Calculating Allele Frequencies
## alleleFrequency calculates freq for the A allele, over all SNPs
Afreqs <- alleleFrequency(genoData) 
## See how it calculates for M and F, and for Both!
head(Afreqs)


## Once we are done, we can close the GDS file:
require(gdsfmt)
showfile.gds(close = TRUE)
# we can also use close(gds)

# Concluding remarks
## Individuals in HCHS/SOL can be related
#### Both household related or genetic related 
## So, the simulated dataset also have related individuals!
## Those were created through sampling the real matrices
### 1. Household matrices (HH)
### 2. Genetic Relatedness Matrices (GRM).

## Household Matrix
HH.mat <- getobj(file.path(dir,
                           "SISG_houshold_matrix.RData"))
### The HH matrix have 1 in the i, j entry 
### If the i, j individuals live in the same household.
HH.mat[1:5, 1:5]
# According to HH matrix, only 19 individuals live in  
# the same house as other people in the study
sum(rowSums(HH.mat) > 1)

## Genetic Relatedness Matrix
kin.mat <- getobj(file.path(dir,
                            "SISG_relatedness_matrix.RData"))
# There are negative kinship values
# and diagonal values are not exactly 1. This is okay
kin.mat[1:5, 1:5]


# Exercises
## Use the GWASTools manual, your R knowledge, and the commands
## we learned to perform the following tasks and answer the questions:

# Open the GDS file
gds <- GdsGenotypeReader(file.path(dir, "SISG_snp_dosages.gds"))
# Connect it with SNP annot and Sample annot
snpAnnot <- getobj(file.path(dir, "SISG_snp_dosages_snpAnnot.RData"))
scanAnnot <- getobj(file.path(dir, "SISG_phenotypes.RData"))
genoData <- GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot = scanAnnot)


### 1. Compare the variance of the trait “trait” between the UW and the UNC groups
#### Tidyverse Approach:
pData(scanAnnot) |>
  as_tibble() |>
  group_by(group) |>
  summarise(trait_var = var(trait))

####################
#### GWASTools Approach (Need to work on that!):
#### UW and UNC groups are encoded in scanAnnot, in group variable
head(pData(scanAnnot))
#### Separate those who are from UW and from UNC
##### First, extract group variable
group_var <- getVariable(scanAnnot, varname = c("group"))
##### Then, identify the positions of those values in the variable
uw_pos <- which(group_var == "uw")
unc_pos <- which(group_var == "unc")
##### Finally, extract those who are from each group
uw_id <- getScanID(genoData)[uw_pos]
unc_id <- getScanID(genoData)[unc_pos]
# Trait in UW
uw_data <- pData(scanAnnot)[uw_pos[1]:uw_pos[length(uw_pos)], ]
getVariable(uw_data, varname = "trait")
####################

### 2. Plot a graph comparing the effect allele frequencies between the groups
#### Tidyverse Approach:
## Overall Allele Frequency
alleleFrequency(genoData) # Integer vector with IDs of scan to exclude
snpID <- getSnpID(gds)
length(snpID) # There are 7463 SNPs
## UW allele frequency
### getGenotype() extracts number of A alleles 
getGenotype(genoData, scan = uw_pos)
UW_af <- getGenotypeSelection(genoData, scan = uw_pos) |>
  as_tibble() |>
  mutate(snp = row_number()) |>
  pivot_longer(-snp, names_to = "sample", values_to = "count") |>
  group_by(snp) |>
  summarise(sum = sum(count),
            freq = sum / 7463) |>
  mutate(group = "UW")

## UNC allele frequency
getGenotype(genoData, scan = unc_pos)
UNC_af <- getGenotypeSelection(genoData, scan = unc_pos) |>
  as_tibble() |>
  mutate(snp = row_number()) |>
  pivot_longer(-snp, names_to = "sample", values_to = "count") |>
  group_by(snp) |>
  summarise(sum = sum(count),
            freq = sum / 7463) |>
  mutate(group = "UNC")

UW_af |>
  inner_join(UNC_af, by = "snp") |>
  select(snp, freq.x, freq.y) |>
  rename("UW" = "freq.x",
         "UNC" = "freq.y") |>
  pivot_longer(c(UW, UNC), names_to = "group", values_to = "freq") |>
  ggplot(aes(group, freq)) +
  geom_boxplot()


### 3. What is the genomic position of the SNP with the largest EAF difference between the UW and the UNC groups?
#### EAF apparently means effect allele Frequency
#### In order to find the genomic position with the largest difference, 
#### First, we need to find those differences
##### Tidyverse Approach:
diff_af <- UW_af |>
  inner_join(UNC_af, by = "snp") |>
  rename("UW" = "freq.x",
         "UNC" = "freq.y") |>
  mutate(diff = UW - UNC)

diff_af |>
  ggplot(aes(snp, diff)) +
  geom_point() # So it is a negative one

diff_af |>
  arrange(diff) # snp 178

# Lets see which SNP is the 178 and where it is located
GWASTools::getPosition(genoData)[178]
GWASTools::getSnpVariable(genoData, "rsID")[178]

### 4. What is the proportion of diseased individuals in males and females? 
### and in the UW and UNC groups?
#### I'm assuming 1 is diseased and 0 is healthy,
#### Since this info is not on metadata
getMetadata(scanAnnot)

##### Tidyverse apporach:
## Between Males and Females
pData(scanAnnot) |>
  as_tibble() |>
  group_by(sex) |>
  summarise(count = n())
## Females = 272
## Males =   228
pData(scanAnnot) |>
  as_tibble() |>
  group_by(sex) |>
  summarise(health = sum(disease)) |>
  ungroup() |>
  mutate(total = c(272, 228),
         prop_health = health/total,
         prop_disease = 1 - prop_health)

## Between UNC and UW
## UNC = 300 ind
## UW = 200 ind
pData(scanAnnot) |>
  as_tibble() |>
  group_by(group) |>
  summarise(health = sum(disease)) |>
  ungroup() |>
  mutate(total = c(300, 200),
         prop_health = health/total,
         prop_disease = 1 - prop_health)
  
  
### 5. Extract the genotypes of rs12033927 and rs17390062.
### What is the LD between them in the combined sample? 
### in the UW group? 
### in the UNC group?
#### Tidyverse approach:
##### Discovering which snps are those
pData(snpAnnot) |>
  as_tibble() |>
  filter(rsID %in% c("rs12033927", "rs17390062")) 

getGenotypeSelection(genoData, snp = c(785, 854))
#### CHP: GWASTools do LD??? I'm not finding how to do LD in GWASTools

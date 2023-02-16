# Admixture Mapping

# Lecture 11. Admixture Mapping 
## Finally!
## Need to finish as well. 
## Important things to keep in mind:

## How to transform Gnomix results into GDS?
## How to generate the admixDataList?

# What is Admixture Mapping?
## It is a type of association analysis that can be performed
## Only on admixed populations

## Human Populations around the world can be considered subpopulations
# They have distinct genomic features
# Different allele frequencies
# Even population-specific variants

# Our population, REDS-III, is an three-way admixed pop
## With contributions from African, European and Native American ancestry

# Since chromosomes are inherited from parents in a process of
# recombination, for each person, segments of the chromosome
# were inherited from a specific ancestry (and ancestor).

# There are softwares that estimates Local Ancestry
### for each person, “say" which ancestry this interval was inherited from


# Let’s look at a file with local ancestry information
## Load libraries
library("GWASTools")
library(tidyverse)
library(GENESIS)

dir <- "/home/chpassos/association_mapping/sisg2017/"

# Sample Annot
scanAnnot <- getobj(file.path(dir,
                              "SISG_phenotypes.RData"))

# This is the file with LAI?
admixAnnot <- getobj(file.path(dir,
                               "SISG_local_ancestry_snpAnnot.RData"))
dim(admixAnnot)
head(pData(admixAnnot))


# We need to use different functions (from before) to read the data
## Because now we have genotype dosages/count per ancestry!
require(gdsfmt)
gds <- openfn.gds(file.path(dir,
                            "SISG_local_ancestry_2.gds"))
gds
# Yeah, I can't open this file
# But I'm guessing is that it is one snp per line,
# and the two ancestries columns there represents counts of each ancestry
## What exactly is this dosage?

# Vector with ancestries
ancestries <- c("CEU", "MEX")
# Empty list
genoDataList <- vector(mode = "list", length = 2)
# Naming each page of the list
names(genoDataList) <- ancestries


for (ancestry in ancestries){
  gds.reader <- GdsGenotypeReader(gds,
                                  genotypeVar=paste0(ancestry, "_dosage"))
  genoDataList[[ancestry]] <- GenotypeData(gds.reader,
                                           scanAnnot=scanAnnot)
} # Don't know exactly what this loop did

# Now we have a list with genotype readers.
# We can apply functions presented before on each of the list components
genoDataList[["CEU"]]
getGenotype(genoDataList$CEU)[1:4,1:4]
getScanID(genoDataList$MEX)[1:4]

# The local ancestry “genotypes” are not actually counts of genotypes.
#### They are counts of ancestry!
# So if person p1 have count of 2 for its CEU dosage in interval 3,
# it means that interval 3 in the p1’s two chromosomes 1 was inherited from a
# CEU ancestor.
## In admixture mapping, we test the local ancestry counts. ##
# Again
# In admixture mapping, we test the local ancestry counts.

# Why is it meaningful? intuition:
### Because if a genotype is more frequent in one population (pop1);
##### And it is also associated with a trait;
###### People with more of the intervals spanning the genotypes
###### inherited from pop1, are more likely to have the trait.

# Mathematically, it can be shown that
## I If there are two ancestral populations;
## There is a single causal variant in an interval;
## The genotype effect size β is the same in the two populations;
## The allele frequency in the two ancestral populations are f1 and f2;
## Then the effect size estimated by testing the LAI counts is:
####: β(f1 − f2)

# This is lower than the effect size of the genotype β.
# If f1 == f2, this “effect” equals zero, and so admixture mapping is not useful.
# If f1 = 0 or f2 = 0, it’ll equal f2β or f1β.
# Regular association mapping may actually give something smaller.

# Suppose we had a specific variant g genotyped on all peopled.
# And we had the local ancestry counts in the interval spanning g.
# Then, testing the genotype association will be more powerful
# than testing the local ancestry.

# However...
## Sometimes we don’t have g genotyped, but we do have
### ancestry counts for an interval spanning g.
## There are many less local ancestry intervals than genotypes
### Reducing the p-value required for significance of LAI associations.
## The causal g may be rare, making genotype associations
### unstable, even if the genotype is in the data.
#### While LAI associations do not suffer from this instability.

# Admixture mapping is useful in these cases #





# Consider the case where we do have the causal genotype.
## Its effect β = 0.06, 
## and the frequency of the local ancestry from pop1 is 0.5.
### Significance of association results: 5 × 10−8
### Significance of admixture mapping results: 5.7 × 10−5
# See figure for difference in power comparing association to
# admixture mapping as a function of ancestry-specific frequencies.

## Admixture Mapping is best used when allele frequencies differences 
## Or even frequency of the local ancestry ??
# is high in admixed pop

# Additional settings in which admixture mapping is useful:
## The effect sizes is different between ancestries (even if the allele frequencies are the same across ancestries!).
## Multiple causal variants in the LAI.

# Take-home message: admixture mapping can sometimes detect
# association regions that association mapping cannot.

# In practice:
gds <- openfn.gds(file.path(dir,
                            "SISG_local_ancestry_2.gds"))
ancestries <- c("CEU", "MEX")
genoDataList <- vector(mode = "list", length = 2)
names(genoDataList) <- ancestries
for (ancestry in ancestries){
  gds.reader <- GdsGenotypeReader(gds,
                                  genotypeVar=paste0(ancestry, "_dosage"))
  genoDataList[[ancestry]] <- GenotypeData(gds.reader,
                                           scanAnnot=scanAnnot)
}

# Separating our covariates
covariates <- c("EV1", "EV2", "sex", "age", "group")
# Separating our outcome
outcome <- "trait"
# Correlation Matrices
## Household matrix
HH.mat <- getobj(file.path(dir,
                           "SISG_houshold_matrix.RData"))
## Genetic Relatedness Matrix
kin.mat <- getobj(file.path(dir,
                            "SISG_relatedness_matrix.RData"))
# Put both on a list
covMatList <- list(HH = HH.mat, kinship = kin.mat)

# Fit the null mixed model as before
nullmod <- fitNullModel(scanAnnot,
                        outcome = outcome,
                        covars = covariates,
                        cov.mat = covMatList,
                        family = "gaussian")


# and test using the GENESIS function admixMapMM.(MM: mixed model)
## Again, as GENESIS has been updated, functions became deprecated and replaced
# For admixture mapping, use instead admixMap

###############################
###############################

# option 1: one GDS file per ancestry
afrfile <- system.file("extdata", "HapMap_ASW_MXL_local_afr.gds", package="GENESIS")
amerfile <- system.file("extdata", "HapMap_ASW_MXL_local_amer.gds", package="GENESIS")
eurfile <- system.file("extdata", "HapMap_ASW_MXL_local_eur.gds", package="GENESIS")
files <- list(afr=afrfile, amer=amerfile, eur=eurfile)
gdsList <- lapply(files, GdsGenotypeReader)

# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(data.frame(
  scanID=getScanID(gdsList[[1]]), stringsAsFactors=FALSE))

# generate a phenotype
set.seed(4)
nsamp <- nrow(scanAnnot)
scanAnnot$pheno <- rnorm(nsamp, mean=0, sd=1)
set.seed(5)
scanAnnot$covar <- sample(0:1, nsamp, replace=TRUE)

genoDataList <- lapply(gdsList, GenotypeData, scanAnnot=scanAnnot)

# iterators
# if we have 3 ancestries total, only 2 should be included in test
genoIterators <- lapply(genoDataList[1:2], GenotypeBlockIterator)

# fit the null mixed model
null.model <- fitNullModel(scanAnnot, outcome="pheno", covars="covar")

# run the association test
myassoc <- admixMap(genoIterators, null.model)
head(myassoc)

lapply(genoDataList, close)


# option 2: create a single file with multiple ancestries
# first, get dosages from all ancestries
library(gdsfmt)
dosages <- lapply(files, function(f) {
  gds <- openfn.gds(f)
  geno <- read.gdsn(index.gdsn(gds, "genotype"))
  closefn.gds(gds)
  geno
})
lapply(dosages, dim)

# create a new file with three dosage matrices, keeping all
# sample and snp nodes from one original file
tmpfile <- tempfile()
file.copy(afrfile, tmpfile)
gds <- openfn.gds(tmpfile, readonly=FALSE)
delete.gdsn(index.gdsn(gds, "genotype"))
add.gdsn(gds, "dosage_afr", dosages[["afr"]])
add.gdsn(gds, "dosage_amer", dosages[["amer"]])
add.gdsn(gds, "dosage_eur", dosages[["eur"]])
closefn.gds(gds)
cleanup.gds(tmpfile)

# read in GDS data, specifying the node for each ancestry
gds <- openfn.gds(tmpfile)
gds
genoDataList <- list()
for (anc in c("afr", "amer", "eur")){
  gdsr <- GdsGenotypeReader(gds, genotypeVar=paste0("dosage_", anc))
  genoDataList[[anc]] <- GenotypeData(gdsr, scanAnnot=scanAnnot)
}

# iterators
genoIterators <- lapply(genoDataList[1:2], GenotypeBlockIterator)

# run the association test
myassoc <- admixMap(genoIterators, null.model)

close(genoDataList[[1]])
unlink(tmpfile)

myassoc
























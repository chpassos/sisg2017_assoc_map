# Admixture Mapping
## I need to finish this lecture

# Lecture 10. Generalization Analysis

# What is "Generalization'?
## It's the replication of an association between a genetic variant and a trait,
## discovered in one population,
## to another population!

# Most association studies were performed in European populations
### These are often detected in very large GWAS (100000 individuals)


# Why perform an generalization analysis?
## To know if the association holds true in other populations
#### This may not always be true

## To gain power by limiting the number of variants tested for
## associations to those already previously reported.

## Because we need to perform replication analysis, but we do not
## have access to an independent study with the same type of
## population and/or the same trait.


# An intuitive approach to generalization analysis:
### Take the list of SNP associations reported in a paper
### Test the same SNPs with the same trait in your data
### Report the significant associations.

#### What should be the p-value threshold to report associations? ####
# It was developed an generalization testing framework that
# originated in the replication analysis literature.
# Combining test results (p-values):
#### from both the discovery study and our study (the follow-up)
#### and calculate an r-value. (for every SNP).
# These r-values take into account multiple testing (of both studies),
# And are used like p-values.
# Since they are already adjusted for multiple testing, an
# association is generalized if the r-value< 0.05

# The generalization framework also takes into account the
# direction of associations.

# If the estimated association is negative in one study, and
# positive in the other, the association will not generalize.

# The generalization R package have an example from the HCHS/SOL platelet count paper.
# I We first load this package

require(generalize)
library(tidyverse)

# The generalization R package has an example data set.
# It has results reported by Geiger et al., 2011, 
## and matched association results from the HCHS/SOL.
# Generalization analysis is done for one study at a time


# load the data set from the package
data("dat")
# look at the column names:
matrix(colnames(dat), ncol = 3)


# The data.frame with the example provides all information we
# need for generalization analysis.

head(dat)
# Using a function from generalize package
## Apparently, it compares the effect allele and other alleles,
## also changes the sign of the effect between studies,
## if it was different.
dat.matched <- matchEffectAllele(dat$rsID,
                                 study2.effect = dat$study2.beta,
                                 study1.alleleA = dat$study1.alleleA,
                                 study2.alleleA = dat$study2.alleleA,
                                 study1.alleleB = dat$study1.alleleB,
                                 study2.alleleB = dat$study2.alleleB)


# Rearranging some columns....
dat$study2.beta <- dat.matched$study2.effect
dat$alleleA <- dat$study1.alleleA
dat$alleleB <- dat$study1.alleleB
dat$study1.alleleA <- dat$study1.alleleB <-
  dat$study2.alleleA <- dat$study2.alleleB <- NULL
head(dat)

# Test for Generalization:
gen.res <- testGeneralization(dat$rsID, dat$study1.pval,
                              dat$study2.pval, dat$study1.n.test[1],
                              study1.effect = dat$study1.beta,
                              study2.effect = dat$study2.beta,
                              directional.control = TRUE,
                              control.measure = "FDR" )
head(gen.res)
require(ggplot2,quietly = TRUE)
require(gridExtra,quietly = TRUE)
require(RColorBrewer,quietly = TRUE)
figure.out <- paste0(getwd(),
                     "/Generalization_example.pdf")
prepareGenResFigure(dat$rsID, dat$study1.beta,
                    dat$study1.se, dat$study2.beta, dat$study2.se,
                    gen.res$generalized, gen.res$gen.rvals,
                    dat$study1.n.test[1],
                    output.file = figure.out,
                    study1.name = "Study1",
                    study2.name = "Study2")

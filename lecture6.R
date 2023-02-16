# Admixture Mapping

# Lecture 7. Dealing with heterogeneity: group-specific variances and stratified analyses


# Admixed populations are complicated

# REDS-III for example is admixed with three ancestral populations
## European, African, Amerindian
## The proportions of admixture differ between background groups

# But these people can also differ in other aspects of daily lives
# Due to lifestyle and other environmental exposure differences

# So, both genetic and environmental differences translate to
# differences in phenotypic variability, disease prevalence.

# Prevalence of diseases also may differ between subgroups in the admixed pop

# These differences may occur in the VARIANCE of the trait!

# Dealing with heterogeneity: HETEROGENEOUS VARIANCES
## Recall the linear mixed model for quantitative traits:
## yi = xTi β + giα + bi + Ei, i = 1, . . . , n,

# Homogeneous Variances
# With homogeneous variances: Var(Eik ) = σ2e for all i = 1,..., n. (The usual model)


# Heterogeneous variances
# If there are well-defined sub-groups, 
# one can fit a model which assigns each group its own residual variance.
# With heterogeneous variances: Var(Eik ) = σ2e,k.
# Where K is each subgroup
# Suppose that each participant i is associated with a subgroup k = 1, ... ,K. Then:
# yik = xTik β + gik α + bik + Eik , i = 1, ... , n, k = 1, ... ,K
# with subscript k added to denote that person i is in group k


# Linear mixed models with heterogeneous variances

# Load libraries
library(GWASTools)
library(GENESIS)
library(tidyverse)

# Define dir object
dir <- paste0("/home/chpassos/association_mapping/sisg2017/")


# Scan Annot object
scanAnnot <- getobj(file.path(dir,
                              "SISG_phenotypes.RData"))
scanAnnot

varLabels(scanAnnot)[1:4]
varLabels(scanAnnot)[5:8]
# Separating our covariates
covariates <- c("EV1", "EV2", "sex", "age", "group")
# Separating our outcome of interest
outcome <- "trait" # A quantitative variable right here
# Loading our correlation matrices
## Household matrix
HH.mat <- getobj(file.path(dir,
                           "SISG_houshold_matrix.RData"))
## Genetic Relatedness matrix
kin.mat <- getobj(file.path(dir,
                            "SISG_relatedness_matrix.RData"))
## Put both in a list
covMatList <- list(HH = HH.mat, kinship = kin.mat)

# Fit the Null Model
## Note the group.var argument below
nullmod <- fitNullModel(scanAnnot,
                        outcome = outcome,
                        covars = covariates,
                        cov.mat = covMatList,
                        group.var = "group",
                        family = "gaussian")
## Results of the Null Model
names(nullmod)
nullmod$varComp
nullmod$fixef
# Unfortunately, function varCompCI(nullmod, prop = TRUE) is
# not supported for heterogeneous variances.

# Association testing as usual
gds <- GdsGenotypeReader(file.path(dir,
                                   "SISG_snp_dosages.gds"))
snpAnnot <- getobj(file.path(dir,
                             "SISG_snp_dosages_snpAnnot.RData"))
genoData <- GenotypeData(gds,
                         snpAnnot=snpAnnot, scanAnnot = scanAnnot)
genoIterator <- GenotypeBlockIterator(genoData, 
                                      snpBlock=10000)

assoc_het <- assocTestSingle(genoIterator,
                null.model = nullmod,
                BPPARAM = BiocParallel::SerialParam())
head(assoc)

## Exercises
# 1. Compare the p-values and effect estimates using the plot()
# command from GWAS with and without heterogeneous variances.
### You can compare -log(p-value,10).
## R: With heterogeneous variances (calculated above)
plot(assoc$Score.pval)

# With homogeneous variances
nullmod <- fitNullModel(scanAnnot,
                        outcome = outcome,
                        covars = covariates,
                        cov.mat = covMatList,
                        family = "gaussian")
assoc_hom <- assocTestSingle(genoIterator,
                             null.model = nullmod,
                             BPPARAM = BiocParallel::SerialParam())

p_het <- assoc_het |>
  as_tibble() |>
  select(variant.id, chr, pos,
         Score.pval) |>
  mutate(neg_logpval_het = -log(Score.pval, 10))
p_hom <- assoc_hom |>
  as_tibble() |>
  select(variant.id, chr, pos,
         Score.pval) |>
  mutate(neg_logpval_hom = -log(Score.pval, 10))

p_het |>
  inner_join(p_hom, by = "variant.id") |>
  ggplot(aes(neg_logpval_het, neg_logpval_hom)) +
  geom_point() # Is that it?

# 2. Can you use heterogeneous variances in logistic regression?
## R: I don't know.


# Dealing with heterogeneity: STRATIFIED ANALYSIS
## Heterogeneous variance model accounts for differences in residual variances
### While other effects remain common to all groups (fixed effects, genetic variances).

## STRATIFIED ANALYSIS allows for a different model, with different parameters, in each stratum.
# One can inspect results in each stratum, and combine results
# across strata in meta-analysis

# Combining results across strata in meta-analysis is not trivial,
# due to genetic relatedness and shared household (environmental correlation) 
# among individuals in different groups
 
# Method to do it: MetaCor.
# Meta-analysis when correlations between strata exist, and are modeled.
##### Mixed-models based approach

# Accounting for correlation may help with
## Inflation control
## Power of analysis

# Lets do it
require(MetaCor)

# Load data again
gds <- GdsGenotypeReader(file.path(dir,
                                   "SISG_snp_dosages.gds"))
snpAnnot <- getobj(file.path(dir,
                             "SISG_snp_dosages_snpAnnot.RData"))
scanAnnot <- getobj(file.path(dir,
                              "SISG_phenotypes.RData"))
genoData <- GenotypeData(gds,
                         snpAnnot=snpAnnot, scanAnnot = scanAnnot)

# Using the MetaCor package is slightly more complicated then using GENESIS
# It takes as arguments actual genotype counts and covariates info
### So we need to extract them first

# We can use a function from the GENESIS package.
covariates <- c("EV1", "EV2", "sex", "age")

## Fuuuck. The following function does not exist in GENESIS anymore
## And I can't find the equivalent in newer versions of the package..
## Guess I'll have to stop this lecture here.


dat <- GENESIS:::createDesignMatrix(
  scanData = pData(scanAnnot),
  outcome = outcome,
  covars = covariates,
  scan.include = scanAnnot$scanID)

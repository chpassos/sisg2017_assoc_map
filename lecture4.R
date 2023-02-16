# Admixture Mapping

# Lecture 4. Quantitaive trait mixed model GWAS

# Load Libraries
library(tidyverse)
library(GWASTools)
library(GENESIS)

# Define our dir object
dir <- "/home/chpassos/association_mapping/sisg2017/"

# Reading the scanAnnotation file
scanAnnot <- getobj(file.path(dir,
                              "SISG_phenotypes.RData"))
scanAnnot
head(pData(scanAnnot))
varLabels(scanAnnot)[1:4] # Fist four variables
varLabels(scanAnnot)[5:8] # Remaining variables

# Which of these variables are covariates?
covariates <- c("EV1", "EV2", "sex", "age", "group")

# Which is out outcome?
outcome <- "trait"


# Reading our correlations matrices
## Household matrix
HH.mat <- getobj(file.path(dir,
                           "SISG_houshold_matrix.RData"))
## Kinship matrix
kin.mat <- getobj(file.path(dir,
                            "SISG_relatedness_matrix.RData"))
## Put both in a list
### List with two pages, each being a matrix
covMatList <- list(HH = HH.mat, kinship = kin.mat)

# Fit the Null Model
## Since 2017 (when this module was offered), there has been some
## changes in the way the packages works.
## Nowadays, we need to use different functions to 
## fit the null model
## Now, we use fitNullModel()
nullmod <- fitNullModel(scanAnnot,
                        outcome = outcome,
                        covars = covariates,
                        cov.mat = covMatList,
                        family = "gaussian")
# Results of the null model
names(nullmod)
nullmod$varComp # V_kinship = 0
varCompCI(nullmod, prop = TRUE) # Proportion kinship = 0 (Not heritable!)
nullmod$fixef
## Although using different functions, at least the two variables
## above are the same!

# Tamar Sofer simulated another dataset
## Apparently, the trait above wasn't so heritable
require(mvtnorm)
n <- nrow(kin.mat)
set.seed(101)
new.trait <-
  nullmod$model.matrix %*% matrix(c(4, 5, -2, 1, 4,2)) +
  matrix(rmvnorm(n = 1, mean = rep(0, n),
                 sigma = diag(rep(120, n)) +
                   80*kin.mat + 40*HH.mat))
scanAnnot$new.trait <- as.numeric(new.trait)
set.seed(101)
nullmod <- fitNullModel(scanAnnot,
                        outcome = "new.trait",
                        covars = covariates,
                        cov.mat = covMatList)
nullmod$varComp # V_kinship = 9
varCompCI(nullmod, prop = TRUE) # Proportion kinship = 4%
## Again, very similar to the results in the module..
# The proportion of variance due to kinship/genetic relatedness is heritability
# The heritability of "trait" is estimated to be 4%, with 95% confidence interval (-96, 104)%.

### Outputs of Null Model
# varComp: the variance component estimates for the random effects
# fixef: a data.frame with point estimates, standard errors, test statistics, and p-values for each of the fixed effect covariates
# fit: a data.frame with the outcome, the fitted values, and various residuals from the model



# After estimating variance components in the “null model"
# they are assumed fixed
## We now use this null model object in association testing.
### Reading our GDS file
gds <- GdsGenotypeReader(file.path(dir,
                                   "SISG_snp_dosages.gds"))
### SNP Annotation file
snpAnnot <- getobj(file.path(dir,
                             "SISG_snp_dosages_snpAnnot.RData"))
### Creating the GenotypeData (GDS + SNP Annot + Scan Annot)
genoData <- GenotypeData(gds,
                         snpAnnot = snpAnnot,
                         scanAnnot = scanAnnot)
## To fit the linear mixed model, we need to use newer functions in GENESIS
# But a new feature has been introduced since 2017
# The function is assocTestSingle()
# But it appears that it doesn't accept GenotypeData as a input
# Instead, it accepts as a input the output of GenotypeBlockIterator() 
# So, we need to run this function first..
## What happens is: In 2017 this function didn't exist! So I'm doing on my own.

# We need to decide how many SNPs we want to read at a time
# This is done via GenotypeBlockIterator object that defines blocks of SNPs
## Default is 10000:
genoIterator <- GenotypeBlockIterator(genoData, 
                                      snpBlock=10000)

# Now we can move to the linear mixed model:
assoc <- assocTestSingle(genoIterator,
                         null.model = nullmod,
                         BPPARAM = BiocParallel::SerialParam())
# By default, the function will perform association tests at all SNPs in the genoData object.

### Outputs of Linear Mixed Model
# The assocTestSingle function will return a data.frame with summary 
# information from the association test for each SNP. 
# Each row corresponds to a different SNP.
head(assoc)
# variant.id: the unique snp ID
# chr: the chromosome
# pos: the position
# n.obs: the number of samples analyzed at that SNP
# freq: the frequency of the tested (“A”) allele
# MAC: the minor allele count
# Score: the value of the score function
# Score.SE: the estimated standard error of the score
# Score.Stat: the score Z test statistic
# Score.pval: the p-value based on the score test statistic
# Est: an approximation of the effect size estimate (beta) for that SNP
# Est.SE: an approximation of the standard error of the effect size estimate
# PVE: an approximation of the proportion of phenotype variance explained

close(gds)

# Inflation in Genome-Wide Association Studies
# Fundamental assumption in GWAS:
## Most genetic variants are not associated with the outcome.
# Test statistics are mostly distributed “under the null"
# λgc =  median(observed test statistics) / median(expected distribution of test statistics)
## Ideally, λ = 1.
## Because the bulk of the association are null.
# We make a QQ plot to examine the results
## q-q plots are used to evaluate inflation


### Exercises
# 1. Use the results from the GWAS that we ran on slide 18, and
# the function qqPlot() from the GWASTools package to make a
# q-q plot of the p-values from the Wald test

# R: The Wald test done in previous versions of GENESIS was actually 
# an approximation to the Wald test based on the Score test, 
# not an exact Wald test. 
# Because the underlying code of the Score and Wald tests was nearly identical, 
# the only meaningful difference was that the Wald test would
# return effect size estimates (Est). 
# The same effect size estimates from this approximation are still 
# available in the new output as the "Est" column.
GWASTools::qqPlot(assoc$Score.pval)



# 2. What is the inflation factor? use the R code 
pchisq(0.5, df = 1, lower.tail = FALSE) 
# to obtain the median of the expected distribution of the test statistics, and
median(assoc$Est, na.rm = TRUE) 
# to obtain the median of the observed test statistics.

# R: The inflation factor is calculated with the above equation
## Remember: GENESIS went through changes over time
## So that it does not provide Wald estimates anymore
## However, according to a maintainer of the package,
## EST column is similar to Wald
median(assoc$Est, na.rm = TRUE)/pchisq(0.5, df = 1, lower.tail = FALSE) 

# 3. Can you evaluate whether the GWAS is too inflated or deflated?
## R: Actually, no.
## I don't know if the above result is even CORRECT. :(

# 4. Use the function manhattanPlot() from the GWASTools
# package to make a Manhattan plots for these p-values.
GWASTools::manhattanPlot(assoc$Score.pval,
                         chromosome = rep(1, nrow(assoc)))

# 5. Set the significance threshold line in the Manhattan plot to be
# 0.05 divided by the number of tested variants (Bonferroni correction).
GWASTools::manhattanPlot(assoc$Score.pval,
                         chromosome = rep(1, nrow(assoc)),
                         signif = 0.05/nrow(assoc))


# 6. Show results in the Manhattan plot only for variants with
# imputation quality ("info") at least 0.8, or genotyped.
## R: Do I need to run another null model and another LMM?
### Identifying these variants
genotyped_variants <- pData(snpAnnot) |>
  as_tibble() |>
  filter(info >= 0.8) |>
  select(snpID) |>
  pull()

scores_genotyped_variants <- assoc |>
  as_tibble() |>
  filter(variant.id %in% genotyped_variants) 

GWASTools::manhattanPlot(scores_genotyped_variants$Score.pval,
                         chromosome = rep(1, nrow(scores_genotyped_variants)),
                         signif = 0.05/nrow(scores_genotyped_variants))


# 7. Which variant is most associated with "trait" among all variants (according to p-value)?
## R:
assoc |>
  as_tibble() |>
  arrange(desc(Score.pval)) |>
  slice(1)

# 8. Use the parameter snp.include of function assocTestMM() to
# test only variants in positions 1029889 - 2136826 on
# chromosome 1. Which variant has the most significant p-value?
## R: To restrict which SNPs I want in the LMM, we do this in
## while defining our genotype block iterator:
## The SNPs to include are:
snpID <- pData(snpAnnot) |>
  as_tibble() |>
  filter(position >= 1029889 & position <= 2136826) |>
  pull(snpID)
genoIterator <- GenotypeBlockIterator(genoData, 
                                      snpInclude = snpID)

## After definig our block of SNPs, run LMM:
assoc <- assocTestSingle(genoIterator,
                         null.model = nullmod,
                         BPPARAM = BiocParallel::SerialParam())
## Of this new Model, the most significant SNP is:
assoc |>
  as_tibble() |>
  arrange(desc(Score.pval)) |>
  slice(1)
GWASTools::manhattanPlot(assoc$Score.pval,
                         chromosome = rep(1, nrow(assoc)),
                         signif = 0.05/nrow(assoc))

# 10. Use the parameter scan.include of function fitNullMM() to
# perform association testing only in people from the UW group.
# Is the most significant variant the same as before?






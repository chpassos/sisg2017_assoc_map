# Admixture Mapping

# Lecture 5. Disease Trait Mixed Model GWAS
library(tidyverse)

# Instead of looking at a quantitative trait, now we are going to look at
# a BINARY outcome
# Modeled as D = 1 (diseased) or D = 0 (disease free)

# The basic logistic regression for a binary outcome is:
# logit[p(Di = 1)] = xTi β + giα, i = 1, . . . , n.

# Di: is the disease status of person i
# xi: is a vector of adjusting covariates (age, sex, etc.)
# β is a vector of their effects.
# gi: is the dosage or count of the genotype allele of interest.
# logit(u) = log[u/(1 − u)], is a function that ensures that
# estimated disease probabilities - u - will always be in the range
# (0, 1) (while logit(u) could be anything).


# Note: there is no “residual”. In linear regression, the residual
# induces the variability. Here, we directly model a probability,
# which induces a variability.


# The basic assumption in this model is that observations are
# "independent and identically distributed" (i.i.d.)
# If this assumption does not hold for your data, 
# Then use MIXED MODELS instead


# The logistic mixed model states that the disease probabilities of
# people who are somehow close or similar to each other, are
# more similar to each other than the disease probabilities of
# people who are not close or similar.


# One way to model this is using random effects.
# For example, if there was one source of such similarities
# between disease probabilities:
## logit[p(Di = 1)] = xTi β + giα + bi, i = 1, . . . , n,
# with bi a random effect, increasing/decreasing the baseline odds of the disease.

# Random effects reflect here similarity in disease odds 
# across individuals using a correlation structure

# As in linear regression, we use matrices to model 
# the correlations between random effects across individuals.

# Logistic mixed models are similar to logistic regression, 
# with the addition of random effects.

# So how are logistic mixed models practically different from linear
# mixed models?
## Variance components are still estimated (but no variance term corresponding to independent errors).
## There is no straight-forward interpretation of heritability based on variance components.
## Computationally, logistic models, and logistic mixed models, take longer to fit (compared to their linear counterparts).
## Logistic mixed models for more than a single correlation matrix are implemented in the software GMMAT and R package GENESIS (the same algorithm).

# We still fit a null model, as in linear regression.
# We use the null model and the "working traits" to test genotype-disease associations

# Take-home message: the “null model” for binary traits takes 4 times
# longer to fit than that for quantitative traits. Afterwards
# computation time is the same.


# In the past, people used linear mixed models instead of logistic
# mixed models.
# I Because it saved a lot of computation time.
# I Is it okay to use linear mixed models?
#   I Sometimes. But better not!
#   I Basic assumption made by linear mixed models: residual
# variance is the same for all people.
# I Basic assumption made by logistic model: if someone has a
# # probability p of disease, the variance of her outcome is p(1 − p).
# Chen et al. (2016, AJHG) showed in “Control for population
# structure and relatedness for binary traits in genetic association
# studies via logistic mixed models" that when
# I MAF differ between sub-populations in the study
# I Disease prevalence differ between these sub-populations
# I LMM test statistics can be either too significant (inflated), or
# too conservative (deflated).
## So basically, don't use linear mixed models in binary outcomes
## just because it's faster..

# Let's try it
# We first load our scanAnnotation object.
library(GWASTools)
library(GENESIS)

dir <- "/home/chpassos/association_mapping/sisg2017/"

scanAnnot <- getobj(file.path(dir,
                              "SISG_phenotypes.RData"))
scanAnnot
varLabels(scanAnnot)[1:4]
varLabels(scanAnnot)[5:8]
# Select covariates
covariates <- c("EV1", "EV2", "sex", "age", "group")
# Select the outcome
outcome <- "disease"
# Load correlation matrices
## Household matrix
HH.mat <- getobj(file.path(dir,
                           "SISG_houshold_matrix.RData"))
## Genetic Relatedness Matrix
kin.mat <- getobj(file.path(dir,
                            "SISG_relatedness_matrix.RData"))
# Put both matrices in a list
covMatList <- list(HH = HH.mat, kinship = kin.mat)



# Fit the null model
## As said in the previous lecture, GENESIS was updated and some functions
## became deprecated. Using newer versions here..
# Ideally, a generalized linear mixed model (GLMM) would be fit for a binary phenotype;
# however, fitting a GLMM is much more computationally demanding than fitting an LMM. 
# To provide a compuationally efficient approach to fitting such a model,
# fitNullModel uses the penalized quasi-likelihood (PQL) approximation to the GLMM (Breslow and Clayton).
nullmod <- fitNullModel(scanAnnot,
                        outcome = outcome, 
                        covars = covariates,
                        cov.mat = covMatList, 
                        family = "binomial")
## Looking at the results.
names(nullmod)
nullmod$varComp
nullmod$fixef
## Similar to results showed in the lecture...

# After estimating variance components in the “null model", they are assumed fixed.
# We now use this null model object in association testing.

# Load GDS
gds <- GdsGenotypeReader(file.path(dir,
                                   "SISG_snp_dosages.gds"))
# Load SNP annotation 
snpAnnot <- getobj(file.path(dir,
                             "SISG_snp_dosages_snpAnnot.RData"))
# Genotype Data: GDS + SNP Annot + Sample Annot
genoData <- GenotypeData(gds,
                         snpAnnot=snpAnnot, scanAnnot = scanAnnot)

# We cannot use a Wald test for logistic mixed models
# Score tests are "under the null", so they are realistic for GWAS based on logistic mixed models.
## Fit the model
# Again, as functions became deprecated, I'll do this on my own
## Generate the genotype block iterator:
genoIterator <- GenotypeBlockIterator(genoData, 
                                      snpBlock=1000)
# Then, fit the model
assoc <- assocTestSingle(genoIterator, null.model = nullmod,
                         BPPARAM = BiocParallel::SerialParam())
## Looking at the results
head(assoc)

close(gds)

### Exercises..

# 1.Use the results from the GWAS we ran on slide 19. Use the
# function qqPlot() from the GWASTools package to make a q-q
# plot figure of the p-values from the Score test.
GWASTools::qqPlot(assoc$Score.pval)


# 2. Use the following approximation between the Score and Wald
# test to obtain log(ORs) and ORs for the SNP effects:
### β = score/var(score)
### SE(β) = 1/ sqrt(var(score))
### OR = exp(β).
# Which SNP has the highest odds ratio (OR)?
### R: 
beta <- assoc$Score/var(assoc$Score)
se_beta <- 1/sqrt(var(assoc$Score))
OR <- exp(beta)
OR |>
  as_tibble() |>
  mutate(row_n = row_number()) |>
  arrange(desc(value)) |>
  slice(1) |>
  pull(row_n) # SNP number 1359 has the highest OR. Which SNP is it?
assoc[assoc$variant.id == 1359, ]
pData(snpAnnot) |>
  as_tibble() |>
  filter(snpID == 1359) # I think that is our SNP.

# 3. Use the function manhattanPlot() from the GWASTools
# package to generate a Manhattan plots for the Score test
# p-values.  Did any of the SNPs achieve genome-wide significance? 
# (p-value= 5 × 10−8) array-wide significance?
GWASTools::manhattanPlot(assoc$Score.pval,
                         chromosome = rep(1, nrow(assoc)))
## R: Apparently, there are three SNPs that achieved genomewide significance!

# 4. Which variant is most associated with "disease" among all variants?
assoc |>
  as_tibble() |>
  arrange(desc(Score.pval)) |>
  slice(1)


# 5. Run linear mixed model GWAS instead of logistic, 
# treating "disease" as a quantitative trait.
## by first fitting a new null model using fitNullMM()
nullmod <- fitNullModel(scanAnnot,
                        outcome = outcome,
                        covars = covariates,
                        cov.mat = covMatList,
                        family = "gaussian")
genoIterator <- GenotypeBlockIterator(genoData, 
                                      snpBlock=10000)
assoc2 <- assocTestSingle(genoIterator,
                         null.model = nullmod,
                         BPPARAM = BiocParallel::SerialParam())
## Compare the p-values obtained by the two methods.
## You can use a scatter plot or a q-q plot.
assoc_binary_tbl <- assoc |>
  as_tibble()
assoc2_quant_tbl <- assoc2 |>
  as_tibble()

assoc_binary_tbl |>
  inner_join(assoc2_quant_tbl, by = "variant.id") |>
  rename("binary_pvalue" = "Score.pval.x",
         "quant_pvalue" = "Score.pval.y") |>
  ggplot(aes(binary_pvalue, quant_pvalue)) +
  geom_point() +
  labs(x = "Pvalues in Binary Mixed Model",
       y = "Pvalues in Linear Mixed Model") +
  theme_bw()

## Do the two GWAS have the same top SNPs?
GWASTools::manhattanPlot(assoc2$Score.pval,
                         chromosome = rep(1, nrow(assoc2)))

assoc |>
  as_tibble() |>
  arrange(desc(Score.pval)) |>
  slice(1:10)
assoc2 |>
  as_tibble() |>
  arrange(desc(Score.pval)) |>
  slice(1:10)
### R: No.

# 6. Use the parameter scan.include of function fitNullMM to fit
# perform association testing only in people from the UNC group.
## Does the this GWAS have the same top SNP as the GWAS that was run on all participants?
  



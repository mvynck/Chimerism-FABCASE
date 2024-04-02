##############################
# LOAD LIBRARIES
##############################

library(poibin)

##############################
# LOAD FUNCTIONS
##############################

source("fabceseFuns.R")

##############################
# EXAMPLE
##############################

# empirical MAFs
# obtained from Vynck et al. 2022
#   see https://github.com/mvynck/Chimerism-nMarkers/tree/main/Data
MAF1 <- c(0.4059633, 0.4357798, 0.3302752, 0.4128440, 0.4564220, 0.4839450, 0.3990826, 0.4403670, 0.4747706,
          0.3348624, 0.4311927, 0.4449541, 0.3922018, 0.3876147, 0.4128440, 0.3876147, 0.4380734, 0.4770642,
          0.4633028, 0.4885321, 0.4380734, 0.3394495, 0.4220183, 0.4380734)

# point estimates for fraction of transplants with at least three informative markers
nInf <- 3
(estimateMUD <- sapply(1:length(MAF1), cdfVec, MAFtoProb(MAF1, donrec="mud", type="i"), nInf - 1))
(estimateSIB <- sapply(1:length(MAF1), cdfVec, MAFtoProb(MAF1, donrec="sib", type="i"), nInf - 1))
(estimatePC <- sapply(1:length(MAF1), cdfVec, MAFtoProb(MAF1, donrec="pc", type="i"), nInf - 1))

# limits
B=999 # number of bootstraps to perform
AF.dep <- 2*(95+56+11) # number of alleles used for MAF estimation
b <- matrix(MAFvec(MAF1, B*length(MAF1), AF.dep), nrow = B, byrow = TRUE) # generate bootstrapped MAFs

# get the 95% CI for MUD transplants, type I markers only
#   note that the getCI function returns a sequence of informativity rates
#   use all 24 markers, i.e., obtain element [24] if interested in the full set
MUDci <- c(getCI(MAF1, b = b, nInf = 3, rel="mud", type="i", quantiles = c(0.025, 0.975))$lower[24],
           getCI(MAF1, b = b, nInf = 3, rel="mud", type="i", quantiles = c(0.025, 0.975))$upper[24])
# get the 95% CI for sibling transplants, type I markers only
SIBci <- c(getCI(MAF1, b = b, nInf = 3, rel="sib", type="i", quantiles = c(0.025, 0.975))$lower[24],
           getCI(MAF1, b = b, nInf = 3, rel="sib", type="i", quantiles = c(0.025, 0.975))$upper[24])
# get the 95% CI for parent-child transplants, type I markers only
PCci <- c(getCI(MAF1, b = b, nInf = 3, rel="pc", type="i", quantiles = c(0.025, 0.975))$lower[24],
           getCI(MAF1, b = b, nInf = 3, rel="pc", type="i", quantiles = c(0.025, 0.975))$upper[24])

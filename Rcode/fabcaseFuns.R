convertMAF <- function(MAF){
  # allele frequency (range 0-1) to minor allele frequency
  return(min(MAF, 0.5)-max(MAF-0.5, 0))
}

MAFvec <- function(MAF, reps, n){
  # resample MAFs with given number of alleles (n)
  # obtain "reps" repeats (number of bootstraps)
  rbinom(reps, n, MAF)/n
}

cdfVec <- function(use, probs, nmark){
  # obtain a probability from the poisson-binomial 
  #   cumulative distribution function
  1-ppoibin(nmark, probs[1:use])
}

MAFtoProb <- function(MAFvec, donrec = c("sib", "mud", "pc"), type=c("i", "i-ii")){
  # convert minor allele frequency to probability of being informative
  #   take into account:
  #     relatedness: sibling (sib), unrelated (mud) or parent-child (pc)
  #     markers types: only type-I markers (recommended), type-I and -II markers
  #       see doi 10.1016/j.cca.2023.117452 for type-I/-II marker rationale
  if(donrec == "mud"){
    if(type=="i"){
      probs <- 2*(MAFvec^2*(1-MAFvec )^2+(MAFvec*(1-MAFvec )^3)+MAFvec^3*(1-MAFvec ))
    }
    if(type=="i-ii"){
      probs <- 2*(MAFvec^2*(1-MAFvec )^2+2*((MAFvec*(1-MAFvec )^3)+MAFvec^3*(1-MAFvec )))
    }
  }
  if(donrec == "sib"){
    if(type=="i"){
      probs <- 1/2*(3*MAFvec^2*(1-MAFvec )^2+2*(MAFvec*(1-MAFvec )^3)+2*MAFvec^3*(1-MAFvec ))
    }
    if(type=="i-ii"){
      probs <- 1/2*(5*MAFvec^2*(1-MAFvec )^2+4*((MAFvec*(1-MAFvec )^3)+MAFvec^3*(1-MAFvec )))
    }
  }
  if(donrec == "pc"){
    if(type=="i"){
      probs <- 2*MAFvec^2*(1-MAFvec )^2+(MAFvec*(1-MAFvec )^3)+MAFvec^3*(1-MAFvec )
    }
    if(type=="i-ii"){
      probs <- 2*(2*MAFvec^2*(1-MAFvec )^2+ (MAFvec*(1-MAFvec )^3)+MAFvec^3*(1-MAFvec ))
    }
  }
  return(probs)
}

getCI <- function(maf.est, b, nInf = 3, rel = c("sib", "mud", "pc"), type = c("i", "i-ii"), quantiles = c(0.025, 0.975)){
  # obtain the percentile confidence interval using estimated MAFs (maf.est)
  #   with b the boostrapped MAFs (with B rows, length(maf.est) columns)
  #   informativity defined as at least nInf informative markers
  #   take into account relatedness (rel) and type (type), see MAFtoProb function
  #   quantiles indicate the bootstrap quantiles to obtain, e.g., symmetric 95% CI by default
  B <- nrow(b)
  d <- matrix(0, nrow = B, ncol = length(maf.est))
  for(iter in 1:B){
    d[iter,] <- sapply(1:length(maf.est), cdfVec, MAFtoProb(b[iter,], donrec=rel, type=type), nInf - 1)
  }
  upper <- vector("numeric", length(maf.est))
  lower <- vector("numeric", length(maf.est))
  # get percentile CIs
  for(mark in 1:length(maf.est)){
    upper[mark] <- sort(d[, mark])[quantiles[2]*(B+1)]
    lower[mark] <- sort(d[, mark])[quantiles[1]*(B+1)]
  }
  return(list(upper=upper,lower=lower))
}
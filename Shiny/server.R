library(rhandsontable)
library(shinydashboard)
library(shinyBS)
library(shiny)
library(ggplot2)
library(rmarkdown)
library(poibin)

#number of alleles
numSamples<-20

##########################
# FUNCTIONS
##########################

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

calcInf <- function(MAF, minInf, AF.dep = 1000, type="i"){
  # calculate informativity for the given MAF vector
  # with informativity defined as observing at least minInf markers
  # and with AF.dep alleles observed
  # and for type "type" markers (type-I = "i" recommended)

  # in the Shiny app, do 999 bootstraps by default
  B = 999

  maf.est <- MAF
  # get bootstrapped MAFs
  b <- matrix(MAFvec(maf.est, B*length(maf.est), AF.dep), nrow = B, byrow = TRUE)
  
  # estimate informativity rates for all relatedness types
  # first MUD, then SIB, then PC

  # MUD
  estimate.mud <- estimate <- sapply(1:length(maf.est), cdfVec, MAFtoProb(maf.est, donrec="mud", type=type), minInf-1)
  d <- matrix(0, nrow = B, ncol = length(maf.est))
  for(iter in 1:B){
    d[iter,] <- sapply(1:length(maf.est), cdfVec, MAFtoProb(b[iter,], donrec="mud", type=type), minInf-1)
  }
  quantiles <- c(0.025, 0.975)
  upper <- vector("numeric", length(maf.est))
  lower <- vector("numeric", length(maf.est))
  for(mark in 1:length(maf.est)){
    upper[mark] <- sort(d[, mark])[quantiles[2]*(B+1)]
    lower[mark] <- sort(d[, mark])[quantiles[1]*(B+1)]
  }
  (upper.final <- upper)
  (lower.final <- lower)
  
  # SIBLING
  estimate.sib <- estimate <- sapply(1:length(maf.est), cdfVec, MAFtoProb(maf.est, donrec = "sib", type = type), minInf-1)
  d <- matrix(0, nrow = B, ncol = length(maf.est))
  for(iter in 1:B){
    d[iter,] <- sapply(1:length(maf.est), cdfVec, MAFtoProb(b[iter,], donrec = "sib", type = type), minInf-1)
  }
  quantiles <- c(0.025, 0.975)
  upper <- vector("numeric", length(maf.est))
  lower <- vector("numeric", length(maf.est))
  for(mark in 1:length(maf.est)){
    upper[mark] <- sort(d[, mark])[quantiles[2]*(B+1)]
    lower[mark] <- sort(d[, mark])[quantiles[1]*(B+1)]
  }
  (upper.final <- upper)
  (lower.final <- lower)
  
  # Parent-child
  estimate.pc <- estimate <- sapply(1:length(maf.est), cdfVec, MAFtoProb(maf.est, donrec = "pc", type), minInf-1)
  d <- matrix(0, nrow = B, ncol = length(maf.est))
  for(iter in 1:B){
    d[iter,] <- sapply(1:length(maf.est), cdfVec, MAFtoProb(b[iter,], donrec = "pc", type), minInf-1)
  }
  quantiles <- c(0.025, 0.975)
  upper <- vector("numeric", length(maf.est))
  lower <- vector("numeric", length(maf.est))
  for(mark in 1:length(maf.est)){
    upper[mark] <- sort(d[, mark])[quantiles[2]*(B+1)]
    lower[mark] <- sort(d[, mark])[quantiles[1]*(B+1)]
  }
  (upper.final <- upper)
  (lower.final <- lower)  
  
  return(list(estimate.mud,
              lower.final.mud,
              upper.final.mud,
              estimate.sib,
              lower.final.sib,
              upper.final.sib,
              estimate.pc,
              lower.final.pc,
              upper.final.pc))
}


adjust.data <- function(df_){
    # adjust the data frame on input change

    nr <- nrow(df_)
    if(nr<input$nsample){
      df_[(nr+1):input$nsample,1:(ncol(df_)-1)]<-0.0
      df_$Include[(nr+1):input$nsample]<-TRUE
      row.names(df_)<-1:nrow(df_)
    }else{
      df_<- df_[1:input$nsample,]
    }
    df_$MAF<-as.numeric(df_$MAF)
    df_$Include<-as.logical(df_$Include)  
    # get the informativity rates and CIs
    restab <- calcInf(df_$MAF[df_$Include == TRUE], 
                      input$ninf,
                      input$nallele,
                      input$eligible)
    # organize the informativity rates and CIs
    restab <- data.frame(Markers=1:length(restab[[1]]),
                         Informativity.Unrelated=restab[[1]],
                         Unrel.95LL=restab[[2]],
                         Unrel.95UL=restab[[3]],
                         Informativity.Sibling=restab[[4]],
                         Sib.95LL=restab[[5]],
                         Sib.95UL=restab[[6]],
                         Informativity.PC=restab[[7]],
                         PC.95LL=restab[[8]],
                         PC.95UL=restab[[9]])
    # prettify the table header
    colnames(restab) <- c("Markers", "Informativity (unrelated)",
                          "95% CI LL", "95% CI UL",
                          "Informativity (sibling)",
                          "95% CI LL", "95% CI UL",
                          "Informativity (parent-child)",
                          "95% CI LL", "95% CI UL")
    # show highest informativity rates first, most informative
    restab <- restab[order(restab$Markers, decreasing = TRUE),]
    # render the table nicely
    output$tableInf <- renderTable(restab, spacing="xs", digits = 3, striped=TRUE)
    return(df_)
}


##################################
# SERVER BODY
##################################

server<-function(input, output,session)({

values <- reactiveValues(data=data.frame(MAF = c(rep(0.000,(numSamples))), Include = TRUE))

observe({
	if(!is.null(input$hot)){
		values$data <- hot_to_r(input$hot)
		values$data<-adjust.data(isolate(values$data))
	}
  
})

observeEvent(input$nsample,({
	
	if(!is.null(input$hot)){
		values$data <- hot_to_r(input$hot)
	}
  values$data<-adjust.data(values$data)
  	
})
)


observeEvent(input$exampleLoad,{
	data.temp<-read.csv("exampledata.csv")
	updateSliderInput(session,"nsample",value=nrow(data.temp))
	values$data<-data.temp
})


output$hot <- renderRHandsontable({
        rhandsontable(values$data,readOnly=F, useTypes= TRUE) %>%
      		hot_table(highlightCol = TRUE, highlightRow = FALSE) %>%
          hot_col(col="MAF", format = "0.000")
})   

})

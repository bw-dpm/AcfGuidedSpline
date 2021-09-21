library(mgcv)

###################
# Validate formula
###################
validate.formula <- function(form, DEBUG=FALSE){
  valid          <- TRUE
  fmla.out       <- NULL
  fmla.spl.type  <- NULL
  fmla.spl.term  <- NULL
  fmla.expl      <- NULL 
  fmla.covs      <- NULL
  
  f.str <- as.character(form)
  outcome <- f.str[2]
  other <- unlist(strsplit(gsub(" ", "", f.str[3], fixed = TRUE),"\\+"))
  
  w.covs <- FALSE
  if(length(other) > 1){w.covs <- TRUE}
  
  # Verify there is only one spline term in the formula ("s" or "bs")
  bss <- length(grep("^bs\\(", other, value=T))
  ss  <- length(grep("^s\\(", other, value=T))
  
  if(DEBUG){
    print(paste("bss=",bss))
    print(paste("ss=",ss))
    print(paste("valid=", valid))
  }
  
  # Issue error if formula is not valid
  if(bss == 0 & ss == 0){
    print("Error: no spline terms found in formula")
    valid <- FALSE
  }else if((bss != 0 & bss > 1) | (ss != 0 & ss > 1)){
    print("Error: only 1 spline term is allowed in formula"); 
    valid <- FALSE
  }else if(bss + ss > 1){
    print("Error: only 1 spline term is allowed in formula"); 
    valid <- FALSE
  }else{
    # Formula is ok
    if(bss == 1){
      fmla.spl.type <- "b.spline"
    }else{
      fmla.spl.type <- "s.spline"
    }
    
    # Spline term
    fmla.spl.term <- grep("s\\(", other, value=T)
    
    # Expl. variable
    fmla.expl <- strsplit(strsplit(fmla.spl.term,"\\(")[[1]][2], ")")[[1]][1]
    
    # Covariates
    if(w.covs){
      fmla.covs <- other[!other %in% fmla.spl.term]
    }
  }
  
  if(DEBUG){
    print(paste("valid (post-check):", valid))
    print(paste("fmla.spl.type:", fmla.spl.type))
    print(paste("fmla.spl.term:", fmla.spl.term))
    print(paste("fmla.expl:", fmla.expl))
    print(paste("w.covs:", w.covs))
    print(paste("fmla.covs:", fmla.covs))
  }
  
  return(list(fmla.valid=valid, fmla.out=outcome, fmla.spl.type=fmla.spl.type,
              fmla.spl.term=fmla.spl.term, fmla.expl=fmla.expl, fmla.covs=fmla.covs))  
}


####################
# atf.mgcv.spline
####################
#' @title Autocorrelation Guided Spline Regression
#'
#' @description  Spline optimization by residual autocorrelation analysis
#' @param formula  A regression formula with a spline term:
#' \itemize{
#' \item "bs": non-penalized
#' \item "s": penalized
#' }
#' @param dframe   Data frame
#' @param df.range The range of df parameter to be optimized
#' @param boxlag   The lag value for the Ljung Box test (default=20)
#' @param spl.type The spline basis type: "cr" cubic, "tp" thin plate (default="cr")
#' @param debug    Debug mode: prints debug information when set to TRUE (default=FALSE).
#' @return A list:
#' \itemize{
#' \item{"best.par"  : }{ Optimized parameter value}
#' \item{"best.pval" : }{ Ljung-Box pvalue corresponding to best.par}
#' \item{"box.pvals" : }{ Vector of p-values relative to the range of parameters}
#' \item{"best.fit"  : }{ GAM regression object corresponding to the best.par} 
#' }
#' @keywords atf.mgcv.spline
#' @export
#' @examples
#' # Create a "iid + drift" signal
#' N <- 1000
#' drift <- sin(29.6*((1:N/N) + 0.2))/((1:N)/N + 0.2)
#' iid <- runif(n=N, min=-1.0, max=1.0)
#' y <- iid + drift
#' mydf <- data.frame(y=y,x=1:N)
#'
#' # Add "random" covariates
#' mydf$cov1 <- sample(N)
#' mydf$cov2 <- sample(N)
#' 
#' # Run atf.mgcv.spline
#' fit <-   atf.mgcv.spline(y ~ bs(x) + cov1 + cov2, dframe=mydf, 3:100, DEBUG=TRUE)
#'
#' # plot the "true" drift
#' lines(drift, col="blue")
#'
#' # Plot LB p-values
#' plot(fit$box.pvals, type="l")
########################################################################
atf.mgcv.spline <- function(formula, dframe, df.range, spl.type="cr", boxlag=20, DEBUG=FALSE){
  
  # Print function call
  print(match.call())
  
  # Validate formula
  fmla <- validate.formula(formula, DEBUG)
  
  # Exit if formula is not valid
  if(!fmla$fmla.valid){
    return(0)
  }
  
  # Get the expl. variable
  expl.str <- fmla$fmla.expl
  
  # Sort data frame by expl. variable 
  dframe <- dframe[order(dframe[[expl.str]]), ]
  
  # Extract explanatory variable  
  expl <- dframe[[expl.str]] 
  
  # Extract covariates
  covs <- ""
  str(fmla$fmla.covs)
  if(!is.null(fmla$fmla.covs)){
    covs <- paste(" + ", paste(fmla$fmla.covs, collapse=" + "), sep="")
    if(DEBUG){
      print(paste("covs:", covs))
    }
  }
  
  # Spline basis type
  if(spl.type == "cr"){
    bs <- "bs='cr'"
  }else if(spl.type == "tp"){
    bs <- "bs='tp'"
  }
  
  # Get the fx term (fx=TRUE for unpenalized spline)
  fx <- TRUE
  if(fmla$fmla.spl.type == "s.spline"){
    fx <- FALSE
  }
  
  
  # Loop through the dfs
  box.pvals <- c(); best.fit <- NULL; best.k <- c();
  resp <- c(); spl.term <- c()
  
  for(k in df.range){
    spl <- paste("s(", expl.str,", ",bs,", k=",k, ", fx=",fx,")" , covs, sep="")
    new.f <- paste(fmla$fmla.out, "~", spl, sep=" ")
    
    if(DEBUG){
      print(paste("Formula: ", new.f))
    }
    
    # Run the regression
    fit <- mgcv::gam(as.formula(new.f), data=dframe)
  
    # Residuals    
    res <- residuals(fit)
    
    # Ljung-Box test of autocorrelations on the residuals 
    pvalue <- Box.test(res, lag=boxlag, type = "Ljung-Box")$p.value
    box.pvals <- c(box.pvals,pvalue)
    
    # Is this the "best" fit?
    if(max(box.pvals) == pvalue)
    {
      best.fit <- fit
      best.k   <- k
    }
  }
  
  # Add names to vector of pvalues
  names(box.pvals) <- df.range
  
  if(DEBUG){
    print(paste("Best param:", best.k))
    plot(expl, dframe[[fmla$fmla.out]], col="grey", xlab=expl.str, ylab=fmla$fmla.out)
    lines(expl, best.fit$fitted.values, col="red")
  }
  
  if(max(box.pvals) < 0.05){warning("P-value relative to best df < 0.05")}
  return(list(best.par=best.k, best.pval=max(box.pvals), box.pvals=box.pvals, best.fit=best.fit))
}

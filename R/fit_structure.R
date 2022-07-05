# VBGfit package
# 5 july 2022

# Description -------------------------------------------------------------
#' Fit a population structure to data.
#'
#' This function fits a timeseries of population structures to a dataset with size at age data.
#' To do so, it is assumed that the size at a given age at a given point in time is normally distributed and individuals follow a Von Bertalanffy growth curve.
#' The fitted asymptotic size is either constant or is allowed to vary between intervals.
#'
#' The model is fitted by maximizing the log likelyhood using the \code{link{nloptr}} package.
#'
#' @import nloptr
#' @import stats
#'
#' @param samples A discrete (integer) value indicating the timepoint of the observation. Timepoints should have a constant interval and consecutive timepoints should have consecutive values. Or a dataset with 3 or 4 colums containing timepoints, ages, sizes and (optionally) weights as described.
#' @param ages A discrete (integer) value indicating the age of the observed individual. Ages should be on the same timescale as the timepoints. Unused if samples is a dataframe.
#' @param sizes A numeric value for the size of the observed individual. Unused if samples is a dataframe.
#' @param weights (optional) The weight of an observation. For example the number of times a given age-size combination is observed. If weights are not provided all weights are set to 1. Unused if samples is a dataframe.
#'
#' @param initpars A vector with initial parameters. Either an 5 element vector containing the elements and order as described below or  an elaborate vector with all parameters as created by \code{\link{makepoppars}}.
#' \itemize{
#' \item{the mean asymptotic size}
#' \item{the variance in asymptotic size}
#' \item{the mean length at the youngest age in the model}
#' \item{the variance in the size at the youngest age in the model}
#' \item{The growth rate of the Von Bertalanffy growth curve}
#' }
#'
#' @param model A specific character expression specifying whether the asymptotic size is allowed to vary between timepoints.
#' \describe{
#' \item{\code{CONS}}{The asymptotic size does not vary between years}
#' \item{\code{VARY}}{The asymptotic size is allowed vary between years}
#' \item{\code{BOTH}}{Fit both models}
#' }
#'
#' @param feedinggroups (optional) Indicate whether multiple feeding levels need to be estimated for different size- or agegroups.
#' \describe{
#' \item{\code{AGE}}{Feedinggroups are clustered by age}
#' \item{\code{SIZE}}{Feedingroups are clustered by size. note that the average size of a cohort is used and individual variation within a cohort is not considered when assigning the feedinggroup}
#' }
#'
#' @param feedingbounds (Only if feedingroups are used) Vector with the boundaries of the feeding groups. The start of the fist feedinggroup and the end of the last feedinggroups do not need to be indicated.
#'
#' @param logscale Boolean to indication whether the parameter space is explored on a log scale.
#' @param shrink Boolean to indicating whether the average size is allowed to shrink or should kept constant when the size exceeds the asymptotic size.
#' @param opts options passed to the \code{\link{nloptr}} function.
#'
#' @return list with detailes o the fitted models with a constant asymptotic size (CONS) and a varying asymptotic size (VARY).
#' For every fitted model the lists contain the following elements:
#' \describe{
#' \item{samples}{A dataframe with the observations used to fit the model.}
#' \item{globpars}{A vector with global parameters used for the modelfit containing
#' number of parameters,
#' number of timepoints,
#' lowest timepoint,
#' number of ages,
#' lowest age,
#' number of samples,
#' Shrinking allowance indication,
#' Log transformation indication,
#' feedinggroup indication
#' and number of feedinggroups}
#' \item{poppars}{Fitted model parameters in a vector as constructed by \code{\link{makepoppars}}.}
#' \item{asympsize}{Estimated mean and and variance of the asymptotic size per feedinggroup per year.}
#' \item{popstruct}{A list with the predicted mean and variance in size for every timepoint and age as constructed by \code{\link{predict_structure}}.}
#' \item{minLL}{The negative log likelyhood of the fitted model.}
#' \item{niter}{The number of iterations needed for the optimization.}
#' \item{initpars}{The initial parameters used for the optimization in a vector as constructed by \code{\link{makepoppars}}.}
#' \item{feedingbounds}{The boundaries of the feedinggroups used for the model.}
#' \item{opts}{The options used by \code{\link{nloptr}}.}
#' }
#'
#' @export
#'

fit_structure <- function(samples, ages=NA, sizes=NA, weights=NA, initpars, model = "BOTH", feedinggroups = NA, feedingbounds = NA, logscale = TRUE, shrink =TRUE,  opts = list(algorithm = "NLOPT_LN_SBPLX", xtolrel = 1E-6, maxeval = 1E6)){

# Check and process input -------------------------------------------------

  #### Check samples

  # check wether samples consists of single entries or a dataframe
  if(is.vector(samples)){
    # check ages
    if(is.na(ages) || !is.vector(ages) || length(ages) != length(samples)){
      stop("ages is missing or not a vector with the same length as samples")
    }
    #check sizes
    if(is.na(sizes) || !is.vector(sizes) || length(sizes) != length(samples)){
      stop("sizes is missing or not a vector with the same length as samples")
    }
    #check weights

    #set weights to 1 if not provided
    if(length(weights) == 1 && sum(is.na(weights))==1){
      weights <- 1
      warning("Equal weights assumed")
    }
    else if(  !is.vector(weights) || length(weights) != length(samples) ){
      stop("weights not a vector with the same length as samples, Equal weights assumed")
      weights <- 1
    }

    #combine into dataset
    samples = data.frame(timepoints = samples, ages = ages, sizes=sizes, weights = weights)
    }
  else if(is.data.frame(samples) && ncol(samples) == 3){
    names(samples) <- c("timepoints","ages","sizes")
    samples$weights <- 1
    warning("Equal weights assumed")
  }
  else if(!is.data.frame(samples) || ncol(samples) != 4){
    stop("Samples is not a vector with timepoints or a dataframe with 3 or 4 four columns containing timepoints, ages, sizes and (optional) weights")
    samples = data.frame(timepoints = NA, ages = NA, sizes=NA, weights = NA)
  }

  # Convert samples to correct datatype
  samples[,1] <- as.integer(samples[,1])
  samples[,2] <- as.integer(samples[,2])
  samples[,3] <- as.numeric(samples[,3])
  samples[,4] <- as.numeric(samples[,4])

  # Select complete cases
  if(sum(complete.cases(samples)) == 0 ){
    stop(paste("Incomplete observations removed.", as.character(sum(complete.cases(samples))), "observations remaining."))
    samples <- samples[complete.cases(samples),]
  }
  else if(sum(complete.cases(samples)) != nrow(samples) ){
    warning(paste("Incomplete observations removed.", as.character(sum(complete.cases(samples))), "observations remaining."))
    samples <- samples[complete.cases(samples),]
  }

  # obtain glob parameters from samples
  ntime <- max(samples[,1])-min(samples[,1])+1
  mintime <- min(samples[,1])
  nages <- max(samples[,2])-min(samples[,2])+1
  minage <- min(samples[,2])
  nobs <-  nrow(samples)

  if(ntime <= 1){
    stop("The dataset should contain more than 1timepoint with complete data")
  }

  if(nages <= 1){
    stop("The dataset should contain more than 1 ageclass with complete data")
  }

  #### Check other input

  # check modeltype
  if(model != "CONS" && model != "VARY" && model != "BOTH"){
    warning("model is not one of 'CONS', 'VARY' or 'BOTH', and is reset to default: 'BOTH'")
    model = "BOTH"
  }

  # check shrink
  if( shrink !=TRUE && shrink !=FALSE ){
    warning("shrink should be TRUE or FALSE, set to TRUE")
    shrink = TRUE
  }
  # convert shrink to integer
  if(!shrink) shrinkpar = 0
  else shrinkpar = 1

  # check logscale
  if( logscale !=TRUE && logscale !=FALSE ){
    warning("Logscale should be TRUE or FALSE, set to TRUE")
    logscale = TRUE
  }
  #convert logscale to integer
  if(!logscale) logpar = 0
  else logpar = 1

  # check feedinggroup
  if(is.na(feedinggroups)) grouppar = 0
  else if( feedinggroups == "AGE") grouppar = 1
  else if( feedinggroups == "SIZE") grouppar = 2
  else{
    warning("Feedinggroups is not one of 'AGE' or 'SIZE' and is ignored")
    grouppar = 0
  }

  # check feedingbounds
  if(grouppar > 0){
    feedingbounds <- unique(feedingbounds[!is.na(feedingbounds)])
    ngroups <- length(feedingbounds)+1
    if(ngroups == 1){
      warning("No correct boundaries for feedinggroups found. Feedinggroups ignored")
      grouppar = 0
      feedingbounds = 0
    }
    else{
      message(paste(as.character(ngroups-1),"unique boundaries of feedinggroups found, continuing with",as.character(ngroups),"feedinggroups"))
    }
  }
  else{
    ngroups = 1;
    feedingbounds = 0;
  }


  # check initpars
  npars <- length(initpars)

  if(anyNA(initpars)){
    stop("NA values are not allowed in initpars")
  }
  else if(sum(initpars <= 0) && logpar==1){
    stop(paste("initpars may not contain values <= 0 if parameters are logtrasformed (logscale=TRUE)"))
  }
  else if(npars == 2*ngroups + 2*ntime + 2*(nages-1) + 1  && model != "CONS"){
    warning("Model set to 'CONS' based on length of initpars")
    model = "CONS"
  }
  else if(npars == 2*(ntime-1)*ngroups + 2*ntime + 2*(nages-1) + 1 && model != "VARY"){
    warning("Model set to 'VARY' based on length of initpars")
    model = "VARY"
    }
  else if( npars != 5 ){
    stop("Size of initpars is oncorrect. Initpars should be a vector of length 5 or a vector as produced by 'makepoppars()'")
  }


# Fit model to data -------------------------------------------------

  #### make output list
  parout <- list(CONS = NA, VARY = NA)

  #### preform fit with constant feedinglevel
  if(model == "CONS" || model == "BOTH"){

    # create vector with all startpars
    if(npars == 5){
      initpars_long <- makepoppars(fRlm_avg = initpars[1], fRlm_var = initpars[2], len0_avg = initpars[3], len0_var = initpars[4], rB = initpars[5], time0_avg=NULL, time0_var=NULL, ntime = ntime, nages = nages,  model="CONS", ngroups = ngroups)
    }
    else initpars_long <- initpars

    # logtransform if needed
    if(logpar == 1) initpars_long2 <- log(as.numeric(initpars_long))
    else initpars_long2 <- as.numeric(initpars_long)

    # create vector with global variables
    globpars <- as.integer(c(length(initpars_long), ntime, mintime, nages, minage, nobs, shrinkpar, logpar, grouppar, ngroups))

    # preform actual fit
    print("fitting model with constant asymptotic size")

    consfitout <- nloptr(x0 = initpars_long2,
                         eval_f = minLL_structure_ ,
                         opts = opts,
                         globpars = globpars,
                         feedingbounds = feedingbounds,
                         samples = samples,
                         modeltype = as.integer(0))

    print("modelfit with constant asymptotic size succeeded")

    # extract backtransform parameter solution
    if(logpar == 1) consparsout <- exp(consfitout$solution)
    else consparsout <- consfitout$solution

    # make dataframe with asymptotic sizes
    asympsize <- data.frame(year = seq(from = mintime, to = mintime+ntime - 2, by = 1))

    for(i in 1:ngroups){
      groupasympsize <- data.frame(avg = rep(consparsout[i],ntime-1), var=rep(consparsout[ngroups + i],ntime-1))
      names(groupasympsize) <- c(paste0("avg_group",as.character(i)), paste0("var_group",as.character(i)) )
      asympsize <- cbind(asympsize, groupasympsize)
    }

    # save needed output
    consout <- list(samples = samples,
                    globpars = globpars,
                    poppars =  consparsout,
                    asympsize = asympsize,
                    popstruct = predict_structure_(consfitout$solution, globpars, as.integer(0), feedingbounds),
                    minLL = consfitout$objective,
                    niter = consfitout$iterations,
                    initpars = initpars_long,
                    feedingbounds = feedingbounds,
                    opts = consfitout$options,
                    allout =  consfitout)

    parout$CONS <- consout
  }


  #### preform fit with varying feedinglevel
  if(model == "VARY" || model == "BOTH"){

    if(npars == 5){
      initpars_long <- makepoppars(fRlm_avg = initpars[1], fRlm_var = initpars[2], len0_avg = initpars[3], len0_var = initpars[4], rB = initpars[5], time0_avg=NULL, time0_var=NULL, ntime = ntime, nages = nages,  model="VARY", ngroups = ngroups)
    }
    else initpars_long <- initpars

    # logtransform if needed
    if(logpar == 1) initpars_long2 <- log(as.numeric(initpars_long))
    else initpars_long2 <- as.numeric(initpars_long)

    # create vector with global variables
    globpars <- as.integer(c(length(initpars_long), ntime, mintime, nages, minage, nobs, shrinkpar, logpar, grouppar, ngroups))

    # preform actual fit
    print("fitting model with varying asymptotic size")

    varyfitout <- nloptr(x0 = initpars_long2,
                         eval_f = minLL_structure_ ,
                         opts = opts,
                         globpars = globpars,
                         feedingbounds = feedingbounds,
                         samples = samples,
                         modeltype = as.integer(1))

    print("modelfit with varying asymptotic size succeeded")

    # extract backtransform parameter solution
    if(logpar == 1) varyparsout <- exp(varyfitout$solution)
    else varyparsout <- varyfitout$solution

    # make dataframe with asymptotic sizes
    asympsize <- data.frame(year = seq(from = mintime, to = mintime+ntime - 2, by = 1))

    for(i in 1:ngroups){
      groupasympsize <- data.frame(avg =varyparsout[((i-1)*(ntime-1)+1):(i*(ntime-1))], var =varyparsout[(ngroups*(ntime-1)+(i-1)*(ntime-1)+1):(ngroups*(ntime-1)+i*(ntime-1))])
      names(groupasympsize) <- c(paste0("avg_group",as.character(i)), paste0("var_group",as.character(i)) )
      asympsize <- cbind(asympsize, groupasympsize)
    }

    # save needed output
    varyout <- list(samples = samples,
                    globpars = globpars,
                    poppars =  varyparsout,
                    asympsize = asympsize,
                    popstruct = predict_structure_(varyfitout$solution, globpars, as.integer(1), feedingbounds),
                    minLL = varyfitout$objective,
                    niter = varyfitout$iterations,
                    initpars = initpars_long,
                    feedingbounds = feedingbounds,
                    opts = varyfitout$options,
                    allout =  varyfitout)

    parout$VARY <- varyout
  }

  # return output
  return(parout)

}

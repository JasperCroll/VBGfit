# VBGfit package
# 5 july 2022

#' Calculate the negative log-likelyhood for a model with given parameters for the given dataset.
#'
#' @param samples A discrete (integer) value indicating the timepoint of the observation. Timepoints should have a constant interval and consecutive timepoints should have consecutive values. Or a dataset with 3 or 4 columns containing timepoints, ages, sizes and (optionally) weights as described.
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
#'
#' @param feedingbounds (Only if feedingroups are used) Vector with the boundaries of the feeding groups. The start of the fist feedinggroup and the end of the last feedinggroups do not need to be indicated.
#'
#' @param logscale Boolean to indication whether the parameter space is explored on a log scale.
#' @param shrink Boolean to indicating whether the average size is allowed to shrink or should kept constant when the size exceeds the asymptotic size.
#'
#' @return The negative log likelyhood of the model defined by the model parameters given the supplied samples.
#'
#' @export
#' @useDynLib VBGfit minLL_growthcurve

minLL_structure <- function(samples, ages=NA, sizes=NA, weights=NA, initpars, model = "BOTH", feedinggroups = NA, feedingbounds = NA, logscale = TRUE, shrink = TRUE){

# Check and process input -------------------------------------------------

  #### Check samples

  # check whether input is provided seperately or as dataframe
  if(is.vector(samples)){
    # check ages
    if(is.na(ages) || !is.vector(ages) || length(ages) != length(samples)){
      stop("ages is missing or not a vector with the same length as samples")
      ages = NA
    }
    #check sizes
    if(is.na(sizes) || !is.vector(sizes) || length(sizes) != length(samples)){
      stop("sizes is missing or not a vector with the same length as samples")
      sizes = NA
    }
    #check weights
    #set weights to 1 if not provided
    if(length(weights) == 1 && sum(is.na(weights))==1){
      weights <- 1
      warning("Equal weights assumed")
    }
    else if( !is.vector(weights) || length(weights) != length(samples) ){
      warning("weights not a vector with the same length as samples, all weights set to 1")
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
    stop("The dataset should contain more than 1 timepoint with complete data")
  }

  if(nages <= 1){
    stop( "The dataset should contain more than 1 ageclass with complete data")
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
    warning( "Feedinggroups is not one of 'AGE' or 'SIZE' and is ignored")
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
    stop( "NA values are not allowed in initpars")
  }
  else if(sum(initpars <= 0) && logpar==1){
    stop( paste("initpars may not contain values <= 0 if parameters are logtrasformed (logscale=TRUE)"))
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
    stop(
      msg = "Size of initpars is incorrect. Initpars should be a vector of length 5 or a vector as produced by 'makepoppars()'")
  }


  initpars_short <- initpars

  # Calculate negative log likelyhood -------------------------------------------------

  #### make output list
  LLout <- list(CONS = NA, VARY = NA)

  #### preform calculation with constant feedinglevel
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

    # do actual calculation
    consout <- .Call('minLL_growthcurve', initpars, globpars, feedingbounds, samples, as.integer(0))

    # save output
    LLout$CONS <- consout
  }

  if(model == "VARY" || model == "BOTH"){
    # create vector with all startpars
    if(npars == 5){
      initpars_long <- makepoppars(fRlm_avg = initpars[1], fRlm_var = initpars[2], len0_avg = initpars[3], len0_var = initpars[4], rB = initpars[5], time0_avg=NULL, time0_var=NULL, ntime = ntime, nages = nages,  model="VARY", ngroups = ngroups)
    }
    else initpars_long <- initpars

    # logtransform if needed
    if(logpar == 1) initpars_long2 <- log(as.numeric(initpars_long))
    else initpars_long2 <- as.numeric(initpars_long)

    # create vector with global variables
    globpars <- as.integer(c(length(initpars_long), ntime, mintime, nages, minage, nobs, shrinkpar, logpar, grouppar, ngroups))

    # do actual calculation
    varyout <- .Call('minLL_growthcurve', initpars, globpars, feedingbounds, samples, as.integer(1))

    # save output
    LLout$VARY <- varyout
  }
  return(LLout)
}



#' Internal! Calculate the negative log Likelyhood given a sampleset
#'
#' This function is for internal use only.
#' This function is a wrapper for a C function but does not check input arguments and therefore easely craches R.
#' Do use \code{\link{minLL_structure}} instead.
#'
#' This function calculates the negative log likelyhood of a population structure given a sampleset
#'
#' @param poppars A numeric vector with log transformation of the population parameters as calculated by \code{\link{makepoppars}}.
#' @param globpars A integer vector containing the following global parameters:
#' \itemize{
#' \item{The number of parameters}
#' \item{The number of timepoints}
#' \item{The lowest timepoint in the model}
#' \item{The number of ages}
#' \item{The lowest age in the model}
#' \item{Number of observations in sampleset (might be NA)}
#' }
#' @param samples A dataset of observations containing age-size measurements at various timepoints. The three columns of the dataset should correspond to the following information int the given order:
#' \itemize{
#' \item{A descrete (integer) value for the timepoint}
#' \item{A descrete (integer) value for the age of the individual}
#' \item{A numericvalue of the size of an individual}
#' }
#' @param modeltype An integer 0 for a model with constant asymptotic size or integer 1 for a model with a asymptotic size varying between timesteps.
#'
#'
#' @return The negative log likelyhood of a population structure given a samplest
#'
#' @noRd
#'
#' @useDynLib VBGfit minLL_growthcurve
minLL_structure_ <- function(poppars, globpars, feedingbounds, samples, modeltype){
  minLL <- .Call('minLL_growthcurve', poppars, globpars, feedingbounds, samples, modeltype)
  return(minLL)
}

# VBGfit package
# 15 march 2022
# latest version on https://bitbucket.org/JCroll/vbgfit/


#'  Calculate the population structure for a given parameterset
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
#' @param ntime number of timepoints in the model.
#' @param nages number of ages in the model.
#' @param mintime first timepoint in the model.
#' @param minage  Lowest age in the model.
#'
#' @param model A specific character expression specifying whether the asymptotic size is allowed to vary between timepoints.
#' \describe{
#' \item{\code{CONS}}{The asyptotic size does not vary between years}
#' \item{\code{VARY}}{The asyptotic size is allowed vary between years}
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
#' @param shrink Boolean to indicating whether the average size is allowed to shrink or should kept constant when the size exceeds the asymptotic size.

#' @return a list with the population structures predicted by the selected models. The population structures are saved as a matrix with the average and variance in asymptotic size in which the columns indicate the timepoints and the rows the different ages in the model.
#'
#' @export
#'
#' @useDynLib VBGfit predict_growthcurve
#'

predict_structure <- function(initpars, ntime, nages, mintime = 0, minage = 0, model = "BOTH", feedinggroups = NA, feedingbounds = NA, shrink = TRUE){

# Check and process input -------------------------------------------------

  #### Check input
  inputCheck <- ArgumentCheck::newArgCheck()

  # check modeltype
  if(model != "CONS" && model != "VARY" && model != "BOTH"){
    ArgumentCheck::addWarning(
      msg = "model is not one of 'CONS', 'VARY' or 'BOTH', and is reset to default: 'BOTH'",
      argcheck = inputCheck)
    model = "BOTH"
  }

  # check shrink
  if( shrink !=TRUE && shrink !=FALSE ){
    ArgumentCheck::addWarning(
      msg = "shrink should be TRUE or FALSE, set to TRUE",
      argcheck = inputCheck)
    shrink = TRUE
  }
  # convert shrink to integer
  if(!shrink) shrinkpar = 0
  else shrinkpar = 1

  # check feedinggroup
  if(is.na(feedinggroups)) grouppar = 0
  else if( feedinggroups == "AGE") grouppar = 1
  else if( feedinggroups == "SIZE") grouppar = 2
  else{
    ArgumentCheck::addWarning(
      msg = "Feedinggroups is not one of 'AGE' or 'SIZE' and is ignored",
      argcheck = inputCheck)
    grouppar = 0
  }

  # check feedingbounds
  if(grouppar > 0){
    feedingbounds <- unique(feedingbounds[!is.na(feedingbounds)])
    ngroups <- length(feedingbounds)+1
    if(ngroups == 1){
      ArgumentCheck::addWarning(
        msg = "No correct boundaries for feedinggroups found. Feedinggroups ignored",
        argcheck = inputCheck)
      grouppar = 0
      feedingbounds = 0
    }
    else{
      ArgumentCheck::addWarning(
        msg = paste(as.character(ngroups-1),"unique boundaries of feedinggroups found, continuing with",as.character(ngroups),"feedinggroups"),
        argcheck = inputCheck)
    }
  }
  else{
    ngroups = 1;
    feedingbounds = 0;
  }

  # check initpars
  npars <- length(initpars)

  if(anyNA(initpars)){
    ArgumentCheck::addError(
      msg = "NA values are not allowed in initpars",
      argcheck = inputCheck)
  }
  else if( (npars == 2*ngroups + 2*ntime + 2*(nages-1) + 1 ) ){
    ArgumentCheck::addWarning(
      msg = "Model set to 'CONS' based on length of initpars",
      argcheck = inputCheck)
    model = "CONS"
    modeltype <- as.integer(0)
  }
  else if(npars == 2*(ntime-1)*ngroups + 2*ntime + 2*(nages-1) + 1 ){
    ArgumentCheck::addWarning(
      msg = "Model set to 'VARY' based on length of initpars",
      argcheck = inputCheck)
    model = "VARY"
    modeltype <- as.integer(1)
  }
  else if(npars == 5){
    ArgumentCheck::addWarning(
      msg = "Model set to 'CONS' based on length of initpars",
      argcheck = inputCheck)
    model = "CONS"
    modeltype <- as.integer(0)
    initpars <- makepoppars(fRlm_avg = initpars[1], fRlm_var = initpars[2], len0_avg = initpars[3], len0_var = initpars[4], rB = initpars[5], time0_avg=NULL, time0_var=NULL, ntime = ntime, nages = nages,  model="CONS", ngroups = ngroups)
    npars = length(initpars)
  }
  else if( npars != 5 ){
    ArgumentCheck::addError(
      msg = "Size of initpars is incorrect. Initpars should be a vector of length 5 or a vector as produced by 'makepoppars()'",
      argcheck = inputCheck)
  }

  ArgumentCheck::finishArgCheck(inputCheck)

  # Calculate sturcture -------------------------------------------------

  globpars <- as.integer(c(npars, ntime, mintime, nages, minage, 0, shrinkpar, 0, grouppar, ngroups))

  initpars <- as.numeric(initpars)

  # run calculation

  structure <- .Call('predict_growthcurve', initpars, globpars, modeltype, feedingbounds)
  names(structure) <- c("avg", "var")
  structure$avg <- data.frame(structure[[1]])
  structure$var <- data.frame(structure[[2]])
  colnames(structure$avg) <- seq(from = globpars[3], to = globpars[3] + globpars[2] - 1, by = 1)
  colnames(structure$var) <- seq(from = globpars[3], to = globpars[3] + globpars[2] - 1, by = 1)
  rownames(structure$avg) <- seq(from = globpars[5], to = globpars[5] + globpars[4] - 1, by = 1)
  rownames(structure$var) <- seq(from = globpars[5], to = globpars[5] + globpars[4] - 1, by = 1)

  return(structure)


}



#' Internal! Calculate the population structure
#'
#' This function is for internal use only.
#' This function is a wrapper for a C function but does not check input arguments and therefore easely craches R.
#' Do use \code{\link{predict_structure}} instead.
#'
#' This function calculates a population structure based on a vector of population parameters.
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
#' @param modeltype An integer 0 for a model with constant asymptotic size or integer 1 for a model with a asymptotic size varying between timesteps.
#'
#' @return a list containing a matrix with the predicted average size and variance in size for every time and age combination in the model.
#'
#'@noRd
#'
#' @useDynLib VBGfit predict_growthcurve
#'
predict_structure_ <- function(poppars, globpars, modeltype, feedingbounds){
  structure <- .Call('predict_growthcurve', poppars, globpars, modeltype, feedingbounds)
  names(structure) <- c("avg", "var")
  structure$avg <- data.frame(structure[[1]])
  structure$var <- data.frame(structure[[2]])
  colnames(structure$avg) <- seq(from = globpars[3], to = globpars[3] + globpars[2] - 1, by = 1)
  colnames(structure$var) <- seq(from = globpars[3], to = globpars[3] + globpars[2] - 1, by = 1)
  rownames(structure$avg) <- seq(from = globpars[5], to = globpars[5] + globpars[4] - 1, by = 1)
  rownames(structure$var) <- seq(from = globpars[5], to = globpars[5] + globpars[4] - 1, by = 1)

  return(structure)
}




# VBGfit package
# 15 march 2022
# latest version on https://bitbucket.org/JCroll/vbgfit/

#' Make parameter vector defining the population stucture
#'
#' @param fRlm_avg A single value for the mean asymptotic size,
#' or a vector with \code{ngroups} elements with a mean asymptotic size per feedinggroup
#' or if the model is \code{VARY} a vector with size \code{ntime-1} with the mean asymptotic size between every timestep
#' or if the model is \code{VARY} a vector with size \code{ngroups*(ntime-1)} with the mean asymptotic size per year for every feedinggroup (clustered per feedinggroup)
#'
#' @param fRlm_var A single value for the variance in asymptotic size,
#' or a vector with \code{ngroups} elements with the variance in asymptotic size per feedinggroup
#' or if the model is \code{VARY} a vector with size \code{ntime-1} with the variance in asymptotic size between every timestep
#' or if the model is \code{VARY} a vector with size \code{ngroups*(ntime-1)} with the variance in asymptotic per year for every feedinggroup (clustered per feedinggroup)
#'
#' @param len0_avg Vector of length \code{ntime} containing the average size at the youngest age in the model for every timepoint. If a single value is provided this is repeated for all timepoints.
#' @param len0_var Vector of length \code{ntime} containing the variance in size at the youngest age in the model for every timepoint. If a single value is provided this is repeaded for all timepoints.
#'
#' @param rB single value for the Von Bertalanffy growth rate of the growth curve.
#'
#' @param time0_avg (optional) vector of length \code{nage-1} containing the mean size at age at the first timepoint except for the mean size at the youngest age. If not given, this is calculated with the parameters for the first timepoint.
#' @param time0_var (optional) vector of length \code{nage-1} containing the variance in size at age at the first timepoint except for the variance in size at the youngest age.  If not given, this is calculated with the parameters for the first timepoint.
#'
#' @param ntime number of timepoints in the model
#' @param nages number of ages in the model
#'
#' @param model A specific character expression specifying whether the asymptotic size is allowed to vary between timepoints.
#' \describe{
#' \item{\code{CONS}}{The asymptotic size does not vary between years}
#' \item{\code{VARY}}{The asymptotic size is allowed vary between years}
#' \item{\code{BOTH}}{Fit both models}
#' }
#'
#' @param ngroups Number of feeding groups in the model
#'
#' @return A vector with all model variables that are optimized by \link{fit_structure}.
#' If model is \code{CONS} the vector contains the following elements in order:
#' \itemize{
#' \item{average asymptotic size (1 element if model is \code{CONS}, \code{ngroups*(ntime-1)} elements if model is \code{VARY})}
#' \item{variance in asymptotic size (1 element if model is \code{CONS}, \code{ngroups*(ntime-1)} elements if model is \code{VARY})}
#' \item{average length at the youngest age in the model (\code{ntime} elements)}
#' \item{variance in length at the youngest age in the model (\code{ntime} elements)}
#' \item{average length at age at the first timepoint in the model except for the youngest age (\code{nages} elements)}
#' \item{variance in length at age at the first timepoint in the model except for the youngest age (\code{nages} elements)}
#' \item{The growth rate of the Von Bertalanffy growth curve (1 element)}
#' }
#'
#' @export
#'
makepoppars <- function(fRlm_avg, fRlm_var, len0_avg, len0_var, rB, time0_avg=NULL, time0_var=NULL, ntime = NULL, nages = NULL,  model="BOTH", ngroups = 1){

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

  # check ntime
  if(is.null(ntime) && length(fRlm_avg) == 0 &&  length(fRlm_var) == 0 &&  length(len0_avg) == 1 &&  length(len0_var) == 1 ){
    ArgumentCheck::addError(
      msg = "ntime is missing and number of timepoints could not be derived. ",
      argcheck = inputCheck)
  }
  else if(is.null(ntime)){
    ntimeoptions <- c(length(fRlm_avg)/ngroups+1, length(fRlm_var)/ngroups+1, length(fRlm_avg)+1, length(fRlm_var)+1, length(len0_avg), length(len0_var))
    if(sum((ntimeoptions>=2 && ntimeoptions%%1==0)) == 0){
      ArgumentCheck::addError(
        msg = "ntime is missing and number of timepoints could not be derived. ",
        argcheck = inputCheck)
    }
    else{
    ntime <- min( ntimeoptions[which(ntimeoptions>=2 && ntimeoptions%%1==0)], na.rm = TRUE )
    ArgumentCheck::addWarning(
      msg = paste("ntime is set to",as.character(ntime),"based on the shortest length of the vectors with fRlm_avg, fRlm_var, time0_avg and time0_var"),
      argcheck = inputCheck)
    }
  }

  #check nages
  if(is.null(nages) && is.null(time0_avg) &&  is.null(time0_var)  ){
    ArgumentCheck::addError(
      msg = "nages is missing and number of ages could not be derived. ",
      argcheck = inputCheck)
    }
  else if(is.null(nages)){
    nagesoptions <- c(length(time0_avg)+1, length(time0_var)+1)
    nages <- min( nagesoptions[which(nagesoptions>=2)] )
    ArgumentCheck::addWarning(
      msg = paste("nages is set to",as.character(nages),"based on the shortest length of the vectors with time0_avg and time0_var"),
      argcheck = inputCheck)
  }

  #check len0_avg
  if(length(len0_avg) > ntime){
    ArgumentCheck::addWarning(
      msg = paste("too many elements in len0_avg, only the first", as.character(ntime), "elements used"),
      argcheck = inputCheck)
    len0_avg <- len0_avg[1:(ntime)]
  }

  #check len0_var
  if(length(len0_var) > ntime){
    ArgumentCheck::addWarning(
      msg = paste("too many elements in len0_var, only the first", as.character(ntime), "elements used"),
      argcheck = inputCheck
    )
    len0_var <- len0_var[1:(ntime)]
  }

  # check time0_avg
  if(length(time0_avg) > ntime-1){
    ArgumentCheck::addWarning(
      msg = paste("too many elements in time0_avg, only the first", as.character(nages-1), "elements used"),
      argcheck = inputCheck
    )
    time0_avg <- time0_avg[1:(nages-1)]
  }

  # check time0_var
  if(length(time0_var) > ntime-1){
    ArgumentCheck::addWarning(
      msg = paste("too many elements in time0_var, only the first", as.character(nages-1), "elements used"),
      argcheck = inputCheck)
    time0_var <- time0_var[1:(nages-1)]
  }

  ArgumentCheck::finishArgCheck(inputCheck)


# Calculate further arguments -------------------------------------------------

  processCheck <- ArgumentCheck::newArgCheck()

  # Calculate estimate of size at age at time 0
  if(is.null(time0_avg)){
    time0_avg <- len0_avg[1] * exp(-rB*seq(1, nages-1, 1 )) + fRlm_avg[1]*(1-exp(-rB*seq(1, nages-1, 1 )))
    if(ngroups > 1){
      ArgumentCheck::addWarning(
        msg = "not accounted for feeding groups when calculating the size at age at the first timepoint",
        argcheck = processCheck)
    }
  }

  # Calculate variance of size at age at time 0
  if(is.null(time0_var)){
    time0_var <- len0_var[1] * exp(-2*rB*seq(1, nages-1, 1 )) + fRlm_var[1]*(1-exp(-rB*seq(1, nages-1, 1 )))^2
    if(ngroups > 1){
      ArgumentCheck::addWarning(
        msg = "not accounted for feeding groups when calculating the variance in size at age at the first timepoint",
        argcheck = processCheck)
    }
  }

  # expand length at age 0
  if(length(len0_avg)==1) len0_avg <- rep(len0_avg, ntime)

  if(length(len0_var)==1) len0_var <- rep(len0_var, ntime)

  # further calculations
  if(model == "CONS" || model == "BOTH"){
    if(length(fRlm_avg) > ngroups){
      ArgumentCheck::addWarning(
        msg = paste("too many elements in fRlm_avg for constant feedinglevel, only the first", as.character(ngroups), "elements used"),
        argcheck = processCheck)
      fRlm_avg_CONS <- fRlm_avg[1:ngroups]
    }
    else if(length(fRlm_avg) == ngroups){
      fRlm_avg_CONS <- fRlm_avg
    }
    else if(length(fRlm_avg) < ngroups && length(fRlm_avg) > 1){
      ArgumentCheck::addWarning(
        msg = "too many elements in fRlm_avg for constant feedinglevel, only the first elements is used",
        argcheck = processCheck)
      fRlm_avg_CONS <- rep(fRlm_avg[1], ngroups)
    }
    else{
      fRlm_avg_CONS <- rep(fRlm_avg[1], ngroups)
    }

    if(length(fRlm_var) > ngroups){
      ArgumentCheck::addWarning(
        msg = paste("too many elements in fRlm_var for constant feedinglevel, only the first", as.character(ngroups), "elements used"),
        argcheck = processCheck)
      fRlm_var_CONS <- fRlm_var[1:ngroups]
    }
    else if(length(fRlm_var) == ngroups) fRlm_var_CONS <- fRlm_var
    else if(length(fRlm_var) < ngroups && length(fRlm_var) > 1){
      ArgumentCheck::addWarning(
        msg = "too many elements in fRlm_var for constant feedinglevel, only the first elements is used",
        argcheck = processCheck)
      fRlm_var_CONS <- rep(fRlm_var[1], ngroups)
    }
    else fRlm_var_CONS <- rep(fRlm_var[1], ngroups)

    poppars_CONS = c(fRlm_avg_CONS, fRlm_var_CONS, len0_avg, len0_var, time0_avg, time0_var, rB)
  }

  if(model == "VARY" || model == "BOTH"){
    if(length(fRlm_avg) > ngroups*(ntime - 1)){
      ArgumentCheck::addWarning(
        msg = paste("too many elements in fRlm_avg for varying feedinglevel, only the first", as.character(ngroups*(ntime -1)), "elements used"),
        argcheck = processCheck)
      fRlm_avg_VARY <- fRlm_avg[1:(ngroups*(ntime - 1))]
    }
    else if(length(fRlm_avg) == ngroups*(ntime - 1)) fRlm_avg_VARY <- fRlm_avg
    else if(length(fRlm_avg) == ntime - 1 ) fRlm_avg_VARY <- rep(fRlm_avg, ngroups)
    else if(length(fRlm_avg) == ngroups ){
      fRlm_avg_VARY <- c()
      for(i in 1:(ngroups)) fRlm_avg_VARY <- c(fRlm_avg_VARY, rep(fRlm_avg[i], (ntime-1)))
    }
    else if(length(fRlm_avg) > 1){
      ArgumentCheck::addWarning(
        msg = "Number of elements in fRlm_avg does not match the number of needed elements for varying feeding level, only the first element is used.",
        argcheck = processCheck)
      fRlm_avg_VARY <- rep(fRlm_avg[1],ngroups*(ntime-1))
    }
    else fRlm_avg_VARY <- rep(fRlm_avg[1],ngroups*(ntime-1))

    if(length(fRlm_var) > ngroups*(ntime - 1)){
      ArgumentCheck::addWarning(
        msg = paste("too many elements in fRlm_var for varying feedinglevel, only the first", as.character(ngroups*(ntime -1)), "elements used"),
        argcheck = processCheck)
      fRlm_var_VARY <- fRlm_var[1:(ngroups*(ntime - 1))]
    }
    else if(length(fRlm_var) == ngroups*(ntime - 1)) fRlm_var_VARY <- fRlm_var
    else if(length(fRlm_var) == ntime - 1 ) fRlm_var_VARY <- rep(fRlm_var, ngroups)
    else if(length(fRlm_var) == ngroups ){
      fRlm_var_VARY <- c()
      for(i in 1:(ngroups)) fRlm_var_VARY <- c(fRlm_var_VARY, rep(fRlm_var[i], (ntime-1)))
    }
    else if(length(fRlm_var) > 1){
      ArgumentCheck::addWarning(
        msg = "Number of elements in fRlm_var does not match the number of needed elements for varying feeding level, only the first element is used.",
        argcheck = processCheck)
      fRlm_var_VARY <- rep(fRlm_var[1],ngroups*(ntime-1))
    }
    else fRlm_var_VARY <- rep(fRlm_var[1],ngroups*(ntime-1))

    poppars_VARY = c(fRlm_avg_VARY, fRlm_var_VARY, len0_avg, len0_var, time0_avg, time0_var, rB)
  }

  ArgumentCheck::finishArgCheck(processCheck)

  # return values

  if(model == "CONS") return(poppars_CONS)
  else if(model == "VARY") return(poppars_VARY)
  else return( list(CONS= poppars_CONS, VARY = poppars_VARY) )

}

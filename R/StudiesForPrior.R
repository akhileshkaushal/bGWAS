###### Functions for study selection for the prior ######



#' Get list of studies that can be used to build the prior
#'
#' Get list of studies with available GWAS summary statistics
#' @param IDs numeric,
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @param verbose boolean, default = FALSE
#' @return Data.Frame of details for all available studies
#' @export

list_priorGWASs <- function(IDs=NULL, Z_matrices = "~/ZMatrices/", verbose=F) {
  if(!is.null(IDs) && !is.numeric(IDs)) stop("ID : should be numeric")
  #Studies = data.table::fread(system.file("Data/AvailableStudies.tsv", package="bGWAS"), showProgress = FALSE)
  Studies = data.table::fread(paste0(Z_matrices, "/AvailableStudies.tsv"), showProgress = FALSE, data.table = F)
  if(!is.null(IDs)){
    Studies=Studies[IDs, ]
  }
  return(Studies)
}


# #' Get prior GWASs names of studies that can be used to build the prior
# #'
# #' Get names of studies with available GWAS summary statistics
# #' @param IDs The IDs of the prior GWASs to be listed, by default=NULL, returns all prior GWASs names (numeric)
# #' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
# #'        (character)
# #' @param verbose boolean, default = FALSE
# #' @return Names of the prior GWASs

list_names <- function(IDs=NULL, Z_matrices = "~/ZMatrices/", verbose=F) {
  if(is.null(IDs)) {
    Names = data.table::fread(paste0(Z_matrices, "/AvailableStudies.tsv"), select = "Name", showProgress = FALSE, data.table=F)
  } else {
    # check that the IDs exists
    Studies = list_priorGWASs(Z_matrices=Z_matrices)
    if(!all(IDs %in% Studies$ID)) print("Please check the IDs, some of them do not match")
    Names = subset(Studies, ID %in% IDs)
  }
  return(Names$Name)
}


# #' Get traits from studies that can be used to build the prior
# #'
# #' Get list of traits with available GWAS summary statistics
# #' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
# #'        (character)
# #' @param verbose boolean, default = FALSE
# #' @return List of traits

list_traits <- function(Z_matrices = "~/ZMatrices/", verbose=F) {
  Traits = data.table::fread(paste0(Z_matrices, "/AvailableStudies.tsv"), select = "Trait", showProgress = FALSE, data.table = F)$Trait
  return(unique(Traits))
}


# #' Get consortia from studies that can be used to build the prior
# #'
# #' Get list of consortia with available GWAS summary statistics
# #' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
# #'        (character)
# #' @param verbose boolean, default = FALSE
# #' @return List of consortia


list_consortia <- function(Z_matrices = "~/ZMatrices/", verbose=F) {
  Consortia = data.table::fread(paste0(Z_matrices, "/AvailableStudies.tsv"), select = "Consortium", showProgress = FALSE, data.table=F)$Consortium
  return(unique(Consortia))
}


#' select_priorGWASs
#'
#' Allow the quick selection of a subset of studies for prior based on 3 criteria. First, include all the prior GWASs
#' specified (if all including parameters are NULL, include all studies), and then remove all the prior GWASs specified
#' (if all excluding parameters are NULL, keep all studies included at the step before)
#' @param include_names list of prior GWASs names
#' @param include_traits list of trait
#' @param include_consortium vector, list of consortium ### TO BE DONE?
#' @param exclude_names list of file names
#' @param exclude_traits list of trait
#' @param exclude_consortium vector, list of consortium ### TO BE DONE?
#' @param Z_matrices The path to the folder containing Z-Matrices, \code{default="~/ZMatrices/"}
#'        (character)
#' @param verbose boolean, default = FALSE
#' @return IDs of studies that meet the criteria
#' @examples
#'   AllStudies = list_priorGWASs()
#'   list_names()
#'   MyStudies = select_priorGWASs(include_traits=c("Heart Rate", "Body mass index", "Smoking"))
#'   #AllStudies[AllStudies$ID %in% MyStudies, ]
#' @export


select_priorGWASs <- function(include_names=NULL, include_traits=NULL, includeConsortia=NULL,
                          exclude_names=NULL, exclude_traits=NULL, excludeConsortia=NULL, Z_matrices = "~/ZMatrices/", verbose=F) {
  # Check parameters
  if(is.null(c(include_names,include_traits,includeConsortia, exclude_names,exclude_traits,excludeConsortia))) stop("You did not specify any criteria for the selection.")
  Studies = list_priorGWASs(Z_matrices = Z_matrices)

  if(!all(include_names %in% list_names(Z_matrices = Z_matrices))) stop("Some names specified in include_names are not correct")
  if(!all(exclude_names %in% list_names(Z_matrices = Z_matrices))) stop("Some names specified in exclude_names are not correct")
  if(!all(include_traits %in% list_traits(Z_matrices))) stop("Some trait specified in include_traits are not correct")
  if(!all(exclude_traits %in% list_traits(Z_matrices))) stop("Some traits specified in exclude_traits are not correct")
  if(!all(includeConsortia %in% list_consortia(Z_matrices))) stop("Some consortia specified in include_consortia are not correct")
  if(!all(excludeConsortia %in% list_consortia(Z_matrices))) stop("Some consortia specified in exclude_consortia are not correct")

### TO BE DONE
  includeConsortia = NULL
  excludeConsortia = NULL

  n = nrow(Studies)
  Crit = c()

  # if inclusion criteria
  if(!is.null(includeConsortia) | !is.null(include_traits) | !is.null(include_names)){
    AllStudies = Studies
    Studies = Studies[0,] # create empty data table
    n = nrow(Studies)

    ## A : do all inclusion
    # 1st : selection based on Study Name
    if(!is.null(include_names)){
      Studies <- rbind(Studies, AllStudies[AllStudies$Name %in% include_names,])
      Crit = c(Crit, c("inclusion of names"))
      nS = nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nS, " prior GWAS(s) have been included when selection using include_names"))
    }
    # 2nd : selection based on Trait
    if(!is.null(include_traits)){
      Studies <- rbind(Studies, AllStudies[AllStudies$Trait %in% include_traits,])
      Crit = c(Crit, c("inclusions of traits"))
      nT =  nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nT, " prior GWAS(s) have been included when selection using include_traits"))
    }
    # 3rd : selection based on Consortium
    if(!is.null(includeConsortia)){
      Studies <- rbind(Studies, AllStudies[AllStudies$Consortium %in% includeConsortia,])
      Crit = c(Crit, c("inclusion of consortia"))
      nC = nrow(Studies) - n
      n = nrow(Studies)
      if(verbose) print(paste0(nC, " prior GWAS(s) have been included when selection using includeConsortia"))
    }
  }

  ## B : do all exclusion
  # 1st : selection based on Study Name
  if(!is.null(exclude_names)){
    Studies = subset(Studies, !(Name %in% exclude_names))
    Crit = c(Crit, c("exclusion of names"))
    nS = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nS, " prior GWAS(s) have been removed when selection using exclude_names"))
  }
  # 2nd : selection based on Trait
  if(!is.null(exclude_traits)){
    Studies = subset(Studies, !(Trait %in% exclude_traits))
    Crit = c(Crit, c("exclusion of traits"))
    nT = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nT, " prior GWAS(s) have been removed when selection using exclude_traits"))
  }
  # 3rd : selection based on Consortium
  if(!is.null(excludeConsortia)){
    Studies = subset(Studies, !(Consortium %in% excludeConsortia))
    Crit = c(Crit, c("exclusion of consortia"))
    nC = n - nrow(Studies)
    n = nrow(Studies)
    if(verbose) print(paste0(nC, " prior GWAS(s) have been removed when selection using excludeConsortia"))
  }
  if(verbose) print(paste0(n, " studies left after selection on ", paste(Crit, collapse=", ")))
  return(Studies$ID)
}



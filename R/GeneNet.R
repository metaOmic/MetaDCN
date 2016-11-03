##' Generate correction and adjacency matrice for data and permutation
##'
##' This function will generate correction and adjacency matrix for data and application
##'
##' @title GeneNet
##' @param data a list of matrice (studies) with row names as genes and column 
##' names as samples.
##' @param labels a list of vectors to specify the case/control labels for 
##' each subject corresponidng to the columns in each data matrix
##' @param caseName a character string denoting case label
##' @param controlName a character string denoting control label
##' @param meanFilter a number to specify qunatile cutoff for mean
##' @param SDFilter a number to specify qunatile cutoff for SD
##' @param edgeCutoff a number to specify qunatile cutoff for correlation to 
##' to define an edge
##' @param permutationTimes a number to specify how many permutations to be 
##' used
##' @param CPUNumbers a number to specify how many CPUs are used for parallel,
##' must be less than permutationTimes
##' @param folder a character string folder name to store results
##' @param pathwayDatabase a list with each element as a vector of 
##' pathway genes 
##' @param silent TRUE/FALSE to specify if suppress screen output
##' @return a list and several RData files saved in the current working directory
##' \item{GeneNet }{A list containing important parameters to pass to SearchBM and MetaDCN functions}
##' \item{AdjacencyMatrice.RData }{A list of adjacency matrice for case and control in each study in the order of case studies and control studies}
##' \item{CorrelationMatrice.RData }{A list of correlation matrice for case and control each study in the order of case studies and control studies}
##' \item{CorrelationMatrice.RData }{A list of correlation matrice for case and control each study in the order of case studies and control studies}
##' \item{AdjacencyMatricePermutationP.RData }{A list of correlation matrice for case and control in each study in the order of case studies and control studies for permutation P}
##' @author Li Zhu
##' @import snow
##' @import snowfall
##' @import igraph
##' @import genefilter
##' @export
##' @examples 
##' data(pathwayDatabase)

GeneNet <- function(data, labels, caseName, controlName, meanFilter=0.2, 
  SDFilter=0.2, edgeCutoff=0.004, permutationTimes=10, CPUNumbers=1, 
  folder="MetaDCN", pathwayDatabase, silent=FALSE){
  
  GeneNetRes <- list()
  GeneNetRes$caseName <- caseName
  GeneNetRes$controlName <- controlName
  GeneNetRes$permutationTimes <- permutationTimes
  GeneNetRes$folder <- folder
  GeneNetRes$pathwayDatabase <- pathwayDatabase
  GeneNetRes$CPUNumbers <- CPUNumbers

  if(CPUNumbers > 1){
    if(CPUNumbers > permutationTimes){
      stop("CPUNumbers should be smaller than or equal to permutationTimes.")
    }else if(CPUNumbers < 2){
      stop("CPUNumbers should be greater than or equal to 2.")
    }
  }
  
  data2 <- sapply(1:length(data), function(x) data[[x]][, 
    which(labels[[x]] %in% c(caseName,controlName))])
  data <- data2

  labels2 <- sapply(1:length(data), function(x) labels[[x]][ 
    which(labels[[x]] %in% c(caseName,controlName))])
  labels <- labels2
  
  caseIndex <- sapply(1:length(data), function(x) 
    which(labels[[x]] == caseName))
  controlIndex <- sapply(1:length(data), function(x) 
    which(labels[[x]] == controlName))

  caseStudyIndex <- seq(1:length(data))
  controlStudyIndex <- seq(length(data)+1, 2*length(data))

  ### generate network (permute_index is only needed for permutation parallel)
  NetworkGeneration(data, caseIndex, controlIndex, caseStudyIndex, 
    controlStudyIndex, meanFilter, SDFilter, edgeCutoff, choose="real",  
    permuteIndex=NA, folder, silent=FALSE) 

  ### generate network and search modules for permutations
  if (CPUNumbers > 1){
    paraRep <- round(permutationTimes/CPUNumbers)
    CPUNumbers2 <- permutationTimes%%CPUNumbers

    ## generate network for permutations
    for (nr in 1:paraRep){
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers) 
      snowfall::sfExport("data", "caseIndex", "controlIndex", "caseStudyIndex",
        "controlStudyIndex", "meanFilter", "SDFilter", "edgeCutoff",
        "permutationTimes", "paraRep", "CPUNumbers2", "CPUNumbers", "nr", 
        "folder")
      snowfall::sfClusterApply(seq(1, CPUNumbers), 
        function(x) NetworkGeneration(data, caseIndex, controlIndex, 
          caseStudyIndex, controlStudyIndex, meanFilter, SDFilter, 
          edgeCutoff, choose="permute", permuteIndex=((nr-1)*CPUNumbers+x), 
        folder, silent=FALSE)) 
      snowfall::sfStop()
    }
    if(CPUNumbers2>0){
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("data", "caseIndex", "controlIndex", 
        "caseStudyIndex", "controlStudyIndex", "meanFilter", "SDFilter", 
        "edgeCutoff", "permutationTimes", "paraRep", "CPUNumbers2", 
        "CPUNumbers", "nr","folder")
      snowfall::sfClusterApply(seq(1,CPUNumbers2), 
        function(x) NetworkGeneration(data, caseIndex, controlIndex, 
              caseStudyIndex, controlStudyIndex, choose="permute", meanFilter, 
              SDFilter, edgeCutoff, permuteIndex=(paraRep*CPUNumbers+x), 
              folder, silent=FALSE)) 
      snowfall::sfStop()
    }
    
  }else {   ### without parallel
    ### generate permutation network
    sapply(seq(1,permutationTimes), function(x) 
      NetworkGeneration(data, caseIndex, controlIndex, 
                  caseStudyIndex, controlStudyIndex, meanFilter, 
                  SDFilter, edgeCutoff, choose="permute", 
                  permuteIndex=x, folder, silent=FALSE))
  }
  return(GeneNetRes)
}


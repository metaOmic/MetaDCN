##' Generate correction and adjacency matrice for data and permutation
##'
##' This function will generate correction and adjacency matrix for data and 
##' permutations.
##'
##' @title GeneNet
##' @param data a list of matrice (studies) with rows as features and columns 
##' as samples.
##' @param labels a list of vectors to specify the case/control labels for 
##' each subject corresponidng to the columns in each data matrix.
##' @param caseName a character string denoting case label.
##' @param controlName a character string denoting control label.
##' @param meanFilter a number between 0 and 1. Features with mean below this
##' cutoff will be filtered out.
##' @param SDFilter a number between 0 and 1. Features with standard deviation
##' (SD) will be filtered out.
##' @param edgeCutoff a number between 0 and 1, denoting the proportion of 
##' edges to kept.
##' @param permutationTimes a number to specify how many permutations to be 
##' used. Large number of permutations will be very time-consuming.
##' @param CPUNumbers a number to specify how many CPUs are used for parallel,
##' if multicores exists. Must be less than or equal to the permutationTimes.
##' @param folder folder path to store results.
##' @param pathwayDatabase a list with each element as a vector of 
##' pathway genes, and names as pathway names.
##' @param silent TRUE/FALSE to specify if suppress screen output.
##' @return GeneNet returns a list of information which will be used for 
##' SearchBM and MetaDCN function, and several RData files stored in folder 
##' path. 
##' @return List of basic informations for SearchBM and MetaDCN input:
##' \item{caseName }{case names from arguments}
##' \item{controlName }{control names from arguments}
##' \item{permutationTimes }{permutation times from arguments}
##' \item{folder }{folder path from arguments}
##' \item{pathwayDatabase }{a list of pathways from arguments}
##' \item{CPUNumbers }{CPU numbers from arguments}
##' @return AdjacencyMatrice.RData is a list of adjacency matrice for case and 
##' control in each study in the order of case studies and control studies.
##' @return CorrelationMatrice.RData is a list of correlation matrice for case 
##' and control in each study in the order of case studies and control studies.
##' @return AdjacencyMatricePermutationP.RData is a list of correlation 
##' matrice for case and control in each study in the order of case studies 
##' and control studies for permutation P.
##' @author Li Zhu (liz86@pitt.edu)
##' @import snow
##' @import snowfall
##' @import igraph
##' @import genefilter
##' @export
##' @examples 
##' data(pathwayDatabase)
##' data(example)
##' GeneNetRes <- GeneNet(data_2l, labels_2l, caseName="inv(16)", controlName="t(8;21)", meanFilter=0.8, SDFilter=0.8, edgeCutoff=0.1, permutationTimes=4, CPUNumbers=1, folder="MetaDCN", pathwayDatabase, silent=FALSE)
##' SearchBMRes <- SearchBM(GeneNetRes, MCSteps=500, jaccardCutoff=0.8, repeatTimes=3, outputFigure=TRUE, silent=FALSE)
##' MetaDNCRes <- MetaDCN(GeneNetRes, SearchBMRes, FDRCutoff=0.3, w1=NULL, silent=FALSE)

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
      NoOutput <- snowfall::sfClusterApply(seq(1, CPUNumbers), 
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
      NoOutput <- snowfall::sfClusterApply(seq(1,CPUNumbers2), 
        function(x) NetworkGeneration(data, caseIndex, controlIndex, 
              caseStudyIndex, controlStudyIndex, choose="permute", meanFilter, 
              SDFilter, edgeCutoff, permuteIndex=(paraRep*CPUNumbers+x), 
              folder, silent=FALSE)) 
      snowfall::sfStop()
    }
    
  }else {   ### without parallel
    ### generate permutation network
    NoOutput <- sapply(seq(1,permutationTimes), function(x) 
      NetworkGeneration(data, caseIndex, controlIndex, 
                  caseStudyIndex, controlStudyIndex, meanFilter, 
                  SDFilter, edgeCutoff, choose="permute", 
                  permuteIndex=x, folder, silent=FALSE))
  }
  return(GeneNetRes)
}


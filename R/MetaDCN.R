##' Assemble basic modules into supermodules
##'
##' This function will assemble basic modules detected from SearchBM into 
##' supermodules.
##' 
##' @title MetaDCN
##' @param GeneNetRes A list from GeneNet function 
##' @param SearchBMRes A list from SearchBM function 
##' @param FDRCutoff A number to specify FDR cutoff for basic modules
##' @param w1 A number chosen from (100, 200, ..., 700) to specify the weight1 
##' used in objective function (optional). If not specified, w1 from SearchBM 
##' function will be used (recommended).
##' @param silent TRUE/FALSE to specify if suppress screen output
##' @return MetaDNC will return a list and files for Cytoscape.
##' @return list containing:
##' \item{w1}{w1 used}  
##' \item{BMInCase }{Summary of basic modules in case controling FDR}
##' \item{BMInConrtol }{Summary of basic modules in control controling FDR}
##' \item{Supermodule }{Summary of supermodules}
##' @return CytoscapeFiles is a folder containing files for Cytoscape.
##' @author Li Zhu (liz86@pitt.edu)
##' @import snow
##' @import snowfall
##' @import igraph
##' @import genefilter
##' @export

MetaDCN <- function(GeneNetRes, SearchBMRes, FDRCutoff, w1=NULL, silent=FALSE){
  
  caseName <- GeneNetRes$caseName
  controlName <- GeneNetRes$controlName
  permutationTimes <- GeneNetRes$permutationTimes
  folder <- GeneNetRes$folder
  pathwayDatabase <- GeneNetRes$pathwayDatabase

  MetaDCNRes <- list()

  if(is.null(w1)){
    w1 <- SearchBMRes$w1
  }
  MetaDCNRes$w1 <- w1

  BMInCase <- read.csv(paste(folder, 
    "/basic_modules_summary_forward_weight_", MetaDCNRes$w1, ".csv", sep=""), 
    header=TRUE)
  MetaDCNRes$BMInCase <- BMInCase[which(BMInCase[,"FDR"] < 
    FDRCutoff),]
  rownames(MetaDCNRes$BMInCase)<-NULL

  BMInControl <- read.csv(paste(folder, 
    "/basic_modules_summary_backward_weight_", MetaDCNRes$w1, ".csv", sep=""),
     header=TRUE)
  MetaDCNRes$BMInControl <- BMInControl[which(BMInControl[,"FDR"] <
    FDRCutoff),]
  rownames(MetaDCNRes$BMInControl) <- NULL

  forwardNum <- nrow(MetaDCNRes$BMInCase)
  backwardNum <- nrow(MetaDCNRes$BMInControl)
  
  if(forwardNum == 0 & backwardNum ==0 ){
    stop(paste("No basic module has FDR < ", FDRCutoff, sep=""))
  }

  if(silent == FALSE){
    cat(paste(forwardNum, "modules are generated in forward direction\n"))
    cat(paste(backwardNum, "modules are generated in backward direction\n"))
  }

  ### use the parameters to do module assembly
  MetaDCNRes$supermodule <- ModuleAssembly(MetaDCNRes$w1, FDRCutoff,
   caseName, controlName, pathwayDatabase, permutationTimes, folder)

  return(MetaDCNRes)
}


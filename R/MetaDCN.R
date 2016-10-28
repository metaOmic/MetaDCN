##' Assemble basic modules into supermodules
##'
##' This function will assemble basic modules found from SearchBM into 
##' supermodules
##' 
##' @title MetaDCN
##' @param GeneNetRes A list from GeneNet function 
##' @param FDRCutoff A number to specify FDR cutoff for basic modules
##' @param w1 A number chosen from (100, 200, ..., 700) to specify the weight1 
##' used in objective function (optional). If not specified, the w1 with the 
##' most modules will be chosen (recommended).
##' @param silent TRUE/FALSE to specify if suppress screen output
##' @return a list and one zip file for Cytoscape
##' \item{w1}{w1 used}  
##' \item{ModuleInCase }{Summary of basic modules in case}
##' \item{ModuleInConrtol }{Summary of basic modules in control}
##' \item{Supermodule }{Summary of supermodules}
##' \item{module_assembly_edge_node_list.zip}{Files for Cytoscape}
##' @author Li Zhu
##' @import snow
##' @import snowfall
##' @import igraph
##' @import genefilter
##' @export
##' @examples 
##' data(pathwayDatabase)

MetaDCN <- function(GeneNetRes, FDRCutoff, w1=NULL, silent=FALSE){
  
  caseName <- GeneNetRes$caseName
  controlName <- GeneNetRes$controlName
  permutationTimes <- GeneNetRes$permutationTimes
  outputPrefix <- GeneNetRes$outputPrefix
  pathwayDatabase <- GeneNetRes$pathwayDatabase

  MetaDCNRes <- list()

  ### selecte parameters
  forwardList <- read.csv(
    paste(outputPrefix, "_threshold_forward.csv", sep=""))[, 2:5]
  backwardList <- read.csv(
    paste(outputPrefix, "_threshold_backward.csv", sep=""))[, 2:5]

  if(is.null(w1)){
    indexMax <- which.max(rowSums(forwardList+backwardList))
    weightList <- seq(from=100, to=700, by=100)
    w1 <- weightList[indexMax]
    if(silent == FALSE){
      cat(paste("w1 parameter chosen is:", w1, "\n"))
    }
  }else{
    if(!(w1 %in% seq(100,700,by=100))){
      stop("Please choose w1 from (100, 200, ..., 700)")
    }
  }

  MetaDCNRes$w1 <- w1

  weightTempIndex <- which(colnames(forwardList) == paste("FDR_", 
    round(FDRCutoff*100), sep=""))
  forwardNum <- forwardList[w1/100, weightTempIndex]
  backwardNum <- backwardList[w1/100, weightTempIndex]
  
  if(forwardNum == 0 & backwardNum ==0 ){
    stop(paste("No basic module has FDR < ", FDRCutoff, sep=""))
  }

  if(silent == FALSE){
    cat(paste(forwardNum, "modules are generated in forward direction\n"))
    cat(paste(backwardNum, "modules are generated in backward direction\n"))
  }

  ModuleInCase <- read.csv(paste(outputPrefix, 
    "_basic_modules_summary_forward_weight_", MetaDCNRes$w1, ".csv", sep=""), 
    header=TRUE)
  MetaDCNRes$ModuleInCase <- ModuleInCase[which(ModuleInCase[,"FDR"] < 
    FDRCutoff),]
  rownames(MetaDCNRes$ModuleInCase)<-NULL

  ModuleInControl <- read.csv(paste(outputPrefix, 
    "_basic_modules_summary_backward_weight_", MetaDCNRes$w1, ".csv", sep=""),
     header=TRUE)
  MetaDCNRes$ModuleInControl <- ModuleInControl[which(ModuleInControl[,"FDR"] <
    FDRCutoff),]
  rownames(MetaDCNRes$ModuleInControl) <- NULL

  ### use the parameters to do module assembly
  MetaDCNRes$supermodule <- ModuleAssembly(MetaDCNRes$w1, FDRCutoff,
   caseName, controlName, pathwayDatabase, permutationTimes, outputPrefix)

  return(MetaDCNRes)
}


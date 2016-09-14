##' Detect meta-analysis differential co-expression network (DCN)
##'
##' This function will generate basic modules and supermodules
##'
##' @title MetaDCN
##' @param data a list of matrix (study) with row names as genes and column 
##' names as samples.
##' @param labels a list of vectors to specify the case/control labels for 
##' each subject corresponidng to the columns in each data matrix
##' @param caseName is the character string denoting case label
##' @param controlName is the character string denoting control label
##' @param meanFilter a number to specify qunatile cutoff for mean
##' @param SDFilter a number to specify qunatile cutoff for SD
##' @param FDRCutoff a number to specify FDR cutoff for basic modules
##' @param edgeCutoff a number to specify qunatile cutoff for correlation to 
##' to define an edge
##' @param MCSteps a number to specify the maximum number of MC steps for 
##' simulated annealing
##' @param jaccardCutoff a number to specify the maximum jaccard index allowed 
##' @param permutationTimes a number to specify how many permutations are used
##' @param repeatTimes a number to specify how many repeats of different seed 
##' modules to be used
##' @param outputFigure TRUE/FALSE to specify if figures are generated
##' @param outputPrefix a character string for output file prefix
##' @param parallel TRUE/FALSE to specify if parallel computing to be used
##' @param CPUNumbers a number to specify how many CPUs are used for parallel,
##' must be less than permutationTimes
##' @param pathwayDatabase a list with each element as a vector of 
##' pathway genes 
##' @param silent TRUE/FALSE to specify if suppress screen output
##' @return several csv files and zip files for figures
##' @author Li Zhu
##' @import snow
##' @import snowfall
##' @import igraph
##' @import genefilter
##' @export
##' @examples 
##' data(study.eg)

MetaDCN <- function(data, labels, caseName, controlName, meanFilter=0.2, 
  SDFilter=0.2, FDRCutoff=0.3, edgeCutoff=0.004, MCSteps=500, 
  jaccardCutoff=0.8, permutationTimes=10, repeatTimes=10, outputFigure=TRUE, 
  outputPrefix="MetaDCN", parallel=FALSE, CPUNumbers=10, pathwayDatabase, 
  silent=FALSE){
  
  if(parallel == TRUE){
    if(CPUNumbers > permutationTimes){
      stop("CPUNumbers should be smaller than or equal to permutationTimes!")
    }else if(CPUNumbers < 2){
      stop("CPUNumbers should be greater than or equal to 2!")
    }
  }
  
  caseIndex <- sapply(1:length(data), function(x) 
    which(labels[[x]] == caseName))
  controlIndex <- sapply(1:length(data), function(x) 
    which(labels[[x]] == controlName))

  caseStudyIndex <- seq(1:length(data))
  controlStudyIndex <- seq(length(data)+1, 2*length(data))

  ### generate network (permute_index is only needed for permutation parallel)
  NetworkGeneration(data, caseIndex, controlIndex, caseStudyIndex, 
    controlStudyIndex, meanFilter, SDFilter, edgeCutoff, choose="real",  
    permuteIndex=NA, silent=FALSE) 

  ### generate network and search modules for permutations
  if (parallel == TRUE){
    paraRep <- round(permutationTimes/CPUNumbers)
    CPUNumbers2 <- permutationTimes%%CPUNumbers

    ## generate network for permutations
    for (nr in 1:paraRep){
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers) 
      snowfall::sfExport("data", "caseIndex", "controlIndex", "caseStudyIndex",
        "controlStudyIndex", "meanFilter", "SDFilter", "edgeCutoff",
        "permutationTimes", "paraRep", "CPUNumbers2", "CPUNumbers", "nr")
      snowfall::sfClusterApply(seq(1, CPUNumbers), 
        function(x) NetworkGeneration(data, caseIndex, controlIndex, 
          caseStudyIndex, controlStudyIndex, meanFilter, SDFilter, 
          edgeCutoff, choose="permute", permuteIndex=((nr-1)*CPUNumbers+x), 
        silent=FALSE)) 
      snowfall::sfStop()
    }
    if(CPUNumbers2>0){
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("data", "caseIndex", "controlIndex", 
        "caseStudyIndex", "controlStudyIndex", "meanFilter", "SDFilter", 
        "edgeCutoff", "permutationTimes", "paraRep", "CPUNumbers2", 
        "CPUNumbers", "nr")
      snowfall::sfClusterApply(seq(1,CPUNumbers2), 
        function(x) NetworkGeneration(data, caseIndex, controlIndex, 
              caseStudyIndex, controlStudyIndex, choose="permute", meanFilter, 
              SDFilter, edgeCutoff, permuteIndex=(paraRep*CPUNumbers+x), 
              silent=FALSE)) 
      snowfall::sfStop()
    }
    
    ### generate energy list for permuted networks 
    for(nr in 1:paraRep){
      ## forward 
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2")
      snowfall::sfClusterApply(seq(1,CPUNumbers),function(x) 
        ModuleSearchPermutation(direction="forward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=(nr-1)*CPUNumbers+x))
      snowfall::sfStop()

      ## backward
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2")
      snowfall::sfClusterApply(seq(1,CPUNumbers),function(x) 
        ModuleSearchPermutation(direction="backward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=(nr-1)*CPUNumbers+x))
      snowfall::sfStop()
    }
    if(CPUNumbers2>0){
      ## forward 
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2")
      snowfall::sfClusterApply(seq(1,CPUNumbers2),function(x) 
        ModuleSearchPermutation(direction="forward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=paraRep*CPUNumbers+x))
      snowfall::sfStop()

      ## backward
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2")
      snowfall::sfClusterApply(seq(1,CPUNumbers2),function(x) 
        ModuleSearchPermutation(direction="backward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=paraRep*CPUNumbers+x))
      snowfall::sfStop()
    }

    ### generate true modules 
    directionTwo<-c("forward","backward")
    snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=length(directionTwo)) 
    snowfall::sfExport("pathwayDatabase", "MCSteps", "repeatTimes", 
      "permutationTimes", "caseName", "controlName", "outputFigure", 
      "outputPrefix", "jaccardCutoff", "directionTwo")
    snowfall::sfClusterApply(seq(1, length(directionTwo)), function(x) 
      ModuleSearch(direction=directionTwo[x], MCSteps, permutationTimes,  
      repeatTimes, jaccardCutoff, caseName, controlName, outputFigure, 
      outputPrefix, pathwayDatabase))
    snowfall::sfStop()

  }else {   ### without parallel
    ### generate permutation network
    sapply(seq(1,permutationTimes), function(x) 
      NetworkGeneration(data, caseIndex, controlIndex, 
                  caseStudyIndex, controlStudyIndex, meanFilter, 
                  SDFilter, edgeCutoff, choose="permute", 
                  permuteIndex=x, silent=FALSE))
    ### generate permutation energy list (forward)
    sapply(seq(1,permutationTimes),function(x) 
      ModuleSearchPermutation(direction="forward", MCSteps,
              permutationTimes, repeatTimes, jaccardCutoff, 
              permuteIndex=x))
    ### generate permutation energy list (backward)
    sapply(seq(1,permutationTimes),function(x) 
            ModuleSearchPermutation(direction="backward", MCSteps,
                    permutationTimes, repeatTimes, jaccardCutoff, 
                    permuteIndex=x))
    ### generate true modules
    directionTwo<-c("forward","backward")
    sapply(seq(1, length(directionTwo)), function(x) ModuleSearch(
      direction=directionTwo[x], MCSteps, permutationTimes,  
      repeatTimes, jaccardCutoff, caseName, controlName, outputFigure, 
      outputPrefix, pathwayDatabase))
  }

  ### selecte parameters
  forwardList <- read.csv(
    paste(outputPrefix, "_threshold_list_forward.csv", sep=""))[, 2:5]
  backwardList <- read.csv(
    paste(outputPrefix, "_threshold_list_backward.csv", sep=""))[, 2:5]
  indexMax <- which.max(rowSums(forwardList+backwardList))
  weightList <- seq(from=100, to=700, by=100)
  if(silent == FALSE){
    cat(paste("w1 parameter chosen is:", weightList[indexMax], "\n"))
  }

  weightTempIndex <- which(colnames(forwardList) == paste("FDR_", 
    round(FDRCutoff*100), sep=""))
  forwardNum <- forwardList[indexMax, weightTempIndex]
  backwardNum <- backwardList[indexMax, weightTempIndex]
  
  if(forwardNum == 0 & backwardNum ==0 ){
    stop("No basic module has FDR < 0.3!")
  }

  if(silent == FALSE){
    cat(paste(forwardNum, "modules are generated in forward direction\n"))
    cat(paste(backwardNum, "modules are generated in backward direction\n"))
  }
  
  ### use the parameters to do module assembly
  ModuleAssembly(weightList[indexMax], FDRCutoff, caseName, controlName, 
    pathwayDatabase, permutationTimes, outputPrefix)
}


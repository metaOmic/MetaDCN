##' Search for differentially co-expressed basic modules
##'
##' This function will search for basic modules differentially co-expressed 
##' between case and control
##' 
##' @title SearchBM
##' @param GeneNetRes a list from GeneNet function. 
##' @param MCSteps a number to specify the maximum number of MC steps for 
##' simulated annealing. Large number will be time-consuming.
##' @param jaccardCutoff a number to specify the maximum jaccard index allowed 
##' between basic modules.
##' @param repeatTimes a number to specify how many repeats of different seeds 
##' to be used. Large number will be time-consuming.
##' @param outputFigure TRUE/FALSE to specify if figures are generated.
##' @param silent TRUE/FALSE to specify if suppress screen output.
##' @return SearchBM will return a list and several Rdata, csv and png files
##' saved in the folder path specified in GeneNet inputs.
##' @return List of basic module information:
##' \item{w1 }{w1 weight with the most basic modules detected}
##' \item{BMInCase}{data matrix listing the information of basic modules
##' higher correlated in case}
##' \item{BMInControl}{data matrix listing the information of basic 
##' modules higher correlated in control}
##' @return permutation_energy_direction_p.Rdata is a list of energies of 
##' basic modules from permutation p.
##' @return basic_modules_summary_direction_weight_w.csv is a summary of basic
##' modules detected using weight w in forward/backward search.
##' @return threshold_direction.csv is a table listing number of basic modules 
##' detected under different FDRs in forward/backward search.
##' @return figure_basic_module_c_repeat_r_direction_weight_w.png is a plot of
##' basic module from component c repeat r using weight w in 
##' forward/backward search.
##' @author Li Zhu (liz86@pitt.edu)
##' @import snow
##' @import snowfall
##' @import igraph
##' @import genefilter
##' @export

SearchBM <- function(GeneNetRes, MCSteps=500, jaccardCutoff=0.8, 
  repeatTimes=10, outputFigure=TRUE, silent=FALSE){
  
  caseName <- GeneNetRes$caseName
  controlName <- GeneNetRes$controlName
  permutationTimes <- GeneNetRes$permutationTimes
  folder <- GeneNetRes$folder
  pathwayDatabase <- GeneNetRes$pathwayDatabase
  CPUNumbers <- GeneNetRes$CPUNumbers

  if(CPUNumbers > 1){
    if(CPUNumbers > permutationTimes){
      stop("CPUNumbers should be smaller than or equal to permutationTimes.")
    }else if(CPUNumbers < 2){
      stop("CPUNumbers should be greater than or equal to 2.")
    }
  }
  
  ### generate network and search modules for permutations
  if (CPUNumbers > 1){
    paraRep <- round(permutationTimes/CPUNumbers)
    CPUNumbers2 <- permutationTimes%%CPUNumbers

    ### generate energy list for permuted networks 
    for(nr in 1:paraRep){
      ## forward 
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2",
        "folder")
      NoOutput <- snowfall::sfClusterApply(seq(1,CPUNumbers),function(x) 
        ModuleSearchPermutation(direction="forward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=(nr-1)*CPUNumbers+x, folder))
      snowfall::sfStop()

      ## backward
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2",
        "folder")
      NoOutput <- snowfall::sfClusterApply(seq(1,CPUNumbers),function(x) 
        ModuleSearchPermutation(direction="backward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=(nr-1)*CPUNumbers+x, folder))
      snowfall::sfStop()
    }
    if(CPUNumbers2>0){
      ## forward 
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2",
        "folder")
      NoOutput <- snowfall::sfClusterApply(seq(1,CPUNumbers2),function(x) 
        ModuleSearchPermutation(direction="forward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=paraRep*CPUNumbers+x, folder))
      snowfall::sfStop()

      ## backward
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2",
        "folder")
      NoOutput <- snowfall::sfClusterApply(seq(1,CPUNumbers2),function(x) 
        ModuleSearchPermutation(direction="backward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=paraRep*CPUNumbers+x, folder))
      snowfall::sfStop()
    }

    ### generate true modules 
    directionTwo<-c("forward","backward")
    snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=length(directionTwo)) 
    snowfall::sfExport("pathwayDatabase", "MCSteps", "repeatTimes", 
      "permutationTimes", "caseName", "controlName", "outputFigure", 
      "folder", "jaccardCutoff", "directionTwo")
    NoOutput <- snowfall::sfClusterApply(seq(1, length(directionTwo)), function(x) 
      ModuleSearch(direction=directionTwo[x], MCSteps, permutationTimes,  
      repeatTimes, jaccardCutoff, caseName, controlName, outputFigure, 
      folder, pathwayDatabase))
    snowfall::sfStop()

  }else {   ### without parallel
    ### generate permutation energy list (forward)
    NoOutput <- sapply(seq(1,permutationTimes),function(x) 
      ModuleSearchPermutation(direction="forward", MCSteps,
              permutationTimes, repeatTimes, jaccardCutoff, 
              permuteIndex=x, folder))
    ### generate permutation energy list (backward)
    NoOutput <- sapply(seq(1,permutationTimes),function(x) 
            ModuleSearchPermutation(direction="backward", MCSteps,
                    permutationTimes, repeatTimes, jaccardCutoff, 
                    permuteIndex=x, folder))
    ### generate true modules
    directionTwo<-c("forward","backward")
    NoOutput <- sapply(seq(1, length(directionTwo)), function(x) ModuleSearch(
      direction=directionTwo[x], MCSteps, permutationTimes,  
      repeatTimes, jaccardCutoff, caseName, controlName, outputFigure, 
      folder, pathwayDatabase))
  }

  ### selecte parameters
  forwardList <- read.csv(
    paste(folder, "/threshold_forward.csv", sep=""))[, 2:5]
  backwardList <- read.csv(
    paste(folder, "/threshold_backward.csv", sep=""))[, 2:5]

  indexMax <- which.max(rowSums(forwardList + backwardList))
  weightList <- seq(from=100, to=700, by=100)
  w1 <- weightList[indexMax]
  if(silent == FALSE){
    cat(paste("w1 parameter chosen is:", w1, "\n"))
  }
  SearchBMRes <- list()
  SearchBMRes$w1 <- w1

  ### output basic modules
  SearchBMRes$BMInCase <- read.csv(paste(folder, 
    "/basic_modules_summary_forward_weight_", w1, ".csv", sep=""), 
    header=TRUE)
  rownames(SearchBMRes$BMInCase)<-NULL

  SearchBMRes$BMInControl <- read.csv(paste(folder, 
    "/basic_modules_summary_backward_weight_", w1, ".csv", sep=""),
     header=TRUE)
  rownames(SearchBMRes$BMControl) <- NULL

  return(SearchBMRes)
}


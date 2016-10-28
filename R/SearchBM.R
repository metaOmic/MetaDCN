##' Search for differentially co-expressed basic modules
##'
##' This function will search for basic modules differentially co-expressed 
##' between case and control
##' 
##' @title SearchBM
##' @param GeneNetRes a list from GeneNet function 
##' @param MCSteps a number to specify the maximum number of MC steps for 
##' simulated annealing
##' @param jaccardCutoff a number to specify the maximum jaccard index allowed 
##' between basic modules
##' @param repeatTimes a number to specify how many repeats of different seeds 
##' to be used
##' @param outputFigure TRUE/FALSE to specify if figures are generated
##' @param silent TRUE/FALSE to specify if suppress screen output
##' @return Several Rdata, csv and png files saved in the current working directory
##' \item{outputPrefix_permutation_energy_direction_p.Rdata}{A list of energies of basic modules from permutation p}
##' \item{outputPrefix_basic_modules_summary_direction_weight_w.csv }{A summary of basic modules found using weight w in forward/backward search}
##' \item{outputPrefix_threshold_direction.csv }{A table listing number of basic modules found under different FDRs in forward/backward search}
##' \item{outputPrefix_figure_basic_module_c_repeat_r_direction_weight_w.png }{Plot of basic module from component c from repeat r using weight w in forward/backward search}
##' @author Li Zhu
##' @import snow
##' @import snowfall
##' @import igraph
##' @import genefilter
##' @export
##' @examples 
##' data(pathwayDatabase)

SearchBM <- function(GeneNetRes, MCSteps=500, jaccardCutoff=0.8, 
  repeatTimes=10, outputFigure=TRUE, silent=FALSE){
  
  if(CPUNumbers > 1){
    if(CPUNumbers > permutationTimes){
      stop("CPUNumbers should be smaller than or equal to permutationTimes.")
    }else if(CPUNumbers < 2){
      stop("CPUNumbers should be greater than or equal to 2.")
    }
  }
  
  caseName <- GeneNetRes$caseName
  controlName <- GeneNetRes$controlName
  permutationTimes <- GeneNetRes$permutationTimes
  outputPrefix <- GeneNetRes$outputPrefix
  pathwayDatabase <- GeneNetRes$pathwayDatabase
  CPUNumbers <- GeneNetRes$CPUNumbers


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
        "outputPrefix")
      snowfall::sfClusterApply(seq(1,CPUNumbers),function(x) 
        ModuleSearchPermutation(direction="forward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=(nr-1)*CPUNumbers+x, outputPrefix))
      snowfall::sfStop()

      ## backward
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2",
        "outputPrefix")
      snowfall::sfClusterApply(seq(1,CPUNumbers),function(x) 
        ModuleSearchPermutation(direction="backward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=(nr-1)*CPUNumbers+x, outputPrefix))
      snowfall::sfStop()
    }
    if(CPUNumbers2>0){
      ## forward 
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2",
        "outputPrefix")
      snowfall::sfClusterApply(seq(1,CPUNumbers2),function(x) 
        ModuleSearchPermutation(direction="forward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=paraRep*CPUNumbers+x, outputPrefix))
      snowfall::sfStop()

      ## backward
      snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=CPUNumbers2) 
      snowfall::sfExport("jaccardCutoff", "MCSteps", "repeatTimes", 
        "permutationTimes", "paraRep", "CPUNumbers", "nr", "CPUNumbers2",
        "outputPrefix")
      snowfall::sfClusterApply(seq(1,CPUNumbers2),function(x) 
        ModuleSearchPermutation(direction="backward", MCSteps,
        permutationTimes, repeatTimes, jaccardCutoff, 
        permuteIndex=paraRep*CPUNumbers+x, outputPrefix))
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
    ### generate permutation energy list (forward)
    sapply(seq(1,permutationTimes),function(x) 
      ModuleSearchPermutation(direction="forward", MCSteps,
              permutationTimes, repeatTimes, jaccardCutoff, 
              permuteIndex=x, outputPrefix))
    ### generate permutation energy list (backward)
    sapply(seq(1,permutationTimes),function(x) 
            ModuleSearchPermutation(direction="backward", MCSteps,
                    permutationTimes, repeatTimes, jaccardCutoff, 
                    permuteIndex=x, outputPrefix))
    ### generate true modules
    directionTwo<-c("forward","backward")
    sapply(seq(1, length(directionTwo)), function(x) ModuleSearch(
      direction=directionTwo[x], MCSteps, permutationTimes,  
      repeatTimes, jaccardCutoff, caseName, controlName, outputFigure, 
      outputPrefix, pathwayDatabase))
  }

}


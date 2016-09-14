## This function will search basic modules in permutation studies
##
## @title ModuleSearchPermutation
## @param direction a character string, either "forward" or "backward"
## @param MCsteps a number to specify the maximum number of MC steps for 
## simulated annealing
## @param permutationTimes a number to specify how many permutations are used
## @param repeatTimes a number to specify how many repeats to be used
## @param jaccardCutoff a number to specify the maximum jaccard index allowed 
## @param permuteIndex a number to indicate which permute is working on, 
## only necessary when parallel=TRUE
## @return one RData file containing a list of module energies and module names
## @import igraph
## @author Li Zhu
ModuleSearchPermutation <- function(direction, MCSteps, permutationTimes,repeatTimes, jaccardCutoff, permuteIndex){
  
  options(stringsAsFactors = FALSE)
  set.seed(1234)
  
  load(paste("adjAllPermutation", permuteIndex, ".Rdata", sep=""))

  data <- adjAll
  geneNum <- length(genes)
  if (direction == "forward"){
    highStudyIndex <- caseStudyIndex
  }else{
    highStudyIndex <- controlStudyIndex
  }
  lowStudyIndex <- setdiff(seq(1,length(data)), highStudyIndex)

  ## construct the second order map and generate connected component
  edgeMap <- matrix(0, nrow=choose(geneNum,2), length(data))
  for (i in 1:length(data)) {
    edgeMap[, i] <- data[[i]][upper.tri(data[[i]])]
  }

  ## find the p value that occur on the disease side
  weight <- (seq(1,length(data)) %in% caseStudyIndex)*1
  pvalueTmp <- genefilter::rowttests(edgeMap, as.factor(weight))$p.value
  foldChange <- abs(rowMeans(edgeMap[, caseStudyIndex]) - 
    rowMeans(edgeMap[, controlStudyIndex]))
  indexTmp <- which(pvalueTmp < 0.1 & foldChange > 0.1)

  if(length(indexTmp) > 0){ # if such component exsits
    ## fill-in initial edges
    graphInitial <- matrix(0, nrow=geneNum, ncol=geneNum)
    graphInitial[upper.tri(graphInitial)][indexTmp] <- 1
    graphInitial <- matrix(graphInitial, nrow=geneNum)
    graphInitial <- makeSymm(graphInitial)
    rownames(graphInitial) <- 1:geneNum
    colnames(graphInitial) <- 1:geneNum
    gInitial <- igraph::graph.adjacency(graphInitial, mode="undirected", 
      add.colnames=TRUE)
    connectedComponent <- igraph::clusters(gInitial)
    
    ## find initial modules
    idConnectComponents <- which(connectedComponent$csize >= 3)
    if(length(idConnectComponents) > 0){  # if such component exsits
      componentMember <- list()
      for (i in 1:length(idConnectComponents)) {
        memberIndex <- which(connectedComponent$membership == 
          idConnectComponents[i])
        componentMember[[i]] <- memberIndex
      }
      
      ## search for differential modules
      weight1List <- seq(from=100, to=700, by=100)
      permutationEnergyList <- list()
      moduleNames <- list()
      countWeight <- 0
      for(weightTmp in weight1List){
        moduleNamesTmp <- NULL
        countWeight <- countWeight+1
        energyTmp <- NULL
        set.seed(1234)
        for (ccc in 1:length(componentMember)) {
          oldSetRepeat <- geneNameRepeat <- oldEnergyRepeat <- 
          moduleNamesRepeat <- list()
          for (rrr in 1:repeatTimes) {
            ##########################
            ###set search space bound 
            ##########################
            upperBound <- 30
            lowerBound <- 3
            ###############################
            ###initial simlation parameters
            ###############################
            if (length(componentMember[[ccc]]) <= 3) {
              next
            } 
            if (length(componentMember[[ccc]]) > 30) {
              memberIndex <- sample(componentMember[[ccc]], 15)
            } else {
              memberIndex <- componentMember[[ccc]]
            }
            oldSet <- index1Original <- memberIndex
            temp <- updateSearchSpace(data, oldSet, highStudyIndex)
            index2Original <- temp$set
            indexAll <- unique(c(index1Original, index2Original))
            
            ####configure data1 make it less memory intensive
            data1 <- lapply(data, function(x) x[indexAll, indexAll])
            gene1 <- genes[indexAll]
            oldSet <- 1:length(index1Original)
            oldTrialSet <- setdiff(1:length(indexAll), oldSet)
            
            #### Simulated annealing 
            accept <- 0
            trace <- NULL
            updateSteps <- 400
            oldEnergy <- energyEvaluation(data1, oldSet, highStudyIndex, 
              lowStudyIndex, w1=weightTmp)[4]
            KT <- oldEnergy/300 
            initialSet <- oldSet
            initalTrialSet <- oldTrialSet
            initialEnergy <- oldEnergy
            for (i in 1:MCSteps) {
              pRemove <- 0.5
              #####random number to decide transition (add or delete node)
              rand <- runif(1)
              if (length(oldSet) >= upperBound) {
                ###if module more than upper_bound, then delete
                temp <- removeNode(oldSet, oldTrialSet)    
              } else if (length(oldSet) <= lowerBound) {
                ###if module less than lower_bound, then add
                temp<- addNode(oldSet, oldTrialSet) 
              } else {  
                if (rand < pRemove) {    
                  temp <- removeNode(oldSet, oldTrialSet)
                } else {  
                  temp <- addNode(oldSet, oldTrialSet) 
                }
              }
              newSet <- temp$x
              newTrialSet <- temp$trialSet
              newEnergy <- energyEvaluation(data1, newSet, highStudyIndex, 
                lowStudyIndex, w1=weightTmp)[4]  
              if (i%%updateSteps == 0) {    
                KT <- KT*0.95
                ratio <- accept/updateSteps
                if (ratio > 0.5) {
                  KT <- KT*0.5
                }
                accept <- 0
                if (ratio < 0.02) {
                  break
                }
              }
              
              if (runif(1) < exp(-(newEnergy - oldEnergy)/KT)){
                oldEnergy <- newEnergy
                oldSet <- newSet
                oldTrialSet <- newTrialSet
                accept <- accept + 1
              }
              trace <- c(trace, oldEnergy)
            }
            geneName<-gene1[oldSet]
            
            if (rrr == 1) {
              geneNameRepeat <- c(geneNameRepeat, list(geneName))
              oldEnergyRepeat <- c(oldEnergyRepeat, list(oldEnergy))
              moduleNamesRepeat <- c(moduleNamesRepeat, list(paste(geneName,collapse="//")))
            } else {
              jaccard <- sapply(1:length(geneNameRepeat), function(x) 
                getJaccard(geneName, geneNameRepeat[[x]]))
              jaccardIndex <- jaccard >= jaccardCutoff # index if too much overlap
              energyIndex <- oldEnergy < unlist(oldEnergyRepeat) # index if energy is lower
              jaccardEnergyIndex <- which(jaccardIndex & energyIndex)
              
              if (sum(jaccardIndex) == 0) {     # no big overlap
                geneNameRepeat <- c(geneNameRepeat, list(geneName))
                oldEnergyRepeat <- c(oldEnergyRepeat, list(oldEnergy))
                moduleNamesRepeat <- c(moduleNamesRepeat, 
                  list(paste(geneName, collapse="//")))
              } else {                        # exist big overlap
                if (length(jaccardEnergyIndex) == 1){
                  geneNameRepeat[jaccardEnergyIndex] <- list(geneName)
                  oldEnergyRepeat[jaccardEnergyIndex] <- list(oldEnergy)
                  moduleNamesRepeat[jaccardEnergyIndex] <- 
                    list(paste(geneName, collapse="//"))
                } else if (sum(jaccardEnergyIndex) > 1) {
                  geneNameRepeat[jaccardEnergyIndex[1]] <- list(geneName)
                  geneNameRepeat <- geneNameRepeat[-jaccardEnergyIndex[-1]]
                  oldEnergyRepeat[jaccardEnergyIndex[1]] <- list(oldEnergy)
                  oldEnergyRepeat <- oldEnergyRepeat[-jaccardEnergyIndex[-1]]
                  moduleNamesRepeat[jaccardEnergyIndex[1]] <- list(paste(geneName,
                    collapse="//"))
                  moduleNamesRepeat <- moduleNamesRepeat[-jaccardEnergyIndex[-1]]
                  } 
                }
              }  
            }
          energyTmp <- c(energyTmp, oldEnergyRepeat)
          moduleNamesTmp <- c(moduleNamesTmp, moduleNamesRepeat)
        }
        permutationEnergyList[[countWeight]] <- energyTmp 
        moduleNames[[countWeight]] <- moduleNamesTmp
      }
    }
  }
  
  save(permutationEnergyList, moduleNames, 
    file=paste("permutation_network_energy_",direction, "_", permuteIndex, 
      ".Rdata",sep=""))
}  
  

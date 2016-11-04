## This function will search basic modules
##
## @title ModuleSearch
## @param direction a character string, either "forward" or "backward"
## @param MCsteps a number to specify the maximum number of MC steps for 
## simulated annealing
## @param permutationTimes a number to specify how many permutations are used
## @param repeatTimes a number to specify how many repeats to be used
## @param jaccardCutoff a number to specify the maximum jaccard index allowed
## @param caseName is the character string denoting case label
## @param controlName is the character string denoting control 
## @param outputFigure TRUE/FALSE to specify if figures are generated
## @param folder a character string for the folder name to store output file 
## prefix
## @param pathwayDatabase a list with each element as a vector of 
## pathway genes 
## @return one csv file for modules, one csv file for pathways, one csv file 
## for threshold, one zip file for basic module plots
## @import igraph
## @author Li Zhu
ModuleSearch<-function(direction, MCSteps, permutationTimes, repeatTimes,
  jaccardCutoff, caseName, controlName, outputFigure, folder, 
  pathwayDatabase){
  
  options(stringsAsFactors = FALSE) 
  set.seed(1234)
  
  ##############################
  ###1. use the real datasets
  ##############################
  pathwayDatabase <- pathwayDatabase[sapply(pathwayDatabase,length)<250]
  load(paste(folder, "/AdjacencyMatrice.Rdata", sep=""))
  data <- adjAll
  studyNum <- length(data)/2
  studyName <- c(paste(caseName, 1:studyNum), paste(controlName, 1:studyNum))
  if(direction == "forward"){
    highStudyIndex <- caseStudyIndex
  }else{
    highStudyIndex <- controlStudyIndex
  }
  lowStudyIndex <- setdiff(seq(1, length(data)), highStudyIndex)
  numGene <- length(genes)

  ##########################################################################
  ####step 2 construct the second order map and generate connected component
  ##########################################################################
  edgeMap <- matrix(0, nrow=choose(numGene,2), length(data))
  for(i in 1:length(data)){
    edgeMap[,i]=data[[i]][upper.tri(data[[i]])]
  }
  colnames(edgeMap) <- studyName

  ###find the p value that occur on the disease side
  weight <- (seq(1,length(data)) %in% caseStudyIndex)*1
  pvalueTmp <- genefilter::rowttests(edgeMap, as.factor(weight))$p.value
  foldChange <- abs(rowMeans(edgeMap[, caseStudyIndex]) - 
    rowMeans(edgeMap[,controlStudyIndex]))
  indexTmp <- which(pvalueTmp < 0.1 & foldChange > 0.1)
  if(length(indexTmp) == 0){
    stop("No edges found! Please reduce edgeCutoff, meanFilter, or SDFilter")
  }

  ####put the seeds in a network and try to get the connected component as the seeds
  graphInitial <- matrix(0, nrow=numGene, ncol=numGene)
  graphInitial[upper.tri(graphInitial)][indexTmp] <- 1
  graphInitial <- matrix(graphInitial, nrow=numGene)
  graphInitial <- makeSymm(graphInitial)
  rownames(graphInitial) <- 1:numGene
  colnames(graphInitial) <- 1:numGene
  gInitial <- igraph::graph.adjacency(graphInitial, mode="undirected", 
    add.colnames=TRUE)
  connectedComponent <- igraph::clusters(gInitial)
  
  ####only consider connected component more than size 3
  idConnectComponents <- which(connectedComponent$csize >= 3)
  if(length(idConnectComponents) == 0){
    stop("No connected componenets found! Please reduce edgeCutoff, meanFilter, or SDFilter")
  }
  componentMember <- list(length(idConnectComponents))
  for (i in 1:length(idConnectComponents)) {
    memberIndex <- which(connectedComponent$membership == 
      idConnectComponents[i])
    componentMember[[i]] <- memberIndex
  }
  
  #####################################################################
  ####step 3 sample differential modules
  #####################################################################
  ###to tune to parameter
  weight1List <- seq(from=100, to=700, by=100)
  countWeight <- 0
  countWeight1 <- 0
  thresholdList <- matrix(0, length(weight1List), 4)
  rownames(thresholdList) <- weight1List
  colnames(thresholdList) <- c("FDR_10","FDR_20","FDR_30","FDR_40")
  for(weightTmp in weight1List){
    #system(paste("mkdir /",weightTmp,"BasicModulePlot/",sep=""))
    countWeight <- countWeight+1;
    countWeight1 <- countWeight1+1;
    set.seed(1234)
    print(paste("use weight1=", weightTmp))
    summary <- matrix(0, nrow=length(componentMember)*repeatTimes, ncol=13)
    colnames(summary) <- c("Component Number", "Repeat Index", 
      "Gene Set", "Size", "Module Energy", "diff_mean","diff_var", 
      "Density in Case", "Density in Control", "pathway name", 
      "Pathway p-value", "Pathway q-value", "Matched genes")
    pathwayResult <- NULL
    count <- 0
    for (ccc in 1:length(componentMember)) {
      if (length(componentMember[[ccc]]) <= 3) {
        next
      }
      oldSetRepeat <- geneNameRepeat <- oldEnergyRepeat <- list()
      for(rrr in 1:repeatTimes){
        ##########################
        ###set search space bound 
        ##########################
        upperBound <- 30
        lowerBound <- 3
        
        ###############################
        ###initial simlation parameters
        ###############################
        if (length(componentMember[[ccc]]) > 30) {
          memberIndex  <-  sample(componentMember[[ccc]], 10)  
        } else {
          memberIndex  <-  componentMember[[ccc]]
        }
        oldSet  <-  index1Original  <-  memberIndex
        temp <- updateSearchSpace(data, oldSet, highStudyIndex)
        index2Original <- temp$set
        indexAll <- unique(c(index1Original, index2Original))
        
        ####configure data1 make it less memory intensive
        data1 <- lapply(data, function(x) x[indexAll,indexAll])
        genes1 <- genes[indexAll]
        oldSet <- 1:length(index1Original)
        oldTrialSet <- setdiff(1:length(indexAll), oldSet)
    
        accept <- 0
        trace <- NULL
        updateSteps <- 400
        oldEnergy <- energyEvaluation(data1, oldSet, highStudyIndex, 
          lowStudyIndex, w1=weightTmp)[4]
        KT <- oldEnergy/300
        initialSet <- oldSet
        initalTrialSet <- oldTrialSet
        initialEnergy <- oldEnergy
        for(i in 1:MCSteps){
          pRemove <- 0.5
          #####random number to decide transition (add or delete node)
          rand <- runif(1)
          if (length(oldSet)>=upperBound) {
            ###if module more than upper_bound, then delete
            temp <- removeNode(oldSet, oldTrialSet)
            pNewToOld <- 1/(length(temp$trialSet))
            pOldToNew <- 1/(length(temp$x))
          } else if (length(oldSet) <= lowerBound){
            ###if module less than lower_bound, then add
            temp <- addNode(oldSet, oldTrialSet) 
            pOldToNew <- 1/(length(temp$trialSet))
            pNewToOld <- 1/(length(temp$x))
          }else{  
            if(rand<pRemove){    
              temp <- removeNode(oldSet, oldTrialSet)
              pNewToOld <- 1/(length(temp$trialSet))
              pOldToNew <- 1/(length(temp$x))
            }else{  
              temp <- addNode(oldSet, oldTrialSet) 
              pOldToNew <- 1/(length(temp$trialSet))
              pNewToOld <- 1/(length(temp$x))
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
          
          if(runif(1) < exp(-(newEnergy - oldEnergy)/KT)*pNewToOld/pOldToNew) {
            oldEnergy <- newEnergy
            oldSet <- newSet
            oldTrialSet <- newTrialSet
            accept <- accept + 1
          }
          trace <- c(trace, oldEnergy)
        }
        geneName <- genes1[oldSet]
        
        if(rrr==1){
          geneNameRepeat <- c(geneNameRepeat, list(geneName))
          oldEnergyRepeat <- c(oldEnergyRepeat, list(oldEnergy))
        }else{
          jaccard <- sapply(1:length(geneNameRepeat), function(x) 
            getJaccard(geneName, geneNameRepeat[[x]]))
          jaccardIndex <- jaccard >= jaccardCutoff # index if too much overlap
          energyIndex <- oldEnergy < unlist(oldEnergyRepeat) # index if energy is lower
          jaccardEnergyIndex <- which(jaccardIndex & energyIndex)
          
          if (sum(jaccardIndex) == 0) { # no big overlap
            geneNameRepeat <- c(geneNameRepeat, list(geneName))
            oldEnergyRepeat <- c(oldEnergyRepeat, list(oldEnergy))
          } else { # exist big overlap
            if (length(jaccardEnergyIndex) == 1) {
              geneNameRepeat[jaccardEnergyIndex] <- list(geneName)
              oldEnergyRepeat[jaccardEnergyIndex] <- list(oldEnergy)
            }else if (sum(jaccardEnergyIndex) > 1) {
              geneNameRepeat[jaccardEnergyIndex[1]] <- list(geneName)
              geneNameRepeat <- geneNameRepeat[-jaccardEnergyIndex[-1]]
              oldEnergyRepeat[jaccardEnergyIndex[1]] <- list(oldEnergy)
              oldEnergyRepeat <- oldEnergyRepeat[-jaccardEnergyIndex[-1]]
            }
          } 
        }
      }
      
      for(rrr in 1:length(geneNameRepeat)){
        ###print the final configuration
        pathwayInfo <- gsaFisher(geneNameRepeat[[rrr]], genes, pathwayDatabase,topNum=3, sort=TRUE)
        if (outputFigure == TRUE) {
          png(file=paste(folder, "/Basic_modules_figures_weight_", weightTmp, 
            "/Basic_module_component_", ccc, 
            "_repeat_", rrr, "_weight_", weightTmp, "_", direction, ".png", 
            sep=""), 
          width = 600, 
          height = 480, units = "px", pointsize = 18)
          printNetworks(data, geneNameRepeat[[rrr]], studyName, 
            nodeName=genes1, a=2, b=studyNum)
          dev.off()
        }
        
        summary[count+rrr,] <- c(ccc, rrr, paste(geneNameRepeat[[rrr]], 
          collapse="/"), length(geneNameRepeat[[rrr]]), oldEnergyRepeat[[rrr]],
          mean(getDensity(data, geneNameRepeat[[rrr]], caseStudyIndex) - getDensity(data, geneNameRepeat[[rrr]], controlStudyIndex)), 
          sd(getDensity(data, geneNameRepeat[[rrr]], caseStudyIndex) - getDensity(data, geneNameRepeat[[rrr]], controlStudyIndex)), 
          paste(format(getDensity(data, geneNameRepeat[[rrr]], caseStudyIndex),
            digits=2), collapse="//"), 
          paste(format(getDensity(data, geneNameRepeat[[rrr]], controlStudyIndex), digits=2), collapse="//"), 
          paste(rownames(pathwayInfo), collapse="///"), 
            paste(pathwayInfo,sep="///"))
        pathwayInfo <- gsaFisher(geneNameRepeat[[rrr]], genes, pathwayDatabase,
          topNum=length(pathwayDatabase), sort=FALSE)
        pathwayResult <- cbind(pathwayResult, pathwayInfo[,1])
      }
      count <- count+length(geneNameRepeat)
    }
    
    ####try to calculate FDR for different cutoff
    summary <- summary[1:count,]
    permutationEnergys <- list()
    for (mmm in 1:permutationTimes) {
      load(paste(folder, "/permutation_energy_",direction, "_", mmm, 
        ".Rdata",sep=""))
      permutationEnergys[[mmm]] <- permutationEnergyList
    }
    summaryFDR <- matrix(0, dim(summary)[1], 2)
    colnames(summaryFDR) <- c("p_value", "FDR")
    for(FDRIndex in 1:dim(summary)[1]){ 
      energyTrue <- as.numeric(summary[FDRIndex,5])
      temp1 <- sapply(1:permutationTimes,function(x) sum(permutationEnergys[[x]][[countWeight]] <= energyTrue))
      temp2 <- sapply(1:permutationTimes,function(x) length(permutationEnergys[[x]][[countWeight]]))
      summaryFDR[FDRIndex, 1] <- (sum(temp1)+1)/(sum(temp2)+1)
    }
    summaryFDR[,2] <- p.adjust(summaryFDR[,1], method="BH")
    summary <- cbind(summary, summaryFDR)
    pathwayResult <- apply(pathwayResult, 2, as.numeric)
    rownames(pathwayResult) <- names(pathwayDatabase)
        
    #########criterions including FDR and size of the module    
    ## consider all repeats
    thresholdList[countWeight1, 1]  <-  sum(summaryFDR[, 2] < 0.1)
    thresholdList[countWeight1, 2]  <-  sum(summaryFDR[, 2] < 0.2)
    thresholdList[countWeight1, 3]  <-  sum(summaryFDR[, 2] < 0.3)
    thresholdList[countWeight1, 4]  <-  sum(summaryFDR[, 2] < 0.4)
    
    write.csv(summary, file=paste(folder, "/basic_modules_summary_", direction, "_weight_", weightTmp, ".csv", sep=""),row.names=FALSE)
  }
  write.csv(thresholdList, file=paste(folder, "/threshold_", 
    direction, ".csv", sep=""))
}


## This function will assemble modules
##
## @title ModuleAssembly
## @param weightChosen a number indicating the weight chosen
## @param FDRCutoff a number to FDR cutoff for basic modules
## @param caseName is the character string denoting case label
## @param controlName is the character string denoting control
## @param pathwayDatabase a list with each element as a vector of 
## pathway genes 
## @param permutationTimes a number to specify how many permutations are used
## @return one csv file for supermodule summary, a zip file for Cytoscape
## @author Li Zhu
ModuleAssembly <- function(weightChosen, FDRCutoff, caseName, controlName, 
  pathwayDatabase, permutationTimes, folder){
  options(stringsAsFactors = FALSE)
  set.seed(1234)

  pathwayDatabase <- pathwayDatabase[sapply(pathwayDatabase,length) < 250]
  pathLength <- sapply(1:length(pathwayDatabase), function(x) 
    length(pathwayDatabase[[x]]))
  names(pathLength) <- names(pathwayDatabase)
  summaryFDRModules <- read.csv(paste(folder, 
    "/basic_modules_summary_forward_weight_", weightChosen, ".csv", sep=""))
  summaryFDRModules1 <- read.csv(paste(folder, 
    "/basic_modules_summary_backward_weight_", weightChosen, ".csv", sep=""))
  
  # apply FDR cutoff
  summaryFDRModules <- summaryFDRModules[which(
    summaryFDRModules[,"FDR"] <= FDRCutoff), ]
  summaryFDRModules1 <- summaryFDRModules1[
  which(summaryFDRModules1[, "FDR"] <= FDRCutoff),]
  forwardIndex <- seq(1:nrow(summaryFDRModules))
  forwardIndex <- sapply(1:length(forwardIndex), 
    function(x) paste("H",forwardIndex[x], sep=""))
  backwardIndex <- seq(1:nrow(summaryFDRModules1))
  backwardIndex <- sapply(1:length(backwardIndex), 
    function(x) paste("L", backwardIndex[x], sep=""))
  
  forBackIndex <- c(forwardIndex, backwardIndex)
  
  moduleList <- lapply(summaryFDRModules[, 3], function(x) 
    strsplit(x, split="/")[[1]])
  moduleList1 <- lapply(summaryFDRModules1[, 3], function(x) 
    strsplit(x, split="/")[[1]])
  moduleListAll <- c(moduleList, moduleList1)  # combine forward and backward
  
  load(paste(folder, "/AdjacencyMatrice.Rdata", sep=""))
  load(paste(folder, "/CorrelationMatrice.Rdata", sep=""))
  data <- adjAll
  studyNum <- length(data)/2
  studyName <- c(paste(caseName, 1:studyNum), paste(controlName, 1:studyNum))
  
  background <- unique(unlist(pathwayDatabase))
  summaryPathwayEnrichment <- sapply(moduleListAll, function(x) 
    as.numeric(gsaFisher(x, background, pathway=pathwayDatabase, sort=FALSE)[,1])) # pathway*modules matrix with p-values
  summaryPathwayEnrichment <- apply(summaryPathwayEnrichment, 2, as.numeric)
  rownames(summaryPathwayEnrichment) <- names(pathwayDatabase)
  
  combineStat <- rowSums(2*log(summaryPathwayEnrichment))
  summaryPathwayEnrichmentSort <- summaryPathwayEnrichment[order(combineStat),]
  
  pathwayThreshold <- 0.05
  topModuleAssembly <- 150
  moduleAssemblySummary <- matrix(0, nrow=topModuleAssembly, ncol=13)
  colnames(moduleAssemblySummary) <- c("pathway_name", "pathway_size",
    "p_value", "size", "gene_in_set", "q_value", "module_num", 
    "module_index","density1","density2", "mean_diff", "sd_diff", "gene_in_path")
  count <- 0
  for (i in 1:topModuleAssembly) {
    pathwayName <- rownames(summaryPathwayEnrichmentSort)[i]
    indexModuleEnrichPathway <- which(as.numeric(
      summaryPathwayEnrichmentSort[i,]) < pathwayThreshold)
    if (length(indexModuleEnrichPathway) >= 2) {
      groupGenes <- moduleListAll[indexModuleEnrichPathway]
      tempLength <- sapply(groupGenes, length)      
      pathwayIndex <- which(names(pathwayDatabase) == pathwayName)
      
      if (length(groupGenes) >= 3) {
        temp1 <- lapply(apply(matrix(combn(1:length(groupGenes), 2),nrow=2), 2, list), unlist)
        temp2 <- lapply(apply(matrix(combn(1:length(groupGenes), 3),nrow=3), 2, list), unlist)
        allCombination <- c(temp1,temp2)
        }else{
          temp1 <- lapply(apply(matrix(combn(1:length(groupGenes), 2),nrow=2), 2, list), unlist)
          allCombination <- temp1
        }
        pvalueTmp <- array(0,dim=length(allCombination))
        for (j in 1:length(allCombination)) {
          groupTmp <- unique(unlist(groupGenes[allCombination[[j]]]))
          pathwayInfo <- gsaFisher(groupTmp, genes, 
            pathwayDatabase[pathwayIndex], topNum=1, sort=FALSE)
          pvalueTmp[j] <- pathwayInfo$pvalue
        }
        indexMin <- which.min(pvalueTmp)
        indexComb <- allCombination[[indexMin]]
        
        groupGenesUnique <- unique(unlist(groupGenes[indexComb]))
        pathwayInfo <- gsaFisher(groupGenesUnique, genes, pathwayDatabase[pathwayIndex], topNum=1, sort=FALSE)
        enrichPvalue <- pathwayInfo$pvalue
        matchPathwayGenes <- strsplit(pathwayInfo$matchedGene, split="/")[[1]]
        
        nodeList <- matrix(0, nrow=length(groupGenesUnique), ncol=2)
        nodeList[,1] <- groupGenesUnique
        indexTmp <- match(matchPathwayGenes, groupGenesUnique)
        nodeList[indexTmp, 2] <- TRUE
        nodeList[-indexTmp, 2] <- FALSE
        
        indexCombIndex <- t(sapply(groupGenesUnique, function(x) 
          sapply(indexComb, function(y) sum(which(groupGenes[[y]] == x))) != 0))
        nodeList <- cbind(nodeList, indexCombIndex, 
          apply(indexCombIndex, 1, function(x) paste(x, collapse="_")))
        colnames(nodeList) <- c("Gene", "Hit", 
          paste("Module", 1:dim(indexCombIndex)[2], sep=""), "Class")

        oldSet <- match(groupGenesUnique,genes)
        density1 <- getDensity(data, oldSet, caseStudyIndex) 
        density2 <- getDensity(data, oldSet, controlStudyIndex)
        meanDiff <- mean(density1 - density2)
        sdDiff <- sd(density1 - density2)
        medianSD <- median(apply(corAll[[1]], 1, sd))
        corMatrix <- lapply(corAll, function(x) x[oldSet, oldSet])
        
        genesNew <- genes[oldSet]
        pairDiffNetwork <- lapply(1:length(caseStudyIndex), function(x) 
         corMatrix[[caseStudyIndex[x]]] - corMatrix[[controlStudyIndex[x]]])
        meanDiffNetwork <- Reduce("+", pairDiffNetwork)/length(caseStudyIndex)
        sdDiffNetwork <- lapply(1:length(caseStudyIndex), function(x) 
          (pairDiffNetwork[[x]] - meanDiffNetwork)^2)
        sdDiffNetwork1 <- sqrt(Reduce("+", sdDiffNetwork)/length(caseStudyIndex))
        
        edgeName <- outer(oldSet, oldSet, paste)
        edgeName1 <- edgeName[upper.tri(edgeName)]
        edgeNameColumn <- t(sapply(edgeName1, function(x) 
          strsplit(x, split=" ")[[1]]))
        edgeNameColumn1 <- apply(edgeNameColumn, 2, as.numeric)
        weightValue <- meanDiffNetwork/(sdDiffNetwork1 + medianSD)
        weightValue1 <- weightValue[upper.tri(weightValue)]
        
        weightValue1Abs <- abs(weightValue1)
        
        pvalueList <- array(0, dim=length(weightValue1Abs))
        for(ttt in 1:permutationTimes){
          groupPermute <- sample(c(caseStudyIndex, controlStudyIndex))
          index1 <- 1:length(caseStudyIndex)
          groupPermute1 <- groupPermute[index1]
          groupPermute2 <- groupPermute[-index1]
          pairDiffNetworkPermute <- lapply(1:length(caseStudyIndex), function(x) 
            corMatrix[[groupPermute1[x]]] - corMatrix[[groupPermute2[x]]])
          meanDiffNetworkPermute <- Reduce("+", pairDiffNetworkPermute)/length(caseStudyIndex)
          sdDiffNetworkPermute <- lapply(1:length(caseStudyIndex),function(x) (pairDiffNetworkPermute[[x]]-meanDiffNetworkPermute)^2)
          sdDiffNetwork1Permute <- sqrt(Reduce("+",sdDiffNetworkPermute)/length(caseStudyIndex))
          weightValuePermute <- meanDiffNetworkPermute/(sdDiffNetwork1Permute+medianSD)
          weightValue1Permute <- weightValuePermute[upper.tri(weightValuePermute)]
          rank1 <- rank(-weightValue1Abs);
          rank2 <- rank(c(-weightValue1Abs,-abs(weightValue1Permute)))[1:length(weightValue1Abs)]
          pvalueList <- pvalueList + rank2 - rank1
        }      
        pvalueList <- pvalueList/(length(weightValue1Abs)*permutationTimes)
        
        color <- array(0, dim=length(weightValue1Abs))
        color[weightValue1 >= 0] <- "red"
        color[weightValue1 < 0] <- "blue"
        
        edgeList <- cbind(genes[as.numeric(edgeNameColumn1[,1])], genes[as.numeric(edgeNameColumn1[,2])], weightValue1, pvalueList)
        colnames(edgeList) <- c("GeneA", "GeneB", "Zscore", "Pvalue")

        write.table(edgeList, file=paste(folder, "/CytoscapeFiles/", pathwayName, 
          "_edge_list.txt", 
          sep=""), row.names=FALSE, 
        quote=FALSE, sep="\t")
        
        write.table(nodeList, file=paste(folder, "/CytoscapeFiles/", pathwayName, 
          "_node_list.txt", sep=""), row.names=FALSE, 
        quote=FALSE, sep="\t")      
        if (length(indexComb) >= 1) {
          count <- count+1
          moduleAssemblySummary[count, 1] <- pathwayName
          moduleAssemblySummary[count, 2] <- pathLength[pathwayName]
          moduleAssemblySummary[count, 3] <- enrichPvalue
          moduleAssemblySummary[count, 4] <- length(groupGenesUnique)
          moduleAssemblySummary[count, 5] <- length(matchPathwayGenes)
          moduleAssemblySummary[count, 7] <- length(indexComb)
          moduleAssemblySummary[count, 8] <- paste(forBackIndex[indexModuleEnrichPathway[indexComb]],collapse=",", sep="")
          moduleAssemblySummary[count, 9] <- paste(format(density1, digits=2),
            collapse="//")
          moduleAssemblySummary[count, 10] <- paste(format(density2, digits=2),
            collapse="//")
          moduleAssemblySummary[count, 11] <- meanDiff
          moduleAssemblySummary[count, 12] <- sdDiff
          indexPath <- which(names(pathwayDatabase) == pathwayName)
          intersectGene <- intersect(toupper(pathwayDatabase[[indexPath]]),
            toupper(groupGenesUnique))
          moduleAssemblySummary[count, 13] <- paste(intersectGene, 
            collapse=",")
        }
      }else{
        next
      }
    }    
    moduleAssemblySummary <- moduleAssemblySummary[1:count, ]
    
    moduleAssemblySummary[,6] <- p.adjust(as.numeric(moduleAssemblySummary[,2])
      ,method="BH")
    moduleAssemblySummary <- moduleAssemblySummary[order(as.numeric(
      moduleAssemblySummary[,6])),]
    
    write.csv(moduleAssemblySummary, 
      file=paste(folder, "/module_assembly_summary_weight_", 
        weightChosen, ".csv", sep=""))

    return(moduleAssemblySummary)      
  }






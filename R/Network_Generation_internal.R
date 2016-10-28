## This function will generate corelation and adjacency matrices
##
## @title NetworkGeneration
## @param data a list of matrix with length as 2 times the number of studies
## each matrix is a study of case/control with rows as genes and columns as 
## samples
## @param caseIndex a list of vectors indicating the index for cases in each 
## study
## @param controlIndex a list of vectors indicating the index for cases in 
## each study
## @param caseStudyIndex a vector indicating the study index for cases
## @param controlStudyIndex a vector indicating the study index for controls
## @param meanFilter a number to specify qunatile cutoff for mean
## @param SDFilter a number to specify qunatile cutoff for SD
## @param edgeCutoff a number to specify qunatile cutoff for correlation to 
## to define an edge
## @param choose is either "real" for study, or "permute" for permutated
## studies
## @param permuteIndex a number to indicate which permute is working on, 
## only necessary when parallel=TRUE
## @param silent TRUE/FALSE to specify if suppress screen output
## @return two RData files are saved in working directory, with each 
## containing a data list of correlation matrices or aadjacency matrices, a 
## vector of gene names, caseStudyIndex, and controlStudyIndex.
## @author Li Zhu

NetworkGeneration <- function(data, caseIndex, controlIndex, caseStudyIndex,
  controlStudyIndex, meanFilter, SDFilter, edgeCutoff, choose, permuteIndex, 
  outputPrefix, silent=FALSE) {

  options(stringsAsFactors = FALSE)
  set.seed(1234)
  
  ##some filtering based on rank sum
  meanRank <- apply(sapply(data, function(x) rowMeans(x)), 2, rank)
  SDRank <- apply(sapply(data, function(x) apply(x, 1, sd)), 2, rank)
  meanMean <- rowSums(meanRank)
  SDMean <- rowSums(SDRank)
  meanCutoff <- quantile(meanMean, meanFilter, na.rm=TRUE)
  filterIndex <- which(meanMean >= meanCutoff)
  SDCutoff <- quantile(SDMean[filterIndex], SDFilter, na.rm=TRUE)
  filterIndex1 <- which(SDMean[filterIndex] >= SDCutoff)
  filterIndex <- filterIndex[filterIndex1]
  if(silent == FALSE){
    cat(paste("total gene: ", dim(data[[1]])[1], "\n", sep=""))
  }
  
  for (i in 1:length(data)){
    data[[i]] <- data[[i]][filterIndex,]
  }
  genes <- rownames(data[[1]])
  if(silent == FALSE){
    cat(paste("after filtering gene: ", dim(data[[1]])[1], "\n", sep=""))
  }
  
  #############################################
  ###get adjacency matrix for each network
  #############################################
  geneNum <- length(genes)
  dataTmp <- list()
  
  ## Permute or not
  for (i in 1:length(data)) {
    dataTmp[[i]] <- data[[i]][1:geneNum,]
    if (choose == "permute") {
      tmp <- sample(1:dim(dataTmp[[i]])[2])
      dataTmp[[i]] <- dataTmp[[i]][,tmp]
    }
  }
    
  corDisease <- list()
  corControl <- list()
  for (i in 1:length(data)) {
    if(silent == FALSE){
      cat(paste("network pair ", i, " constructed\n"))
    }  
    corDisease[[i]] <- cor(t(dataTmp[[i]][, caseIndex[[i]]]),
      use="pairwise.complete.obs", method="spearman")
    ###set all NAN to correlation 0
    corDisease[[i]][is.na(corDisease[[i]])] <- 0
    ###set the diagonal term to 0
    diag(corDisease[[i]]) = 0 
    corControl[[i]]<-cor(t(dataTmp[[i]][, controlIndex[[i]]]), 
      use="pairwise.complete.obs", method="spearman")
    ###set all NAN to correlation 0
    corControl[[i]][is.na(corControl[[i]])] <- 0
    ###set the diagonal term to 0
    diag(corControl[[i]]) = 0
  }  
  
  ## cutoff for making edges
  adjDisease <- list()
  adjControl <- list()
  for (i in 1:length(data)) {
    adjDisease[[i]] <- abs(corDisease[[i]])
    adjCutoff <- quantile(abs(corDisease[[i]]), 1 - edgeCutoff)
    adjDisease[[i]][adjDisease[[i]] >= adjCutoff] <- 1
    adjDisease[[i]][adjDisease[[i]] < adjCutoff] <- 0
    
    ## calculate the property for control network
    adjControl[[i]] <- abs(corControl[[i]])
    adjCutoff <- quantile(abs(corControl[[i]]), 1 - edgeCutoff)
    adjControl[[i]][adjControl[[i]] >= adjCutoff] <- 1
    adjControl[[i]][adjControl[[i]] < adjCutoff] <- 0
  }
  
  ## remove isolated nodes
  removeCaseIndex <- which(rowSums(sapply(adjDisease, rowSums)) == 0)
  removeControlIndex <- which(rowSums(sapply(adjControl, rowSums)) == 0)
  removeIndex <- intersect(removeCaseIndex, removeControlIndex)
  
  keepIndex <- setdiff(1:dim(adjDisease[[1]])[1], removeIndex)
  genes <- genes[keepIndex]

  ## create adjacent matrix and correlation matrix
  adjAll=list()
  corAll=list()
  
  for (i in caseStudyIndex) {
    adjAll[[i]] <- adjDisease[[i]][keepIndex, keepIndex]
  }
  for (i in controlStudyIndex) {
    adjAll[[i]] <- adjControl[[i-length(data)]][keepIndex, keepIndex]
  }
  
  for (i in caseStudyIndex) {
    corAll[[i]] <- corDisease[[i]][keepIndex, keepIndex]
  }
  for (i in controlStudyIndex) {
    corAll[[i]] <- corControl[[i-length(data)]][keepIndex, keepIndex]
  }
    
  ## save matrix
  if (choose == "real") {
    save(adjAll, genes, caseStudyIndex, controlStudyIndex, 
      file=paste(outputPrefix, "_AdjacencyMatrice.Rdata", sep=""))
    save(corAll, genes, caseStudyIndex, controlStudyIndex, 
      file=paste(outputPrefix, "_CorrelationMatrice.Rdata", sep=""))
  }else if (choose == "permute") {
    save(adjAll, genes, caseStudyIndex, controlStudyIndex, 
      file=paste(outputPrefix, "_AdjacencyMatricePermutation", 
        permuteIndex, ".Rdata", sep=""))
  }  
    
}





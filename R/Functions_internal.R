## Some internal functions
## @author Li Zhu

########################
# pathway analysis
########################
gsaFisher <- function(x, background, pathway, topNum=length(pathway), 
  sort=TRUE) {
  countTable <- matrix(0, 2, 2)
  x <- toupper(x)
  background <- toupper(background)
  index <- which(toupper(background) %in% toupper(x) == FALSE)
  backgroundNonGeneList <- background[index]
  x <- toupper(x)
  pathway <- lapply(pathway, function(x) intersect(toupper(background), 
    toupper(x)))
  getFisher <- function(path) {
    res <- NA
    ####in the gene list and in the pathway
    countTable[1, 1] <- sum(x %in% path)
    ####in the gene list but not in the pathway
    countTable[1, 2] <- length(x) - countTable[1,1]
    ####not in the gene list but in the pathway
    countTable[2, 1] <- sum(backgroundNonGeneList %in% path)
    ####not in the gene list and not in the pathway
    countTable[2, 2] <- length(backgroundNonGeneList) - countTable[2, 1]       
    matchedGene <- x[x %in% path]
    matchNum <- length(matchedGene)
    overlapInfo <- array(0, dim=4)
    names(overlapInfo) <- c("DE in Geneset", "DE not in Genese", 
      "NonDE in Geneset", "NonDE out of Geneset")
    overlapInfo[1] <- countTable[1, 1]
    overlapInfo[2] <- countTable[1, 2]
    overlapInfo[3] <- countTable[2, 1]
    overlapInfo[4] <- countTable[2, 2]
    if (length(countTable) == 4) {
      res <- fisher.test(countTable, alternative="greater")$p}
    return(list(pvalue=res, matchedGene=matchedGene, matchNum=matchNum,
                fisherTable=overlapInfo))
  }
  pval <- array(0, dim=length(pathway))
  
  matchedGeneList <- list(length(pathway))
  num1 <- array(0, dim=length(pathway))
  num2 <- matrix(0, nrow=length(pathway), ncol=4)
  colnames(num2) <- c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
  for (i in 1:length(pathway)) {
    result <- getFisher(pathway[[i]])
    pval[i] <- result$pvalue
    matchedGeneList[[i]] <- result$matchedGene
    num1[i] <- result$matchNum
    num2[i,] <- result$fisherTable
  }
  names(pval) <- names(pathway)
  qval <- p.adjust(pval, "BH")
  
  sigIndex<-sort(pval, decreasing=FALSE, index.return=TRUE, na.last=NA)$ix  
  if(sort==TRUE){
    sigIndexTop <- sigIndex[1:topNum]
  }else{
    sigIndexTop <- 1:length(pathway)
  }
  
  matchedGeneFraction <- array(0, dim=topNum)
  matchedGene <- array(0, dim=topNum)
  for (i in 1:topNum) {  
    matchedGene[i] <- paste(matchedGeneList[[sigIndexTop[i]]], collapse="/")    
  }
  
  summary<-data.frame(pvalue = pval[sigIndexTop], 
    qvalue = qval[sigIndexTop], 
    matchedGene = matchedGene
  )
  a <- format(summary, digits=3)      
  return(a)
}

######################
# Network functions
######################
makeSymm <- function(a){
  ## make upper tri to low tri
  ind <- lower.tri(a) 
  a[ind] <- t(a)[ind] 
  return(a) 
}
## generate the studies
edgeConserve <- function(x){
  S1 <- 5
  S2 <- 9
  if(x == 1){
    p <- rbeta(1, S1, 1)
    t <- rbern(1, p)
  }else if(x == 0){
    p <- rbeta(1, 1, S2)
    t <- rbern(1, p)
  }
  return(t)
}
edgeUnconserve <- function(x, pNull){
  t <- rbern(1, pNull)
  return(t)
}

ttest <- function(x,y,...) {
  if(all(is.na(x)) == 1 | all(is.na(y)) == 1) {
    return(NA)
  } else {
    m <- x[which(is.na(x) == 0)]
    n <- y[which(is.na(y) == 0)]
    obj <- try(t.test(m, n,...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) 
    else return(obj$p.value)    
  }
}

printNetworks <- function(data, memberIndex, studyName, nodeName, a, b){
  par(mfrow=c(a,b), oma=c(0,0,0,0), mar=c(0,0,0,0), mgp=c(0,0,0), cex.main=1)
  studyNum <- length(data)
  half1 <- 1:(studyNum/2)
  half2 <- (studyNum/2+1):studyNum
  col <- array("grey",length(memberIndex))
  densityTemp <- lapply(1:length(data), function(x) igraph::graph.adjacency(
    data[[x]][memberIndex, memberIndex], mode="undirected", add.colnames=NA, 
    add.rownames=NA))
  densityAll <- sapply(densityTemp, igraph::graph.density)
  indexDensity <- which.max(abs(densityAll[half1] - densityAll[half2]))
  if (densityAll[half1[indexDensity]] < densityAll[half2[indexDensity]]){
    indexDensity <- half2[indexDensity]
  }
  tempLayout <- igraph::layout.circle(densityTemp[[indexDensity]])
  for (i in 1:(length(data))) {
    gTemp<-igraph::graph.adjacency(data[[i]][memberIndex, memberIndex], 
      mode="undirected", add.colnames=NA, add.rownames=NA)    
    igraph::V(gTemp)$label.cex <- 1
    plot(gTemp, layout=tempLayout, vertex.label=memberIndex, edge.width=2, 
        vertex.color=col, margin=c(-1,0,-1.6,0))
    title(paste(studyName[i],"\n", format(igraph::graph.density(gTemp), digits=3)), line=-3)
  }  
}

getDensity <- function(data, oldSet, group){
  data1 <- lapply(data[group], function(x) x[oldSet,oldSet])
  data1Graph <- lapply(data1, function(x) igraph::graph.adjacency(x,mode="undirected"))
  densityAll <- sapply(data1Graph, igraph::graph.density)
  return(densityAll)
}

##############################
# Simulated annealing function
##############################

## energy function for getting dense differential modules
energyEvaluation <- function(data, memberNew, group1Index, group2Index, 
  w1=200){
  sizeTemp <- length(memberNew)
  data1 <- lapply(data, function(x) x[memberNew, memberNew])
  densityAll <- sapply(data1, function(x) sum(rowSums(x))/(dim(x)[1]*(dim(x)[1]-1)))
  x1 <- sizeTemp  # greater the better
  x2 <- mean(densityAll[group1Index] - densityAll[group2Index])  
  x3 <- 1-sd(densityAll[group1Index] - densityAll[group2Index])  
  energy1 <- exp(-10/30 * (x1))
  energy2 <- exp(-5 * (x2))
  energy3 <- exp(-5 * (x3))
  w2 <- w3 <- (1000-w1)/2
  energy <- (w1*energy1+w2*energy2+w3*energy3)
  energyAll <- c(w1*energy1, w2*energy2, w3*energy3, energy)
  return(energyAll)
}

## add or removing nodes
addNode <- function(x, trialSet) {
  if (length(trialSet) == 0) {
    trialSet <- numeric()
  }
  if(length(trialSet) == 0) {
    temp <- NULL
    temp$x <- x
    temp$trialSet <- numeric()
    return(temp)
  }
  if (length(trialSet) > 1) {
    a <- sample(c(trialSet), size=1)
  }else{
    a <- trialSet
  }
  temp <- NULL
  temp$x <- c(x, a)
  temp$trialSet <- setdiff(trialSet, a)
  return(temp)
}

removeNode <- function(x, trialSet){  
  if (length(x) <= 1) {
    return(x)
  }  
  setNum <- length(x)
  setToChoose <- x
  a <- sample(1:setNum, 1)  
  trialSet <- c(trialSet, setToChoose[a])
  x <- setdiff(x, setToChoose[a])
  temp <- NULL
  temp$x <- x
  temp$trialSet <- trialSet
  return(temp)
}

updateSearchSpace <- function(data, currentState, caseGroupIndex) {
  temp <- Reduce("+", data[caseGroupIndex])
  temp1 <- matrix(0, dim(temp)[1], dim(temp)[2])
  ###at least has connection to 1 of the studies
  temp1[temp >= 2] <- 1
  result <- NULL
  gTemp <- igraph::graph.adjacency(temp1, mode="undirected")
  setTemp <- igraph::neighborhood(gTemp, 1, currentState, "all")
  setTemp1 <- setdiff(unique(unlist(setTemp)), currentState)
  result$set <- setTemp1
  return(result)
}

getJaccard <- function(x, y){
  return(length(intersect(x, y))/length(union(x, y)))
}


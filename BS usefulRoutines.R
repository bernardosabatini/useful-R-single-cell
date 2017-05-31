#Useful routines from Bernardo


#____________________________
# ROUTINES to interact with FindAllMarkers and FindAllICs

# topNGenes - extract top genes for a particular cluster from an object produced by FindAllMarkers
#   markers : the object created by FindAllMarkers 
#   clustNum : the number of the cluster to analyze
#   N : how many genes to return.  Defacults to 20
topNGenes <- function(markers, clustNum, N=20) {
  head(subset(markers, cluster==clustNum), N)
}

# topTGenes - extract top genes from a FindAllMarkers object with avg_diff set by threshold
#   markers : the object created by FindAllMarkers 
#   clustNum : the number of the cluster to analyze
#   thresh : the threshold to use.  Interpreted as requesting avg_diff>=0 if thresh>=0 and avg_diff<thresh otherwise
#   N : limit of number of genes to be return but defaults to NULL which means return all genes
topTGenes <- function(markers, clustNum, thresh=0.1, N=NULL) {
  if (thresh>=0) {
    if (is.null(N))
      subset(markers, (cluster==clustNum)&(avg_diff>thresh))
    else
      head(subset(markers, (cluster==clustNum)&(avg_diff>thresh)), N)
  } else {
    if (is.null(N))
      subset(markers, (cluster==clustNum)&(avg_diff<thresh))
    else
      head(subset(markers, (cluster==clustNum)&(avg_diff<thresh)), N)
  }
}

# cycleAllTop - 
#   sObj : a Seurat object
#   markers : results of FindAllMarkers
#   clusterIds : a list of the cluster numbers to cycle through
#   nGenes : how many genes for each cluster to display
#   allAtOnce : if F, one FeaturePlot per gene.  if T, put them all in one panel.  T only suitable for upto 4 genes or so
#   
#   Cycles through the top nGenes in each of the specified clusters and does the FeaturePlot for each.
#   recommended that you first open a pdf output file (pdf('filename.pdf)), run cycleAllTop, and then close pdf (dev.off())
cycleAllTop <- function(sObj, markers, clusterIds, nGenes, allAtOnce=F) {
  for (counter in clusterIds) {
    print(paste("Cluster", counter))
    tt=topTGenes(markers, counter, 0, nGenes)
    cGenes=tt$gene
    print(paste("Top", nGenes, "enriched genes"))
    print(cGenes)
    
    textplot(paste("Cluster", counter, "Top", nGenes, "enriched genes"),cex=2)
    textplot(cGenes,cex=2)
    
    if (allAtOnce) FeaturePlot(sObj, cGenes)  
    else for (gene in cGenes) FeaturePlot(sObj, gene)
  }
}



#______________________________________________
# ROUTINES to get info on ICs
icComps <- function(sObj, icName, howMany=5) {
  data=sObj@ica.x[,icName]
  absSort=sort.int(abs(data), decreasing = T, index.return = T)
  returnData=data[1:howMany]
  names(returnData) <- rownames(sObj@ica.x)[absSort$ix[1:howMany]]
  return(returnData)
}


#______________________________________________
# ROUTINES to manipulate genes and gene names

# gList - returns gene names selected by grep 
#   sObj : the Seurat object
#   preGene : string with to use for matching in grep
gList <- function(sObj, preGene) rownames(sObj@data)[grep(preGene, rownames(sObj@data))]



#______________________________________________
# ROUTINES to get cell names from clusters

# cCells - return cell names of all cells in a cluster
#   sObj : a Seurat object
#   clustNum : the number of the cluster
cCells <- function(sObj, clustNum) {
  if (length(clustNum)==1) colnames(sObj@data)[which(sObj@ident %in% clustNum)]
  else colnames(sObj@data)[which(sObj@ident %in% clustNum)]
}

# randomCells - return the names of a random subsets of the cells an Seurat object
randomCells <- function(sObj, fract=0.1) {
  nCells=dim(sObj@data)[2]
  colnames(sObj@data)[sample(1:nCells, floor(nCells*fract), replace = F)]
}


#__________________________________________
# ROUTINES to display data

# smallTSNEPlot - calls TSNEPlot for random subset of cells
#   sObj : the Seurat object
#   fract : optional argument to specify the fraction of cells to use.  Default 0.1
#   ... : anything you want to pass to TNSEPlot (like do.label=T)
smallTSNEPlot <- function(sObj, fract=0.1, ...) TSNEPlot(sObj, cells.use = randomCells(sObj, fract), ...)



#__________________________________
# ROUTINES for gene levels and correlations
# probably not ready for general use yet

#undoExp - reverse the normalized log(1+counts) operation
#   a : the log normalized data (sObj@raw.data or sObj@data)
# It is slow!
# useage - rawCounts = undoExp(sObj@data)
undoExp <- function(a) {
  aa = (exp(a)-1)/10000
  aax<-sapply(1:ncol(aa),function(y) {
    txptsum<-1/(sort(unique(aa[,y]))[2])
    txptsum*aa[,y]
  })
  aax<-round(aax)
  colnames(aax)<-colnames(aa)
  rownames(aax)<-rownames(aa)
  return(aax)
}

normGeneData <- function(gData, totalExp=10000) apply(gData, 2, function(x)(totalExp*x/sum(x)))

# mycor - simple wrapper for cor returning whole number percentages
mycor <- function (c1, c2) {
  round(100*cor(c1, c2))
}

# myp - plots data as a function of number of entries
myp <- function(tempD,labS=NULL) {
  if (is.null(labS)) labS=deparse(substitute(tempD))
  plot(1:length(tempD), tempD, xlab=labS)
}


# myps - plots sorted data as a function of number of entries
myps <- function(tempD) {myp(sort(tempD), deparse(substitute(tempD)))}

# cg - caluculate gene correlations
cg <- function(eData, genes, cells=NULL) {
  if (is.null(cells)) {
    if (length(genes)==1) {
      mycor(eData[genes, ], t(eData))
    }
    else {
      mycor(colMeans(eData[genes, ]), t(eData))
    }
  }
  else {
    if (length(genes)==1) {
      mycor(eData[genes, cells], t(eData[,cells]))
    }
    else {
      mycor(colMeans(eData[genes, cells]), t(eData[,cells]))
    }
  }
}


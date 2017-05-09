#Useful routines from Bernardo

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


# randomCells - return the names of a random subsets of the cells an Seurat object
randomCells <- function(sObj, fract=0.1) {
  nCells=dim(sObj@data)[2]
  colnames(sObj@data)[sample(1:nCells, floor(nCells*fract), replace = F)]
}

# smallTSNEPlot - calls TSNEPlot for random subset of cells
#   sObj : the Seurat object
#   fract : optional argument to specify the fraction of cells to use.  Default 0.1
#   ... : anything you want to pass to TNSEPlot (like do.label=T)
smallTSNEPlot <- function(sObj, fract=0.1, ...) TSNEPlot(sObj, cells.use = randomCells(sObj, fract), ...)

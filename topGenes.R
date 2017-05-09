topNGenes <- function(markers, clustNum, N=20) {
  head(subset(markers, cluster==clustNum), N)
}

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

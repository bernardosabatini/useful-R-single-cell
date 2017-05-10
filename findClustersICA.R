FindClustersICA <- function(object, genes.use = NULL, ic.use = NULL, pc.use = NULL, k.param = 30,
                            k.scale = 25, plot.SNN = FALSE, prune.SNN = 1/15,
                            save.SNN = FALSE, reuse.SNN = FALSE, do.sparse = FALSE, 
                            modularity.fxn = 1, resolution = 0.8, algorithm = 1,
                            n.start = 100, n.iter = 10, random.seed = 0, print.output = TRUE,
                            temp.file.location = NULL){
  
  # for older objects without the snn.k slot
  if(typeof(validObject(object, test = T)) == "character"){
    object@snn.k <- numeric()
  }
  
  snn.built <- FALSE
  if (.hasSlot(object, "snn.dense")) {
    if (length(object@snn.dense) > 1) {
      snn.built <- TRUE
    }
  }
  if (.hasSlot(object, "snn.sparse")) {
    if (length(object@snn.sparse) > 1) {
      snn.built <- TRUE
    }
  }
  
  if((missing(genes.use) && missing(ic.use) && missing(pc.use) && missing(k.param) && missing(k.scale) && 
      missing(prune.SNN) && snn.built) || reuse.SNN){
    save.SNN <- TRUE
    if (reuse.SNN && !snn.built){
      stop("No SNN stored to reuse.")
    }
    if (reuse.SNN && (!missing(genes.use) || !missing(pc.use) || !missing(k.param) || 
                      !missing(k.scale) || !missing(prune.SNN))){
      warning("SNN was not be rebuilt with new parameters. Continued with stored SNN. To suppress this
                      warning, remove all SNN building parameters.")
    }
  }
  # if any SNN building parameters are provided or it hasn't been built, build a new SNN
  else{
    print("calling build")
    object <- BuildSNNICA(object, genes.use, ic.use, pc.use, k.param, k.scale,
                          plot.SNN, prune.SNN, do.sparse, print.output)
  }
  print("passed line 42")
  # deal with sparse SNNs
  if (length(object@snn.sparse) > 1) {
    SNN.use <- object@snn.sparse
  } else {
    SNN.use <- object@snn.dense
  }
  for (r in resolution) {
    object <- RunModularityClustering(object, SNN.use, modularity.fxn, r,
                                      algorithm, n.start, n.iter, random.seed,
                                      print.output, temp.file.location)
    object <- GroupSingletons(object, SNN.use)
    name <- paste("res.", r, sep = "")
    object <- StashIdent(object, name)
  }
  
  if (!save.SNN) {
    object@snn.sparse <- sparseMatrix(1, 1, x = 1)
    object@snn.dense <- matrix()
    object@snn.k <- integer()
  }
  return(object)
}

RunModularityClustering <- function(object, SNN = matrix(), modularity = 1,
                                    resolution = 0.8, algorithm = 1,
                                    n.start = 100, n.iter = 10, random.seed = 0,
                                    print.output = TRUE, temp.file.location = NULL){
  
  seurat.dir <- system.file(package="Seurat")
  
  ModularityJarFile <- paste(seurat.dir,
                             "/java/ModularityOptimizer.jar", sep = "")
  
  seurat.dir <- paste0(strsplit(seurat.dir, "/")[[1]][0:(length(strsplit(seurat.dir, "/")[[1]])-1)], collapse = "/")
  seurat.dir <- paste0(seurat.dir, "/", sep = "")
  
  diag(SNN) <- 0
  if (is.object(SNN)) {
    SNN <- as(SNN, "dgTMatrix")
    edge <- cbind(i = SNN@j, j = SNN@i, x = SNN@x)
  } else {
    swap <- which(SNN != 0, arr.ind = TRUE) - 1
    temp <- swap[, 1]
    swap[, 1] <- swap[, 2]
    swap[, 2] <- temp
    edge <- cbind(swap, SNN[which(SNN != 0, arr.ind = TRUE)])
  }
  rownames(edge) <- NULL
  colnames(edge) <- NULL
  
  edge <- edge[!duplicated(edge[, 1:2]), ]
  
#  temp.file.location <- set.ifnull(temp.file.location, seurat.dir)
  print("FindClusterICA lines 95 commmented out")
  
  unique_ID <- sample(10000 : 99999, 1)
  edge_file <- paste(temp.file.location, "edge_", unique_ID, ".txt", sep = "")
  output_file <- paste(temp.file.location, "output_", unique_ID, ".txt", sep = "")
  while (file.exists(edge_file)) {
    unique_ID <- sample(10000 : 99999, 1)
    edge_file <- paste(temp.file.location, "edge_", unique_ID, ".txt", sep = "")
    output_file <- paste(temp.file.location, "output", unique_ID, ".txt", sep = "")
  }
  if (print.output) {
    print.output <- 1
  }
  else {
    print.output <- 0
  }
  
  write.table(x = edge, file = edge_file, sep = "\t", row.names = FALSE,
              col.names = FALSE)
  if (modularity == 2 && resolution > 1){
    stop("error: resolution<1 for alternative modularity")
  }
  command <- paste("java -jar", shQuote(ModularityJarFile), shQuote(edge_file), shQuote(output_file),
                   modularity, resolution, algorithm, n.start, n.iter,
                   random.seed, print.output, sep = " ")
  system(command, wait = TRUE)
  ident.use <- read.table(file = output_file, header = FALSE, sep = "\t")[, 1]
  
  object <- SetIdent(object, object@cell.names, ident.use)
  file.remove(edge_file)
  file.remove(output_file)
  return (object)
}


GroupSingletons <- function(object, SNN){
  # identify singletons
  singletons <- c()
  for (cluster in unique(object@ident)) {
    if (length(WhichCells(object, cluster)) == 1) {
      singletons <- append(singletons, cluster)
    }
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- unique(object@ident)
  cluster_names <- setdiff(cluster_names, singletons)
  connectivity <- vector(mode="numeric", length = length(cluster_names))
  names(connectivity) <- cluster_names
  for (i in singletons) {
    for (j in cluster_names) {
      subSNN = SNN[WhichCells(object, i), match(WhichCells(object, j), colnames(SNN))]
      if (is.object(subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
      } else {
        connectivity[j] <- mean(subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(names(connectivity[mi]), 1)
    object <- SetIdent(object, cells.use = WhichCells(object,i), 
                       ident.use = closest_cluster)
    
  }
  if (length(singletons) > 0){
    print(paste(length(singletons), "singletons identified.", length(unique(object@ident)), "final clusters."))
  }
  return(object)
}
NormalizeData <- function(object, total.expr=1e4, use.raw=0) {
  genes.use=rownames(object@data)
  cells.use=colnames(object@data)
  
  bin.size <- 1000
  max.bin <- floor(length(cells.use)/bin.size) + 1
  cat("Scaling data ", file = stderr())
  pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())


  if (use.raw) data.use=object@raw.data[genes.use, cells.use]
  else data.use=object@data
  
  
  for(i in 1:max.bin) {
    my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
    my.inds <- my.inds[my.inds <= length(cells.use)]
    #print(my.inds)
    geneSumByCell <- colSums(data.use[,cells.use[my.inds]])
    
    new.data <- sweep(data.use[,cells.use[my.inds]], 2, geneSumByCell/total.expr, "/")    
    object@data[,cells.use[my.inds]] <- new.data
    setTxtProgressBar(pb, i)  
  } 
  
  return(object)
}
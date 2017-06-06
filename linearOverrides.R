AverageExpression <-  function(object,genes.use=NULL,return.seurat=F,add.ident=NULL,use.linear=F,...) {
  genes.use=set.ifnull(genes.use,rownames(object@data))
  genes.use=unique(ainb(genes.use,rownames(object@data)))
  ident.use=object@ident
  if (!is.null(add.ident)) {
    new.data=FetchData(object,add.ident)
    new.ident=paste(object@ident[rownames(new.data)],new.data[,1],sep="_")
    object=SetIdent(object,rownames(new.data),new.ident)
  }
  data.all=data.frame(row.names = genes.use)
  for(i in levels(object@ident)) {
    temp.cells=WhichCells(object,i)
    if (length(temp.cells)==1) data.temp=(object@data[genes.use,temp.cells])
    if (length(temp.cells)>1) 
      if (use.linear) data.temp=apply(object@data[genes.use,temp.cells],1,mean)
      else data.temp=apply(object@data[genes.use,temp.cells],1,expMean)
      data.all=cbind(data.all,data.temp)
      colnames(data.all)[ncol(data.all)]=i
  }
  colnames(data.all)=levels(object@ident)
  if (return.seurat) {
    toRet=new("seurat",raw.data=data.all)
    toRet=Setup(toRet,project = "Average",min.cells = 0,min.genes = 0,is.expr = 0,...)
    return(toRet)
  }
  return(data.all)
}

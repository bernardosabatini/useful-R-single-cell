ncg <- function(gene) length(which(big20@data[gene,]>0))

ncc <- function(gene) colnames(big20@data)[which(big20@data[gene,]>0)]

gl <- function(preGene, gGenes) {gGenes[grep(preGene, gGenes)]}

glo <- function(preGene, sObj) {
  gGenes=rownames(sObj@data)
  gGenes[grep(preGene, gGenes)]
  }


cCells <- function(clust, ddd) {
  if (length(clust)==1) colnames(ddd@data)[which(ddd@ident %in% clust)]
  else colnames(ddd@data)[which(ddd@ident %in% clust)]
  }



mycor <- function (c1, c2) {
  round(100*cor(c1, c2))
}


myp <- function(tempD,labS=NULL) {
  if (is.null(labS)) labS=deparse(substitute(tempD))
  plot(1:length(tempD), tempD, xlab=labS)
}

myps <- function(tempD) {myp(sort(tempD), deparse(substitute(tempD)))}


cg <- function(genes, cells=NULL, objSnn) {
  if (is.null(cells)) {
    if (length(genes)==1) {
      mycor(objSnn[genes, ], t(objSnn))
    }
    else {
      mycor(colMeans(objSnn[genes, ]), t(objSnn))
    }
  }
  else {
    if (length(genes)==1) {
      mycor(objSnn[genes, cells], t(objSnn[,cells]))
    }
    else {
      mycor(colMeans(objSnn[genes, cells]), t(objSnn[,cells]))
    }
  }
}

cgB <- function(genes, cells=NULL) {
  if (is.null(cells)) {
    if (length(genes)==1) {
      mycor(bigDNLscale[genes, ], t(bigDNLscale))
    }
    else {
      mycor(colMeans(bigDNLscale[genes, ]), t(bigDNLscale))
    }
  }
  else {
    if (length(genes)==1) {
      mycor(bigDNLscale[genes, cells], t(bigDNLscale[,cells]))
    }
    else {
      mycor(colMeans(bigDNLscale[genes, cells]), t(bigDNLscale[,cells]))
    }
  }
}

cgN <- function(genes, cells=NULL) {
  if (is.null(cells)) {
    if (length(genes)==1) {
      mycor(bigNei[genes, ], t(bigNei))
    }
    else {
      mycor(colMeans(bigNei[genes, ]), t(bigNei))
    }
  }
  else {
    if (length(genes)==1) {
      mycor(bigNei[genes, cells], t(bigNei[,cells]))
    }
    else {
      mycor(colMeans(bigNei[genes, cells]), t(bigNei[,cells]))
    }
  }
}  

cgS <- function(genes, cells=NULL) {
  if (is.null(cells)) {
    if (length(genes)==1) {
      mycor(bigSnn[genes, ], t(bigSnn))
    }
    else {
      mycor(colMeans(bigSnn[genes, ]), t(bigSnn))
    }
  }
  else {
    if (length(genes)==1) {
      mycor(bigSnn[genes, cells], t(bigSnn[,cells]))
    }
    else {
      mycor(colMeans(bigSnn[genes, cells]), t(bigSnn[,cells]))
    }
  }
}  



mt <- function (cp, th) {
  if (th>0) {
    colnames(cp)[which(cp>th)]
  }
  else {
    colnames(cp)[which(cp<th)]
  }
}

mtv2 <- function (cp, th) {
  if (th>0) {
    cp[, sort(colnames(cp)[which(cp>th)])]
  }
  else {
    cp[, sort(colnames(cp)[which(cp<th)])]
  }
}

mtv <- function (cp, th) {
  if (th>0) {
    sort(cp[, colnames(cp)[which(cp>th)]])
  }
  else {
    sort(cp[, colnames(cp)[which(cp<th)]])
  }
}

scG <- function (genes) {
  for (gene in genes) {
       print(gene)
       print(resData[gene,1])
       }
}

topC <- function(cp, N=20) {
  cp[which(is.na(cp))] <- 0
  print(length(cp))
  nEnt = length(cp)
  tts = sort(cp[,colnames(cp)])
  if (N < 0)
    tts[1:-N]
  else
    tts[(nEnt - N):nEnt]
}


undoExp = function(a) {
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

cellClass=FetchData(big20, vars.all = "ident")
big20cells=colnames(big20@data)
resData=big20@data
for (counter in 0:20) {
  print(counter)
  cellIndex=which(cellClass==counter)
  resData[,cellIndex]=(exp(big20@data[,cellIndex])-1)-(exp(bigAvgExp[,as.character(counter)])-1)
}

bigAvgExpNew=bigAvgExp
for (counter in 0:20) {
  print(counter)
  cellIndex=which(cellClass==counter)
  bigAvgExpNew[,as.character(counter)]=rowMeans(bigDNLscale[,cellIndex])
}

# residuals but try first normalizing total expression per cell
bigD=big20@data[goodGenes,]
#bigDNL=exp(bigD-1)
bigDNL=undoExp(bigD)
bigDNL[is.na(bigDNL)] <- 0
bigDNLs=colSums(bigDNL, na.rm=T)
bigDNLscale=scale(bigDNL, center=F, scale=bigDNLs)

bigAvgExp=bigAvgExp[goodGenes,]
bigAvgExpN=exp(bigAvgExp-1)
bigAvgExpN[is.na(bigAvgExpNL)] <- 0
bigAEs=colSums(bigAvgExpN, na.rm=T)
bigAvgExpNL=scale(bigAvgExpN, center=F, scale=bigAEs)

cellClass=FetchData(big20, vars.all = "ident")
big20cells=colnames(big20DNL)

resData=bigDNLscale

for (counter in 0:20) {
  print(counter)
  cellIndex=which(cellClass==counter)
#  resData[,cellIndex]=bigDNLscale[,cellIndex] - bigAvgExpNL[,as.character(counter)]
  resData[,cellIndex]=bigDNLscale[,cellIndex] - bigAvgExpNew[,as.character(counter)]
}


neiAvg <- function(cellI) {
  rowMeans(bigDNLscale[,kn$nn.index[cellI,]])
}

snnAvg <- function(cellI, scaleM=bigDNLscale) {
  rowMeans(scaleM[,snnSort[cellI,1:10]])
}

snnMedian <- function(geneI, cellI) {
  median(as.numeric(big20@data[geneI,snnSort[cellI,1:30]]))
}

nNei=30
nCells=15780
allCellNames=rownames(big20@snn.dense)
snnSort=matrix(data=NA, nCells, nNei)
for (counter in 1:15780) {
  s3=sort.int(big20@snn.dense[counter,],decreasing = T)[2:(nNei+1)]
  ii=which(allCellNames %in% names(s3))
  snnSort[counter,]=ii
}

bigSnn=bigDNLscale
for (ind in 1:nCells) bigSnn[,ind]=bigDNLscale[,ind]-snnAvg(ind)

calcScaleM <- function(dataM) {
  # pass in the @data matrix with genes you care about like big20@data[goodGenes,]
  print("start")
  
  dataMRaw=undoExp(dataM)
  dataMRaw[is.na(dataMRaw)] <- 0 # fill in NA
  dataMRawS=colSums(dataMRaw, na.rm=T)
  return(scale(dataMRaw, center=F, scale=dataMRawS))
}

calcSnnRes <- function(sObj) { # like big20
  nNei=30 # use nearest 30 neighbors # try 10 too
  nCells=ncol(sObj@data)
  allCellNames=rownames(sObj@snn.dense)
  snnSort=matrix(data=NA, nCells, nNei)
  print("start")
  for (counter in 1:nCells) {
    if (counter/10==floor(counter/10)) print(counter)
    
    s3=sort.int(sObj@snn.dense[counter,],decreasing = T)[2:(nNei+1)]
    ii=which(allCellNames %in% names(s3))
    snnSort[counter,]=ii
  }
  return(snnSort)
}

finishSnnRes <- function(dScale) {
  dSnn=dScale
  nCells=ncol(dScale)
  for (ind in 1:nCells) dSnn[,ind]=dScale[,ind]-snnAvg(ind)
}

randomCells <- function(sObj, fract=0.1) {
  nCells=dim(sObj@data)[2]
  colnames(sObj@data)[sample(1:nCells, floor(nCells*fract), replace = F)]
}



FindICs <- function(object, ident.1, ident.2 = NULL, ics.use = NULL, sigma.use = 2, 
                        test.use = "bimod", min.pct = 0.1, min.diff.pct = -Inf, print.bar = FALSE,
                        only.pos = FALSE, max.cells.per.ident = Inf, random.seed = 1, 
                        latent.vars = "nUMI", min.cells = 3) {
  
  ics.use=set.ifnull(ics.use, 1:ncol(object@ica.x)) # get the numbers of ICs to use
  icNames.use=paste("IC", ics.use, sep="") # turns the numbers into names
  
  if (max.cells.per.ident < Inf) object=SubsetData(object,max.cells.per.ident = max.cells.per.ident,random.seed = random.seed)
  # in case the user passed in cells instead of idedimntity classes
  if (length(as.vector(ident.1) > 1) && any(as.character(ident.1) %in% object@cell.names)) {
    cells.1=ainb(ident.1,object@cell.names)
  } else {
    cells.1=WhichCells(object,ident.1)
  }
  
  # if NULL for ident.2, use all other cells
  if (length(as.vector(ident.2) > 1) && any(as.character(ident.2) %in% object@cell.names)) {
    cells.2=ainb(ident.2,object@cell.names)
  } else {
    if (is.null(ident.2)) {
      cells.2=object@cell.names
    }
    else {
      cells.2=WhichCells(object, ident = ident.2)
    }
  }
  cells.2=anotinb(cells.2,cells.1)
  
  #error checking
  if (length(cells.1) == 0) {
    print(paste("Cell group 1 is empty - no cells with identity class", ident.1))
    return(NULL)
  }
  if (length(cells.2) == 0) {
    print(paste("Cell group 2 is empty - no cells with identity class", ident.2))
    return(NULL)
  }
  
  if (length(cells.1) < min.cells) {
    print(paste("Cell group 1 has fewer than", as.character(min.cells), "cells in identity class", ident.1))
    return(NULL)
  }
  if (length(cells.2) < min.cells) {
    print(paste("Cell group 2 has fewer than", as.character(min.cells), " cells in identity class", ident.2))
    return(NULL)
  }
  
  #gene selection (based on percent expressed)
  # thresh.min=object@is.expr
  # data.temp1=round(apply(object@ic.rot[genes.use, cells.1, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  # data.temp2=round(apply(object@data[genes.use, cells.2, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  # data.alpha=cbind(data.temp1,data.temp2); colnames(data.alpha)=c("pct.1","pct.2")
  # alpha.min=apply(data.alpha,1,max); names(alpha.min)=rownames(data.alpha); genes.use=names(which(alpha.min>min.pct))
  # alpha.diff=alpha.min-apply(data.alpha,1,min); 
  # genes.use=names(which(alpha.min>min.pct&alpha.diff>min.diff.pct))
  zData.1=t(scale(t(object@ica.rot[cells.1, icNames.use]), center = T, scale = T));
  zData.2=t(scale(t(object@ica.rot[cells.2, icNames.use]), center = T, scale = T));

    print(dim(zData.1))
  print(dim(zData.2))
  print(mean(zData.1[,"IC45"]))
  print(mean(zData.2[,"IC45"]))
  print(mean(object@ica.rot[cells.1, "IC45"]))
  print(mean(object@ica.rot[cells.2, "IC45"]))
  
  data.temp1=round(apply(zData.1,2,function(x)return(length(which(abs(x)>=sigma.use))/length(x))),3)
  data.temp2=round(apply(zData.2,2,function(x)return(length(which(abs(x)>=sigma.use))/length(x))),3)
  data.alpha=cbind(data.temp1,data.temp2); colnames(data.alpha)=c("pct.1","pct.2")
  
  #IC selection (based on average difference)
  data.1=apply(object@ica.rot[cells.1, icNames.use, drop = F],2,mean)
  data.2=apply(object@ica.rot[cells.2, icNames.use, drop = F],2,mean)
  total.diff=round((data.1-data.2),3)
  
  ics.diff = icNames.use # we can inmplement an IC threshold but unclear why #names(which(abs(total.diff)>thresh.use))
  #  genes.use=ainb(genes.use,genes.diff)
  
  #perform DR
  if (test.use=="bimod") to.return=DiffExpTestIC(object,cells.1,cells.2,icNames.use,print.bar) #ONLY ONE IMPLEMENTED
  if (test.use=="roc") to.return=MarkerTest(object,cells.1,cells.2,genes.use,print.bar)
  if (test.use=="t") to.return=DiffTTest(object,cells.1,cells.2,genes.use,print.bar)
  if (test.use=="tobit") to.return=TobitTest(object,cells.1,cells.2,genes.use,print.bar)
  if (test.use=="negbinom") to.return=NegBinomDETest(object,cells.1,cells.2,genes.use,latent.vars,print.bar, min.cells)
  if (test.use=="poisson") to.return=PoissonDETest(object,cells.1,cells.2,genes.use,latent.vars,print.bar, min.cells)
  
  #return results
  to.return[,"avg_diff"]=total.diff[rownames(to.return)]
  to.return=cbind(to.return,data.alpha[rownames(to.return),])
  if (test.use=="roc") {
    to.return=to.return[order(-to.return$power,-to.return$avg_diff),]
  } else to.return=to.return[order(to.return$p_val,-to.return$avg_diff),]
  
  
  if(only.pos) to.return=subset(to.return,avg_diff>0)
  return(to.return)
}

DiffExpTestIC <- function(object, cells.1,cells.2,icNames.use,print.bar=FALSE) {
            iterate.fxn=lapply; if (print.bar) iterate.fxn=pblapply
            p_val=unlist(iterate.fxn(icNames.use,function(x)diffLRT(as.numeric(object@ica.rot[cells.1,x]),as.numeric(object@ica.rot[cells.2,x]))))
            to.return=data.frame(p_val,row.names = icNames.use)
            return(to.return)
          }

diffLRT = function(x,y,xmin=1) {
  lrtX=bimodLikData(x)
  lrtY=bimodLikData(y)
  lrtZ=bimodLikData(c(x,y))
  lrt_diff=2*(lrtX+lrtY-lrtZ)
  return(pchisq(lrt_diff,3,lower.tail = F))
}

bimodLikData=function(x,xmin=0) {
  x1=x[x<=xmin]
  x2=x[x>xmin]
  xal=minmax(length(x2)/length(x),min=1e-5,max=(1-1e-5))
  likA=length(x1)*log(1-xal)
  mysd=sd(x2)
  if(length(x2)<2) {
    mysd=1
  }
  likB=length(x2)*log(xal)+sum(dnorm(x2,mean(x2),mysd,log=TRUE))
  return(likA+likB)
}

set.ifnull <- function(a,b) if (is.null(a)) return(b) else return(a)

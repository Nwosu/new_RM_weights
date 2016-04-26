PerformEquating <-
function(listcountry, tol, tol2, mincommon, step){
  n = length(listcountry)
  items = colnames(listcountry[[1]]$mat.res)
  m = length(items)
  a = NULL
  b = NULL
  se.b = NULL
  names = NULL
  eqset = NULL
  tab.resp = sapply(1:length(listcountry), function(i) 
    colSums(listcountry[[i]]$XX, na.rm = T))
  colnames(tab.resp) = sapply(1:length(listcountry), function(i) listcountry[[i]]$country)
  v3  = matrix(NA, length(listcountry), 1)
  for(i in 1:length(listcountry)){
    v3[i,] = min(tab.resp[,i])
  }
  rownames(v3) = sapply(1:length(listcountry), function(i) listcountry[[i]]$country)
  sub.names = rownames(v3)[v3<10] 
  warning(paste("There are too few positive answers", sub.names,"."))
  # Case 1: 5-9 affirmative responses
  yes1 = which(v3<10 & v3>4)
  country.yes1 = rownames(v3)[yes1]
  which.yes1 = which.yes2 = list()
  if(length(yes1)!=0){
  for(i in 1:length(yes1)){
    which.yes1[[i]] = which(tab.resp[,yes1[i]]<10)
    if(sum(tab.resp[,yes1[i]]<=4)) warning("YES1 contains some YES2 responses")
  }
  }
  # Case 2: 1-4 affirmative responses
  yes2 = which(v3<=4)
  country.yes2 = rownames(v3)[yes2]
  if(length(yes2)!=0){
  for(i in 1:length(yes2)){
    which.yes2[[i]] = which(tab.resp[,yes2[i]]<10)
    which.yes2.ser = which(tab.resp[,yes2[i]]<=4)
  }
  } else{
    which.yes2=NULL
    which.yes2.ser=NULL
  }
  # In this second case, person parameters have to be estimated using the global standard
  for (i in 1:n){
    a = rbind(a,listcountry[[i]]$a)
    b = rbind(b,listcountry[[i]]$b)
    se.b = rbind(se.b,listcountry[[i]]$se.b)
    names = rbind(names, listcountry[[i]]$country)
    eqset = c(eqset,i)
  }  
  rownames(a) = names
  rownames(b) = names
  rownames(se.b) = names
  colnames(b) = items
  colnames(se.b) = items
  b.equated = b
  #scaling item parameters to have sd = 1 
  sd.b = apply(b, 1, sd)*sqrt((m-1)/m) 
  b.scaled = NULL
  for (i in 1:m){
    b.scaled = cbind(b.scaled, b[,i]/sd.b)
  }
  #initializing equating cycle
  # determine which items are common to the 'global' by only considering the equalization set of countries
  b.tot = apply(b.scaled[eqset,], 2, median)
  diff.mat.scaled = NULL
  for (i in eqset) {
    diff.mat.scaled = rbind(diff.mat.scaled, b.scaled[i,]-b.tot)
  }  
  common=matrix(0,nrow(b),ncol(b))
  common[abs(diff.mat.scaled) < tol2] = 1
  id.comm = matrix(as.logical(common),nrow(common))
  rownames(id.comm) = names
  colnames(id.comm) = items
  if(length(yes1)!=0){
  for(i in 1:length(yes1)){
    id.comm[yes1[i],which.yes1[[i]]] = FALSE
  }
  }
  if(length(yes2)!=0){
  for(i in 1:length(yes2)){
    id.comm[yes2[i],which.yes2[[i]]] = FALSE
  }
  }
  # define adjustment parameters needed to rescale on mean and sd of common items only
  aa = NULL
  bb = NULL
  # define parameters to scale initial estimates of a and b
  aaa = NULL
  bbb = NULL
  discr = tol+1
  while (discr > tol) {
    for(i in 1:n) {
      bb[i] = sd(b.tot[id.comm[i,]])/sd(b.scaled[i,id.comm[i,]])
      aa[i] = mean(b.tot[id.comm[i,]])-mean(b.scaled[i,id.comm[i,]])*bb[i]
    }
    #new scaled matrix of all item parameters
    for(i in 1:m){
      b.scaled[,i] = aa + b.scaled[,i]*bb
    }
    # redetermining the totals based on new scaled values
    # only of the common ones in the equalization set
    b.common = b.scaled[eqset,]
    b.common[id.comm[eqset,]==FALSE] = NA
    b.tot = NULL
    b.tot = apply(b.common, 2, function(x) {median(x, na.rm = TRUE)})
    # re-normalizing the total to have mean zero and unit standard deviation
    b.tot = b.tot - mean(b.tot)
    b.tot = b.tot/(sd(b.tot)*sqrt((m-1)/m))
    # redetermining the matrix of distances
    for (i in 1:n) {
      diff.mat.scaled[i,] = b.scaled[i,]-b.tot
    }  
    olddiscr = discr
    discr = max(tol, max(abs(diff.mat.scaled[id.comm])))
    number = apply(id.comm[eqset,],1,sum)
    if (min(number) < mincommon) {
      warning(paste("There are too few items left in", names[number< mincommon],
                  ". Therefore, taking it out of the equalization set"))
      eqset = eqset[number > mincommon]
    }    
    else if (min(number) >= mincommon) {
      common=matrix(0,nrow(b),ncol(b))
      common[abs(diff.mat.scaled) < tol2] = 1
      new.id.comm = matrix(as.logical(common),nrow(b),ncol(b))
      tol2 = min(tol2 - step, discr)
      id.comm = new.id.comm
      if(length(yes1)!=0){
        for(i in 1:length(yes1)){
          id.comm[yes1[i],which.yes1[[i]]] = FALSE
        }
      }
      if(length(yes2)!=0){
        for(i in 1:length(yes2)){
          id.comm[yes2[i],which.yes2[[i]]] = FALSE
        }
      }
      diff = discr
    }
    #recalculate aaa and bbb as we have redifend b.tot   
    for (i in 1:n)  {
      bbb[i] = sd(b.tot[id.comm[i,]])/sd(b[i,id.comm[i,]])
      aaa[i] = mean(b.tot[id.comm[i,]])-mean(b[i,id.comm[i,]])*bbb[i]
    }
  }
  colnames(b.scaled) = items
#   write.table(rbind(b.scaled,b.tot), "Standardized estimates.txt")
  rownames(id.comm) = names
  colnames(id.comm) = items
  which.yes2.notser = which.yes2
  which.yes2 = which.yes2.ser
  return(list(aa = aa, bb = bb, aaa=aaa, bbb=bbb, id.comm = id.comm, 
              b = b, b.scaled = b.scaled, b.tot = b.tot, diff.mat = diff.mat.scaled, 
              items = items, names = names, diff = diff, eqset = eqset, number = number,
              yes1 = yes1, yes2 = yes2, which.yes1 = which.yes1, which.yes2 = which.yes2,
              which.yes22 = which.yes2.notser))
}

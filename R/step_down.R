.UnweightedStepdownControl<-function(p.value, n.discovery, q.spend, tau = 0.3){ ### Proposed FDR method
  len.p.value = length(p.value)
  alpha.i = (n.discovery  + (1:len.p.value))/(len.p.value+1-(1:len.p.value))*q.spend
  alpha.i = pmin(alpha.i/(1+alpha.i),tau)   # thresholding cutoff values by tau
  index.x = which(sort(p.value)/alpha.i > 1)
  if(length(index.x) > 0){
    cut.off = alpha.i[min(index.x)]
    output = (p.value <= cut.off)
  }else{
    cut.off = tau
    output = rep(TRUE,len.p.value)
  }
  return(list(cut.off = cut.off, output = output))
}



.WeightedStepdownControl<-function(p.value, weight, n.discovery, q.spend, tau = 0.3){ ### Proposed FDR method
  len.p.value = length(p.value)
  alpha.i = array(NA, dim = len.p.value)
  for(i in 1:len.p.value){
    alpha.i[i] = (n.discovery + sum(weight[1:i]))/(sum(weight[i:len.p.value]))*q.spend
  }
  alpha.i = pmin(alpha.i/(1+alpha.i),tau)   # thresholding cutoff values by tau
  index.x = which(sort(p.value)/alpha.i > 1)
  if(length(index.x) > 0){
    cut.off = alpha.i[min(index.x)]
    output = (p.value <= cut.off)
  }else{
    cut.off = tau
    output = rep(TRUE,len.p.value)
  }
  return(list(cut.off = cut.off, output = output))
}



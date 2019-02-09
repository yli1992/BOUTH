setClass("Tree", representation = list(node = "data.frame",
                        edge = "list", test = "list", info = "list"))


setMethod("initialize", "Tree", function(.Object, ...) {
  value <- callNextMethod()
  value
})


.build.tree <-function(anno.table, na.symbol = 'unknown'){
  if(is.null(anno.table)) stop("Incorrect number of parameters.")
  n.rank = ncol(anno.table)
  num = id = xparent = parent=level= NULL
  count = 0

  #check if there is any ambiguous symbol (identical symbol for different nodes)
  tot.nodes = NULL
  for(i in n.rank:1){
    tot.nodes = c(tot.nodes, setdiff(unique(anno.table[,i]),na.symbol))
  }
  if(length(tot.nodes)!=length(unique(tot.nodes))){
    stop("At least two nodes from different levels have the same name.")
  }

  #Build a "Tree" object
  for(i in n.rank:1){
    this.level.node = setdiff(unique(anno.table[,i]),na.symbol)
    this.level.len = length(this.level.node)
    if(this.level.len > 0){
      id = c(id, this.level.node)
      num = c(num, c((count+1): (count+this.level.len)))
      count = count + this.level.len
      level = c(level, rep(n.rank+1-i, this.level.len))
      for(j in 1:this.level.len){
        t = which(anno.table[,i] == this.level.node[j])
        k = i
        while(k > 1){
          ts = unique(anno.table[t,k-1])
          if(length(ts) > 1){
            stop("This is not a tree.")
          }else{
            if(ts != na.symbol){
              break
            }else{
              k = k-1
            }
          }
        }
        if(i > 1){
          xparent = c(xparent, ts)
        }else{
          xparent = c(xparent, "head")
        }
      }
    }
  }
  names(num) = id
  for(i in 1:count-1){
    parent[i] = as.double(num[xparent[i]])
  }
  parent = c(parent, count+1)


  leaf.inner.cutoff = min(parent)
  edge.list =  vector("list", count)
  for(i in 1:count){
    if(i < leaf.inner.cutoff){
      edge.list[[i]] = integer(0)
    }else{
      edge.list[[i]] = which(parent==i)
    }
  }

  names(num) = NULL
  dec.level = rep(NA, length(level))

  tree = new(Class = "Tree", node = data.frame(num,id, parent,
              level,dec.level, stringsAsFactors=FALSE), edge = edge.list,
             test= list(reject=logical(0),pvalue=numeric(0), num= numeric(0)),
             info = list(alpha = numeric(0), nl = numeric(0), ql = numeric(0))
             )
  return(tree)

}



setMethod("summary","Tree", function(object){
  node = slot(object, "node")
  node.id = node$id
  test = slot(object, "test")
  cat("Number of hypotheses:", length(node.id), "\n")
  if(length(test$reject)==0 || length(test$pvalue)==0){
    cat("None hypothesis testing has been done.", "\n")
  }else{


    result.to.print = data.frame(rev(node.id[test$reject]), formatC(rev(test$pvalue[test$reject]),format = "e",digits =2))
    colnames(result.to.print) = c("node","p.value")
    cat("Number of tree discoveries:", sum(test$reject), "\n")
    cat("FDR/FAR controlled by:", object@test$num, "\n")
    n.to.print <- min(nrow(result.to.print), 10)
    print(result.to.print[1:n.to.print,], row.names=FALSE)
    if(n.to.print < nrow(result.to.print)) {
      cat('[only 10 top-level significant hypotheses shown]', '\n')
    }
  }
})



setGeneric(".sig", function(object) standardGeneric(".sig"))


setMethod(".sig","Tree", function(object){
  node = slot(object, "node")
  test = slot(object, "test")

  if(length(test$reject)==0){
    return("None hypothesis testing has been done")
  }
  return(test$reject)
})




setGeneric(".top", function(object) standardGeneric(".top"))


setMethod(".top","Tree", function(object){
  node = slot(object, "node")
  test = slot(object, "test")

  if(length(test$reject)==0){
    return("None hypothesis testing has been done")
  }
  nn = length(node$id)
  discoveries = which(test$reject)
  nd = length(discoveries)
  p.values = signif(test$pvalue[test$reject],3)
  top.discoveries = array(TRUE, dim = nd)

  for(i in 1:nd){
    t = discoveries[i]
    while(t < nn){
      t = node$parent[t]
      if(test$reject[t]){
        top.discoveries[i] = FALSE
        break
      }
    }
  }
  result = array(FALSE, nn)
  result[discoveries[top.discoveries]]=TRUE
  return(result)
})




setGeneric(".query", function(object, x, na.symbol) standardGeneric(".query"))


setMethod(".query","Tree", function(object, x, na.symbol){
  node = slot(object, "node")
  test = slot(object, "test")

  if(length(test$reject)==0){
    return("None hypothesis testing has been done")
  }
  nn = length(node$id)
  nq = length(x); curr = NULL;
  name.vec = array(NA, nq)
  name.mat = matrix(NA, nrow = nq, ncol = max(node$level))
  for(i in 1:nq){
    curr = x[i]
    Path.v = Path.m= NULL
    while(curr <= nn){
      if(is.null(Path.v)){
        Path.v = node$id[curr]
        Path.m = c(node$id[curr], rep(" ", node$level[curr]-1))
      }else{
        Path.v = paste(node$id[curr], Path.v, sep=',')
        Path.m = c(node$id[curr], Path.m)
      }
      curr.next = node$parent[curr]
      if(curr < nn){
        level.curr = node$level[curr]; level.next =node$level[curr.next]
        if(level.next > level.curr + 1){
          Path.v = paste(paste0(rep(na.symbol,level.next-level.curr-1),collapse=','), Path.v, sep = ',')
          Path.m = c(rep(na.symbol, level.next-level.curr-1), Path.m)
        }
      }
      curr = curr.next
    }

    name.vec[i] =  Path.v
    name.mat[i,] = Path.m

  }
  return(list(name.vec = name.vec, name.mat = name.mat))
})


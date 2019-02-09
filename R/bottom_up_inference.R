.StoufferCombination<-function(pvalue){
  nan.pvalue = pvalue[!is.na(pvalue)]
  n = sum(nan.pvalue<=1)
  if(n>0){
    z.stat = sum(qnorm(1-nan.pvalue))/sqrt(n)
    return(1-pnorm(z.stat))
  }else{
    return(0)
  }
}


.bottom.up.unweighted<-function(tree, pvalue, far = 0.10, tau = 0.3){## pvalue is a vector
  ## if a node A is detected later than another node B, A is less significant than B.
  CompareTwoNodes<-function(a,b){
    a = as.vector(a)
    b = as.vector(b)
    if(a[1] > b[1]){return(TRUE)
    }
    if(a[1] < b[1]){return(FALSE)
    }
    if(a[1] == b[1]){
      if(a[2] > b[2]){return(TRUE)
      }else{return(FALSE)
      }
    }
  }


  FindLeastSigNode<-function(x){ #x is an n*2 matrix
    pos.min = 1
    n = nrow(x)
    if(n>1){
      for(i in 2:n){
        if(CompareTwoNodes(x[i,],x[pos.min,])){pos.min = i
        }
      }
    }
    return(pos.min)
  }


  n.node = nrow(tree@node)
  max.depth = max(tree@node$level)
  if(length(pvalue)!=sum(tree@node$level==1)){
    stop("Incorrect dimension of p-values!")
  }
  if(!is.null(names(pvalue))){
    level.one.id = tree@node$id[tree@node$level==1]
    if(setequal(level.one.id, names(pvalue))){
      pvalue = pvalue[level.one.id]
    }else{
      stop("Mis-matched symbols between p-values and tree")
    }
  }

  trim.stat = matrix(NA, nrow = n.node, ncol = 2) # when a node is detected, keep the record
  cum.rej = 0; last.cut.off = 0
  sig.in = array(TRUE,n.node); sig.out = array(FALSE,n.node)
  pvalue.in = pvalue.out = dec.level = array(NA, n.node); pvalue.in[1:length(pvalue)]= pvalue
  q.spent = array(0, max.depth)
  alpha.out = array(NA, max.depth)
  nl.out = array(NA, max.depth)

  for(i in 1:max.depth){

    this.level.set = which(tree@node$level == i)
    nl.out[i] = length(this.level.set)
    q.this.level = far/n.node*nl.out[i] ### q-spending by proportion
    q.spent[i] = q.this.level
    this.level.tree = tree@edge[this.level.set]


    # removing significant children from previous levels
    if(i>1){
      pre.sig.nodes = which(sig.in[1:(min(this.level.set)-1)])
      this.level.tree = lapply(this.level.tree,function(x)return(setdiff(x,pre.sig.nodes)))
      indices = (sapply(this.level.tree,length) > 0)
      ## swap with its least significant child node
      na.nodes = this.level.set[!indices]
      for(j in na.nodes){
        j.child = tree@edge[[j]]
        least.sig.child = j.child[FindLeastSigNode(trim.stat[j.child, ,drop=FALSE])]
        sig.out[least.sig.child] = FALSE
        dec.level[j] = dec.level[least.sig.child]
        dec.level[least.sig.child] = NA
        sig.out[j] = TRUE
        pvalue.out[j] = pvalue.out[least.sig.child]
        trim.stat[j,] = trim.stat[least.sig.child, ]
      }

      ## finish the trimming step
      this.level.set = this.level.set[indices]
      this.level.tree = this.level.tree[indices]
    }


    if(length(this.level.set)>0){
      if(i>1){
        this.level.p = sapply(this.level.tree,function(x)return(.StoufferCombination(pvalue.in[x])))
      }else{
        this.level.p = pvalue
      }
      this.level.test.summary = .UnweightedStepdownControl(this.level.p, cum.rej, q.this.level, tau)
      this.level.detect = this.level.test.summary$output
      sig.in[this.level.set] = sig.out[this.level.set] = this.level.detect
      pvalue.out[this.level.set] = this.level.p
      dec.level[this.level.set] = i
      trim.stat[this.level.set[this.level.detect],1] = i #detected at the i-th level
      trim.stat[this.level.set,2] =this.level.p          #p value when it is detected
      cum.rej = cum.rej + sum(this.level.detect)
      last.cut.off = this.level.test.summary$cut.off
      alpha.out[i] = last.cut.off

      if(sum(!this.level.detect)>0){
        pvalue.in[this.level.set[!this.level.detect]]=
          (this.level.p[!this.level.detect]-last.cut.off)/(1-last.cut.off)
      }

    }
  }
  tree@info$alpha = alpha.out
  tree@info$nl = nl.out
  tree@info$ql = q.spent
  tree@test$reject = sig.out; tree@test$pvalue = pvalue.out; tree@node$dec.level = dec.level
  tree@test$num = far
  return(tree)
}






.bottom.up.weighted<-function(tree, pvalue, far = 0.10, tau = 0.3, mode = 'one.stage'){

  get_weight<-function(hash.id, tree){## calculate the least favorable weights

    n.node = nrow(tree@node)
    weight = array(1, dim = n.node)
    hash.xid = hash.id
    subtree<-vector("list",n.node)
    this = hash.id[1]

    while(this <= n.node){
      parent = tree@node$parent[this]
      if(!parent %in% hash.id){
        hash.id =  c(hash.id, parent)
      }
      hash.id = hash.id[-1]
      if(parent <= n.node){
        subtree[[parent]] = c(subtree[[parent]], this)
        this = hash.id[1]
      }else{
        break;
      }
    }

    ## compute least favorable weights
    LFW<-vector("list",n.node); LFW[hash.xid] = 1
    for(i in 1:n.node){
      if(length(subtree[[i]]) > 0){
        LFW[[i]] = sort(unlist(LFW[subtree[[i]]]), decreasing = TRUE)
        LFW[[i]][1] = LFW[[i]][1] + 1
      }
    }
    return(sort(LFW[[n.node]]))
  }


  n.node = nrow(tree@node)
  max.depth = max(tree@node$level)
  if(length(pvalue)!=sum(tree@node$level==1)){
    stop("Incorrect dimension of p-values!")
  }
  if(!is.null(names(pvalue))){
    level.one.id = tree@node$id[tree@node$level==1]
    if(setequal(level.one.id, names(pvalue))){
      pvalue = pvalue[level.one.id]
    }else{
      stop("Mis-matched symbols between p-values and tree")
    }
  }
  q.spent = array(0, max.depth)
  alpha.out = array(NA, max.depth)
  nl.out = array(NA, max.depth)

  if(mode == 'one.stage'){
    for(i in 1:max.depth){
      nl.out[i] = sum(tree@node$level == i)
      q.spent[i] = far/n.node*nl.out[i]
    }
  }
  if(mode == 'stageI.in.two.stage'){
    q.spent[1] = far
    nl.out = sapply(c(1:max.depth),function(x)return(sum(tree@node$level==x)))
  }
  if(mode == 'stageII.in.two.stage'){
    if(max.depth<2){stop('Two stage analysis requires at least 2 levels.')}
    n.node1 = n.node - length(which(tree@node$level == 1))
    for(i in 2:max.depth){
      q.spent[i] = far/n.node1*length(which(tree@node$level == i))
    }
  }

  cum.rej = 0; last.cut.off = 0
  sig.out = array(FALSE,n.node)
  pvalue.in = pvalue.out = dec.level = array(NA, n.node); pvalue.in[1:length(pvalue)]= pvalue


  for(i in 1:max.depth){
    this.level.set = which(tree@node$level == i)
    q.this.level = q.spent[i] ### by default q-spending by proportion
    this.level.tree = tree@edge[this.level.set]

    # removing significant children from previous levels
    if(i>1){

      ## nodes that are significant before testing this level
      pre.sig.nodes = which(sig.out[1:(min(this.level.set)-1)])
      ## nodes to be tested now
      this.level.tree = lapply(this.level.tree,function(x)return(setdiff(x,pre.sig.nodes)))
      ind = (sapply(this.level.tree,length) > 0)
      above.level.set = which(tree@node$level >= i)
      xind = (sapply(lapply(tree@edge[above.level.set],function(x)return(setdiff(x,pre.sig.nodes))),length)>0)
      xind = (xind | sig.out[above.level.set])
      na.nodes = above.level.set[!xind]
      sig.out[na.nodes] = TRUE
      pvalue.out[na.nodes] = NA


      if(length(na.nodes)>0){
        # find its least significant offspring on the i-1 th level
        last.level.set = intersect(which(tree@node$level == (i-1)),which(sig.out))
        for(j in 1:length(na.nodes)){
          jind =  (tree@node$parent[last.level.set]==na.nodes[j])
          pvalue.out[na.nodes[j]] = max(pvalue.out[last.level.set[jind]])
        }
        dec.level[na.nodes] = i-1
      }
      if(mode == 'stageI.in.two.stage') next;

      cum.rej = cum.rej + sum(!xind)
      this.level.set = this.level.set[ind]
      this.level.tree = this.level.tree[ind]
    }


    if(length(this.level.set)>0){
      if(i>1){
        this.level.p = sapply(this.level.tree,function(x)return(.StoufferCombination(pvalue.in[x])))
      }else{
        this.level.p = pvalue
      }

      weight = get_weight(this.level.set, tree)
      this.level.test.summary = .WeightedStepdownControl(this.level.p, weight, cum.rej, q.this.level, tau)
      this.level.detect = this.level.test.summary$output
      sig.out[this.level.set] = this.level.detect
      pvalue.out[this.level.set] = this.level.p
      dec.level[this.level.set] = i


      cum.rej = cum.rej + sum(this.level.detect)
      last.cut.off = this.level.test.summary$cut.off
      alpha.out[i] = last.cut.off


      if(sum(!this.level.detect)>0){
        pvalue.in[this.level.set[!this.level.detect]]=
          (this.level.p[!this.level.detect]-last.cut.off)/(1-last.cut.off)
      }
    }

  }

  tree@info$alpha = alpha.out
  tree@info$nl = nl.out
  tree@info$ql = q.spent

  tree@test$reject = sig.out; tree@test$pvalue = pvalue.out; tree@node$dec.level = dec.level
  if(mode=='stageI.in.two.stage'){#1001
    tree@test$num = last.cut.off
  }else{
    tree@test$num = far
  }
  return(tree)
}


#' Bottom-up Tree Hypothesis Tests
#'
#' This function implements the (one-stage and two-stage) bottom-up approach to testing hypotheses that 
#' have a branching tree dependence structure, with false discovery rate control. 
#' Our motivating example comes from testing the association between a trait of 
#' interest and groups of microbes that have been organized into operational 
#' taxonomic units (OTUs) or amplicon sequence variants (ASVs).  Given p-values 
#' from association tests for each individual OTU or ASV, we would like to know 
#' if we can declare that a certain species, genus, or higher taxonomic grouping 
#' can be considered to be associated with the trait. If a large proportion of 
#' species from that genus influence the trait, we should conclude the genus 
#' influences the trait.  Conversely if only a few of the species from a genus 
#' are non-null, then a better description of the microbes that influence occurrence 
#' of the trait is a list of associated species.  Finding taxa that can be said 
#' to influence a trait in this sense is the first goal of our approach.  
#' The second goal is to locate the highest taxa in the tree for which we 
#' can conclude many taxa below, but not any ancestors above, influence risk; 
#' we refer to such taxa as driver taxa. 
#' @param anno.table an \code{n} by \code{l} data.frame that specifies the
#'   annotation (i.e., grouping) of \code{n} leaf nodes at \code{l} levels. The
#'   \code{l} levels are ordered with the highest level (the root node) first 
#'   and the lowest level (the leaf nodes) last, so that each row represents a
#'   path from the root node to a leaf node.
#' @param pvalue.leaves a vector of p-values at leaf nodes, the names of which 
#'   should match the values in the last column of \code{anno.table}.
#' @param na.symbol a character string which is to be inpretreted as \code{NA}
#'   values. The default is `unknown'.
#' @param far nominal false assignment rate (FAR), the error rate in analogy
#'   with the false discovery rate. If \code{far} takes one value, it
#'   corresponds to the one-stage bottom-up test with overal FAR control at the
#'   specified level. If \code{far} takes a vector of two values, it corresponds
#'   to the two-stage bottom-up test with the test of level-1 nodes controlling
#'   FAR at the level of the first value of \code{far} and the test of the
#'   remaining levels controlling FAR at the level of the second value of
#'   \code{far}. The default is 0.1.
#' @param is.weighted a logical value indicating whether the weighted or
#'   unweighted test is performed. The default is `TRUE'.
#' @param tau a pre-specified constant to prevent nodes with large p-values from
#'   being detected if a large number (say, \code{m}) of null hypotheses can be
#'   easily rejected because of very low p-values, in which case \code{q} x
#'   \code{m} nodes with large p-values can be said to be detected, while still
#'   controlling the overall error rate at level \code{q}. The default is 0.3.
#' @return A list consisting of 
#' \item{results.by.level}{a data frame that gives the number of detected nodes at each level along with information on the test for that level.} 
#'   \item{results.by.node}{a data frame that gives detailed results at each node (i.e., leaf and inner nodes), including the (derived) p-value, the indicator of a detection or not, and the indicator of a driver node or not.} 
#'   \item{tree}{A tree structure used in the function \code{graphlan}.}
#'   
#' @references Li, Y., Satten, GA., Hu Y.-J. "A bottom-up approach to testing hypotheses that have a
#' branching tree dependence structure, with false discovery rate control" XXX(2018)
#' 
#' @examples
#' data(IBD)
#' 
#' ## one-stage weighted bottom-up test on the IBD data
#' test.1 = bouth(anno.table = IBD$tax.table, pvalue.leaves = IBD$pvalue.otus,
#' na.symbol = "unknown", far = 0.1, is.weighted = TRUE)
#'
#' ## extract all detected nodes
#' test.1$results.by.node[test.1$results.by.node$is.detected, ]
#'
#' ## extract all detected driver nodes
#' test.1$results.by.node[test.1$results.by.node$is.driver,]
#'
#'
#' ## two-stage (weighted) bottom-up test on the IBD data
#' test.2 = bouth(anno.table = IBD$tax.table, pvalue.leaves = IBD$pvalue.otus,
#' na.symbol = "unknown", far = c(0.05, 0.05))
#'
#'
#' @export
bouth<-function(anno.table, pvalue.leaves, na.symbol ='unknown', far = 0.10, tau = 0.3, is.weighted = TRUE){
  # this version supports one-stage & two-stage tests

  if(is.null(anno.table) || is.null(pvalue.leaves)) stop("Incorrect number of parameters")
  if(!is.null(colnames(anno.table))){
    level.name = rev(colnames(anno.table))
  }else{
    level.name = rev(paste("level", c(1:ncol(anno.table)), sep='-'))
  }
  if(is.null(names(pvalue.leaves))){
    names(pvalue.leaves) = anno.table[,ncol(anno.table)]
  }
  # build a tree object based on the input table
  tree = .build.tree(anno.table = anno.table, na.symbol = na.symbol)
  if(max(pvalue.leaves)==1){
    warning("There are p-values that equal exactly 1.")
    pvalue.leaves[pvalue.leaves==1] = 0.5*(1+max(c(0,pvalue.leaves[pvalue.leaves<1])))
  }

  # one-stage approach
  if(length(far)==1){
    if(is.weighted){
      results = .bottom.up.weighted(tree, pvalue.leaves, far, tau)
    }else{
      results = .bottom.up.unweighted(tree, pvalue.leaves, far, tau)
    }
  }
  # two-stage approach
  if(length(far)==2){
    results = .two.stage.test(anno.table, pvalue.leaves, na.symbol, far.leaves = far[1], far.inner = far[2],
                              tau)

  }


  options(digits=3)

  ## create a table for all nodes tested in a bottom-up procedure
  n.node = length(results@node$level); max.level = max(results@node$level)
  show.node = .query(results, c(1:n.node), na.symbol); colnames(show.node$name.mat) = rev(level.name)
  show.pvalue = results@test$pvalue
  show.detected = .sig(results)
  show.driver =  .top(results)
  results.by.node = data.frame(show.node$name.mat, show.pvalue, show.detected, show.driver,stringsAsFactors = FALSE)
  colnames(results.by.node) = c(rev(level.name),"pvalue", "is.detected", "is.driver")
  results.by.node = results.by.node[order(show.node$name.vec),]
  rownames(results.by.node) = NULL

  ## create a table to save characteristics for the bottom-up procedure
  detections.per.level = sapply(c(1:max.level),function(x)return(sum(results@node$level[results@test$reject]==x)))
  results.by.level = data.frame(c(1:max.level),level.name,
                         results@info$ql, results@info$alpha,
                         results@info$nl, detections.per.level,
                         stringsAsFactors=FALSE)
  colnames(results.by.level) = c("level","level.name", "q_l", "pvalue.cutoff", "all.nodes", "detected.nodes")


  return(list(tree = results, results.by.node = results.by.node, results.by.level = results.by.level))
}

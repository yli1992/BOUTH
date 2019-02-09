.two.stage.test<-function(anno.table, pvalue,
                         na.symbol = 'unknown',
                         far.leaves = 0.1, far.inner = 0.1, tau = 0.3){
  if(is.null(anno.table) || is.null(pvalue)) stop("Incorrect number of parameters")
  n.leaves = nrow(anno.table)
  n.level = ncol(anno.table)
  if(is.null(names(pvalue))){
    stop('The p-values on leaf nodes should be labeled')
  }

  if(!(n.leaves== length(pvalue)&(setequal(anno.table[,n.level],names(pvalue))))){
    stop('Mismatch between the table and the names of p-values')
  }

  ## the first stage
  table.stage1 = anno.table; colnames(table.stage1)[n.level] = 'leaves'
  tree.stage1 = .build.tree(anno.table = table.stage1, na.symbol = na.symbol)
  test.tree.stage1 = .bottom.up.weighted(tree = tree.stage1, pvalue = pvalue,
                                            far = far.leaves, tau = tau, mode = 'stageI.in.two.stage')
  sig.leaves = test.tree.stage1@node$id[test.tree.stage1@test$reject & (test.tree.stage1@node$level == 1)]
  p.leaves = test.tree.stage1@test$pvalue[test.tree.stage1@test$reject & (test.tree.stage1@node$level == 1)]
  sig.inner.stage1 = test.tree.stage1@node$id[test.tree.stage1@test$reject & (test.tree.stage1@node$level > 1)]
  p.inner.stage1 = test.tree.stage1@test$pvalue[test.tree.stage1@test$reject & (test.tree.stage1@node$level > 1)]


  ## the second stage
  table.stage2 = table.stage1[!table.stage1$leaves%in%sig.leaves, ]
  pvalue.stage2 = (pvalue[!names(pvalue)%in%sig.leaves] -  test.tree.stage1@test$num)/(1- test.tree.stage1@test$num)
  tree.stage2 = .build.tree(anno.table = table.stage2, na.symbol = na.symbol)
  test.tree.stage2 = .bottom.up.weighted(tree = tree.stage2, pvalue = pvalue.stage2,
                                            far= far.inner, tau = tau, mode = 'stageII.in.two.stage')

  # summarize for the output tree
  sig.inner.stage2 = test.tree.stage2@node$id[test.tree.stage2@test$reject]
  p.inner.stage2 = test.tree.stage2@test$pvalue[test.tree.stage2@test$reject]

  reject.out = test.tree.stage1@test$reject; names(reject.out) = test.tree.stage1@node$id
  pvalue.out = test.tree.stage1@test$pvalue; names(pvalue.out) = test.tree.stage1@node$id
  dec.level.out = test.tree.stage1@node$dec.level; names(dec.level.out) = test.tree.stage1@node$id

  reject.out[c(sig.inner.stage1, sig.inner.stage2)] = TRUE; names(reject.out)= NULL
  pvalue.out[sig.inner.stage1] = p.inner.stage1; pvalue.out[test.tree.stage2@node$id] = test.tree.stage2@test$pvalue
  names(pvalue.out) = NULL
  test.tree.stage1@test$reject = reject.out
  test.tree.stage1@test$pvalue = pvalue.out
  test.tree.stage1@test$num = c(far.leaves, far.inner)
  dec.level.out[test.tree.stage2@node$id]=test.tree.stage2@node$dec.level
  names(dec.level.out) = NULL
  test.tree.stage1@node$dec.level = dec.level.out


  test.tree.stage1@info$alpha = c(test.tree.stage1@info$alpha[1], test.tree.stage2@info$alpha[-1])
  test.tree.stage1@info$ql = c(test.tree.stage1@info$ql[1], test.tree.stage2@info$ql[-1])
  return(test.tree.stage1)

}

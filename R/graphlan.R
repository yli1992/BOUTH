#' Visualizing the results of \code{bouth} with GraPhlAn
#' 
#' Create a tree figure for visualizing the results of \code{bouth}
#' 
#' @param bouth.out an output object from the \code{bouth} function.
#' @param show.leaf a logical value indicating whether or not the leaf nodes are
#'   shown in the output figure.
#' @param output.dir the directory for storing the output figure. The default is
#'   the current working directory.
#' @param graphlan.dir the directory where the GraPhlAn package is located.
#'   
#' @note This function is dependent on the python package GraPhlAn. Parameters
#'   \code{clade.separation}, \code{branch.thickness},
#'   \code{branch.bracket.depth}, \code{branch.bracket.width}, and
#'   \code{clade.marker.size} can be set the same as in the python package. Details can be found
#'   at
#'   (\href{https://bitbucket.org/nsegata/graphlan/overview}{https://bitbucket.org/nsegata/graphlan/overview}).
#'   
#'   
#' @references Asnicar, F, et al. "Compact graphical representation of
#'   phylogenetic data and metadata with GraPhlAn." PeerJ 3 (2015): e1029.
#'   
#' @examples
#' data(IBD)
#' 
#' ## performing the one-stage, weighted bottom-up test on the IBD data
#' test.1 = bouth(anno.table = IBD$tax.table, pvalue.leaves = IBD$pvalue.otus,
#' na.symbol = "unknown", far = 0.1, is.weighted = TRUE)
#' 
#' 
#' ## suppose the GraPhlAn package is located at graphlan_directory/
#' graphlan(bouth.out = test.1, graph.dir = graphlan_directory)
#'
#' @export
graphlan<-function(bouth.out, show.leaf = FALSE, output.dir = getwd(), graphlan.dir = NULL,
                   clade.separation = NULL, branch.thickness = NULL, branch.bracket.depth = NULL, branch.bracket.width= NULL,
                   clade.marker.size = NULL){

  if(is.null(graphlan.dir)){stop("Unspecified directory for GraPhlAn.")}
  if(!"graphlan"%in%dir(graphlan.dir)){stop("The GraphlAn hasn't been installed at this directory.")}
  if(is.null(clade.separation)) clade.separation = 0.5
  if(is.null(branch.thickness)) branch.thickness = 0.5
  if(is.null(branch.bracket.depth)) branch.bracket.depth = 0.5
  if(is.null(branch.bracket.width)) branch.bracket.width = 0.25
  if(is.null(clade.marker.size)) clade.marker.size = 15
  basics = c(paste("clade_separation",clade.separation,sep='\t'),
             paste("branch_thickness",branch.thickness,sep='\t'),
             paste("branch_bracket_depth",branch.bracket.depth, sep='\t'),
             paste("branch_bracket_width", branch.bracket.width, sep = '\t'),
             paste("clade_marker_size", clade.marker.size, sep='\t'),
             paste("clade_marker_edge_color", "#555555", sep='\t'),
             paste("clade_marker_edge_width", 0.5, sep='\t'),
             paste("branch_color_from_ancestor", 0, sep='\t'))
  gamma_prop1 = 0.02 # A subtree will be colored if it contains more than 2% of all leaf nodes
  gamma_prop2 = 0.1  # A subtree will have a label if it contains more than 10% of all leaf nodes
  tree = bouth.out$tree
  nn = nrow(tree@node)
  nl = sum(tree@node$level==1)
  nr = max(tree@node$level)
  if(!show.leaf){
    nr = nr - 1
  }
  ## create styles and obtain the mapping table
  sub.nodes = tree@node$num[tree@node$level == max(tree@node$level)-1]
  sub.count = array(0, dim = length(sub.nodes))
  out0 = NULL

  # the mapping table
  table = matrix(NA, nrow = nl, ncol = nr)
  for(i in 1:nl){
    stack = NULL; len.stack = 0
    curr = i
    while(curr <= nn){
      stack = c(curr, stack)
      len.stack = len.stack+1
      curr = tree@node$parent[curr]
    }
    if(!show.leaf & stack[len.stack] <= nl){
      stack = stack[-len.stack]
      len.stack = len.stack-1
    }
    table[i,1:len.stack] = tree@node$id[stack]
    sub.count[sub.nodes%in%stack] = sub.count[sub.nodes%in%stack] + 1
  }
  # colors for subtrees
  sub.nodes = sub.nodes[sub.count >= nl*gamma_prop1]
  sub.count = sub.count[sub.count >= nl*gamma_prop1]
  sub.nodes = tree@node$id[sub.nodes]

  mypalette = RColorBrewer::brewer.pal(length(sub.nodes),"Accent")
  for(i in 1: length(sub.nodes)){
    out0 = c(out0, paste(sub.nodes[i], "annotation_background_color", mypalette[i], sep = '\t'))
    if(sub.count[i] > gamma_prop2){
      out0 = c(out0, paste(sub.nodes[i], "annotation", sub.nodes[i], sep = '\t'))
    }
  }

  table = data.frame(table, stringsAsFactors = FALSE)
  colnames(table) = paste("attri",c(1:nr),sep='_')
  # unique lineages
  table = dplyr::distinct_(table)
  # prepare for the annotation I
  out1 = array(NA, nrow(table))
  for(i in 1:nrow(table)){
    out1[i] = paste(table[i,!is.na(table[i,])], collapse ='.')
  }

  ## label detected nodes
  sig.id = which(tree@test$reject)
  if(!show.leaf){
    sig.id = sig.id[sig.id > nl]
  }
  sig.id = tree@node$id[sig.id]
  out2 = array(NA, length(sig.id))
  for(i in 1:length(sig.id)){
    out2[i] = paste(sig.id[i], "clade_marker_color", "r", sep='\t')
  }

  ## output all files for GraPhlAn
  anno1 = paste(output.dir,"anno1.txt", sep = '/')
  write.table(c(basics,out0, out2), anno1, col.names = F, row.names = F, quote = F)
  anno2 = paste(output.dir,"anno2.txt", sep = '/')
  write.table(out1, anno2, col.names = F, row.names = F, quote = F)

  command = c("#!/bin/sh", paste("export PATH=", graphlan.dir ,"/graphlan/:$PATH",sep=''),
              "graphlan_annotate.py anno2.txt plot.xml --annot anno1.txt",
              "graphlan.py plot.xml plot.png --dpi 150 --size 4 --pad 0.2")

  script = paste(output.dir,"run.sh",sep='/')
  write.table(command, script, col.names = F, row.names = F, quote = F)
  system2("bash", args = "run.sh")
  file.remove(anno1)
  file.remove(anno2)
  file.remove(script)
}

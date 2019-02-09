


#' An example dataset
#'
#' A real data example from a study on human gut microbiome and
#' inflammatory bowel disease (IBD). The preprocessing was performed on
#' the online platform QIITA (\href{https://qiita.ucsd.edu}{https://qiita.ucsd.edu}).
#' @references Halfvarson, J, et al. "Dynamics of the human gut microbiome
#' in inflammatory bowel disease." Nature microbiology 2.5 (2017): 17004.
#' @format
#' \describe{
#'   \item{tax.table}{A 2360 by 8 data frame that gives taxomonic assignment to
#'   2360 operational taxonomic units (OTUs) from 8 levels, which from the
#'   top to bottom are kingdom, phylum, class, order, family, genus, species, and OTUs.
#'   The value `unknown' indicates missing assignment at some levels.}
#'   \item{pvalue.otus}{A vector of p-values for testing the differential
#'   abundance of each of the 2360 OTUs between the ulcerative colitis (UC)
#'   and control groups.}
#'
#' }
#' @examples
#' data(IBD)
"IBD"

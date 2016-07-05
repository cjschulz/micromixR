# Functions to normalize microbial count data

#### Begin clr_samples_rows ####
#' Convert OTU table to a proportion table
#' 
#' This function takes an OTU counts table as input returns a proportion table.
#' The default input is samples as rows (Margin=1), but samples can be set to
#' columns by setting Margin to 2.
#' 
#' @export
#' @examples
#' ## example otu table with samples as rows
#' otu_tab1 <- matrix( 
#' c(2, 4, 3, 1, 5, 7), 
#' nrow=2, ncol=3)  
#' rownames(otu_tab1) <- c("sample1","sample2")
#' colnames(otu_tab1) <- c("otu1","otu2","otu3")
#' 
#' prop_table(otu_tab1)


prop_table <- function(x, Margin=1) {
  if (Margin == 2) {
    x / colSums(x)[col(x)]
  } else {
    x / rowSums(x)[row(x)]
  }
}

#### End clr_samples_rows ####

#### Begin clr_samples_rows ####
#' Calculate the centered log transform (clr) of an otu table.
#' 
#' This function takes an OTU counts table as input returns the clr transform.
#' Samples must be rows, and columns as otus. The default log base is 2.
#' Log base can also be set to 10 or natural log
#' 
#' @references 
#' Aitchison, J. (1986) The Statistical Analysis of Compositional Data,
#' Monographs on Statistics and Applied Probability. Chapman & Hall Ltd.,
#' London (UK). 416p.
#' @export
#' @examples
#' ## example
#' otu_tab1 <- matrix( 
#' c(2, 4, 3, 1, 5, 7), 
#' nrow=2, ncol=3)  
#' rownames(otu_tab1) <- c("sample1","sample2")
#' colnames(otu_tab1) <- c("otu1","otu2","otu3")
#' 
#' clr_samples_rows(otu_tab1)

clr_samples_rows <- function(x, constant=0.5, log=2) {
  if (log == 2) {
    print("clr transorm log base 2")
    log.t <- (log2(x + constant))
  } else if (log == 10) {
    print("clr transorm log base 20")
    log.t <- (log10(x + constant))
  } else if (log == "e") {
    print("clr transorm natural log base")
    log.t <- (log(x + constant))
  }  else 
    print("clr transorm log base 2")
  log.t <- (log2(x + constant))
  
  log.t <- apply(log.t, 2, as.numeric)
  clr.t <- (log.t)-rowMeans(log.t)
  row.names(clr.t) <- row.names(x)
  return(clr.t)
}

#### End clr_samples_rows ####

#### Begin rarefy_samples ####
#' Rarefy OTU table to a specificied sample depth
#' 
#' This function takes an OTU counts table as input returns a subsampled table.
#' The default input is samples as rows (Margin=1), but samples can be set to
#' columns by setting Margin to 2. The rarified table is a multinomial sample 
#' with probabilities set to observed counts witihn each sample. The depth is
#' set by user, often is the min observed depth (i.e. min(rowSums(otu_table))). 
#'  
#' @export
#' @examples
#' ## example otu table with samples as rows
#' otu_tab1 <- matrix( 
#' c(2, 4, 3, 1, 5, 7), 
#' nrow=2, ncol=3)  
#' rownames(otu_tab1) <- c("sample1","sample2")
#' colnames(otu_tab1) <- c("otu1","otu2","otu3")
#' 
#' rarefy_table(otu_tab1)

rarefy_table <- function(x, depth = 1000, Margin = 1) {
  if (Margin == 2) {
    print("Samples as columns.")
    x <- t(x)
  } else {
    print("Samples as rows.")
  }
    otus.r <- apply(x, 1, function(x) rmultinom(1, size=depth, prob = x))
    otus.r <- t(otus.r)
    colnames(otus.r) <- colnames(x)
    return(otus.r)
}

#### End rarefy_table ####
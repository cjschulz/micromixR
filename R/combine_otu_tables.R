#' A function to combine otu tables
#'
#' This function takes two or more otu tables as a list of matrices, and combines them.
#' Must keep rows as samples and taxa as columns (or vice versa, must be consistent)
#' Must provide a list of otu matrices as input (i.e., not a list of dataframes)
#' 
#' @param otu_tabs a list of otu matrices to be combined
#' @export
#' @examples 
#' otu_tab1 <- matrix( 
#' c(2, 4, 3, 1, 5, 7), 
#' nrow=3, ncol=2) 
#' otu_tab2 <- matrix(
#'  c(3, 5, 1, 3, 5, 8),
#'  nrow=2, ncol=3)
#' rownames(otu_tab1) <- c("otu1","otu2","otu3")
#' colnames(otu_tab1) <- c("sample1","sample2")
#' rownames(otu_tab2) <- c("otu3","otu4")
#' colnames(otu_tab2) <- c("sample3","sample4","sample5")
#'
#' otu_tabs <- list(otu_tab1, otu_tab2)
#' combine_otu_tables(otu_tabs)


combine_otu_tables <- function(otu_tabs){
    sample_list <- lapply(otu_tabs, rownames) # next create a list of rownames
    taxa_list <- lapply(otu_tabs, colnames) # create a list of colnames (taxa here)
    all_samp_names <- Reduce(union, sample_list) #then get the union of the list, samples as rows here
    all_tax_names <- Reduce(union, taxa_list)
    all_samples <- matrix(0, nrow=length(all_samp_names), ncol=length(all_tax_names),
                          dimnames=list(all_samp_names, all_tax_names)) #create matrix of correct dimensions
    
    all_samples[rownames(otu_tabs[[1]]), colnames(otu_tabs[[1]])] <- otu_tabs[[1]] #add in first otu table
    #loop through the rest of the otu tables, adding each
    for (i in 2:length(otu_tabs)) { 
        all_samples[rownames(otu_tabs[[i]]), colnames(otu_tabs[[i]])] <- otu_tabs[[i]]
    } 
    
    return(all_samples)
    
}

# need more inforamtion to describe otu table format
# otu_tabs <- list(otu_table1, otu_table2)
# combined_otu_table <- combine_otu_tables(otu_tabs)



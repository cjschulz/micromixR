## need to read new biom 2.0 hfd5 format ##
#setwd("D:/Simulating Microbial Communities/micro data to simulate")

#source("http://bioconductor.org/biocLite.R")
#biocLite(c("rhdf5","biom"))
library(rhdf5)
library(biom)



#' @export

# This generates the matrix columns-wise
generate_matrix <- function(x){
  require(rhdf5)
  require(biom)
  indptr  = x$sample$matrix$indptr+1
  indices = x$sample$matrix$indices+1
  data    = x$sample$matrix$data
  nr = length(x$observation$ids)
  
  counts = sapply(2:length(indptr),function(i){
    x = rep(0,nr)
    seq = indptr[i-1]:(indptr[i]-1)
    x[indices[seq]] = data[seq]
    x
  })
  rownames(counts) = x$observation$ids
  colnames(counts) = x$sample$ids
  # I wish this next line wasn't necessary
  lapply(1:nrow(counts),function(i){
    counts[i,]
  })
}


generate_metadata <- function(x){
  metadata = x$metadata
  metadata = lapply(1:length(x$ids),function(i){
    id_metadata = lapply(metadata,function(j){
      if(length(dim(j))>1){ as.vector(j[,i,drop=FALSE]) }
      else{ j[i] }
    })
    list(id = x$ids[i],metadata=id_metadata)
  })
  return(metadata)
}

namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}


read_hdf5_biom<-function(file_input){
  x = h5read(file_input,"/",read.attributes = TRUE)
  data = generate_matrix(x)
  rows = generate_metadata(x$observation)
  columns = generate_metadata(x$sample)
  shape = c(length(data),length(data[[1]])) # dim(data)
  
  # Experimental -- need to actually load these from file
  id = attr(x,"id")
  vs = attr(x,"format-version")
  format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
  format_url = attr(x,"format-url")
  type = "OTU table"
  #type=attr(x,"type")
  generated_by = attr(x,"generated-by")
  date = attr(x,"creation-date")
  matrix_type = "dense"
  matrix_element_type = "int"
  
  namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
            rows,columns,shape,data)
}


########################
#' @export
# i=input .biom file, o=output .txt file
new_biom_to_otu <- function(i,o) {
  x = read_hdf5_biom(i)
  y <- biom(x)
  z <- biom_data(y)
  otus <- as.matrix(z)
  otus <- as.data.frame(t(otus))  
  write.table(otus, o, sep="\t", col.names=NA, row.names = T)
}

# example
# OBLT <- new_biom_to_otu("Obese_lean_twins_Qiita_77.biom","OBLT_otus.txt")

#########################






#y = read_hdf5_biom("Obese_lean_twins_Qiita_77.biom")
# This won't work now if the type is not acceptable.
#y <- biom(y)
#z <- biom_data(y)
#otus <- as.matrix(z) # This works
#otus <- as.data.frame(t(otus))

#write.table(otus, "Obese_lean_twins_Qiita_77_otus.txt", sep="\t")


#y = read_hdf5_biom("GLM_Qiita_1041.biom")
# This won't work now if the type is not acceptable.
#y <- biom(y)
#z <- biom_data(y)
#otus <- as.matrix(z) # This works
#otus <- as.data.frame(t(otus))

#write.table(otus, "GLM_Qiita_1041_otus.txt", sep="\t")


#y = read_hdf5_biom("NY_CP_soil_Qiita_2104.biom")
# This won't work now if the type is not acceptable.
#y <- biom(y)
#z <- biom_data(y)
#otus <- as.matrix(z) # This works
#otus <- as.data.frame(t(otus))

#write.table(otus, "NY_CP_soil_Qiita_2104_otus.txt", sep="\t")
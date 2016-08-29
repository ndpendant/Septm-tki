#how to install----
#install.packages("devtools")
#library("devtools")
#devtools::install_github("klutometis/roxygen")
#library(roxygen2)
#setwd("~/scripts/R_library/paul")
#document()
#setwd("..")
#install("paul")

#functions----

#' t_test_cols
#'
#' This function performs a t-test on specified columns by row of an entire data frame or matrix and returns the p-values as a vector. Returns an NA value if t.test() encounters an error instead of a p-value.
#' @param data_in Data frame or matrix
#' @param colsA The columns for the first group in the t-test
#' @param colsB The columns for the second group in the t-test
#' @param var_equal_in Are the variances of the two groups equal? Defaults to FALSE
#' @param paired_in Are the two groups paired? Defaults to FALSE
#' @keywords t-test, t.test
#' @export
#' @examples
#' t.test.cols()

t_test_cols = function(data_in,colsA,colsB,var_equal_in=FALSE,paired_in=FALSE) {    
    pvalue_vector = apply(data_in,1,function(x) tryCatch(t.test(as.numeric(x[colsA]),as.numeric(x[colsB]),var_equal = var_equal_in,paired = paired_in)$p.value,error=function(x) NA))
    return(pvalue_vector)
}


#' uniprot_regex
#'
#' This function returns a regular expression matching Uniprot accessions with and/or without isoforms.
#' @param isoforms Should the regex only include isoforms ("only"), not include isoforms ("none"), or include both ("both")? Defaults to "both" 
#' @keywords uniprot
#' @export
#' @examples
#' uniprot_regex()

uniprot_regex = function(isoforms="both") {
    if (isoforms == "both" ) {
        return("[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-[0-9]+)?|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
    }
    else if (isoforms == "only") {
        return("[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}(-[0-9]+)?")
    }
    else if (isoforms =="none") {
        return("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
    }
}

#' default_margins
#'
#' This function sets the default margins to c(5.1,4.1,4.1,2.1).
#' @param No parameters
#' @keywords margins
#' @export
#' @examples
#' default_margins()

default_margins = function() {
    par(mar=c(5.1,4.1,4.1,2.1))
}
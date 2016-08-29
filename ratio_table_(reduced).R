

#This function was created to add simplify the addition of columns in the ratio table.
#' @param  ratiotable data frame/matrix
#' @param  ratename is ratio data subset of the cell_line_meta (ex. erlo_ratio_names). # of elements = 1
#' @param  newcols is the string or vector of strings of the new titles for each column. # of elements = x
#' @param cutoff is a character or vector of strings related to adjustment of data #number of elements = #number of newcols
#' (increase1.5 = '1.5i', decrease1.5 = '1.5d', increase2 = '2i', decrease2 = '2d', # rows with data == '!na.is' )
#' 
# ******NUMBER OF ELEMENTS ASSIGNED TO NEWCOLS MUST EQUAL NUMBER OF ELEMENTS ASSIGNED TO CUTOFF********


genaddform <-function(ratiotable,ratename,newcols,cutoff)
{
  
  for(i in 1:length(newcols))
  {
   x = cutoff[i]
      
      if(x == '1.5i')
      {
        ratiotable[, newcols[i]] = apply(ratiotable[, ratename, drop=F], 1, function(x) {  sum(na.omit(x) > log2(1.5))})
      }
    
      else if(x == '1.5d')
      {
        ratiotable[, newcols[i]] = apply(ratiotable[, ratename, drop=F], 1, function(x){sum(na.omit(x) < log2(0.667))})
      }
    
      else if(x == '2i')
      {
        ratiotable[, newcols[i]] = apply(ratiotable[, ratename, drop=F], 1, function(x){sum(na.omit(x) > 1)})
      }
      else if(x == '2d')
      {
        ratiotable[, newcols[i]] = apply(ratiotable[, ratename, drop=F], 1, function(x){sum(na.omit(x) < -1)})
      }
      else if(x == '!na.is')
      {
        ratiotable[, newcols[i]] = apply(ratiotable[, ratename, drop=F], 1, function(x){sum(!is.na(x))})
      }
  }
  return(ratiotable)
}
      
  


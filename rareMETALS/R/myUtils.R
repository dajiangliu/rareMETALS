#' Safe sum
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
Sum <- function(a)
  {
    return(sum(a,na.rm=TRUE));
  }

#' Safe mean
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
Mean <- function(a)
  {
    return(mean(a,na.rm=TRUE));
  }
    
#' Safe colSums
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
ColSums<- function(a)
  {
    return(colSums(a,na.rm=TRUE));
  }
    
#' Safe rowSums
#'
#' @param a Input vector, matrix etc.
#' @return Return the sum of the elements in a;
#' @export
RowSums<- function(a)
  {
    return(rowSums(a,na.rm=TRUE));
  }

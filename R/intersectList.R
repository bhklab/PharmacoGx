#' Utility to find the intersection between a list of more than two vectors or
#' lists
#' 
#' This function extends the native intersect function to work on two or more
#' arguments.
#' 
#' @examples 
#' list1 <- list('a', 'b', 'c')
#' list2 <- list('a', 'c')
#' list3 <- list('a', 'c', 'd')
#' listAll <- intersectList(list1, list2, list3)
#' listAll
#' 
#' @param ... A list of or any number of vector like objects of the same mode,
#'   which could also be operated on by the native R set operations
#' @return A vector like object of the same mode as the first argument,
#'   containing only the intersection common to all arguments to the function
#' @export


intersectList <-
function(...) {
   args <- list(...)
   nargs <- length(args)
   if(nargs == 0) {
     return(args)
   }
   if (nargs == 1) {
     if (nargs == 1 && is.list(args[[1]])) {
       do.call("intersectList", args[[1]])
     } else {
       return (args[[1]])
     }
   } else if (nargs == 2) {
     return (intersect(args[[1]], args[[2]]))
   } else {
     return (intersect(args[[1]], intersectList(args[-1])))
   }
}

## End


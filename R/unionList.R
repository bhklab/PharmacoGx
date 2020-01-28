#' Utility to find the union between a list of more than two vectors or
#' lists
#' 
#' This function extends the native union function to work on two or more
#' arguments.
#' 
#' @examples 
#' list1 <- list('a', 'b')
#' list2 <- list('a', 'c')
#' list3 <- list('c', 'd')
#' listAll <- unionList(list1, list2, list3)
#' listAll
#' 
#' @param ... A list of or any number of vector like objects of the same mode,
#'   which could also be operated on by the native R set operations
#' @return A vector like object of the same mode as the first argument,
#'   containing all the elements of all arguments passed to the function
#' @export 



# unionList <-
# function(...) {
#    args <- list(...)
#    nargs <- length(args)
#    if(nargs == 0) {
#      return(unlist(args))
#    }
#    if (nargs == 1) {
#      if (nargs == 1 && is.list(args[[1]])) {
#        do.call("unionList", args[[1]])
#      } else {
#        return(unique(args[[1]]))
#      }
#    } else if (nargs == 2) {
#      return(union(args[[1]], args[[2]]))
#    } else {
#      return(union(args[[1]], unionList(args[-1])))
#    }
# }

unionList <-
function(...) {
  args <- list(...)
  nargs <- length(args)
  return(unique(unlist(do.call(c, args))))
}

#' Utility to find the symmetric set difference of a list of two or more 
#' vectors or lists
#' 
#' The function finds the symmetric set differnces between all the arguments, defined as
#' Union(args)-Intersection(args)
#' 
#' @examples 
#' list1 <- list('a', 'b', 'c')
#' list2 <- list('a', 'c')
#' list3 <- list('a', 'c', 'd')
#' listAll <- symSetDiffList(list1, list2, list3)
#' listAll
#' 
#' @param ... A list of or any number of vector like objects of the same mode,
#'   which could also be operated on by the native R set operations
#' @return A vector like object of the same mode as the first argument,
#'   containing only the symmetric set difference
#' @export
symSetDiffList <- function(...){
	return(setdiff(unionList(...), intersectList(...)))
}

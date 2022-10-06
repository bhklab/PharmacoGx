# ==== PharmacoSet Class


#'`[`
#'
#' @examples
#' data(CCLEsmall)
#' CCLEsmall["WM1799", "Sorafenib"]
#'
#' @param x object
#' @param i Cell lines to keep in object
#' @param j Drugs to keep in object
#' @param ... further arguments
#' @param drop A boolean flag of whether to drop single dimensions or not
#'
#'@return Returns the subsetted object
#'
#' @export
setMethod(`[`, 'PharmacoSet', function(x, i, j, ..., drop = FALSE){
  if(is.character(i)&&is.character(j)){
    return(subsetTo(x, cells=i, drugs=j,  molecular.data.cells=i))
  }
  else if(is.numeric(i) && is.numeric(j) && all(as.integer(i)==i) && all(as.integer(j)==j)){
    return(subsetTo(x, cells=sampleNames(x)[i], drugs=treatmentNames(x)[j],  molecular.data.cells=sampleNames(x)[i]))
  }
})
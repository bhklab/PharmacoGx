#' Cosine Permuations
#' 
#' @inherit CoreGx::cosinePerm
#' @inheritParams CoreGx::cosinePerm
#' 
#' @export
cosinePerm <- function(x, y, nperm=1000, 
                       alternative=c("two.sided", "less", "greater"), 
                       include.perm=FALSE, nthread=1)
{
  CoreGx::cosinePerm(x, y, nperm, alternative, include.perm, nthread)
}
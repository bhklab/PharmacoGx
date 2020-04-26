#' Drug sensitivity calling using waterfall plots
#'
#' @inherit CoreGx::callingWaterfall
#' @inheritParams CoreGx::callingWaterfall
#' 
callingWaterfall <-
function(x, type=c("IC50", "AUC", "AMAX"), intermediate.fold=c(4, 1.2, 1.2),
         cor.min.linear=0.95, name="Drug", plot=FALSE)
{
  CoreGx::callingWaterfall(x, type, intermediate.fold, cor.min.linear, name, 
                           plot)
}
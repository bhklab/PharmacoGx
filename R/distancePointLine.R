##========================================================
##  Credits:
##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
##  With grateful thanks for answering our needs!
##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
##========================================================
distancePointLine <-
function(x, #x-coordinate of point
                       y, #y-coordinate of point
                       a, #coefficient in line equation ax + by + c = 0
                       b, #coefficient in line equation ax + by + c = 0
                       c) { #coefficient in line equation ax + by + c = 0
  
  #Function calculates shortest distance between point and line in R^2.
  
  if (!(all(is.finite(c(x, y, a, b, c))))) {
    stop("All inputs to linePtDist must be real numbers.")
  }
  
  return(abs(a * x + b * y + c) / sqrt(a ^ 2 + b ^ 2))
}

  ## end

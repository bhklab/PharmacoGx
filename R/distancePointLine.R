##========================================================
##  Credits:
##  Theory by Paul Bourke http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
##  Based in part on C code by Damian Coventry Tuesday, 16 July 2002
##  Based on VBA code by Brandon Crosby 9-6-05 (2 dimensions)
##  With grateful thanks for answering our needs!
##  This is an R (http://www.r-project.org) implementation by Gregoire Thomas 7/11/08
##========================================================
distancePointLine <-
function (x, y, slope, intercept) {
 ## x, y is the point to test.
 ## slope, intercept is the line to check distance.
 ##
 ## Returns distance from the line.
 ##
 ## Returns 9999 on 0 denominator conditions.
 x1 <- x-10
 x2 <- x+10
 y1 <- x1*slope+intercept
 y2 <- x2*slope+intercept
 dd <- distancePointSegment(x, y, x1, y1, x2, y2)
 return(dd)
}

  ## end

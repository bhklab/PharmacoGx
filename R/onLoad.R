## TODO:: Determine work around for `[` operator that allows dispatch on a single argument
#.onLoad <- function(libname, pkgname) {
#    `[` <- function(x, ...) { if(is(x, 'LongTable')) subset(x, ...) else .Primitive('[')(x, ...) }
#    assign('[', `[`, .GlobalEnv)
#}

.onLoad <- function(libname, pkgname) {
    `[` <- function(x, ...) { if(is(x, 'LongTable')) subset(x, ...) else .Primitive('[')(x, ...) }
    assign('[', `[`, .GlobalEnv)
}

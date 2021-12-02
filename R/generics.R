################################################################################
#' Generic functions for accessing attributes of objects
#'
#' These generic functions are needed to access the objects' attributes.
#' Note that the naming convention \code{getAttribute} was applied, where
#' \code{attribute} is the name of the attribute/slot of the class of the
#' object.
#'
#' @name generics-accessors
#'
#' @param object object from which to get the value
#' @param ... optional parameters; for documentation see the documentation of
#'             the methods to each of the generic.
#'
#' @seealso
#' For an overview on the classes of the framework, and all of their
#' attributes, see the class diagrams in the package description
#' [cf. \code{\link{quantspec-package}}].


## Class-SmoothedPG

#' @name generics-accessors
#' @aliases getBeta
#' @export
setGeneric("getBeta",
           function(object, ...){standardGeneric("getBeta")})

#' @name generics-accessors
#' @aliases getBetaSdNaive
#' @export
setGeneric("getBetaSdNaive", function(object, ...){standardGeneric("getBetaSdNaive")})
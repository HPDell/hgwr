#' @noRd
make.dummy <- function(data) {
    Reduce(cbind, Map(make.dummy.extract, data, names(data)))
}

#' @noRd
#' 
#' @export 
make.dummy.extract <- function(col, name) {
    UseMethod("make.dummy.extract")
}

#' @noRd
#' 
#' @method make.dummy.extract character
make.dummy.extract.character <- function(col, name) {
    lev <- unique(col)
    lev <- lev[-length(lev)]
    dummy <- Reduce(c, lapply(lev, function(x) {
        list(as.numeric(col == x))
    }))
    lev_names <- paste(name, lev, sep = ".")
    names(dummy) <- lev_names
    as.data.frame(dummy)
}

#' @noRd
#' 
#' @method make.dummy.extract factor
make.dummy.extract.factor <- function(col, name) {
    lev <- levels(col)
    lev <- lev[-length(lev)]
    dummy <- Reduce(c, lapply(lev, function(x) {
        list(as.numeric(col == x))
    }))
    lev_names <- paste(name, lev, sep = ".")
    names(dummy) <- lev_names
    as.data.frame(dummy)
}

#' @noRd
#' 
#' @method make.dummy.extract logical
make.dummy.extract.logical <- function(col, name) {
    dummy <- list(as.numeric(col == TRUE))
    names(dummy) <- c(name)
    as.data.frame(dummy)
}

#' @noRd
#' 
#' @method make.dummy.extract default
make.dummy.extract.default <- function(col, name) {
    dummy <- list(as.numeric(col))
    names(dummy) <- c(name)
    as.data.frame(dummy)
}

### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

### 'character_class' must be a string containing a set of character
### specifications like in [0-9\\s] but **without** the square brackets
### e.g. "0-9\\s".
### The function is vectorized with respect to 'x'.
### Returns a 2-col matrix with one row per string in 'x'.
### Sanity check:
###   x <- c("", " ", "a", "2", "  ", " a", " 2", "a ", "aa", "a2", "2 ",
##           "2a", "22", "a23", "a 3", "a3 ", "2 3 ", " a  1 bb 2  33 ")
###   prefix_suffix <- split_numeric_suffix(x, "0-9\\s")
###   reconstructed <- paste0(prefix_suffix[ , 1], prefix_suffix[ , 2])
###   stopifnot(identical(x, reconstructed))
split_numeric_suffix <- function(x, character_class)
{   
    pattern <- sprintf("^.*[^%s]([%s]*)$", character_class, character_class)
    suffix <- sub(pattern, "\\1", x, perl=TRUE)
    prefix <- substr(x, 1L, nchar(x) - nchar(suffix))
    colnames <- c("prefix", "suffix")
    matrix(c(prefix, suffix), ncol=2L, dimnames=list(NULL, colnames))
}


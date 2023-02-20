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
###          "2a", "22", "a23", "a 3", "a3 ", "2 3 ", " a  1 bb 2  33 ")
###   prefix_suffix <- .split_suffix(x, "0-9\\s")
###   reconstructed <- paste0(prefix_suffix[ , 1], prefix_suffix[ , 2])
###   stopifnot(identical(x, reconstructed))
.split_suffix <- function(x, character_class)
{   
    pattern <- sprintf("^.*[^%s]([%s]*)$", character_class, character_class)
    suffix <- sub(pattern, "\\1", x, perl=TRUE)
    prefix <- substr(x, 1L, nchar(x) - nchar(suffix))
    colnames <- c("prefix", "suffix")
    matrix(c(prefix, suffix), ncol=2L, dimnames=list(NULL, colnames))
}

format_organism <- function(organism)
{
    ## Remove leading and trailing whitespaces (like strip() in Python).
    organism <- gsub("^\\s*|\\s*$", "", organism)
    ## Replace multiple whitespaces with single space.
    organism <- gsub("\\s+", " ", organism)
    first_letter <- substring(organism, 1, 1)
    other_letters <- substring(organism, 2)
    paste0(toupper(first_letter), tolower(other_letters))
}

abbreviate_organism_name <- function(organism)
{
    ## If 'organism' has a numeric suffix, we handle it separately.
    ## This is an attempt at doing something sensible with organism names
    ## like "Torque teno virus 1" where we want .abbreviate_organism_name()
    ## to return "Tvirus1".
    ## Note that we allow whitespaces in the numeric suffix but we'll remove
    ## them before adding the suffix back to the abbreviated organism name.
    prefix_suffix <- .split_suffix(organism, "0-9\\s")
    organism <- prefix_suffix[ , "prefix"]

    ## Abbreviate 'organism' e.g.
    ##   "Homo sapiens" -> "Hsapiens"
    ##   "Canis lupus familiaris" -> "Cfamiliaris"
    ##   "Torque teno virus" -> "Tvirus"
    parts <- strsplit(organism, "\\s+")[[1]]
    if (length(parts) <= 1) {
        abbr_organism <- parts
    } else {
        first_letter <- substr(head(parts, 1), 1, 1)
        last_part <- tail(parts, 1)
        abbr_organism <- paste0(first_letter, last_part)
    }

    ## Remove whitespaces from the numeric suffix.
    suffix <- gsub("\\s", "", prefix_suffix[ , "suffix"])

    ## Add numeric suffix back.
    paste0(abbr_organism, suffix)
}

check_pkg_maintainer <- function(pkg_maintainer)
{
    pattern <- "\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>"
    if (!grepl(pattern, pkg_maintainer, ignore.case=TRUE))
        stop(wmsg("please provide a valid email address"))
}

organism2biocview <- function(organism)
{
    parts <- strsplit(organism, " +")[[1]]
    first_part <- head(parts, 1)
    last_part <- tail(parts, 1)
    paste0(first_part, "_", last_part)
}

build_Rexpr_as_string <- function(seqnames)
{
    if (length(seqnames) == 0)
        return("character(0)")
    paste0('c', '(', paste0('"', seqnames, '"', collapse=","), ')')
}


# https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
is_singular <- function( L )
    inherits( try( solve( L ), silent = TRUE ), "try-error" )

## Function to standardise a set of affy data, whereby the mean foreach gene
## becomes zero, and the standard deviation becomes 1.
##
## Mark Cowley,21 Feb 2005
##
affy.standardise <- function(data) {
    means <- apply(data, 1, mean)
    sd <- apply(data, 1, sd)
    z <- (data - means) / sd

    dimnames(z) <- dimnames(data)

    return(z)
}

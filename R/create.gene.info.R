## This function requires that the package has been downloaded and installed eg:
## > library(reposTools)
## > install.packages2("mgu74av2")
## > install.packages2("hgu133a")
## If that fails:
## > source("http://www.bioconductor.org/getBioC.R")
## > getBioC("install.packages2")
##
## WARNING:
## It is not recommended that you change the environments as the colnames will not be correct!
##
## NB: tested on mgu74av2 and hgu133a
## Mark Cowley, 25 Feb 2005
##
create.gene.info <- function(chip="mgu74av2", environments=c("SYMBOL", "GENENAME", "LOCUSID", "UNIGENE",
                                                             "REFSEQ", "ACCNUM", "CHR", "CHRLOC")) {
    require(chip, character.only = TRUE)
    require(mjc)

    psid <- names(as.list(get(paste(chip,"ACCNUM", sep=""))))

    environments <- paste(chip, environments, sep="")

    gene.info <- psid
    for(data in environments) {
        ## 1) eg: data = "mgu74av2SYMBOL", so get(data) gets the environment
        ## 2) convert environment to a list. some of the elements may have more than one value
        ## 3) unlist.getelement only gets the first element from the list
        tmp <- unlist.getelement( as.list(get(data)), 1 )

        stopifnot(names(tmp) == psid)

        gene.info <- cbind(gene.info, tmp)
    }

    colnames(gene.info) <- c("Probe.Set.ID", "Symbol", "Name",
                             "LocusLink", "UniGene", "RefSeq",
                             "GenBank", "Gene.Chr", "Gene.Pos")

    ## Change factors to characters, then the chars to numeric for the 2 cols that need it:
    gene.info <- as.character.data.frame(gene.info)
    gene.info$LocusLink <- as.numeric(gene.info$LocusLink)
    gene.info$Gene.Pos <- as.numeric(gene.info$Gene.Pos)

    detach(pos = match(paste("package", chip, sep=":"), search()))
    cat("Unloading required package:", chip, "\n")

    return( gene.info )
}

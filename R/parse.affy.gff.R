#' parse.affy.gff
#' 
#' Affy GFF files are "nearly" GTF files in that the 9th/last column contains
#' key/value pairs. nearly, because GTF files must begin with the two keys of
#' "gene_id" and "transcript_id". The pairs are "; " seperated, last field may
#' not have the ";"
#' 
#' The 9th field in a gff table contains lots of semicolon seperated key/value
#' pairs. If you'd like to retrieve the value, associated with a particular
#' key,
#' then this method is the goods.
#' eg: "number_independent_probes 3; probeset_id 4027908; exon_cluster_id
#' 1056646;"
#' 
#' @param gff A gff object with at least a \code{$group} attribute
#' @param key Unknown
#' @return Undocumented return value
#' @author Mark Cowley, 2/1/08
#' @export
parse.affy.gff <- function(gff, key) {
    stopifnot(is.character(key))
    #
    # which rows contain the key, in the group (9th/last) column
    #
    rows <- grep(key, gff$group)
    values <- gff$group[rows]
    values <- sub(p("^.*", key, " "), "", values)
    values <- trim(sub(";.*", "", values)) # last entry is not followed by semicolon, but may have trailing space.
    
    res <- gff
    res$group <- NA
    res$group[rows] <- values
    
    return(res)
}


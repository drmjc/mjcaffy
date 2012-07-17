#' Convert Affymetrix ST array to CHIP file
#' 
#' Build a GSEA .chip file on the Affymetrix Gene & Exon Arrays.
#' The GSEA .chip
#' file has these headers: \dQuote{Probe Set ID}, \dQuote{Gene Symbol}, \dQuote{Gene Title}
#' and are tab delimited
#' 
#' Given either a \code{gene.info} object (arg2), or a \code{gene.info.Rda.gz} file (arg 1),
#' create a GSEA chip file, saved to \code{outfile} (arg 3).
#' For mice/rats, you can also provide a homologene \code{data.frame}, or file name,
#' and
#' for those gene symbols that have human homologs, they will be mapped
#' accordingly,
#' else the rat/mouse symbol will be made all uppercase.
#' 
#' Re the homologene file:
#' Warren has written \dQuote{\code{/pwbc/data/Homologene/parseHomologen.py}}
#' which takes the output from homologene FTP site (eg \code{homologeneBuild63.data})
#' and parses it. The required output looks like this:\cr
#' HID\\tHuman.Gene.ID\\tHuman.Gene.Symbol\\tMouse.Gene.ID\\tMouse.Gene.Symbol\\tRat.Gene.ID\\tRat.Gene.Symbol\cr
#' 3\\t34\\tACADM\\t11364\\tAcadm\\t24158\\tAcadm\cr
#' 5\\t37\\tACADVL\\t11370\\tAcadvl\\t25363\\tAcadvl\cr
#' 6\\t38\\tACAT1\\t110446\\tAcat1\\t25014\\tAcat1\cr
#' 7\\t90\\tACVR1\\t11477\\tAcvr1\\t79558\\tAcvr1\cr
#' ......\cr
#' 
#' @param gene.info.file the path to a gene.info.Rda.gz file. see
#'   \code{\link{parse.affy.ST.array.annotation}}
#' @param gene.info a \code{data.frame} created by \code{\link{parse.affy.ST.array.annotation}}.
#'   Must specify one of: \code{gene.info.file} or \code{gene.info}.
#' @param outfile the chip file name. Hint, should not contain hyphens since
#'   GSEA commandline is java-based and hates hyphens.
#' @param convert2human logical: if \code{TRUE}, then convert the gene symbols to human
#'   Gene Symbols. \code{homologene} MUST be supplied
#' @param homologene either the filename or the homologene data frame created
#'   by /pwbc/data/Homologene/parseHomologen.py. Set to \code{NULL}, or 
#'   \code{convert2human=FALSE} to ignore.
#' @return Creates a chip object with 3 columns: \dQuote{Probe Set ID}, \dQuote{Gene Symbol},
#'   \dQuote{Gene Title}
#' @author Mark Cowley, 2008-07-30
#' @export
make.GSEA.chip.affyST <- function(gene.info.file=NULL, gene.info=NULL, outfile=NULL, convert2human=TRUE,
	homologene="/pwbc/data/Homologene/HuMoRn.homologeneBuild65.txt") {

	if( is.null(gene.info.file) && is.null(gene.info) )
		stop("You can either specify the filename where the gene.info object is, or the object itself")

	if( !is.null(gene.info.file) ) {
		objectName <- load(gene.info.file)
		if ( objectName != "gene.info" )
			assign("gene.info", get(objectName))
	}
	stopifnot( !is.null(gene.info) )

	gene.info <- gene.info[order(gene.info[,1]), ]

	
	if( is.null(outfile) ) {
		if( is.null(gene.info.file) ){
			stop("You haven't specified the outfile nor the input filename")
		}
		else {
			outfile <- sub(".gene.info.[Rr][Dd][aA].gz", ".chip", gene.info.file)
			if(convert2human)
				outfile <- sub("\\.mm[0-9]+.chip|\\.rn[0-9]+\\.chip", ".HUMANSYMBOLS.chip", outfile)
			outfile <- file.path(dirname(outfile), gsub("-", "_", basename(outfile)))
		}
	}
	else {
		if( grepl("-", basename(outfile)) ) {
			cat("GSEA hates hyphens in filenames. replacing '-' with '_'.")
			outfile <- file.path(dirname(outfile), gsub("-", "_", basename(outfile)))
		}
	}
	cat("Writing output to", outfile, "\n")
	
	if( convert2human ) {
		# #
		# # Test whether the gene symbols are mostly upper-case, and thus humans,
		# # or whether they are mostly Sentence case, and thus mouse/rats.
		# # in the if statement, these proportions are:
		# # human: 0.939, rat: 0.124, human: 0.121 using the na26 files.
		# # 
		# symbols <- sunique(na.rm(gene.info$Symbol[gene.info$Symbol != "---"]))
		# SYMBOLS <- toupper(symbols)
		# if( sum(symbols == SYMBOLS)/length(symbols) < 0.5 ) {
		if( is.null(homologene) ) {
			warning("It looks like you're making a chip file for a non-human organism,
	I will make these gene symbols upper-case, but recommend that you supply a homologene data object")
			gene.info$Symbol <- toupper(gene.info$Symbol)
		}
		else {
			# removing controls and junk from the gene.info$Symbol column
			gene.info$Symbol[grep("AFFX", gene.info$Symbol)] <- NA
			gene.info$Symbol[gene.info$Symbol==gene.info$ProbeSetID] <- NA
			gene.info$Symbol[grep('^[0-9]+$', gene.info$Symbol)] <- NA # some predicted genes have fully numeric ProbeSetID's, and are VERY unlikely to be Human Gene Symbols.
			
			cat("Processing Homologene\n")
			# do the homologene mapping thang
			if( file.exists(homologene) ) {
				homologene <- read.delim(homologene)
				homologene[homologene=="---"] <- NA
			}
			else if (!is.data.frame(homologene) ) {
				stop("Must specify either the file name of the data.frame of the homologene object.\n")
			}

			if( grepl("Ra|Rn", gene.info.file) ) {
				homologene <- homologene[,c("Rat.Gene.Symbol", "Human.Gene.Symbol")]
			}
			else if( grepl("Mo|Mus|MG", gene.info.file) ) {
				homologene <- homologene[,c("Mouse.Gene.Symbol", "Human.Gene.Symbol")]
			}
			else {
				stop("Can't work out if this is mouse or rat...")
			}
			colnames(homologene)[1] <- "Gene.Symbol"
			homologene <- homologene[rowSums(is.na(homologene)) == 0, ]
			homologene <- homologene[order(homologene$Gene.Symbol), ]

			# 1. look up the orthologs.
			symbols.in <- sunique(na.rm(gene.info$Symbol[gene.info$Symbol != "---"]))
			orthologs <- vecmerge(symbols.in, homologene, by.y="Gene.Symbol", all.x=TRUE, all.y=FALSE)
			cat("Of the", nrow(gene.info), "probesets, representing", nrow(orthologs), "Gene Symbols", 
				nrow(orthologs) - sum(is.na(orthologs[,2])), "had an ortholog from homologene.\n")
			idx <- which(is.na(orthologs[,2]))
			cat("The remaining", length(idx), "were mapped to human symbols by uppercasing them.\n")
			orthologs[idx,2] <- toupper(orthologs[idx,1])
		
			# replace the GeneSymbols
			gene.info <- merge(gene.info, orthologs, by.x="Symbol", by.y="Gene.Symbol", all.x=TRUE, all.y=FALSE, sort=FALSE)
			gene.info$Symbol <- gene.info$Human.Gene.Symbol
		}
		# }
	}

	# make the output object...
	# We can change this to also have an alias column
	# columns <- intersect(c("ProbeSetID", "GeneSymbol", "Description", "GeneSymbol.Redundant"), colnames(gene.info))
	columns <- c("ProbeSetID", "Symbol", "Description")
	chip.data <- gene.info[,columns]
	cn <- c("Probe Set ID", "Gene Symbol", "Gene Title")
	colnames(chip.data) <- cn
	
	# # V1 
	# # export this table to a chip file.
	# header <- paste(cn, collapse="\t")#c("Probe Set ID\tGene Symbol\tGene Title")
	# OUT <- file(outfile, "w")
	# writeLines(header, OUT)
	# write.table(chip.data, OUT, sep="\t", row.names=FALSE, col.names=FALSE, na="---", quote=FALSE)
	# close(OUT)

	# V2 
	# export this table to a chip file.
	write.table(chip.data, outfile, sep="\t", row.names=FALSE, col.names=TRUE, na="---", quote=FALSE)	
}

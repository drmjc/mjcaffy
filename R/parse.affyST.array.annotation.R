#' parse an Affymetrix Gene or Exon ST transcript.csv file
#' 
#' parse an Affymetrix Gene or Exon ST transcript.csv file, and:\cr
#' - create gene.info.raw (essentially an R version of the full csv file)\cr
#' - create gene.info gene.info contains 6 columns:\cr
#' ProbeSetID\cr
#' geneSymbol (the 2nd "//" seperated term in gene_assignment)\cr
#' Description (the 3nd "//" seperated term in gene_assignment)\cr
#' GenBank ID (the 1st "//" seperated term in mrna_assignment, \code{NA}ing out non-GenBanks (ie ENS...))\cr
#' gene_assignment\cr
#' mrna_assignment\cr
#' 
#' Specifying which probes:\cr
#' Sometimes there are missing probes that are in the PGF/CLF file, but not in the
#' CSV file. To be sure, provide a vector of probe ID's for the expected probeset ID's.
#' For those probes that are missing annotation, default annotations will be reported.
#' 
#' @param csv.file path to the transcript.csv file
#' @param mps.file Useful Exon arrays to subset the transcript-level probesets to the core,
#'  extended or full. The mps file specifies which probesets are in each category. 
#'  This should be ignored for Gene Arrays.
#' @param probes a character vector of probe ID's. If \code{NULL} (the default), then 
#'  all probes in the csv file will be returned. if not \code{NULL}, then only those probes
#'  will be returned, in the same order as probes.
#' @param char.limit restrict the length of the text in each column to something manageble
#' @param create.files logical: create the tsv and RDa files?
#' @param verbose logical: verbose output?
#' 
#' @return 
#' For a given an input file, example: \code{csv.file="HuGene-1_0-st-v1.na26.hg18.transcript.csv"}\cr
#'   These files are created:\cr
#'   \dQuote{HuGene-1_0-st-v1.na26.hg18.transcript.Rda.gz}, and\cr
#'   \dQuote{HuGene-1_0-st-v1.na26.hg18.gene.info.Rda.gz}, and\cr
#'   \dQuote{HuGene-1_0-st-v1.na26.hg18.gene.info.tsv}, and\cr
#'   \code{gene.info} is invisibly returned
#' 
#' @author Mark Cowley, 2008-07-24
#' @export
#' @importFrom genomics make.ucsc.string
#' @importFrom excelIO write.xls
parse.affy.ST.array.annotation <- function(csv.file, mps.file=NULL, probes=NULL, char.limit=512, create.files=TRUE, verbose=TRUE) {
	gene.info.rda.file <- NULL
	gene.info.xls.file <- NULL
	
	!grepl("zip$", csv.file) || stop("It looks like you've passed a czv.zip file into the Affy ST annotation parser.\n")

	if( verbose ) cat("Importing CSV file:", basename(csv.file), "\n")
	gene.info.raw <- import.APT(csv.file, keep.first.column=TRUE)
	# gene.info.raw[gene.info.raw == "---"] <- NA
	
	isTranscript <- TRUE
	if( "exon_id" %in% colnames(gene.info.raw) ) {
		gene.info.raw$category <- gene.info.raw$probeset_type
		gene.info.raw$total_probes <- gene.info.raw$probe_count
		isTranscript <- FALSE
	}
	
	if( !is.null(mps.file) ) {
		if( verbose ) cat("Importing MPS file:", basename(mps.file), "\n")
		stopifnot( file.exists(mps.file) )
		mps <- import.mps(mps.file)
		# V1
		gene.info.raw <- gene.info.raw[match(mps$probeset_id, gene.info.raw$probeset_id), ]

		# V2
		# ids <- setdiff(mps$probeset_id, gene.info.raw$probeset_id)
		# gene.info.raw <- vecmerge(mps$probeset_id, gene.info.raw, by.y="probeset_id", all.x=TRUE, all.y=FALSE)
		# rownames(gene.info.raw) <- gene.info.raw$probeset_id

		# fix the output file names.
		if( grepl("core.mps", mps.file, ignore.case=TRUE) )
			csv.file <- sub("csv", "core.csv", csv.file)
		else if( grepl("extended.mps", mps.file, ignore.case=TRUE) )
			csv.file <- sub("csv", "extended.csv", csv.file)
		else if( grepl("full.mps", mps.file, ignore.case=TRUE) )
			csv.file <- sub("csv", "full.csv", csv.file)
		else {
			stop("the mps file doesn't end in [core|extended|full].mps\n")
		}
	}
	

	##############################################
	## genomic locations, UCSC style
	## - symbol, description
	##############################################
	if( verbose ) cat("Getting genomic location of each probeset.\n")
	
	pos <- make.ucsc.string(
		gene.info.raw$seqname,
		gene.info.raw$start,
		gene.info.raw$stop
	)
	idx <- !is.na(pos)
	pos[idx] <- paste(pos[idx], " (",gene.info.raw$strand[idx],")", sep="")
	pos[is.na(pos)] <- ""
	# names(pos) <- rownames(gene.info.raw)
	
	
	##############################################
	## gene-level descriptions
	## - symbol, description, symbols.redundant
	##############################################
	if( verbose ) cat("Getting Gene-level descriptions.\n")
	
	#
	# work out the geneSymbol and geneDescription from the gene_assignment column.
	# ... this is simply the first Symbol and Description that is found.
	#
	tmp <- sub(" /// .*", "", gene.info.raw$gene_assignment)
	tmp <- strsplit(tmp, " // ")
	names(tmp) <- gene.info.raw$probeset_id
	geneSymbol <- sapply(tmp, "[", 2)
	geneDescription <- sapply(tmp, "[", length(tmp[[1]]))
	
	#
	# work out the gene symbol aliases.
	#
	tmp <- strsplit(gene.info.raw$gene_assignment, " /// ")
	tmp2 <- sapply(tmp, function(x) {
		parts <- strsplit(x, " // ")
		parts2 <- sapply(parts, "[", 2)
		res <- paste(unique(parts2), collapse=" /// ")
		res
	})
	tmp2[tmp2=="NA"] <- NA
	# geneSymbolDesc$geneSymbol.Redundant <- tmp2
	geneSymbol.Redundant <- tmp2


	##############################################
	## mRNA-level descriptions
	## - symbol, description, symbols.redundant
	##############################################
	if( verbose ) cat("Getting mRNA-level descriptions.\n")

	#
	# work out the GenBank ID from the mrna_assignment column.
	#
	tmp <- strsplit(gene.info.raw$mrna_assignment, " // ")
	names(tmp) <- gene.info.raw$probeset_id

	mrnaSymbol <- sapply(tmp, "[", 1)
	mrnaDescription <- sapply(tmp, "[", isTranscript*2+1) # ie 1 if !isTranscript, or 3 if isTranscript
	# work out the alternative mRNA ID's
	tmp <- strsplit(gene.info.raw$mrna_assignment, " /// ")
	tmp2 <- sapply(tmp, function(x) {
		parts <- strsplit(x, " // ")
		parts1 <- sapply(parts, "[", 1)
		res <- paste(unique(parts1), collapse=" /// ")
		res
	})
	tmp2[tmp2=="NA"] <- NA
	mrnaSymbol.Redundant <- tmp2
	
	
	##############################################
	## all probesets need a symbol of some sort
	## - symbol, description
	##############################################
	if( verbose ) cat("Getting most curated annotation for each probeset.\n")
	
	
	##############################################
	# Every probeset needs an annotation!!
	# Ideally this is a:
	# gene annotation; if not available, then we need
	# mrna annotation; if not available, then we need
	# control annotation; if not available, then we need
	# "unannotated transcript"
	##############################################

	# 1. use the gene annotation to start off with, and start filling in NA's with lesser quality info.
	bestSymbol <- geneSymbol
	bestDescription <- geneDescription
	# 2. fill in missing info using the mrna info that we just determined.
	idx <- which(is.na(geneSymbol) & !is.na(mrnaSymbol))
	bestSymbol[idx] <- mrnaSymbol[idx]
	bestDescription[idx] <- mrnaDescription[idx]
	# 3. There are some main probes that are completely missing gene_assignment AND mrna_assignment
	idx <- which( gene.info.raw$category == "main" & 
				  is.na(gene.info.raw$gene_assignment) & 
				  is.na(gene.info.raw$mrna_assignment) )
	
	if( length(idx) > 0 ) {
		bestSymbol[idx] <- gene.info.raw$probeset_id[idx]
		bestDescription[idx] <- paste("unannotated transcript: ", pos[idx], sep="")
	}
	# 4. all controls on the chip.
	idx <- which(gene.info.raw$category != "main")
	if( length(idx) > 0 ) {
		bestSymbol[idx] <- gene.info.raw$probeset_id[idx]
		bestDescription[idx] <- gene.info.raw$category[idx]
		# 5. AFFX controls: "--- // --- // AFFX-BioB-3_at, bac_spike // ---..." -> "AFFX-BioB-3_at" and "bac_spike"
		# ... these are in the 2nd field in probeset.csv, or third field for the transcript.csv...
		# therefore just use sub()
		idx <- grep("AFFX-", gene.info.raw$mrna_assignment, fixed=TRUE)
		tmp <- gene.info.raw$mrna_assignment[idx]
		tmp <- sub("^.+AFFX", "AFFX", tmp)
		tmp <- sub(" //.*", "", tmp)
		tmp2 <- strsplit(tmp, ", ")
		bestSymbol[idx] <- sapply(tmp2, "[", 1)
		bestDescription[idx] <- sapply(tmp2, "[", 2)
	}
	
	
	##############################################
	## make gene.info
	##############################################
	if( verbose ) cat("Creating gene.info object.\n")

	gene.info <- data.frame(
		ProbeSetID=rownames(gene.info.raw), 
		transcript_cluster_id=gene.info.raw$probeset_id, # works for probeset.csv and transcript.csv files
		Symbol=bestSymbol,
		Description=bestDescription,
		GeneSymbol=geneSymbol,
		GeneDescription=geneDescription,
		GeneSymbol.Redundant=geneSymbol.Redundant,
		mrnaSymbol,
		mrnaDescription,
		mrnaSymbol.Redundant,
		pos=pos,
		gene_assignment=gene.info.raw$gene_assignment, 
		mrna_assignment=gene.info.raw$mrna_assignment,
		total_probes=gene.info.raw$total_probes,
		category=gene.info.raw$category
	)

	if( "crosshyb_type" %in% colnames(gene.info.raw) ) {
		# then this was probably from a Gene ST array
		gene.info$crosshyb_type <- gene.info.raw$crosshyb_type
		gene.info <- move.column(gene.info, "crosshyb_type", "pos")
	}
	else {
		# then probably an Exon ST array
		gene.info$crosshyb_type <- NA
	}
	##############################################
	

	##############################################
	# limit the really long fields to just the first char.limit characters
	gene.info$gene_assignment <- substr(gene.info$gene_assignment,1,char.limit)
	gene.info$mrna_assignment <- substr(gene.info$mrna_assignment,1,char.limit)
	# ncharGA <- nchar(gene.info$gene_assignment)
	# ncharMA <- nchar(gene.info$mrna_assignment)
	# idx <- which(ncharGA > char.limit)
	# if( length(idx) > 0 ) {
	# 	gene.info$gene_assignment[idx] <- substr(gene.info$gene_assignment[idx], 1, char.limit)
	# }
	# idx <- which(ncharMA > char.limit)
	# if( length(idx) > 0 ) {
	# 	gene.info$mrna_assignment[idx] <- substr(gene.info$mrna_assignment[idx], 1, char.limit)
	# }
	##############################################
	
		
	##############################################
	# reorder rows
	##############################################
	gene.info <- gene.info[order(gene.info$ProbeSetID), ]
	
	##############################################
	# restrict/expand to specified probes.
	##############################################
	if( !is.null(probes) ) {
		!any(duplicated(probes)) || stop("probes must contain unique probeset ID's")
		
		if( length(intersect(gene.info$ProbeSetID, probes)) == 0 ) {
			warning("none of your specified probes were found in the Affymetrix csv annotation file. -- ignoring your probe list.")
		}
		else {
			idx <- probes %in% gene.info$ProbeSetID
			if( any(!idx) ) {
				n <- sum(!idx)
				if( verbose ) cat(sprintf("%d probes were not found in the Affymetrix csv annotation file. -- growing result to include these probes.\n", n))
				gi2 <- gene.info[1:n, ]
				for(i in c("ProbeSetID", "transcript_cluster_id")) gi2[,i] <- probes[!idx]
				for(i in 3:ncol(gi2))                              gi2[,i] <- rep(NA, n)
				for(i in c("gene_assignment", "mrna_assignment"))  gi2[,i] <- rep("unannotated probeset: not found in Affymetrix csv file", n)
				gene.info <- rbind(gene.info, gi2)
			}
			gene.info <- gene.info[match(probes, gene.info[,1]), ]
		}
	}
	##############################################

	##############################################
	## save gene.info
	##############################################
	if( create.files ) {
		if( verbose ) cat("Saving gene.info object.\n")
		gene.info.rda.file <- sub("csv", "gene.info.Rda.gz", csv.file)
		gene.info.xls.file <- sub("csv", "gene.info.xls", csv.file)
		save(gene.info, file=gene.info.rda.file, compress=TRUE)
		write.xls(gene.info, gene.info.xls.file, na="---")
	}
	
	invisible( gene.info )
	##############################################
}
# CHANGELOG:
# 2008-07-24: V1
# 2011-10-05
# - improved roxygen comments
# - added probes argument
# - added test case
# 2011-10-27
# - support added for probeset.csv files

#' test_parse.affy.ST.array.annotation 
#'
#' @param verbose logical
#' @return nothing
#' @author Mark Cowley, 2011-10-05
#' @export
#' @importFrom oligo read.celfiles rma
#' @importFrom Biobase featureNames
test_parse.affy.ST.array.annotation <- function(verbose=TRUE) {

	###############################################################################
	# Mouse
	# NB - the na30 csv was missing some probes wrt featuresNames(rma(read.celfiles("array.CEL"))), 
	# 35519 vs 35556 = 37 probes.
	# This has been fixed in na32: nrow = 35556
	csv.file <- "/misc/FacilityBioinformatics/NfsVeyron/pub/eg.input/MoGene-1_0-st-v1.na30.mm9.transcript.csv"
	gi <- parse.affy.ST.array.annotation(csv.file, probes=NULL, create.files=FALSE, verbose=verbose)
	cat(sprintf("Found %d probes in the transcript.csv file\n", nrow(gi)))
	
	cel.files <- "/misc/FacilityBioinformatics/NfsVeyron/pub/eg.input/MoGene01.CEL"
	cel <- read.celfiles(cel.files, verbose=verbose)
	norm <- rma(cel)
	probes <- featureNames(norm)
	cat(sprintf("Found %d probes in an RMA normalised MoGene CEL file\n", length(probes)))
	gi2 <- parse.affy.ST.array.annotation(csv.file, probes=probes, create.files=FALSE, verbose=verbose)
	cat(sprintf("Found %d probes in the transcript.csv + CEL file file\n", nrow(gi2)))
	###############################################################################

	###############################################################################
	# test probeset csv parsing
	csv.zip.file <- "/Volumes/GRIW/FacilityBioinformatics/NfsVeyron/pub/genepattern/modules/Affymetrix2chip/csv/MoGene-1_1-st-v1.na32.mm9.probeset.csv.zip"
	csv.file <- crossplatform.unzip(csv.zip.file, dest=tempdir())
	csv.file <- csv.file[grep("csv$", csv.file)]
	# csv.file <- "/Users/marcow/src/GenePattern.modules/Affymetrix2chip/dev/HuEx-1_0-st-v2.na32.hg19.probeset.top1k.csv"
	
	gi <- parse.affy.ST.array.annotation(csv.file, probes=NULL, create.files=FALSE, verbose=verbose)
	cat(sprintf("Found %d probes in the probeset.csv file\n", nrow(gi)))
}


# # Slimline version of parse.affy.ST.array.annotation which creates a chip file
# #
# # don't think this will really save much time over parse.affy.ST.array.annotation,
# # and having it here introduces code duplication.
# #
# # currently untested
# netaffx.csv2chip <- function(csv.file, mps.file=NULL, char.limit=512, verbose=TRUE) {
# 	
# 	if( grepT("probeset.csv$", csv.file) ) {
# 		stop("This function currently supports transcript.csv files only. It shouldn't take much more work to get it handle probeset.csv files though.\n")
# 	}
# 	if( verbose ) cat("Importing CSV file:", basename(csv.file), "\n")
# 	gene.info.raw <- import.APT(csv.file, keep.first.column=TRUE)
# 	gene.info.raw[gene.info.raw == "---"] <- NA
# 	
# 	if( grepT("probeset.csv$", csv.file) ) {
# 		gene.info.raw$category <- gene.info.raw$probeset_type
# 	}
# 	
# 	if( !is.null(mps.file) ) {
# 		if( verbose ) cat("Importing MPS file:", basename(mps.file), "\n")
# 		stopifnot( file.exists(mps.file) )
# 		mps <- import.mps(mps.file)
# 		# V1
# 		gene.info.raw <- gene.info.raw[match(mps$probeset_id, gene.info.raw$probeset_id), ]
# 	}
# 	
# 	##############################################
# 	## genomic locations, UCSC style
# 	## - symbol, description
# 	## (useful for when there is no gene symbol)
# 	##############################################
# 	if( verbose ) cat("Getting genomic location of each probeset.\n")
# 	
# 	pos <- rep(NA, nrow(gene.info.raw))
# 	idx <- !is.na(gene.info.raw$seqname)
# 	pos[idx] <- make.ucsc.string(
# 		gene.info.raw$seqname[idx],
# 		gene.info.raw$start[idx],
# 		gene.info.raw$stop[idx])
# 	pos[idx] <- paste(pos[idx], " (",gene.info.raw$strand[idx],")", sep="")
# 	
# 	##############################################
# 	## gene-level descriptions
# 	## - symbol, description, symbols.redundant
# 	##############################################
# 	if( verbose ) cat("Getting Gene-level descriptions.\n")
# 	
# 	#
# 	# work out the geneSymbol and geneDescription from the gene_assignment column.
# 	# ... this is simply the first Symbol and Description that is found.
# 	#
# 	tmp <- strsplit(gene.info.raw$gene_assignment, " // ")
# 	names(tmp) <- gene.info.raw$probeset_id
# 	tmp2 <- lapply(tmp, "[", 2:3)
# 	# geneSymbolDesc <- data.frame(geneSymbol=sapply(tmp2, "[", 1), Description=sapply(tmp2, "[", 2))
# 	geneSymbol <- sapply(tmp2, "[", 1)
# 	geneDescription <- sapply(tmp2, "[", 2)
# 	
# 	##############################################
# 	## mRNA-level descriptions
# 	## - symbol, description, symbols.redundant
# 	##############################################
# 	if( verbose ) cat("Getting mRNA-level descriptions.\n")
# 
# 	#
# 	# work out the GenBank ID from the mrna_assignment column.
# 	#
# 	tmp <- strsplit(gene.info.raw$mrna_assignment, " // ")
# 	names(tmp) <- gene.info.raw$probeset_id
# 
# 	mrnaSymbol <- sapply(tmp, "[", 1)
# 	mrnaDescription <- sapply(tmp, "[", 3)
# 	
# 	##############################################
# 	## all probesets need a symbol of some sort
# 	## - symbol, description
# 	##############################################
# 	if( verbose ) cat("Getting most curated annotation for each probeset.\n")
# 	
# 	#
# 	# Every probeset needs an annotation!!
# 	# Ideally this is a gene annotation; if not available, then we need
# 	# mrna annotation; if not available, then we need
# 	# control annotation; if not available, then we need
# 	# "unannotated transcript"
# 	#
# 	# 1. use the gene annotation to start off with, and start filling in NA's with lesser quality info.
# 	bestSymbol <- geneSymbol
# 	bestDescription <- geneDescription
# 	# 2. fill in missing info using the mrna info that we just determined.
# 	idx <- which(is.na(geneSymbol) & !is.na(mrnaSymbol))
# 	bestSymbol[idx] <- mrnaSymbol[idx]
# 	bestDescription[idx] <- mrnaDescription[idx]
# 	# 3. There are some main probes that are completely missing gene_assignment AND mrna_assignment
# 	idx <- which( gene.info.raw$category == "main" & 
# 				  is.na(gene.info.raw$gene_assignment) & 
# 				  is.na(gene.info.raw$mrna_assignment) )
# 	
# 	if( length(idx) > 0 ) {
# 		bestSymbol[idx] <- gene.info.raw$probeset_id[idx]
# 		bestDescription[idx] <- paste("unannotated transcript: ", pos[idx], sep="")
# 	}
# 	# 4. all controls on the chip.
# 	idx <- which(gene.info.raw$category != "main")
# 	if( length(idx) > 0 ) {
# 		bestSymbol[idx] <- gene.info.raw$probeset_id[idx]
# 		bestDescription[idx] <- gene.info.raw$category[idx]
# 		# 5. AFFX controls: "--- // --- // AFFX-BioB-3_at, bac_spike // ---..." -> "AFFX-BioB-3_at" and "bac_spike"
# 		# ... these are in the 2nd field in probeset.csv, or third field for the transcript.csv...
# 		# therefore just use sub()
# 		idx <- grep("AFFX-", gene.info.raw$mrna_assignment, fixed=TRUE)
# 		tmp <- gene.info.raw$mrna_assignment[idx]
# 		tmp <- sub("^.+AFFX", "AFFX", tmp)
# 		tmp <- sub(" //.*", "", tmp)
# 		tmp2 <- strsplit(tmp, ", ")
# 		bestSymbol[idx] <- sapply(tmp2, "[", 1)
# 		bestDescription[idx] <- sapply(tmp2, "[", 2)
# 	}
# 	
# 	##############################################
# 	## make chip
# 	##############################################
# 	if( verbose ) cat("Creating chip object.\n")
# 	
# 	chip <- data.frame(
# 		Probe.Set.ID=rownames(gene.info.raw),
# 		Gene.Symbol=bestSymbol,
# 		Gene.Description=bestDescription
# 	)
# 
# 	chip <- chip[order(chip$Probe.Set.ID), ]
# 	
# 	return(chip)
# }

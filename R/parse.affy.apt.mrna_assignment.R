# Function to parse the mrna_assignment column out of an ST transcript.csv file
# (an presumably a probeset.csv file, though this hasn't been checked).
#
# These files contain a number of different sources of annotation, including:
#	RefSeq, GenBank, ENSEMBL (of various grades) etc...
#
# Each annotation source found in the annot file will be merged into 2 columns
# in the output, called eg: "ENSEMBL Transcript ID" and "ENSEMBL Transcript DESC",
# or "RefSeq ID" and "RefSeq Desc".
# Thus there will be 1+2x columns in output for the probeset_id, then each of the
#	x annotation sources.
# When multiple ID/Descriptions from the same annotation source are identified for
#	each gene, then they will be seperated by " // ".
#
# More info:
# This needs parsing up into bits:
# annot$mrna_assignment[1]
# [1] "GENSCAN00000069186 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:1:3018707:3044814:1 // chr1 // 100 // 100 // 33 // 33 // 0 /// GENSCAN00000048429 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:11:27072113:27075899:1 // chr1 // 79 // 73 // 19 // 24 // 0 /// GENSCAN00000020837 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:4:8040558:8046315:1 // chr1 // 54 // 73 // 13 // 24 // 0 /// GENSCAN00000035385 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:3:4161125:4188920:-1 // chr1 // 12 // 100 // 4 // 33 // 1"
# This should be 4 elements:
# fields[[1]]
# [1] "GENSCAN00000069186 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:1:3018707:3044814:1 // chr1 // 100 // 100 // 33 // 33 // 0" 
# [2] "GENSCAN00000048429 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:11:27072113:27075899:1 // chr1 // 79 // 73 // 19 // 24 // 0"
# [3] "GENSCAN00000020837 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:4:8040558:8046315:1 // chr1 // 54 // 73 // 13 // 24 // 0"   
# [4] "GENSCAN00000035385 // ENSEMBL Prediction //  cdna:Genscan chromosome:NCBIM37:3:4161125:4188920:-1 // chr1 // 12 // 100 // 4 // 33 // 1"  
# Each field has 9 parts, eg:
## fields[[1]]
# [[1]]
# [1] "GENSCAN00000069186"                                 
# [2] "ENSEMBL Prediction"                                 
# [3] "cdna:Genscan chromosome:NCBIM37:1:3018707:3044814:1"
# [4] "chr1"                                               
# [5] "100"                                                
# [6] "100"                                                
# [7] "33"                                                 
# [8] "33"                                                 
# [9] "0"                                                  
# 
# [[2]]
# [1] "GENSCAN00000048429"                                    
# [2] "ENSEMBL Prediction"                                    
# [3] "cdna:Genscan chromosome:NCBIM37:11:27072113:27075899:1"
# <snip>
#
#
# Parameters:
#	f: an optional file name to the csv annotation file
#	annot: a pre-imported annotation data.frame. Must supply one of f or annot
#		must contain at least mrna_assignment, and probeset_id columns
#	n: if NULL, then process all rows of annot, otherwise, process first n rows.
#
# Value:
#	a data.frame with 'n' rows for each row of annot, and many columns (see description).
#
# Todo:
#	- group columns into NCBI and ENSEMBL.
#	- pick up miRNA's into a miRBase ID/description column
#
# Mark Cowley, 2008-07-21
#
parse.affy.apt.mrna <- function(f=NULL, annot=NULL, n=100) {
	if( !is.null(f) && is.null(annot) )
		annot <- import.APT(f)
	stopifnot( !is.null(annot) )
	
	# annot$mrna_assignment[annot$mrna_assignment == "---"] <- NA
	
	if( !is.null(n) && n < nrow(annot) && n > 0 )
		annot <- annot[1:n,]
	DELIM1 <- " /// "
	DELIM2 <- " // "
	
	fields <- strsplit(annot$mrna_assignment, DELIM1)
	names(fields) <- annot$probeset_id
	fields <- lapply(fields, function(toplevel) { trim(strsplit(toplevel, DELIM2)) })
	
	.annot.sources <- function(ll) {
		sapply(ll, "[", 2)
	}
	annot.sources.counts <- sort(table(unlist(lapply(fields, .annot.sources))), decreasing=TRUE)
	# cat(data.frame(names(annot.sources.counts), annot.sources.counts, stringsAsFactors=FALSE), "\n")
	annot.sources <- names(annot.sources.counts)

	columns <- c("probeset_id", paste(rep(annot.sources, each=2), c("ID", "Desc")))
	res <- as.data.frame(matrix(NA, nrow(annot), length(columns)), stringsAsFactors=FALSE)
	dimnames(res) <- list(annot$probeset_id, columns)
	res$probeset_id <- rownames(res)
	
	for(i in 1:nrow(annot)) {
		assign("i", i, pos=1) # for debugging
		field <- fields[[i]]
		if( length(field) == 0 || (length(field) == 1 && is.na(field)) )
			next
		for(j in 1:length(field))  {
			subfield <- field[[j]]
			cols <- match(paste(subfield[2], c("ID", "Desc")), colnames(res))

			res[i, cols[1] ] <- safepaste(res[i,cols[1] ], subfield[1], sep=DELIM2)
			res[i, cols[2] ] <- safepaste(res[i,cols[2] ], subfield[3], sep=DELIM2)
		}	
	}
	res
}

# Function to parse the gene_assignment column out of an ST transcript.csv file
# (and presumably a probeset.csv file, though this hasn't been checked).
#
# These files contain a number of different sources of annotation, including:
#	RefSeq, GenBank, ENSEMBL (of various grades) etc...
#
# Each annotation source found in the annot file will be merged into 2 columns
# in the output, called eg: "ENSEMBL Transcript ID" and "ENSEMBL Transcript DESC",
# or "RefSeq ID" and "RefSeq Desc".
# Thus there will be 1+2x columns in output for the probeset_id, then each of the
#	x annotation sources.
# When multiple ID/Descriptions from the same annotation source are identified for
#	each gene, then they will be seperated by " // ".
#
# Parameters:
#	f: an optional file name to the csv annotation file
#	annot: a pre-imported annotation data.frame. Must supply one of f or annot
#		must contain at least gene_assignment, and probeset_id columns
#	n: if NULL, then process all rows of annot, otherwise, process first n rows.
#
# Value:
#	a data.frame with 'n' rows for each row of annot, and many columns (see description).
#
# Todo:
#	- group columns into NCBI and ENSEMBL.
#	- pick up miRNA's into a miRBase ID/description column
#
# Mark Cowley, 2008-07-21
#
parse.affy.apt.gene_assignment <- function(f=NULL, annot=NULL, n=100) {
	if( !is.null(f) && is.null(annot) )
		annot <- import.APT(f)
	stopifnot( !is.null(annot) )
	
	# annot$gene_assignment[annot$gene_assignment == "---"] <- NA
	
	if( !is.null(n) && n < nrow(annot) && n > 0 )
		annot <- annot[1:n,]
	DELIM1 <- " /// "
	DELIM2 <- " // "
	fields <- strsplit(annot$gene_assignment, DELIM1)
	names(fields) <- annot$probeset_id
	fields <- lapply(fields, function(toplevel) { trim(strsplit(toplevel, DELIM2)) })
	

	COLS <- c("ID", "Symbol", "Desc", "EntrezGeneID")
	columns <- c("probeset_id", COLS)
	res <- as.data.frame(matrix(NA, nrow(annot), length(columns)), stringsAsFactors=FALSE)
	dimnames(res) <- list(annot$probeset_id, columns)
	res$probeset_id <- rownames(res)
	
	cols <- match(COLS, colnames(res))
	for(i in 1:nrow(annot)) {
		assign("i", i, pos=1) # for debugging
		field <- fields[[i]]
		if( length(field) == 0 || (length(field) == 1 && is.na(field)) )
			next
		for(j in 1:length(field))  {
			subfield <- field[[j]]

			res[i,cols[1] ] <- safepaste(res[i,cols[1] ], subfield[1], sep=DELIM2)
			res[i,cols[2] ] <- safepaste(res[i,cols[2] ], subfield[2], sep=DELIM2)
			res[i,cols[3] ] <- safepaste(res[i,cols[3] ], subfield[3], sep=DELIM2)
			res[i,cols[4] ] <- safepaste(res[i,cols[4] ], subfield[5], sep=DELIM2)
		}	
	}
	# make Entrez Gene ID's numerical.
	res[,5] <- as.numeric( res[,5] )
	res
}

# Allows the pasting of an NA lhs to a sensible rhs.
# todo: allow variable number of arguments that auto trims or changes the NA's
#	- if first or last is NA, then skip them, if middle values are NA, then
#	possibly change them to "" ??
#
# Mark Cowley, 2008-07-21
#
safepasteOriginal <- function(lhs, rhs, sep=", ") {
	if(is.null(lhs) || is.na(lhs))
		res <- rhs
	else
		res <- paste(lhs, rhs, sep=sep)

	return(res)
}
safepaste <- function(..., sep=", ") {
	args <- list(...)
	nonNA <- which( !is.na(args) )
	if( length(nonNA) == 1 )
		args[[nonNA]]
	else if( length(nonNA) > 1 )
		paste(args[nonNA], collapse=sep)
	else
		NA
}

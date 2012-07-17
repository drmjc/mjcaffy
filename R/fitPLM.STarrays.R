#' fitPLM model on Affymetrix ST arrays
#' 
#' Affymetrix ST arrays have some probesets with thousands of probes which makes
#' fitPLM REALLY slow, and prone to running out of memory. Limit the probesets
#' to those with at most max.probes. Most of these probes that are excluded are
#' control probes, but not all.
#' 
#' Using a 4x2.66 GHz Mac Pro with 5 GB RAM, and 6 MoGene arrays, I find
#' max.probes [200,1000]
#' gives similar run time. limiting to
#' 100 probes = 74 s
#' 200 probes = 69 s
#' 500 probes = 67 s
#' 1000 probes = 80 s
#' all probes = > 4min before I get a memory allocation error.
#' My 1.8 GHz MacBook Pro can run probesets with < 1000 probes in 134s which is
#' acceptable.
#' 
#' @param object see fitPLM. My testing assumes that object is an AffyBatch.
#' @param \dots see fitPLM. My testing assumes that object is an AffyBatch.
#' @param subset see fitPLM. My testing assumes that object is an AffyBatch.
#'   subset is over-written by this function.
#' @param max.probes determine the largest number of probes allowed in a
#'   probeset.
#' @param verbose TRUE/FALSE to print a message about how many probes are
#'   removed.
#' @return the outpur from running fitPLM on a subset of probesets
#' @author Mark Cowley, 2009-06-24
#' @export
#' @importFrom affy cleancdfname
#' @importFrom affyPLM fitPLM
fitPLM.STarrays <- function(object, ..., subset=NULL, max.probes=1000, verbose=TRUE) {
	cdfenv <- cleancdfname(cdfName(object))
	require(cdfenv, character.only=TRUE)
	
	# determine which probesets to keep
	nprobes <- unlist(eapply(get(cdfenv), nrow, all.names=TRUE))
	subset <- names(nprobes[nprobes <= max.probes])
	nexcluded <- length(nprobes) - length(subset)
	if( verbose )
		cat("Excluding", nexcluded, "ProbeSets because they contain >", max.probes, "probes.\n")

	plm <- fitPLM(object, ..., subset=subset)
}

# #
# # fitPLM takes AGES on ST arrays, perhaps moGene and RaGene in particular.
# # This is due to some probesets having huge numbers of probes - most of these
# # are control probesets where mrna_assignment is all --- // --- // etc
# #
# # Since we use fitPLM in the QC report purely to get to the residuals and be able
# # to draw NUSE and RLE, can we optimise the fitPLM by subsetting to an optimal
# # number of probes??
# #
# # estimates based on 6 arrays:
# #
# # optimum on largeish machines appears to be 1000, which only eliminates a handful of probes.
# #
# # on veyron                on optimus             on cowbook
# # limiting to              limiting to            limiting to
# #	100 probes  = 92 s     	100 probes  = 74 s    	100 probes  = 117 s
# #	200 probes  = 90 s     	200 probes  = 69 s    	200 probes  = 110 s
# #	500 probes  = 91 s     	500 probes  = 67 s    	500 probes  = 109 s
# #	1000 probes = 125 s    	1000 probes = 80 s    	1000 probes = 135 s
# #	all probes  = > 20min  	all probes  = > 4min  	all probes  = > 8.5 min
# #
# setwd("/pwbc/projects/MelissaMoore")
# load("Robjects/cel.Rda.gz")
# require("affyPLM")
# cdfenv <- "mogene10stv1cdf"
# require( cdfenv, character.only=TRUE )
# 
# nprobes <- unlist(eapply(get(cdfenv), nrow, all.names=TRUE))
# sort(nprobes, decreasing=TRUE)[1:20]
# 	# 10338063 10338056 10338037 10338047 10338067 10338065 10338035 10338060 
# 	#     5326     2500     2466     1218     1189     1160     1022     1020 
# 	# 10338059 10338064 10338003 10338001 10338066 10483871 10455148 10608718 
# 	#     1003      844      723      622      567      317      178      165 
# 	# 10482528 10386495 10509645 10608705 
# 	#      158      155      153      133 
# 
# max.probes <- 100
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
#  #   user  system elapsed 
#  # 91.803   0.275  92.080 
# 
# max.probes <- 200
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
#  #   user  system elapsed 
#  # 90.408   0.097  90.505 
# 
# max.probes <- 500
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
#  #   user  system elapsed 
#  # 91.124   0.186  91.312 
# 
# max.probes <- 1000
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
# #    user  system elapsed 
# # 125.604   0.009 125.613 
# 
# max.probes <- 10000
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
# 	# Error in fitPLM(cel, subset = subset) : 
# 	#   Calloc could not allocate (28355625 of 8) memory
# 	# Timing stopped at: 1223 4.43 1238 
# 
# 
# 
# #
# # on optimus, it's faster than veyron!!
# #
# max.probes <- 100
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
#   #  user  system elapsed 
#   # 74.27    3.27   77.82 
# 
# max.probes <- 200
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
#    # user  system elapsed 
#    # 68.9     3.5    72.5 
# 
# max.probes <- 500
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
#   #  user  system elapsed 
#   # 67.35    3.53   70.70 
# 
# max.probes <- 1000
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
#   #  user  system elapsed 
#   # 80.85    9.36   84.14 
# max.probes <- 10000
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
# 	# (3943,0xa0589720) malloc: *** mmap(size=2273140736) failed (error code=12)
# 	# *** error: can't allocate region
# 	# *** set a breakpoint in malloc_error_break to debug
# 	# Error in fitPLM(cel, subset = subset) : 
# 	#   Realloc could not re-allocate (size -2021830496) memory
# 	# Timing stopped at: 249 53.4 205
# 
# 
# #
# # on cowbook
# #
# max.probes <- 100
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
# 	 #   user  system elapsed 
# 	 # 116.69    5.87  183.38 
# 
# max.probes <- 200
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
# 	 #   user  system elapsed 
# 	 # 109.63    7.29  179.60 
# 
# max.probes <- 500
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
# 	 #   user  system elapsed 
# 	 # 108.79    6.57  168.40 
# 
# max.probes <- 1000
# subset <- names(nprobes[nprobes <= max.probes])
# system.time(plm <- fitPLM(cel, subset=subset))
# 	  #  user  system elapsed 
# 	  # 134.7    12.9   218.4 
# 
# max.probes <- 10000
# subset <- names(nprobes[nprobes <= max.probes])
# 	# system.time(plm <- fitPLM(cel, subset=subset))
# 	# R(7766,0xa05b0720) malloc: *** mmap(size=2273140736) failed (error code=12)
# 	# *** error: can't allocate region
# 	# *** set a breakpoint in malloc_error_break to debug
# 	# Error in fitPLM(cel, subset = subset) : 
# 	#   Realloc could not re-allocate (size -2021830496) memory
# 	# Timing stopped at: 508 49.6 633 

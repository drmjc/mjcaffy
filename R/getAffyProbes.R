#' Function to obtain the probe sequences for an Affy ProbeSet(s)
#' 
#' @author Mark Cowley
#' @export
getAffyProbes <- function(ProbeSetIDs, probePackage="hgu133plus2probe") {
	stopifnot( require(probePackage, character.only=TRUE) )

	if( length(ProbeSetIDs) > 1) {
		res <- list()
		for(i in 1:length(ProbeSetIDs)) {
			res[[i]] <- getAffyProbes(ProbeSetIDs[i], probePackage)
		}
		res
	}
	else {
		res <- as.data.frame( get(probePackage)[get(probePackage)$Probe.Set.Name == ProbeSetIDs,], stringsAsFactors=FALSE)
		res <- res[order(res$Probe.Interrogation.Position),]
		res$ProbeName <- paste(res$Probe.Set.Name, pad(1:nrow(res)), sep="_")
		res <- res[,c("ProbeName", "sequence")]
		res2 <- res$sequence
		names(res2) <- res$ProbeName
		return(res2)
	}
}

#' @importFrom genomics write.fasta
exportAffyProbes.fasta <- function(ProbeSetIDs, probePackage="hgu133plus2probe", outdir=".") {
	if( length(ProbeSetIDs) > 1 ) {
		for(i in 1:length(ProbeSetIDs)) {
			exportAffyProbes.fasta(ProbeSetIDs[1], probePackage, outdir)
		}
	}
	else {
		probes <- getAffyProbes(ProbeSetIDs)
		f <- file.path(outdir, paste(ProbeSetIDs, ".fa", sep=""))
		write.fasta(probes, names(probes), f)
	}
}

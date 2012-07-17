\name{make.GSEA.chip.affyST}
\alias{make.GSEA.chip.affyST}
\title{Convert Affymetrix ST array to CHIP file}
\usage{
  make.GSEA.chip.affyST(gene.info.file = NULL,
    gene.info = NULL, outfile = NULL, convert2human = TRUE,
    homologene = "/pwbc/data/Homologene/HuMoRn.homologeneBuild65.txt")
}
\arguments{
  \item{gene.info.file}{the path to a gene.info.Rda.gz
  file. see \code{\link{parse.affy.ST.array.annotation}}}

  \item{gene.info}{a \code{data.frame} created by
  \code{\link{parse.affy.ST.array.annotation}}.  Must
  specify one of: \code{gene.info.file} or
  \code{gene.info}.}

  \item{outfile}{the chip file name. Hint, should not
  contain hyphens since GSEA commandline is java-based and
  hates hyphens.}

  \item{convert2human}{logical: if \code{TRUE}, then
  convert the gene symbols to human Gene Symbols.
  \code{homologene} MUST be supplied}

  \item{homologene}{either the filename or the homologene
  data frame created by
  /pwbc/data/Homologene/parseHomologen.py. Set to
  \code{NULL}, or \code{convert2human=FALSE} to ignore.}
}
\value{
  Creates a chip object with 3 columns: \dQuote{Probe Set
  ID}, \dQuote{Gene Symbol}, \dQuote{Gene Title}
}
\description{
  Build a GSEA .chip file on the Affymetrix Gene & Exon
  Arrays. The GSEA .chip file has these headers:
  \dQuote{Probe Set ID}, \dQuote{Gene Symbol}, \dQuote{Gene
  Title} and are tab delimited
}
\details{
  Given either a \code{gene.info} object (arg2), or a
  \code{gene.info.Rda.gz} file (arg 1), create a GSEA chip
  file, saved to \code{outfile} (arg 3). For mice/rats, you
  can also provide a homologene \code{data.frame}, or file
  name, and for those gene symbols that have human
  homologs, they will be mapped accordingly, else the
  rat/mouse symbol will be made all uppercase.

  Re the homologene file: Warren has written
  \dQuote{\code{/pwbc/data/Homologene/parseHomologen.py}}
  which takes the output from homologene FTP site (eg
  \code{homologeneBuild63.data}) and parses it. The
  required output looks like this:\cr
  HID\\tHuman.Gene.ID\\tHuman.Gene.Symbol\\tMouse.Gene.ID\\tMouse.Gene.Symbol\\tRat.Gene.ID\\tRat.Gene.Symbol\cr
  3\\t34\\tACADM\\t11364\\tAcadm\\t24158\\tAcadm\cr
  5\\t37\\tACADVL\\t11370\\tAcadvl\\t25363\\tAcadvl\cr
  6\\t38\\tACAT1\\t110446\\tAcat1\\t25014\\tAcat1\cr
  7\\t90\\tACVR1\\t11477\\tAcvr1\\t79558\\tAcvr1\cr
  ......\cr
}
\author{
  Mark Cowley, 2008-07-30
}


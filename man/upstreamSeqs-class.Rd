\name{upstreamSeqs-class}
\docType{class}
\alias{upstreamSeqs-class}
\alias{Nmers,numeric,character,upstreamSeqs-method}
\alias{keys,upstreamSeqs-method}
\alias{organism,upstreamSeqs-method}
\alias{seqs,upstreamSeqs-method}
\alias{show,upstreamSeqs-method}
\alias{Nmers}
\alias{getUpstream} 
\alias{keys} 
\alias{organism} 
\alias{seqs}

\title{Class "upstreamSeqs" }
\description{Container for a collection of upstream sequences }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("upstreamSeqs", ...)}.
Environments are used to store collections of DNAstrings.
}
\section{Slots}{
	 \describe{
    \item{\code{seqs}:}{Object of class \code{"environment"} ~~ }
    \item{\code{chrom}:}{Object of class \code{"environment"} ~~ }
    \item{\code{revComp}:}{Object of class \code{"environment"} ~~ }
    \item{\code{type}:}{Object of class \code{"environment"} ~~ }
    \item{\code{organism}:}{Object of class \code{"character"} ~~ }
    \item{\code{provenance}:}{Object of class \code{"ANY"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{Nmers}{\code{signature(n = "numeric", orf = "character", usobj = "upstreamSeqs")}: obtain
all subsequences of length n as view elements of a DNA string }
    \item{keys}{\code{signature(x = "upstreamSeqs")}: ... }
    \item{organism}{\code{signature(x = "upstreamSeqs")}: ... }
    \item{seqs}{\code{signature(x = "upstreamSeqs")}: ... }
    \item{show}{\code{signature(object = "upstreamSeqs")}: ... }
	 }
}
\author{Vince Carey \email{stvjc@channing.harvard.edu}}
\examples{
showClass("upstreamSeqs")
data(sceUpstr)
sceUpstr
keys(sceUpstr)[1:5]
}
\keyword{classes}

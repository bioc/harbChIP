\name{harbChIP}
\docType{data}
\alias{harbChIP}
\title{Experimental Data Package: harbChIP}
\description{binding ratios and intergenic region data from
a ChIP-chip experiment in yeast}
\usage{
data(harbChIP)
}
\format{
  The format is:
  An \code{ExpressionSetObject} with covariates:
\itemize{
\item \code{txFac}: transcription factor symbol from Harbison website CSV file columnnames
}
}

\note{
derived from web site \url{jura.wi.mit.edu/young_public/regulatory_code/GWLD.html},
binding ratios
}

\examples{
  data(harbChIP)
  harbChIP
  experimentData(harbChIP)
  exprs(harbChIP)[1:6,1:7]
  dim(exprs(harbChIP))
  pData(featureData(harbChIP))[1:6,]
}
\keyword{datasets}


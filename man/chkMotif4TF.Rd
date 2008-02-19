\name{chkMotif4TF}
\alias{chkMotif4TF}
\alias{chkAllUS}
\title{analyze relationship between motif frequency
and binding intensity for selected motif and TF}
\description{analyze relationship between motif frequency
and binding intensity for selected motif and TF
}
\usage{
chkMotif4TF(motif, TF, chset, upstr, bthresh=2, countthresh=0)
}
\arguments{
  \item{motif}{ character string in alphabet known to Biostrings}
  \item{TF}{ name of a TF (sample name in the ChIP-chip data
structure \code{chset}}
  \item{chset}{ an ExpressionSet instance harboring ChIP-chip data
}
  \item{upstr}{ an instance of upstreamSeqs}
  \item{bthresh}{ threshold for binding intensity results to
declare TF 'bound' to the upstream region}
  \item{countthresh}{ threshold for motif count to be considered
'present' in upstream region}
}
\details{
Uses \code{\link[Biostrings]{countPattern}} to perform 
motif count.
}
\value{a list with elements \code{call}, \code{table},
and \code{test}, the latter providing the result of
\code{\link{fisher.test}}
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
# slow
\dontrun{
data(harbChIP)
data(sceUpstr)
chkMotif4TF("CGGCCG", "RDS1", harbChIP, sceUpstr)
}
}
\keyword{ models }

\name{allhex}
\alias{allhex}
\title{utility function: get all hexamers in upstream sequence for an ORF}
\description{utility function: get all hexamers in upstream sequence for an ORF}
\usage{
allhex(orf, usobj)
}
\arguments{
  \item{orf}{ character string, ORF name}
  \item{usobj}{ upstreamSeqs object}
}
\details{
computes Biostrings views
}
\value{
computes Biostrings views
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(sceUpstr)
allhex("YAL001C", sceUpstr)
}
\keyword{ models }

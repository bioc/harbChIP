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
computes Biostrings Views
}
\value{
computes Biostrings Views
}
%\references{  }
\author{Vince Carey \email{stvjc@channing.harvard.edu}}
%\note{  }


%\seealso{  }

\examples{
data(sceUpstr)
allhex("YAL001C", sceUpstr)
}
\keyword{ models }

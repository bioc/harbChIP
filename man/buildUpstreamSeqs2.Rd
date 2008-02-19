\name{buildUpstreamSeqs2}
\alias{buildUpstreamSeqs2}
\title{workflow component -- build an upstreamSeqs instance from a FASTA read}
\description{workflow component -- build an upstreamSeqs instance from a FASTA read
}
\usage{
buildUpstreamSeqs2(fastaRead, organism="sce", provenance="harmen")
}
\arguments{
  \item{fastaRead}{ results of a readFASTA from Biostrings}
  \item{organism}{ string naming organism}
  \item{provenance}{ string or structure describing provenance}
}
\details{
generates an instance of upstreamSeqs
}
\value{
generates an instance of upstreamSeqs
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
# x = readFASTA(...)
# y = buildUpstreamSeqs2(x)
}
\keyword{ models }

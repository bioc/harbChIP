\name{sceUpstr}
\alias{sceUpstr}
\title{Biostrings representations of S. cerevisiae upstream regions}
\description{Biostrings representations of S. cerevisiae upstream regions}
\usage{
data(sceUpstr)
}
\details{
environment-based S4 object with DNAstring elements
}
\value{
environment-based S4 object with DNAstring elements
}

\author{Vince Carey \email{stvjc@channing.harvard.edu}}

\examples{
data(sceUpstr)
sceUpstr
getUpstream("YAL001C", sceUpstr)
}
\keyword{ models }

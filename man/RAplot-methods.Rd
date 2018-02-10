\name{RAplot-methods}
\docType{methods}
\alias{RAplot-methods}
\alias{RAplot,RA-method}
\alias{RAplot,FLUSH-method}
\alias{RAplot}
\title{Array vs RESidual standard deviation plot }
\description{
   Plot the square root of the Array effect as a function of log-Residual standard deviation 
}
\usage{
\S4method{RAplot}{RA}(object, ...)
\S4method{RAplot}{FLUSH}(object, ...)
}
\arguments{
  \item{object}{Object of class "RA" or "FLUSH"}
  \item{\dots}{}
}
\details{

}
\value{
None. Just plot.
}
\references{Stefano Calza, Wolfgang Raffelsberger, Alexander Ploner, Jose Sahel, Thierry Leveillard
and Yudi Pawitan. Filtering genes to improve sensitivity in oligonucleotide
microarray data analysis - 2007 NAR}
\author{Stefano Calza <calza@med.unibs.it> }
\note{}
\seealso{\code{\link{fitRA}}}
\examples{

\dontrun{
demo(FLUSH.tour)

RAplot(choe.fSet,add.rq=T)
}
}
\keyword{methods}
\keyword{hplot}

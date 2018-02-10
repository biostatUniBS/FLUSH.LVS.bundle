\name{Flush-methods}
\docType{methods}
\alias{Flush-methods}
\alias{Flush,ExpressionSet-method}
\alias{Flush,matrix-method}
\alias{Flush}
\title{Flush genes}
\description{
Perform filetring algorithm on an object of class exprSet or matrix,
returning a subsetted object.
}
\usage{
\S4method{Flush}{matrix}(object, RA, check = TRUE, ...)
\S4method{Flush}{ExpressionSet}(object, RA, check = TRUE, onlyExprs =
TRUE, which = 1)
\S4method{Flushit}{ExpressionSet}(object,RA, proportion=0.6, df=10, delta=NULL, lambda=NULL,
                              check = TRUE, onlyExprs=TRUE, which=1,...)

}
\arguments{
  \item{object}{An object of class exprSet or matrix}
  \item{RA}{An object of class RA as from compute.RA and optionally FlushSet}
  \item{proportion}{Proportion of genes to filter out, i.e. the quantile
    to estimate in the quantile regression. Used only if RA doesn't come
    from FlushSet}
  \item{delta}{Number. Weighting parameter. Used only if RA doesn't come
    from FlushSet}
  \item{lambda}{Number. Weighting parameter. Used only if RA doesn't come
    from FlushSet}
  \item{df}{Integer. Degrees of freedom of the B-spline. Used only if RA doesn't come
    from FlushSet}
  \item{check}{Logical. If TRUE some additional check on probes name
    consistency are performed.}
  \item{onlyExprs}{Logical. If TRUE returns an object of class
    ExpressionSet with the selected features. If FALSE an object of
    class FLUSH is returned with the RA object included (for
    plotting)}
  \item{which}{Integer. If more than one proportion is supplied, which one to use
    for filtering. Defaults to 1.}
  \item{...}{additional arguments.Not used.}
}
\details{
}
\value{
  An object of class FLUSH or ExpressionSet (if onlyExprs=TRUE).
}
\references{Stefano Calza, Wolfgang Raffelsberger, Alexander Ploner, Jose Sahel, Thierry Leveillard
and Yudi Pawitan. Filtering genes to improve sensitivity in oligonucleotide
microarray data analysis - 2007 NAR}
\author{Stefano Calza <calza@med.unibs.it> }
\note{
}
\seealso{\code{\link{RAplot}}, \code{\link{fitRA}}, \code{\link{FlushSet}}}
\examples{

\dontrun{
data(FLUSH.RData)

choe.RA <- compute.RA(choe)
choe.fSet <- FlushSet(choe.RA)
choe.Flush <- Flush(choe.MAS5,choe.fSet, onlyExprs = FALSE)
choe.Flush <- Flush(choe.MAS5,choe.choe.RA, proportion=0.6)
}
}
\keyword{methods}
\keyword{manip}


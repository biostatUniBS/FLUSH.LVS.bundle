\name{normalize.lvs}
\alias{normalize.lvs}
\alias{normalize.AffyBatch.lvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LVS normalization function(s) }
\description{
  Performs the normalization step in LVS procedure.
}
\usage{
normalize.lvs(object, ref.fun = c("median", "mean"), lvs.id, use.loess = FALSE, ...)
normalize.AffyBatch.lvs(eset, ref.fun = c("median", "mean"), lvs.id, use.loess = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ a matrix containing unnormalized expression values.}
  \item{eset} { an ExpressionSet or exprSet object.}
  \item{ref.fun}{ function used to compute the reference distribution. Default is \code{median}.}
  \item{lvs.id}{ The Least variant ID Probes selected after quantilic regression.}
  \item{use.loess}{ fit the normalizing curve using a smooth.spline (default) or a loess.}
  \item{\dots}{ \code{\dots} }
}
%\details{
%}
\value{
  For \code{normalize.lvs} a matrix containing LVS normalized expression values. 
  For \code{normalize.AffyBatch.lvs} a normalized ExpressionSet.
  }
 
\references{ S. Calza et al. 'Normalization of oligonucleotide arrays based on the least variant set of genes' (2008, BMCBioinformatics).}
\author{ Stefano Calza <calza@med.unibs.it>, Davide Valentini and Yudi Pawitan.}
% \note{ ~~further notes~~ 
% }
\seealso{ \code{\link{lvs}}, \code{\link{compute.RA}}, \code{\link{lvs.fit}}}
\examples{
\dontrun{
# Starting from an Affybatch object called aBatch 
data.RA <- compute.RA(aBatch)
lvs.id <- lvs.fit(data.RA, proportion=0.6)
lvs.prep <- expresso(aBatch, normalize=FALSE, bgcorrect.method="mas",pmcorrect.method="mas",summary.method="mas")
normalize.AffyBatch.lvs(lvs.prep,lvs.id=lvs.id)
 }}
\keyword{ normalization }
\keyword{ lvs }

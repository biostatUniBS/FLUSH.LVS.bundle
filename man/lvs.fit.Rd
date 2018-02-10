\name{lvs.fit}
\alias{lvs.fit}

\title{ LVS selection function }
\description{
  Selects the Least Variant Set of genes, accordingly to the chosen proportion of genes expected not to vary between arrays. 
}
\usage{
lvs.fit(object, proportion = 0.6, DF=10,...)
}
\arguments{
  \item{object}{an object of class RA or a matrix with residual standard errors and array effects.}
  \item{proportion}{ the proportion below which genes are expected not to vary between samples. Default is set to 0.6. }
  \item{DF}{ the degrees of freedom used for the bi-spline in the quantile regression.}
  \item{\dots}{ \code{\dots} }
}
%\details{
%}
\value{
  a TRUE/FALSE vector with the seleted ID probes from the quantile regression.
   }
 
\references{ S. Calza et al. 'Normalization of oligonucleotide arrays
  based on the least variant set of genes' (2008, BMCBioinformatics).}
\author{ Stefano Calza <calza@med.unibs.it>, Davide Valentini and Yudi Pawitan. }
\seealso{ \code{\link{lvs}}, \code{\link{normalize.lvs}},
  \code{\link{compute.RA}},\code{\link{RA-class}},
  \code{\link{normalize.AffyBatch.lvs}} }
\examples{
\dontrun{
# Starting from an Affibatch object called aBatch 
data.RA <- compute.RA(aBatch)
lvs.id <- lvs.fit(data.RA, proportion=0.6)
lvs.prep <- expresso(aBatch, normalize=FALSE,
bgcorrect.method="mas",pmcorrect.method="mas",
summary.method="mas")
normalize.AffyBatch.lvs(lvs.prep,lvs.id=lvs.id)
 }}
\keyword{ normalization }
\keyword{ lvs }

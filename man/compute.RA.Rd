\name{compute.RA}
\alias{compute.RA}
\alias{fitRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LVS function for the residual standard errors and the array effects }
\description{
  Given an \code{AffyBatch} object fits a robust linear model at probe level and return the residual standard errors and the array effects.
}
\usage{
compute.RA(abatch,bg.RA=c("nome","imm","rma"))
}
%- maybe also 'usage' for other objects documented here.

\value{
 an object of class RA with a 'data' slot containing a two columns data.frame with IDs as rows and RSE and Array effects as columns.
}
\references{ S. Calza et al. 'Normalization of oligonucleotide arrays based on the least variant set of genes' (2008, BMCBioinformatics).}
\author{ Stefano Calza <calza@med.unibs.it>, Davide Valentini and Yudi Pawitan.}

\seealso{ \code{\link{lvs}}, \code{\link{normalize.lvs}}, \code{\link{lvs.fit}}, \code{\link{normalize.AffyBatch.lvs}} \code{\link{fitRA}}}
\examples{
\dontrun{
# Starting from an Affibatch object called aBatch 
data.RA <- compute.RA(aBatch)
lvs.id <- lvs.fit(data.RA, proportion=0.6)
lvs.prep <- expresso(aBatch, normalize=FALSE, bgcorrect.method="mas",pmcorrect.method="mas",summary.method="mas")
normalize.AffyBatch.lvs(lvs.prep,lvs.id=lvs.id)
 }}
\keyword{ normalization }
\keyword{ lvs }

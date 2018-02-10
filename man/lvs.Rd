\name{lvs}
\alias{lvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ LVS Normalization }
\description{
  Performs a normalization of oligonucleotide arrays based on the Least Variant Set of genes.
}
\usage{
lvs(aBatch, proportion = 0.6, bgcorrect.method = "mas", pmcorrect.method = "mas", summary.method = "mas")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aBatch}{ an \code{AffyBatch} object. }
  \item{proportion}{ the proportion below which genes are expected not to vary between samples. Default is set to 0.6. }
  \item{bgcorrect.method}{ a background correction method among methods available in \code{\link[affy]{expresso}} function. Default is set to mas5. }
  \item{pmcorrect.method}{ a PM correction method among methods available in \code{\link[affy]{expresso}} function. Default is set to mas5. }
  \item{summary.method}{ a summarization method among methods available in \code{\link[affy]{expresso}} function. Default is set to mas5.  }
}
\details{
  This function performs a complete LVS procedure (see references) starting from an AffyBatch object, wrapping \code{compute.RA}, \code{lvs.fit} and \code{normalize.AffyBatch.lvs} functions. As a first step it selects Least Variant genes, then it preprocess and summarize the affybatch, using mas5 as default, and finally it normalize the data according to LVS procedure (taking into account identified Least Variant genes as reference Set), reporting expression measures. The default \code{proportion} is set to 0.6, accordingly to the experience that no more than 40\% of genes vary between arrays. 
}
\value{
  an ExpressionSet object, containing the LVS normalized values, with expression measures.
}
\references{ S. Calza et al. 'Normalization of oligonucleotide arrays based on the least variant set of genes' (2008, BMCBioinformatics).}
\author{ Stefano Calza <calza@med.unibs.it>, Davide Valentini and Yudi Pawitan. }
% \note{ ~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{ \code{\link{normalize.lvs}}, \code{\link{compute.RA}}, \code{\link{lvs.fit}}, \code{\link[affy]{expresso}}, \code{\link[affy]{read.affybatch}}} 
\examples{
\dontrun{
# Download and install Choe data package from http://www.meb.ki.se/~yudpaw/
data(choe)
data.norm <- lvs(choe)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lvs }
\keyword{ normalization }

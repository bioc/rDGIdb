\name{plotInteractionsBySource}
\alias{plotInteractionsBySource}
\title{
Interactions by source
}
\description{
Creates a plot to visualize the number of interactions per source.
}
\usage{
plotInteractionsBySource(queryResult, ...)
}
\arguments{
    \item{queryResult}{
A \code{\linkS4class{rDGIdbResult}} object obtained from
\code{\link{queryDGIdb}}.}
    \item{...}{Other arguments passed to barplot.}
}
\value{
A plot is created on the current graphics device.
}
\references{
Wagner AH, Coffman AC, Ainscough BJ, Spies NC, Skidmore ZL, Campbell KM,
Krysiak K, Pan D, McMichael JF, Eldred JM, Walker JR, Wilson RK, Mardis ER,
Griffith M, Griffith OL. DGIdb 2.0: mining clinically relevant drug-gene
interactions. Nucleic Acids Research. 2016 Jan 4;44(D1):D1036-44.
doi:10.1093/nar/gkv1165.
}
\author{
Thomas Thurnherr <thomas.thurnherr@bsse.ethz.ch>
}

\seealso{
\code{\linkS4class{rDGIdbResult}}, \code{\link{queryDGIdb}},
\code{\link{rDGIdbFilters}}
}
\examples{
queryResult <- queryDGIdb('KRAS')
plotInteractionsBySource(queryResult)
}
\keyword{ ~aplot }
\keyword{ ~device }

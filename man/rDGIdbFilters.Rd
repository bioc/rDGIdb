\name{rDGIdbFilters}
\alias{rDGIdbFilters}
\alias{geneCategories}
\alias{interactionTypes}
\alias{sourceDatabases}
\title{
Filter options for rDGIdb queries.
}
\description{
The gene categories, interaction types, and source databases available to
filter queries with \code{\link{queryDGIdb}}.
}
\usage{
geneCategories()

interactionTypes()

sourceDatabases()
}
\value{
\code{geneCategories} returns a character vector of available gene categories
which can be used to filter DGIdb queries with \code{\link{queryDGIdb}}.

\code{interactionTypes} returns a character vector of available interaction
types which can be used to filter DGIdb queries with \code{\link{queryDGIdb}}.

\code{sourceDatabases} returns a character vector of available source databases
which can be used to filter DGIdb queries with \code{\link{queryDGIdb}}.
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
\code{\link{queryDGIdb}}, \code{\linkS4class{rDGIdbResult}},
\code{\link{plotInteractionsBySource}}
}
\examples{
geneCategories()
interactionTypes()
sourceDatabases()
res <- queryDGIdb(genes = "AP1",
            interactionTypes = interactionTypes(),
            geneCategories = geneCategories(),
            sourceDatabases = sourceDatabases())
}
\keyword{ ~interface }
\keyword{ ~database }

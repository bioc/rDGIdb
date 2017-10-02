\name{queryDGIdb}
\alias{queryDGIdb}
\title{
Query DGIdb using R
}
\description{
Provides an interface to query DGIdb from within R using the DGIdb API. The
function is implemented to reflect the as good as possible the web-interface.
}
\usage{
queryDGIdb(genes,
        sourceDatabases = NULL,
        geneCategories = NULL,
        interactionTypes = NULL)
}
\arguments{
    \item{genes}{
A character vector of genes for which drug interactions are queried.
}
    \item{sourceDatabases}{
A character vector of source databases to be queried. To query all
available databases, skip argument or use \code{NULL}. To see available options,
type \code{sourceDatabases()}.
}
    \item{geneCategories}{
A character vector of gene categories to be queried. To query for all
gene categories, skip argument or use \code{NULL}. To see available options,
type \code{geneCategories()}.
}
    \item{interactionTypes}{
A character vector of interaction types to be queried. To query for all
interaction types, skip argument or use \code{NULL}. To see available options,
type \code{interactionTypes()}.
}
}

\details{
Only the first argument (\code{genes}) is required, all other arguments are
optional. When optional arguments are not provided, the query considers all
available possibilities for optional arguments: \code{sourceDatabases},
\code{geneCategories}, and \code{interactionTypes}.
}

\value{
A S4 object of type \code{\linkS4class{rDGIdbResult}}.
}
\references{
Wagner AH, Coffman AC, Ainscough BJ, Spies NC, Skidmore ZL, Campbell KM,
Krysiak K, Pan D, McMichael JF, Eldred JM, Walker JR, Wilson RK, Mardis ER,
Griffith M, Griffith OL. DGIdb 2.0: mining clinically relevant drug-gene
interactions. Nucleic Acids Research. 2016 Jan 4;44(D1):D1036-44.
doi:10.1093/nar/gkv1165.
}
\author{
Thomas Thurnherr <thomas.thurnherr@gmail.com>
}

\seealso{
\code{\linkS4class{rDGIdbResult}}, \code{\link{rDGIdbFilters}},
\code{\link{plotInteractionsBySource}}
}
\examples{
genes <- c("XYZA", "TNF", "IL6", "IL8")
result <- queryDGIdb(genes = genes)

result <- queryDGIdb(genes = genes,
                sourceDatabases = c("MyCancerGenome","TEND","TTD"),
                geneCategories = "CELL SURFACE",
                interactionTypes = c("activator","inhibitor"))
}
\keyword{ ~interface }
\keyword{ ~database }

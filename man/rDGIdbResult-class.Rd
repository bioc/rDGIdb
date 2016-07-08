\name{rDGIdbResult}
\docType{class}
\alias{rDGIdbResult}
\alias{rDGIdbResult-class}
\alias{new,rDGIdbResult-method}
\alias{resultSummary}
\alias{resultSummary,rDGIdbResult-method}
\alias{detailedResults}
\alias{detailedResults,rDGIdbResult-method}
\alias{byGene}
\alias{byGene,rDGIdbResult-method}
\alias{searchTermSummary}
\alias{searchTermSummary,rDGIdbResult-method}

\title{
Result object of a rDGIdb query
}
\description{
S4 helper functions to access DGIdb results, formatted in a similar way as
result tables on the DGIdb web interface.
}
\usage{
\S4method{resultSummary}{rDGIdbResult}(object)

\S4method{detailedResults}{rDGIdbResult}(object)

\S4method{byGene}{rDGIdbResult}(object)

\S4method{searchTermSummary}{rDGIdbResult}(object)
}
%\section{Slots}{
%
%    \code{rDGIdbResult} objects can be created by
%    \code{new(Class = "rDGIdbResult")}. The class contains the following slots:
%
%    \itemize{
%        \item{\code{matchedTerms}:}{A \code{character} vector of input genes
%        for which matches were found in DGIdb.}
%        \item{\code{unmatchedTerms}:}{Drug-gene interactions summarized by the
%        source(s) that reported them.}
%        \item{\code{detailedResults}:}{Search terms matching exactly one gene
%        that has one or more drug interactions.}
%        \item{\code{byGene}:}{Drug interaction count and druggable categories
%        associated with each gene.}
%        \item{\code{searchTermSummary}:}Summary of the attempt to map gene
%        names supplied by the user to gene records in DGIdb.
%        \item{\code{data}:}{Raw \code{data.frame} containing the complete
%        matched term result of the query.}
%    }
%}

\arguments{
    \item{object}{
A \code{rDGIdbResult} object returned by \code{\link{queryDGIdb}}.
    }
}
\value{
\code{resultSummary} returns a \code{data.frame} that summarizes drug-gene
interactions by the source(s) that reported them.

\code{detailedResult} returns a \code{data.frame} with search terms matching
exactly one gene that has one or more drug interactions.

\code{byGene} returns a \code{data.frame} with drug interaction count and
druggable categories associated with each gene.

\code{searchTermSummary} returns a \code{data.frame} that summarizes the
attempt to map gene names supplied by the user to gene records in DGIdb.
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
\code{\link{queryDGIdb}}, \code{\link{rDGIdbFilters}},
\code{\link{plotInteractionsBySource}}
}

\examples{
genes <- c("XYZA", "TNF", "IL6", "IL8")
result <- queryDGIdb(genes = genes)
resultSummary(result)
detailedResults(result)
byGene(result)
searchTermSummary(result)
}

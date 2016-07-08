plotInteractionsBySource <- function(queryResult, ...) {

    if (missing(queryResult)) {
        stop("Input argument 'queryResult' missing.")
    } else if (class(queryResult) != 'rDGIdbResult') {
        stop("Wrong input format, object of class 'rDGIdbResult' required.")
    } else if (nrow(resultSummary(queryResult)) == 0) {
        stop("Input has no interactions to plot.")
    }

    data <- resultSummary(queryResult)
    columnSums <- sort(colSums(data[,3:(ncol(data) - 1)]), decreasing = TRUE)
    op <- par(mar = c(5,13,3,1))
    barplot(columnSums, las = 2, las = 1, horiz = TRUE,
            xlab = "Number of interactions", ... = ...)
    invisible()
}

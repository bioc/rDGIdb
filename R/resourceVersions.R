getWebsite <- function(url) {
    getRequest <- GET(url = url)
    text <- content(getRequest, as = "text", encoding = "ISO-8859-1")
    return(text)
}

resourceVersions <- function() {
    url <- "http://dgidb.genome.wustl.edu/sources/#"
    text <- getWebsite(url)
    array <- strsplit(text, '\n')[[1]]
    interactionId <- grep('class=\'interaction', array)
    name <- sapply(array[interactionId + 4], 
                   function(x) { strsplit(x, '<|>', )[[1]][3] } )
    version <- array[interactionId + 17]
    result <- data.frame(cbind(name, version), stringsAsFactors = FALSE)
    colnames(result) <- c('Name', 'Version')
    rownames(result) <- 1:nrow(result)
    result <- result[order(result$Name),]
    return(result)
}

# drugBankVersion <- function() {
#     url <- "http://www.drugbank.ca/about"
#     text <- getWebsite(url)
#     array <- strsplit(text, '<|>')[[1]]
#     grep('DrugBank Version \\d\\.\\d', array, value = TRUE)
#      
# }


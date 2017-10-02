
queryDGIdb <- function(genes,
                       sourceDatabases = NULL,
                       geneCategories = NULL,
                       interactionTypes = NULL) {
  #,curatedOnly = c(FALSE, TRUE)) {
  
  if (missing(genes)) stop("Need to specify a vector of genes to query.")
  
  if (is.null(genes) || length(genes) == 0 ||
      !is.character(genes) || genes == "") {
    stop("Need to specify a non-empty vector of genes names.")
  }
  
  if (missing(sourceDatabases) | is.null(sourceDatabases)) {
    databases <- NULL
  } else {
    databases <- match.arg(arg = sourceDatabases,
                           choices = sourceDatabases(),
                           several.ok = TRUE)
    databases <- paste(databases, collapse = ",")
  }
  if (missing(geneCategories) | is.null(geneCategories)) {
    categories <- NULL
  } else {
    categories <- match.arg(arg = geneCategories,
                            choices = geneCategories(),
                            several.ok = TRUE)
    categories <- paste(categories, collapse=",")
  }
  if (missing(interactionTypes) | is.null(interactionTypes)) {
    interactions <- NULL
  } else {
    interactions <- match.arg(arg = interactionTypes,
                              choices = interactionTypes(),
                              several.ok = TRUE)
    interactions <- paste(interactions, collapse = ",")
  }
  
  # if (missing(curatedOnly)) {
  #     trustLevel <- NULL
  # } else if (!is.logical(curatedOnly)) {
  #     stop("Argument curatedOnly has to be logical (TRUE or FALSE)")
  # } else if (curatedOnly) {
  #     trustLevel <- "Expert%20cureated"
  # } else {
  #     trustLevel <- NULL
  # }
  
  # Check internet connection
  tryCatch({
    msg <- ""
    r <- GET("http://dgidb.genome.wustl.edu/api/v2/interaction_types.json")
    if (status_code(r) != 200) {
      msg <- "DGIdb service not available."
    }
  }, error = function(err) {
    msg <- "Check internet connection"
  })
  if (msg != "")
    stop(msg)
  
  # Query DGIdb
  cat("Querying DGIDB...")
  queryResult <- queryDgidbPost(genes,
                                interactionSources = databases,
                                geneCategories = categories,
                                interactionTypes = interactions)
  #,trustLevel = trustLevel)
  cat("done!\n")
  
  # Init result class: rDGIdbResult
  result <- new(Class = "rDGIdbResult")
  
  # Set unmatched terms if any
  if (length(queryResult$unmatchedTerms) != 0)
    result <- setUnmatchedTerms(result, queryResult$unmatchedTerms)
  
  # Set matched terms and populate different formats of result tables
  if (!is.null(queryResult$matchedTerms$searchTerm)) {
    
    # Set result data
    result <- setData(result, queryResult$matchedTerms)
    
    # Populate result summary
    result <- setResultSummary(result)
    
    # Populate by gene table
    result <- setByGene(result)
    
    # Populate search term summary
    result <- setSearchTermSummary(result)
    
    #Populate detailed results
    if (nrow(result@resultSummary) > 0)
    result <- setDetailedResults(result)
  }
  
  return(result)
  # End of function queryDGIdb()
}

# Uses httr POST function to query DGIdb. Post instead of get allows
# long list of genes to be queried.
queryDgidbPost <- function(genes, interactionSources, geneCategories,
                           interactionTypes) {
  url <- "http://dgidb.genome.wustl.edu/api/v2/interactions.json"
  body <- list(genes = paste(genes, collapse = ","),
               interaction_sources = interactionSources,
               gene_categories = geneCategories,
               interaction_types = interactionTypes)
  #source_trust_levels = trustLevel)
  body <- body[!sapply(body, is.null)]
  postRequest <- POST(url = url, body = body, encode = 'multipart')
  text <- content(postRequest, as = "text", encoding = "ISO-8859-1")
  if (grepl('error|DOCTYPE', text)) stop("Oops, badly formatted query.")
  if (identical(text, "")) stop("Query response was emtpy.")
  result <- fromJSON(text, simplifyVector = TRUE)
  return(result)
}

sourceDatabases <- function() {
  url <- "http://dgidb.org/api/v2/interaction_sources.json"
  result <- queryAPIget(url)
  return(result)
}

geneCategories <- function() {
  url <- "http://dgidb.org/api/v2/gene_categories.json"
  result <- queryAPIget(url)
  return(result)
}

interactionTypes <- function() {
  url <- "http://dgidb.org/api/v2/interaction_types.json"
  result <- queryAPIget(url)
  return(result)
}

queryAPIget <- function(url) {
  getRequest <- GET(url = url)
  text <- content(getRequest, as = "text", encoding = "ISO-8859-1")
  if (grepl('error|DOCTYPE', text)) stop("Oops, badly formatted query.")
  if (identical(text, "")) stop("Query response was emtpy.")
  result <- fromJSON(text, simplifyVector = TRUE)
  
  return(result)
}


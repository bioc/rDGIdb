
rDGIdbResult <- setClass(
  
  "rDGIdbResult",
  
  slots = c(
    matchedTerms = "character",
    unmatchedTerms = "character",
    resultSummary = "data.frame",
    detailedResults = "data.frame",
    byGene = "data.frame",
    searchTermSummary = "data.frame",
    data = "data.frame" # Matched terms raw data returned from DGIdb
  )
)

setGeneric(name = "setData",
           def = function(theObject, data)
             standardGeneric("setData")
)
setMethod(f = "setData",
          definition = function(theObject, data) {
            theObject@data <- data
            theObject@matchedTerms <- data$searchTerm
            return(theObject)
          })

setGeneric(name = "getData",
           def = function(theObject)
             standardGeneric("getData")
)
setMethod(f = "getData",
          definition = function(theObject) {
            return(theObject@data)
          })

setGeneric(name = "setResultSummary",
           def = function(theObject)
             standardGeneric("setResultSummary")
)
setMethod(f = "setResultSummary", definition = function(theObject) {
    # List of sources for which interactions were found with the genes given
    uniqueSources <- unique(unlist(sapply(
        theObject@data$interactions,
        function(x) {unique(x$source)}, simplify = FALSE)))
    
    # Summarize interactions similar to web interface
    interactionList <- lapply(theObject@data$geneName[
        sapply(theObject@data$interactions, length) != 0],
        getResultSummary, theObject@data, uniqueSources)
    interactionData <- data.frame(do.call(rbind, interactionList),
                                  stringsAsFactors = FALSE)
    if(nrow(interactionData) > 0) {
        # Convert DB and other columns to numeric
        interactionData[,3:(ncol(interactionData) - 2)] <-
            sapply(3:(ncol(interactionData) - 2),
                   function(x) as.numeric(interactionData[,x]))
        interactionData <- interactionData[
            order(as.numeric(interactionData$Score),
                  decreasing = TRUE),]
        rownames(interactionData) <- 1:nrow(interactionData)
    }
    theObject@resultSummary <- interactionData
    return(theObject)
})

setGeneric(name = "setByGene",
           def = function(theObject)
             standardGeneric("setByGene")
)
setMethod(f = "setByGene",
          signature = "rDGIdbResult",
          definition = function(theObject) {
            nrow <- length(theObject@matchedTerms)
            tmp <- data.frame(matrix(nrow = nrow, ncol = 5,
                   dimnames = list(NULL, c('SearchTerm','Gene','GeneName',
                              'DistinctDrugCount','DruggableGeneCategories'))),
                              stringsAsFactors = FALSE)
            tmp[,c('SearchTerm','Gene','GeneName')] <-
              theObject@data[,
                             c('searchTerm','geneName','geneLongName')]
            tmp$DistinctDrugCount <-
              sapply(theObject@data$interactions, function(x) {
                length(unique(x$drugName)) } )
            tmp$DruggableGeneCategories <-
              paste(theObject@data$geneCategories[[1]]$name, collapse = ",")
            tmp <- tmp[order(tmp$SearchTerm),]
            theObject@byGene <- tmp
            return(theObject)
          })

setGeneric(name = "setSearchTermSummary",
           def = function(theObject)
             standardGeneric("setSearchTermSummary")
)
setMethod(f = "setSearchTermSummary",
          signature = "rDGIdbResult",
          definition = function(theObject) {
            
            names <- c('SearchTerm','MatchType','Matches')
            tmp <- data.frame(
              matrix(nrow = 0, ncol = 3, dimnames = list(NULL, names)),
              stringsAsFactors = FALSE)
            if (length(theObject@matchedTerms) > 0) {
              matchedTermSummary <- data.frame(cbind(theObject@matchedTerms,
                   rep('Definite', length(theObject@matchedTerms)),
                   theObject@data$geneName), stringsAsFactors = FALSE)
              colnames(matchedTermSummary) <- names
              tmp <- rbind(tmp, matchedTermSummary)
            }
            # if (length(theObject@unmatchedTerms) > 0) {
            #     noneId <- is.null(theObject@unmatchedTerms$suggestions) |
            #         sapply(theObject@unmatchedTerms$suggestions, length) == 0
            #     unmatchedTermsVector <- unlist(strsplit(
            #       theObject@unmatchedTerms$searchTerm[noneId], split = ', '))
            #     # Unmatched terms without suggestions
            #     unmatchedTermsNone <- cbind(unmatchedTermsVector,
            #         rep('None', length(unmatchedTermsVector)),
            #             rep('None', length(unmatchedTermsVector)))
            #     colnames(unmatchedTermsNone) <- names
            #     tmp <- rbind(tmp, unmatchedTermsNone)
            #     # Unmatched terms with suggestions
            #     unmatchedTermsSuggestions <-
            #         cbind(theObject@unmatchedTerms$searchTerm[!noneId],
            #             rep('Ambiguous', sum(!noneId)),
            #             sapply(theObject@unmatchedTerms$suggestions[!noneId],
            #                 function(x) { paste(x, collapse=',') } ))
            #     colnames(unmatchedTermsSuggestions) <- names
            #     tmp <- rbind(tmp, unmatchedTermsSuggestions)
            # }
            tmp <- tmp[order(tmp$SearchTerm),]
            theObject@searchTermSummary <- tmp
            return(theObject)
          })

setGeneric(name = "setDetailedResults",
           def = function(theObject)
             standardGeneric("setDetailedResults")
)
setMethod(f = "setDetailedResults",
          signature = "rDGIdbResult",
          definition = function(theObject) {
              tmp <- do.call(rbind, apply(theObject@data, 1, function(x) {
                  nrow <- nrow(x$interactions)
                  if (nrow > 0) {
                      data.frame(cbind(rep(x$searchTerm, times = nrow),
                                       rep(x$geneName, times = nrow),
                                       x$interactions[,
                                            c('drugName', 'interactionTypes', 
                                            'sources','pmids')]),
                                 stringsAsFactors = FALSE)
                  }
              }))
              tmp$interactionTypes <- 
                  sapply(tmp$interactionTypes, paste, collapse = ",")
              tmp$sources <- sapply(tmp$sources, paste, collapse = ",")
              tmp$pmids <- sapply(tmp$pmids, paste, collapse = ",")
              colnames(tmp) <-
                  c('SearchTerm', 'Gene', 'Drug', 
                    'InteractionType', 'Source','PMIDs')
              tmp <- tmp[order(tmp$SearchTerm),]
              rownames(tmp) <- 1:nrow(tmp)
              theObject@detailedResults <- tmp
              return(theObject)
          })

setGeneric(name = "setUnmatchedTerms",
           def = function(theObject, unmatchedTerms)
             standardGeneric("setUnmatchedTerms")
)
setMethod(f = "setUnmatchedTerms",
          signature = "rDGIdbResult",
          definition = function(theObject, unmatchedTerms) {
            theObject@unmatchedTerms <- unmatchedTerms
            # data.frame(cbind(unmatchedTerms, 
            #                  rep(NULL, length(unmatchedTerms))), 
            #            stringsAsFactors = FALSE)
            return(theObject)
          })
setGeneric(name = "resultSummary",
           def = function(object)
             standardGeneric("resultSummary")
)
setMethod(f = "resultSummary",
          signature = "rDGIdbResult",
          definition = function(object) {
            return(object@resultSummary)
          })
setGeneric(name = "detailedResults",
           def = function(object)
             standardGeneric("detailedResults")
)
setMethod(f = "detailedResults",
          signature = "rDGIdbResult",
          definition = function(object) {
            return(object@detailedResults)
          })
setGeneric(name = "byGene",
           def = function(object)
             standardGeneric("byGene")
)
setMethod(f = "byGene",
          signature = "rDGIdbResult",
          definition = function(object) {
            return(object@byGene)
          })
setGeneric(name = "searchTermSummary",
           def = function(object)
             standardGeneric("searchTermSummary")
)
setMethod(f = "searchTermSummary",
          signature = "rDGIdbResult",
          definition = function(object) {
            return(object@searchTermSummary)
          })


getResultSummary <- function(gene, output, sources) {
  # Row index with gene interaction information
  idx <- which(output$geneName == gene)
  # Prepare table of interactions (drug vs. interaction DB)
  foundSources <- unique(unlist(output[idx,]$interactions[[1]]$sources))
  result <- #data.frame(cbind(output[idx,]$interactions[[1]]$drugName,
    as.data.frame(t(sapply(output[idx,]$interactions[[1]]$sources,
                           function(x, y) { return(y %in% x) }, foundSources)),
                            stringsAsFactors = FALSE)
  if (length(foundSources) == 1) {result <- t(result)} # result is vector
  dimnames(result) <- 
    list(output[idx,]$interactions[[1]]$drugName, foundSources)
  # Expand matrix to all possible DBs, set multiple occurances of
  # interactions to one, and add gene and drug names
  tmp <- data.frame(matrix(0, nrow = nrow(result), ncol = 4 + length(sources),
            dimnames = list(NULL,c('Gene', 'Drug', sources, 'PMID','Score'))),
                stringsAsFactors = FALSE)
  tmp[,colnames(result)] <- 1*result
  tmp[tmp > 1] <- 1 # Remove double counts
  tmp$Score <- rowSums(tmp) + 
                  sapply(output[idx,]$interactions[[1]]$pmids, length)
  tmp$Gene <- rep(gene, nrow(tmp))
  tmp$Drug <- rownames(result)
  tmp$PMID <- sapply(output[idx,]$interactions[[1]]$pmids, paste, collapse=",")
  # Determine type of interaction
  #resultType <- table(output[idx,]$interactions[[1]]$drugName,
  #                    output[idx,]$interactions[[1]]$interactionType)
  #if (nrow(tmp) == 1) {
  #    tmp$Type <- paste(colnames(resultType), collapse = ",")
  #} else {
  #    listResult <- lapply(split(resultType, seq(nrow(resultType))),
  #                         function(x, names) { names[x>0] },
  #                         colnames(resultType))
  #    tmp$Type <- sapply(listResult, paste, collapse=',')
  #}
  return(as.matrix(tmp))
}

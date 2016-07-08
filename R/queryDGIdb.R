
queryDGIdb <- function(genes,
    sourceDatabases = c("CIViC","CancerCommons","ChEMBL",
        "ClearityFoundationBiomarkers","ClearityFoundationClinicalTrial",
        "DoCM","DrugBank","GuideToPharmacologyInteractions","MyCancerGenome",
        "MyCancerGenomeClinicalTrial","PharmGKB","TALC","TEND","TTD",
        "TdgClinicalTrial"),
    geneCategories = c("abc transporter","b30_2 spry domain",
        "cell surface","clinically actionable","cytochrome p450",
        "dna directed rna polymerase","dna repair","drug metabolism",
        "drug resistance","druggable genome","exchanger",
        "external side of plasma membrane","fibrinogen",
        "g protein coupled receptor","growth factor","histone modification",
        "hormone activity","ion channel","kinase","lipase","lipid kinase",
        "methyl transferase","myotubularin related protein phosphatase",
        "neutral zinc metallopeptidase","nuclear hormone receptor",
        "phosphatidylinositol 3 kinase","phospholipase","protease",
        "protease inhibitor","protein phosphatase","pten family",
        "rna directed dna polymerase","serine threonine kinase",
        "short chain dehydrogenase reductase","thioredoxin",
        "transcription factor binding","transcription factor complex",
        "transporter","tumor suppressor","tyrosine kinase","unknown"),
    interactionTypes = c("activator","adduct","agonist",
        "allosteric modulator","antagonist","antibody","antisense",
        "antisense oligonucleotide","binder","blocker","chaperone",
        "cleavage","cofactor","competitive","immunotherapy","inducer",
        "inhibitor","inhibitory allosteric modulator","inverse agonist",
        "ligand","modulator","multitarget","n/a","negative modulator",
        "other/unknown","partial agonist","partial antagonist",
        "positive allosteric modulator","potentiator","product of",
        "stimulator","suppressor","vaccine")) {
    #,curatedOnly = c(FALSE, TRUE)) {

    if (missing(genes)) stop("Need to specify a vector of genes to query.")

    if (is.null(genes) || length(genes) == 0 ||
        !is.character(genes) || genes == "") {
        stop("Need to specify a non-empty vector of genes names.")
    }

    if (missing(sourceDatabases) |
        all(sourceDatabases() %in% sourceDatabases)) {
        databases <- NULL
    } else {
        databases <- match.arg(arg = sourceDatabases,
                                choices = sourceDatabases(),
                                several.ok = TRUE)
        databases <- paste(databases, collapse = ",")
    }
    if (missing(geneCategories) |
        all(geneCategories() %in% geneCategories)) {
        categories <- NULL
    } else {
        categories <- match.arg(arg = geneCategories,
                                choices = geneCategories(),
                                several.ok = TRUE)
        categories <- paste(categories, collapse=",")
    }
    if (missing(interactionTypes) |
        all(interactionTypes() %in% interactionTypes)) {
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
        r <- GET("http://dgidb.genome.wustl.edu/api/v1/interaction_types.json")
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
    if (!is.null(queryResult$unmatchedTerms$searchTerm))
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
        result <- setDetailedResults(result)
    }

    return(result)
# End of function queryDGIdb()
}

# Uses httr POST function to query DGIdb. Post instead of get allows
# long list of genes to be queried.
queryDgidbPost <- function(genes, interactionSources, geneCategories,
            interactionTypes) {
    url <- "http://dgidb.genome.wustl.edu/api/v1/interactions.json"
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
    sourceDatabases <- c("CIViC","CancerCommons","ChEMBL",
        "ClearityFoundationBiomarkers","ClearityFoundationClinicalTrial","DoCM",
        "DrugBank","GuideToPharmacologyInteractions","MyCancerGenome",
        "MyCancerGenomeClinicalTrial","PharmGKB","TALC","TEND","TTD",
        "TdgClinicalTrial")
    return(sourceDatabases)
}

interactionTypes <- function() {
    interactionTypes <- c("activator","adduct","agonist",
        "allosteric modulator","antagonist","antibody","antisense",
        "antisense oligonucleotide","binder","blocker","chaperone","cleavage",
        "cofactor","competitive","immunotherapy","inducer","inhibitor",
        "inhibitory allosteric modulator","inverse agonist","ligand",
        "modulator","multitarget","n/a","negative modulator","other/unknown",
        "partial agonist","partial antagonist","positive allosteric modulator",
        "potentiator","product of","stimulator","suppressor","vaccine")
    return(interactionTypes)
}

geneCategories <- function() {
    geneCategories <- c("abc transporter","b30_2 spry domain",
        "cell surface","clinically actionable","cytochrome p450",
        "dna directed rna polymerase","dna repair","drug metabolism",
        "drug resistance","druggable genome","exchanger",
        "external side of plasma membrane","fibrinogen",
        "g protein coupled receptor","growth factor","histone modification",
        "hormone activity","ion channel","kinase","lipase","lipid kinase",
        "methyl transferase","myotubularin related protein phosphatase",
        "neutral zinc metallopeptidase","nuclear hormone receptor",
        "phosphatidylinositol 3 kinase","phospholipase","protease",
        "protease inhibitor","protein phosphatase","pten family",
        "rna directed dna polymerase","serine threonine kinase",
        "short chain dehydrogenase reductase","thioredoxin",
        "transcription factor binding","transcription factor complex",
        "transporter","tumor suppressor","tyrosine kinase","unknown")
    return(geneCategories)
}

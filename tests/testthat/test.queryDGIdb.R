library(rDGIdb)

test_that("Genes is a required input", {
    expect_error(queryDGIdb())
    expect_error(queryDGIdb(genes = ''))
    expect_error(queryDGIdb(genes = c(1,2,3)))
})

test_that("Wrong optional arguments", {
    genes <- c("BRAF", "KRAS", "TP53")
    expect_error(queryDGIdb(genes = genes, sourceDatabases = "Alls"))
    expect_error(queryDGIdb(genes = genes, geneCategories = "Alls"))
    expect_error(queryDGIdb(genes = genes, interactionTypes = "Alls"))
})

test_that("Wrong gene names becomes unmatched terms", {
    expect_match(queryDGIdb(genes = c("XYZA", "XYZB"))@unmatchedTerms, c("XYZA, XYZB"))
})

test_that("Query DGIdb and result summary works", {
    result <- queryDGIdb(genes = "BRAF")
    expect_false(is.null(result@resultSummary))
    expect_true(nrow(result@resultSummary) > 0 && ncol(result@resultSummary) > 0)
})

test_that("Returns the right result", {
    result <- queryDgidbPost(genes = c("XYZA", "BRAF"),
        interactionSources = "ChEMBL,MyCancerGenome",
        geneCategories = "clinically actionable",
        interactionTypes = "n/a,inhibitor")
    expect_match(result$unmatchedTerms, "XYZA")
    expect_match(result$matchedTerms$searchTerm, "BRAF")
    expect_is(result$matchedTerms$interactions[[1]], 'data.frame')
})

test_that("Resource versions are returned", {
    versions <- resourceVersions()
    expect_false(is.null(versions))
    expect_true(nrow(versions) > 0 && ncol(versions) == 2)
})

test_that("Helper functions return something", {
    expect_true(is.character(geneCategories()))
    expect_true(is.character(interactionTypes()))
    expect_true(is.character(sourceDatabases()))
})

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
    expect_match(queryDGIdb(genes = c("XYZA", "XYZB"))@unmatchedTerms$searchTerm, "XYZA, XYZB")
})

test_that("Query DGIdb works", {
    result <- queryDGIdb(genes = "BRAF")
    expect_equal(nrow(result@resultSummary), 17)
    expect_equal(result@resultSummary$Score[1], 12)
})

test_that("Returns the right result", {
    result <- queryDgidbPost(genes = c("XYZA", "BRAF"),
        interactionSources = "ChEMBL,MyCancerGenome",
        geneCategories = "clinically actionable",
        interactionTypes = "n/a,inhibitor")
    expect_match(result$unmatchedTerms$searchTerm, "XYZA")
    expect_match(result$matchedTerms$searchTerm, "BRAF")
    expect_is(result$matchedTerms$interactions[[1]], 'data.frame')
})

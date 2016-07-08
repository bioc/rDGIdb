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
        sourceDatabases = c("CIViC","CancerCommons","ChEMBL",
            "ClearityFoundationBiomarkers","ClearityFoundationClinicalTrial",
            "DoCM","DrugBank","GuideToPharmacologyInteractions",
            "MyCancerGenome","MyCancerGenomeClinicalTrial","PharmGKB","TALC",
            "TEND","TTD","TdgClinicalTrial"),
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
            "stimulator","suppressor","vaccine"))
}
\arguments{
    \item{genes}{
A character vector of genes for which drug interactions are queried.
}
    \item{sourceDatabases}{
A character vector of source databases to be queried. To query all
available databases, skip argument or use \code{sourceDatabases()}.
}
    \item{geneCategories}{
A character vector of gene categories to be queried. To query for all
gene categories, skip argument or use \code{geneCategories()}.
}
    \item{interactionTypes}{
A character vector of interaction types to be queried. To query for all
interaction types, skip argument or use \code{interactionTypes()}.
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
Thomas Thurnherr <thomas.thurnherr@bsse.ethz.ch>
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
                geneCategories = "cell surface",
                interactionTypes = c("activator","inhibitor"))
}
\keyword{ ~interface }
\keyword{ ~database }

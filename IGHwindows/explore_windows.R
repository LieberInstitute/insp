library('SummarizedExperiment')
library('recount')

## Load files
projects <- c('sra', 'TCGA', 'SRP012682')
rse <- lapply(projects, function(project) {
    load(paste0('rse_', project, '.Rdata'))
    
    ## Scale counts
    result <- scale_counts(rse, round = FALSE)
    return(result)
})
projects <- names(rse) <- c('SRA', 'TCGA', 'GTEx')

## Define cutoff
cutoff <- 1


presence <- mapply(function(se, project) {
    assays(se)$counts >= cutoff
}, rse, projects, SIMPLIFY = FALSE)

presence_n <- mapply(function(p, project) {
    res <- colSums(p)
    names(res) <- colnames(rse[[project]])
    return(res)
}, presence, projects)
names(presence_n) <- NULL
presence_n <- unlist(presence_n)


## Explore presence by sample
boxplot(presence_n)
head(sort(presence_n, decreasing = TRUE))

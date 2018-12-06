library('SummarizedExperiment')
library('recount')
library('ggplot2')

## Load files
projects <- c('sra', 'TCGA', 'SRP012682')
rse <- lapply(projects, function(project) {
    load(paste0('rse_', project, '.Rdata'))
    
    ## Scale counts
    result <- scale_counts(rse, round = FALSE)
    return(result)
})
projects <- names(rse) <- c('SRA', 'TCGA', 'GTEx')

## To do this automatically, use a for loop around the unique values of
## for(region_of_interest in unique(rowRanges(rse[[1]])$region)) {
##     code below
## }
region_of_interest <- 'IgA-reg1'

## Construct a data.frame that will be easy to use with ggplot2
counts <- mapply(function(se, project) {
    
    ## Subset to region of interest
    se <- subset(se, region == region_of_interest)
    
    ## Assign names so the original code works
    names(se) <- seq_len(nrow(se))
    
    ## Original code, minus width part
    df <- data.frame(
        counts = as.vector(assays(se)$counts),
        region = factor(rep(names(se), each = ncol(se)), levels = names(se)),
        project = rep(project, ncol(se) * nrow(se))
    )
    return(df)
}, rse, projects, SIMPLIFY = FALSE)
counts <- do.call(rbind, counts)
rownames(counts) <- NULL

## Save data for later
save(counts, file = paste0('counts_', region_of_interest, '.Rdata'))

## Save image in a PDF file for the region of interest
## using a height based on the number of boxplots to make
pdf(file = paste0('counts_', region_of_interest, '.pdf'),
    height = round(sum(rowRanges(rse[[1]])$region == region_of_interest) / 10) + 7)
ggplot(counts, aes(y = log10(counts + 1), x = region)) + geom_boxplot() +
    facet_grid(. ~ project) + coord_flip() + theme_bw(base_size = 18)
dev.off()

## This is a big pdf and takes a while to open, alternatively, save to a png (non-vector graphics)
png(file = paste0('counts_', region_of_interest, '.png'),
    height = (round(sum(rowRanges(rse[[1]])$region == region_of_interest) / 10) + 7) * 480 / 7)
ggplot(counts, aes(y = log10(counts + 1), x = region)) + geom_boxplot() +
    facet_grid(. ~ project) + coord_flip() + theme_bw(base_size = 18)
dev.off()

## Exclude 0 counts
png(file = paste0('counts_GT0_', region_of_interest, '.png'),
    height = (round(sum(rowRanges(rse[[1]])$region == region_of_interest) / 10) + 7) * 480 / 7)
ggplot(subset(counts, counts > 0), aes(y = log10(counts + 1), x = region)) + geom_boxplot() +
    facet_grid(. ~ project) + coord_flip() + theme_bw(base_size = 18)
dev.off()

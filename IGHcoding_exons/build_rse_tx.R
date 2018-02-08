## This is an example on how to collapse the data back to the transcripts
## (was meant as an example, but then I ended up collapsing all 3 projects;
## it can still be useful to understand what I did and getting the RPKM
## values)
library('recount')

## Get one of the data sets
load('rse_sra.Rdata')

## We have 1 row per exon and need to group the data by transcript id
rowRanges(rse)

## Custom function for building a 
build_transcript_rse <- function(input_rse) {
    ## Extract counts
    counts <- assays(input_rse)$counts
    ## Group by transcript id
    counts_tx <- split(as.data.frame(counts),
        rowRanges(input_rse)$Transcript_Id)
    ## Sum over the different exons
    counts_final <- do.call(rbind, lapply(counts_tx, colSums))
    
    ## Split exons by transcript id
    exons <- split(rowRanges(input_rse), rowRanges(input_rse)$Transcript_Id)
    txs <- unlist(range(exons))
    txs$Transcript_Id <- names(txs)
    ## Add the length of the exons, for calculating the RPKM
    txs$bp_length <- sum(width(exons))
    output_rse <- SummarizedExperiment(assays = list(counts = counts_final), colData = colData(input_rse), rowRanges = txs)
    return(output_rse)
}

## Build a RSE object with the data grouped by transcript.
## Do this also for GTEx and TCGA loading the appropriate data.
## (see below where I do this for all 3 of them)
rse_tx <- build_transcript_rse(rse)

## compute RPKMs
rpkm_tx <- getRPKM(rse_tx, mapped_var = 'mapped_read_count')

## Or use the following
# getRPKM(rse_tx)
## but you might get some NaN values since some samples have 0 values for all
## 9 transcripts in this data set

## In any way, use whatever you used for your other analysis once you had an
## rse object.



### The following code goes ahead and runs build_transcript_rse() for
## GTEx, sra and TCGA
file_key <- c('sra', 'SRP012682', 'TCGA')
files <- paste0('rse_', file_key, '.Rdata')
for(f in 1:3) {
    message(paste(Sys.time(), 'processing project', file_key[f]))
    load(files[f], verbose = TRUE)
    rse_tx <- build_transcript_rse(rse)
    save(rse_tx, file = paste0('rse_tx_', file_key[f], '.Rdata'))
}

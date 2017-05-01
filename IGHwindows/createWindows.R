## Usage:
# Rscript createWindows.R > createWindows_log.txt 2>&1

library('GenomicRanges')
library('devtools')
library('jaffelab')

coords <- readLines('coords.txt')

IGHranges <- GRanges(
    seqnames = rep('chr14', 25),
    ranges = IRanges(
        start = as.numeric(gsub(',', '', ss(coords, '-', 1))),
        end = as.numeric(gsub(',', '', ss(coords, '-', 2)))
    ),
    locus = c('IgM-coding', 'IgM-reg1', 'IgM-reg2', 'IgG3-coding', 'IgG3-reg1', 'IgG3-reg2', 'IgG1-coding', 'IgG1-reg1', 'IgG1-reg2', 'IgG1-reg3', 'IgA-coding', 'IgA-reg1', 'IgA-reg2', 'IgG2-coding', 'IgG2-reg1', 'IgG2-reg2', 'IgE-coding', 'IgE-reg1', 'IgE-reg2', 'IgG4-coding', 'IgG4-reg1', 'IgG4-reg2', 'IgA2-coding', 'IgA2-reg1', 'IgA2-reg2')
)
names(IGHranges) <- IGHranges$locus
IGHranges$coding <- grepl('coding', IGHranges$locus)

IGHwindows <- slidingWindows(IGHranges, width = 51, step = 25)
nwin <- elementNROWS(IGHwindows)
IGHwindows <- unlist(IGHwindows)
IGHwindows$locus <- rep(IGHranges$locus, nwin)
IGHwindows$coding <- rep(IGHranges$coding, nwin)

## Should all be 51 bp windows
table(width(IGHwindows))

save(IGHwindows, file = 'IGHwindows.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
## Usage:
# Rscript createIGH.R > createIGH_log.txt 2>&1

library('GenomicRanges')
library('devtools')

IGH <- GRanges(
    seqnames = rep('chr14', 10),
    ranges = IRanges(
        start = c(105586890, 105743072, 105771406, 105856219, 105644791, 105626067, 105708666, 105588396, 105601729, 105863197),
        end = c(105935000, 105764000, 105785504, 105863197, 105663500, 105641000, 105718666, 105597400, 105621183, 	
105939755)
    ),
    locus = c('IGH-W', 'IGHG1', 'IGHG3', 'IGHM', 'IGHG2', 'IGHG4', 'IGHA1', 'IGHA2', 'IGHE', 'IGHJ-IGHD'),
    priority = rep(1:4, c(1, 3, 5, 1))
)
names(IGH) <- IGH$locus

save(IGH, file = 'IGH.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
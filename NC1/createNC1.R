## Usage:
# Rscript createNC1.R > createNC1_log.txt 2>&1

library('GenomicRanges')
library('devtools')


nc <- read.table('Regiones-Finales.txt', header = FALSE, stringsAsFactors = FALSE)
colnames(nc) <- c('name', 'chr', 'start', 'end')
nc$chr <- tolower(nc$chr)

NC1 <- GRanges(nc)
names(NC1) <- NC1$name

save(NC1, file = 'NC1.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

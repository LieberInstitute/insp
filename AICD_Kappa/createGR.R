## Usage:
# Rscript createGR.R > createGR_log.txt 2>&1

library('GenomicRanges')
library('sessioninfo')

AICD_dat <- read.csv('AICD_Kappa_regions.csv')
colnames(AICD_dat)[1] <- 'name'

AICD_dat$Start <- as.integer(gsub(',', '', AICD_dat$Start))
AICD_dat$End <- as.integer(gsub(',', '', AICD_dat$End))

AICD <- with(AICD_dat, GRanges(
    seqnames = chromosome,
    ranges = IRanges(
        start = pmin(Start, End),
        end = pmax(Start, End)
    ),
    locus = name)
)
names(AICD) <- AICD$locus

save(AICD, file = 'AICD.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
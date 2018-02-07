## Usage:
# Rscript createIGHcoding_exons.R > createIGHcoding_exons_log.txt 2>&1


## Cargar paquetes necesarios
library('GenomicRanges')
library('devtools')


## Leer las coordenadas que me enviaron y convertirlas en un GRanges object
info <- read.csv('IGH_coding_exons.csv', stringsAsFactors = FALSE)
## Change to integers
info$End <- as.integer(gsub(',', '', info$End))
info$Start <- as.integer(gsub(',', '', info$Start))
## Fix end and starts
s <- mapply(function(s, e) { return(min(c(s, e)))}, info$End, info$Start)
e <- mapply(function(s, e) { return(max(c(s, e)))}, info$End, info$Start)
info$Start <- s
info$End <- e
info$chr <- 'chr14'
info <- GRanges(info)
stopifnot(length(unique(info$Internal_Id)) == length(info))
names(info) <- info$Internal_Id
info

## Rename and save
names(info)
IGHcoding_exons <- info
save(IGHcoding_exons, file = 'IGHcoding_exons.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

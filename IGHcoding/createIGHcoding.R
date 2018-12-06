## Usage:
# Rscript createIGHcoding.R > createIGHcoding_log.txt 2>&1


## Cargar paquetes necesarios
library('GenomicRanges')
library('devtools')


## Leer las coordenadas que me enviaron y convertirlas en un GRanges object
info <- read.table('Coordenadas_IGH_Coding.txt', stringsAsFactors = FALSE)
info$End <- as.integer(gsub('\xca|,', '', info$End))
info$Start <- as.integer(gsub(',', '', info$Start))
info$Chromosome <- 'chr14'
info <- GRanges(info)
info

## Cargar unos datos de recount (solo necesitamos las coordenadas de los disjoint exons)
## DRP000366 solo tiene 1 muestra, por eso lo escogi
library('SummarizedExperiment')
# http://duffel.rail.bio/recount/v2/DRP000366/rse_exon.Rdata
load("/Users/lcollado/Downloads/rse_exon.Rdata")

## Cuantos disjoint exons sobrelapan las regiones que me enviaron
countOverlaps(info, rse_exon)
ov <- findOverlaps(info, rse_exon)

## Calcular el numero de bases codificantes al sumar la longitud de los
## disjoint exons que sobrelapan cada region que me enviaron
width_exons <- sapply(split(subjectHits(ov), queryHits(ov)), function(i) sum(width(rse_exon)[i]))

## Longtiud de las regiones que me mandaron
width_info <- width(info)

## Igualar nombres
names(width_info) <- names(width_exons) <- names(info)

## Longitud de sus regions
width_info
## Longitud de los disjoint exons
width_exons
## Diferencia (ustedes - disjoint exons): numero de bases no codicantes
width_info - width_exons
## Porcentaje no codificante de las regiones que me mandaron
round( (width_info - width_exons ) / width_info * 100, 2)
## Promerio del porcentaje no codificante
mean(round( (width_info - width_exons ) / width_info * 100, 2))

## Cada disjoint exon solo sobrelapa una de las regiones que me enviaron
## (solo checando que si sea asi)
table(table(subjectHits(ov)))

## Obtener unas coordenadas de genes
# http://duffel.rail.bio/recount/v2/DRP000366/rse_gene.Rdata
load("/Users/lcollado/Downloads/rse_gene.Rdata")

## Sus regiones sobrelapan más de 1 gene, por eso tendríamos que usar los
## disjoint exons
countOverlaps(info, rse_gene)

## Guardar los disjoint exons de IGH
IGH_disjoint_exons <- rowRanges(rse_exon)[subjectHits(ov), ]
IGH_disjoint_exons$IGHregion <- names(info)[queryHits(ov)]
length(IGH_disjoint_exons)
IGH_disjoint_exons$gene_id <- names(IGH_disjoint_exons)

names(IGH_disjoint_exons) <- paste0(IGH_disjoint_exons$IGHregion, '_', unlist(sapply(runLength(Rle(queryHits(ov))), seq_len)))

IGH_disjoint_exons
save(IGH_disjoint_exons, file = 'IGH_disjoint_exons.Rdata')

## Tambien guardar en una tabla
df <- data.frame(IGH_disjoint_exons)
df$short_name <- names(IGH_disjoint_exons)
write.table(df, file = 'IGH_disjoint_exons.txt', sep = '\t', row.names = FALSE, quote = FALSE)

## Guardar para despues (de ser necesario)
IGHcoding <- info
save(IGHcoding, file = 'IGHcoding.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
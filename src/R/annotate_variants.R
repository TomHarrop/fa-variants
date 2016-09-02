library(VariantAnnotation)
library(dplyr)

annotation.file <-
  'data/genome/Osativa_323_v7.0.gene_exons.gffread.rRNAremoved.gtf'
variants.file <- 'output/variants/variants_filtered.vcf.gz'


gtf <- rtracklayer::import.gff(annotation.file, format = "gtf",
                               genome = "Osativa_323_v7.0")
txdb <- GenomicFeatures::makeTxDbFromGFF(annotation.file,
                                         format = "gtf")
GenomicFeatures::makeTxDb(gtf)
levels(gtf$type)

vcf <- VariantAnnotation::readVcf(variants.file, genome = "Osativa_323_v7.0")

coding.variants <- locateVariants(vcf, txdb, CodingVariants())

as.data.table(coding.variants) %>%
  group_by(GENEID) %>%
  summarise(variants = n_distinct(QUERYID))



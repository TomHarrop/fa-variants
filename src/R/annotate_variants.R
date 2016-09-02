library(VariantAnnotation)
library(GenomicFeatures) # for transcriptsBy function
library(dplyr)

# genome
genome.file <- "data/genome/Osativa_323_v7.0.fa"
genome <- Rsamtools::FaFile(genome.file)

# annotation
annotation.file <-
  'data/genome/Osativa_323_v7.0.gene_exons.gffread.rRNAremoved.gtf'
gtf <- rtracklayer::import.gff(annotation.file, format = "gtf",
                               genome = "Osativa_323_v7.0")
txdb <- GenomicFeatures::makeTxDbFromGFF(annotation.file,
                                         format = "gtf")

# called variants
variants.file <- 'output/variants/variants_filtered.vcf.gz'
vcf <- VariantAnnotation::readVcf(variants.file, genome = "Osativa_323_v7.0")

# variannts in CDS
coding.variants <- locateVariants(vcf, txdb, CodingVariants())

variant.table <- as.data.frame(coding.variants) %>%
  group_by(GENEID) %>%
  summarise(variants = n_distinct(QUERYID))

variant.table[order(variant.table$variants, decreasing = TRUE),]

# coding changes?
coding <- predictCoding(vcf, txdb, genome)
data.frame(coding, stringsAsFactors = FALSE) %>%
  groupby(GENEID, QUERYID)

nonsyn <- coding[mcols(coding)$CONSEQUENCE == "nonsynonymous"] %>%
  tbl_df() %>%
  group_by(GENEID) %>%
  summarise(n.nonsyn = n_distinct(QUERYID))
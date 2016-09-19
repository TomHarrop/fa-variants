#!/usr/bin/env Rscript

library(VariantAnnotation)
library(GenomicFeatures) # for transcriptsBy function
library(data.table)

command.args <- commandArgs(trailingOnly = TRUE)
# command.args <- c("-v", "output/split_variants/G.variants_filtered.vcf.gz",
#                   "-f", "data/genome/Osativa_323_v7.0.fa",
#                   "-j", "data/genome/Osativa_323_v7.0.gene_exons.gffread.rRNAremoved.gtf",
#                   "-z", "output/cds_variants/G.cds_variants.Rds")
parsed.args <- argparsR::ParseArguments(
  accepted.switches = list(
    `input.vcf` = "-v",
    `input.fa` = "-f",
    `input.gtf` = "-j",
    `other.output` = "-z"),
  command.args)

rutils::GenerateMessage("Find variants in CDS")
rutils::PrintF("input.vcf: %s\n", parsed.args$input.vcf)

# genome
genome <- Rsamtools::FaFile(parsed.args$input.fa)

# annotation
rutils::GenerateMessage("Generating txdb")
txdb <- GenomicFeatures::makeTxDbFromGFF(parsed.args$input.gtf,
                                         format = "gtf")

# called variants
rutils::GenerateMessage("Loading VCF")
vcf <- VariantAnnotation::readVcf(parsed.args$input.vcf,
                                  genome = "Osativa_323_v7.0")

# variannts in CDS
rutils::GenerateMessage("Locating coding variants")
coding.variants <- locateVariants(vcf, txdb, CodingVariants())

# variants per gene
rutils::GenerateMessage("Counting variants per gene")
variant.table <- as.data.table(coding.variants)[, .(
  variants = length(unique(QUERYID))), by = GENEID]

# tidy and sort
setnames(variant.table, "GENEID", "gene")
variant.table[, gene := gsub("\\.MSUv7.*", "", gene)]
setkey(variant.table, gene)

# save output
out.dir <- dirname(parsed.args$other.output)
bn <- gsub("^([[:alpha:]]).*", "\\1", basename(parsed.args$other.output))
log.file <- paste0(out.dir, "/SessionInfo.", bn, ".txt")

rutils::GenerateMessage("Saving output")
rutils::PrintF("  out.dir: %s\n", out.dir)
rutils::PrintF(" log.file: %s\n", log.file)

if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(variant.table, parsed.args$other.output)

s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

rutils::GenerateMessage("Done")

quit(save = "no", status = 0)

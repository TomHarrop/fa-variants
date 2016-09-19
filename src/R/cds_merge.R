#!/usr/bin/env Rscript

library(data.table)

rutils::GenerateMessage("Merge cds_variants results")

command.args <- commandArgs(trailingOnly = TRUE)
# command.args <- c("-y", "test/G.table.Rds",
#                   "-y", "test/B.table.Rds",
#                   "-z", "output/test/cds_variants.Rds")
parsed.args <- argparsR::ParseArguments(
  accepted.switches = list(
    `other.input` = "-y",
    `other.output` = "-z"),
  command.args)

# set up sample names
rutils::GenerateMessage("Generating accession column")
input.names <- sapply(parsed.args$other.input, function(x)
  gsub("^([[:alpha:]]).*", "\\1", basename(x)))
input.names <- plyr::revalue(input.names, c(
  "B" = "Oryza barthii",
  "G" = "Oryza glaberrima",
  "J" = "Oryza sativa japonica",
  "I" = "Oryza sativa indica",
  "R" = "Oryza rufipogon"
), warn_missing = FALSE)
  
# read input
rutils::GenerateMessage("Reading input")
input.tables <- lapply(parsed.args$other.input, readRDS)
names(input.tables) <- input.names

# make long data.table
rutils::GenerateMessage("Binding results")
cds.variants <- rbindlist(input.tables, idcol = "accession")
setkey(cds.variants, gene, accession)

# save output
out.dir <- dirname(parsed.args$other.output)
log.file <- paste0(out.dir, "/SessionInfo.cds_merge.txt")

rutils::GenerateMessage("Saving output")
rutils::PrintF("  out.dir: %s\n", out.dir)
rutils::PrintF(" log.file: %s\n", log.file)

if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

saveRDS(cds.variants, parsed.args$other.output)

s.inf <- rutils::GitSessionInfo()
writeLines(s.inf, log.file)

rutils::GenerateMessage("Done")

quit(save = "no", status = 0)


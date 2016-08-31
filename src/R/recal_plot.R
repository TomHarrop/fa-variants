#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}
printf <- function(...) {
  cat(sprintf(...), file = stderr())
}

GenerateMessage(paste(
  "Plot results of base recalibration to match GATK site",
  "http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr",
  sep = "\n"))

# parse CLI
command.args <- commandArgs(trailingOnly = TRUE)
# command.args <- c("-i", "output/covar_analysis/post_recal_data.table",
#                   "-i", "output/covar_analysis/recal_data.table",
#                   "-o", "test/out.pdf")
parsed.args <- argparsR::ParseArguments(
  accepted.switches = list(
    `output.files` = "-r", `input.files` = "-t"),
  command.args)


# get before and after tables from input.files
after.table.file <- grep(
  "/post_recal_data.table$", parsed.args$input.files, value = TRUE)
before.table.file <- grep(
  "/recal_data.table$", parsed.args$input.files, value = TRUE)
plot.output.file <- grep(
  ".pdf$", parsed.args$output.files, value = TRUE)
log.file <- paste(dirname(plot.output.file),
                  "SessionInfo.recal_plot.txt",
                  sep = "/")

# read tables into R
GenerateMessage("Loading files")
printf("before.table.file: %s\n", before.table.file)
before.table <- gsalib::gsa.read.gatkreport(before.table.file)
printf(" after.table.file: %s\n", after.table.file)
after.table <- gsalib::gsa.read.gatkreport(after.table.file)

# prepare plot data
GenerateMessage("Munging plot data")
events <- c(M = "Mismatches", D = "Deletions", I = "Insertions")
quality.score.table <- rbindlist(list(
  `Original data` = data.table(before.table$RecalTable1),
  `After recalibration` = data.table(after.table$RecalTable1)),
  idcol = "analysis")
quality.score.table[, analysis := factor(analysis, levels = c(
  "Original data", "After recalibration"
))]
quality.score.table[, EventType := factor(
  plyr::revalue(EventType, events), levels = events)]

covariate.table <- rbindlist(list(
  `Original data` = data.table(before.table$RecalTable2),
  `After recalibration` = data.table(after.table$RecalTable2)),
  idcol = "analysis")
covariate.table[, analysis := factor(analysis, levels = c(
  "Original data", "After recalibration"
))]
covariate.table[, EventType := factor(
  plyr::revalue(EventType, events), levels = events)]

cycle.table <- covariate.table[CovariateName == "Cycle"]
cycle.table[, CovariateValue := as.numeric(as.character(CovariateValue))]


nt <- c("A", "T", "C", "G")
di.nt <- paste0(expand.grid(nt, nt)[, 1], expand.grid(nt, nt)[, 2])
dinuc.table <- covariate.table[CovariateName == "Context" &
                                 CovariateValue %in% di.nt]

# set up scale
x <- c(NA, 1e+07, 2e+07, 3e+07, 4e+07, 5e+07, NA)
obs.formatter <- function(x, debug = TRUE) {
  
  # decimal and exponential component
  dc <- as.numeric(gsub("^(.+)e.*", "\\1", x))
  ec <- as.numeric(gsub(".*e\\+0?", "", x))

  # work out label decimals
  dc.c <- as.character(dc)
  dc.c <- dc.c[grep("^[[:digit:]]+\\.[[:digit:]]+", dc.c)]
  dec <- sapply(dc.c, gsub, pattern = "^.*\\.", replacement = "")
  if (length(dec) > 0) {
    dec.no <- max(sapply(dec, length))
  } else {
    dec.no <- 0
  }
  
  # NAs and zeros
  na.lab <- unique(c(which(is.na(dc)), which (is.na(ec))))
  zero.lab <- which(dc == 0)

  # format dc as a string
  dc.str <- format(round(dc, dec.no), nsmall = dec.no, trim = TRUE)
  
  raw.lab <- paste0("\"", dc.str, "\"", "%*%10^", ec)
  raw.lab[na.lab] <- NA
  raw.lab[zero.lab] <- 0

  parse(text = raw.lab)
  }

# 1. Reported quality vs. empirical quality
GenerateMessage("Reported quality vs. empirical quality")
p1 <- ggplot(quality.score.table,
       aes(x = QualityScore, y = EmpiricalQuality,
           size = Observations, colour = EventType)) +
  theme_bw() +
  xlab("Reported quality score") +
  ylab("Empirical quality score") +
  ggtitle("Reported vs. empirical quality") +
  facet_grid(EventType ~ ReadGroup + analysis) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = alpha("black", 0.5)) +
  geom_point(alpha = 0.75) +
  scale_size_area(labels = obs.formatter) +
  scale_color_brewer(palette = "Set1", guide = FALSE)

# 2. Distribution of quality scores
GenerateMessage("Distribution of quality scores")
p2 <- ggplot(quality.score.table,
       aes(x = QualityScore, y = Observations,
           colour = EventType, fill = EventType)) +
  theme_bw() +
  xlab("Reported quality score") +
  ggtitle("Distribution of quality scores") +
  facet_grid(ReadGroup + EventType ~ analysis) +
  geom_density(stat = "identity", alpha = 0.5) +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  scale_fill_brewer(palette = "Set1", guide = FALSE) +
  scale_y_continuous(labels = obs.formatter)

# 3. Residual error by machine cycle
GenerateMessage("Residual error by machine cycle")
p3 <- ggplot(cycle.table,
       aes(x = CovariateValue, y = EmpiricalQuality - QualityScore,
           colour = EventType, size = Observations)) +
  theme_bw() +
  facet_grid(ReadGroup ~ EventType + analysis) +
  xlab("Cycle") +
  ylab(expression("Empirical" - "reported quality")) +
  ggtitle("Residual error by machine cycle") +
  geom_point() +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_size_area(labels = obs.formatter) +
  scale_color_brewer(palette = "Set1", guide = FALSE)

# 4. Residual error by dinucleotide
GenerateMessage("Residual error by dinucleotide")
dinuc.pd <- dinuc.table[, .(
  `Residual error` = mean(EmpiricalQuality - QualityScore)),
  by = .(analysis, ReadGroup, CovariateValue)]
p4 <- ggplot(dinuc.pd,
       aes(x = CovariateValue, y = `Residual error`)) +
  theme_bw() +
  xlab(NULL) + ylab(expression("Empirical" - "reported quality")) +
  ggtitle("Residual error by dinucleotide") +
  facet_grid(ReadGroup ~ analysis) +
  geom_point() +
  geom_point()

# save output
GenerateMessage("Writing plots to file")
printf(" plot.output.file: %s\n", plot.output.file)
outdir <- dirname(plot.output.file)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

pdf(plot.output.file, width = 10, height = 7.5)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# save logs
GenerateMessage("Logging SessionInfo")
printf("         log.file: %s\n", log.file)

sInf <- c(paste("git branch:", system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          paste("date:", date()),
          capture.output(sessionInfo()))
writeLines(sInf, log.file)

GenerateMessage("Done")

quit(save = "no", status = 0)

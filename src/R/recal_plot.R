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

# parse arguments (make this a package later)
ParseCli <- function(command.args) {
  accepted.switches <- c("-i", "-o")
  input.files <- c()
  output.files <- c()
  switch.i <- grep("^-", command.args)
  for (i in switch.i) {
    # check for option starting with "-"
    if (!grepl("^-", command.args[[i]])) {
      break
      
      # check for expected option
    } else if (!command.args[[i]] %in% accepted.switches) {
      stop(paste("Invalid option:", command.args[[i]]))
      
      # check for empty argument
    } else if (i + 1 > length(command.args)) {
      stop(paste("Option", command.args[[i]], "requires an argument"))
    } else if (grepl("^-", command.args[[i + 1]])) {
      stop(paste("Option", command.args[[i]], "requires an argument"))
    } else if (command.args[[i + 1]] == "") {
      stop(paste("Option", command.args[[i]], "requires an argument"))
      
      # do parsing
    } else if (command.args[[i]] == "-i") {
      input.files <- c(input.files, command.args[[i + 1]])
    } else if (command.args[[i]] == "-o") {
      output.files <- c(output.files, command.args[[i + 1]])
    }
  }
  list(input.files = input.files, output.files = output.files)
}

GenerateMessage(paste(
  "Plot results of base recalibration to match GATK site",
  "http://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr",
  sep = "\n"))

# parse CLI
command.args = commandArgs(trailingOnly=TRUE)
parsed.args <- ParseCli(command.args)

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
di.nt <- paste0(expand.grid(nt, nt)[,1], expand.grid(nt, nt)[,2 ])
dinuc.table <- covariate.table[CovariateName == "Context" & 
                                 CovariateValue %in% di.nt]


# 1. Reported quality vs. empirical quality
GenerateMessage("Reported quality vs. empirical quality")
p1 <- ggplot(quality.score.table,
       aes(x = QualityScore, y = EmpiricalQuality,
           size = Observations, colour = EventType)) +
  theme_bw() +
  xlab("Reported quality score") +
  ylab("Empirical quality score") +
  ggtitle("Reported vs. empirical quality") +
  facet_grid(ReadGroup ~ EventType + analysis) +
  geom_abline(slope = 1, intercept = 0, linetype = 2,
              colour = alpha("black", 0.5)) +
  geom_point(alpha = 0.75) +
  scale_size_area() +
  scale_color_brewer(palette = "Set1")

# 2. Distribution of quality scores
GenerateMessage("Distribution of quality scores")
p2 <- ggplot(quality.score.table,
       aes(x = QualityScore, y = Observations,
           colour = EventType, fill = EventType)) +
  theme_bw() +
  xlab("Reported quality score") +
  ggtitle("Distribution of quality scores") +
  facet_grid(ReadGroup ~ EventType + analysis) +
  geom_density(stat = "identity", alpha = 0.5) +
  scale_color_brewer(palette = "Set1", guide = FALSE) +
  scale_fill_brewer(palette = "Set1")

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
  scale_size_area() +
  scale_color_brewer(palette = "Set1")

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
pdf(plot.output.file, width = 10, height= 7.5)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# save logs
GenerateMessage("Logging SessionInfo")
printf("         log.file: %s\n", log.file)

sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
writeLines(sInf, log.file)

GenerateMessage("Done")

quit(save = "no", status = 0)

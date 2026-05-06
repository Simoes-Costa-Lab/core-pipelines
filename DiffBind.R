# ------------------------------------------ #
#   DiffBind Analysis (Reusable Template)    #
# ------------------------------------------ #

## 1. Load packages
library(DiffBind)
library(tidyverse)
library(GenomicRanges)
library(ggplot2)

## 2. Set comparison ID and create output folder
comparison_id <- "H3K27ac"   # Change this for H3K27Lac, etc.
outdir <- file.path("diffbind_DESEQ", comparison_id)
#outdir <- file.path("diffbind_EDGER", comparison_id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Optional: set working directory with BAMs/peaks
setwd("/data/collaborations/joanna_lum/H3K27Lac")

## 3. Helper to read narrowPeak
read_np <- function(path) {
  read_tsv(path, col_names = FALSE) %>%
    dplyr::select(Chr = X1, Start = X2, End = X3) %>%
    distinct() %>%
    GRanges()
}

## 4. Define file paths and sample info
samples <- list(
  list(id="ctrl_rep1",   condition="control", bam="BAM/SI_40576_R1.toGRCh38_nodups.bam",   peak="BAM/SI_40576_R1.toGRCh38_nodups.BAM_peaks.narrowPeak"),
  list(id="ctrl_rep2",   condition="control", bam="BAM/SI_40577_R1.toGRCh38_nodups.bam",   peak="BAM/SI_40577_R1.toGRCh38_nodups.BAM_peaks.narrowPeak"),
  list(id="crispr_rep1", condition="CRISPR",  bam="BAM/SI_40582_R1.toGRCh38_nodups.bam",   peak="BAM/SI_40582_R1.toGRCh38_nodups.BAM_peaks.narrowPeak"),
  list(id="crispr_rep2", condition="CRISPR",  bam="BAM/SI_40583_R1.toGRCh38_nodups.bam",   peak="BAM/SI_40583_R1.toGRCh38_nodups.BAM_peaks.narrowPeak")
)

## 5. Build DBA object
dba_obj <- NULL
for (i in seq_along(samples)) {
  s <- samples[[i]]
  gr <- read_np(s$peak)
  dba_obj <- dba.peakset(
    dba_obj,
    peaks = gr,
    peak.caller = "macs",
    sampID = s$id,
    condition = s$condition,
    replicate = as.integer(gsub("\\D", "", s$id)),
    bamReads = s$bam
  )
}

## 6. Count reads
counts <- dba.count(dba_obj, minOverlap=2, score=DBA_SCORE_TMM_READS_FULL, filter=10)

## 7. Set contrast and analyze
counts <- dba.contrast(counts, group1=counts$masks$CRISPR, group2=counts$masks$control, name1="CRISPR", name2="control")
counts_analyzed <- dba.analyze(counts, method=DBA_DESEQ2)

## 8. Report results
report <- dba.report(counts_analyzed, method=DBA_DESEQ2, th=1, bNormalized=TRUE, DataType=DBA_DATA_FRAME)
write.table(report, file=file.path(outdir, "diffbind_report.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## 9. Save PCA
pdf(file.path(outdir, "PCA.pdf"), width=6, height=6)
dba.plotPCA(counts_analyzed)
dev.off()

## 10. Save MA
df <- as.data.frame(report)
df <- df[, !duplicated(colnames(df))]  # clean up columns

pdf(file.path(outdir, "MA_plot.pdf"), width=6, height=6)
ggplot(df, aes(x = Conc, y = Fold)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Average log2 Concentration",
    y = "log2 Fold Change",
    title = "MA Plot"
  )
dev.off()

## 11. Volcano plot
df$log10FDR <- -log10(df$FDR)
df$Significant <- df$FDR < 0.05 & abs(df$Fold) > 1

pdf(file.path(outdir, "Volcano_plot.pdf"), width=6, height=6)
ggplot(df, aes(x=Fold, y=log10FDR, color=Significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("gray70", "firebrick")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal(base_size=14) +
  labs(title=paste("Volcano Plot:", comparison_id), x="log2 Fold Change", y="-log10(FDR)")
dev.off()

## 12. Extract and save peak sets for GO
depleted <- subset(report, Fold <= -1 & FDR <= 0.05)
enriched <- subset(report, Fold >= 1 & FDR <= 0.05)

write.table(depleted, file=file.path(outdir, "depleted_peaks.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(enriched, file=file.path(outdir, "enriched_peaks.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

saveRDS(depleted, file=file.path(outdir, "depleted_peaks.rds"))
saveRDS(enriched, file=file.path(outdir, "enriched_peaks.rds"))

message("✅ Finished DiffBind for: ", comparison_id)

# Load packages
library(ChIPseeker)
library(GenomicFeatures)
library(clusterProfiler)
library(org.Gg.eg.db)
library(ggplot2)
library(GenomeInfoDb)
library(rtracklayer)
library(ggrepel)
library(stringr)

# Load TxDb
txdb <- loadDb("txdb_galGal6.ncbiRefSeq.sqlite")

# Get narrowPeak files from Peaks/
peak_files <- list.files("Peaks/", pattern = "*_noigg_consensus.bed$", full.names = TRUE)
names(peak_files) <- gsub("_noigg_consensus.bed", "", basename(peak_files))

# Annotate peaks
peak_annot_list <- lapply(peak_files, function(peak_file) {
  # Read peak file
  peaks <- readPeakFile(peak_file)

  # Match seqlevels
  common_seqlevels <- intersect(seqlevels(txdb), seqlevels(peaks))
  seqlevels(txdb, pruning.mode = "coarse") <- common_seqlevels
  seqlevels(peaks, pruning.mode = "coarse") <- common_seqlevels

  # Annotate
  annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
})

# Save annotation tables and plots
for (name in names(peak_annot_list)) {
  peak_anno <- peak_annot_list[[name]]

  # Save annotation table
  write.table(as.data.frame(peak_anno), 
              file = paste0(name, "_annotated_peaks.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
# Genomic annotation summary
anno_stat <- peak_anno@annoStat

# Simplify annotations
anno_stat$Category <- dplyr::case_when(
  anno_stat$Feature %in% c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)") ~ "Promoter (<=3kb)",
  anno_stat$Feature %in% c("1st Intron", "Other Intron") ~ "Intron",
  anno_stat$Feature %in% c("1st Exon", "Other Exon") ~ "Exon",
  anno_stat$Feature %in% c("5' UTR", "3' UTR") ~ "UTR",
  anno_stat$Feature %in% c("Distal Intergenic", "Downstream (<=300)") ~ "Intergenic",
  TRUE ~ "Other"
)

# Summarize simplified categories
# Summarize by Category
df <- anno_stat %>%
  group_by(Category) %>%
  summarise(Frequency = sum(Frequency), .groups = "drop") %>%
  arrange(desc(Category)) %>%
  mutate(
    fraction = Frequency / sum(Frequency),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, -1)),
    label_pos = (ymin + ymax) / 2,
    label = paste0(Category, "\n(", round(100 * fraction), "%)")
  )

# Define custom colors
custom_colors <- c(
  "Promoter (<=3kb)" = "#D36AC2",      # pink/purple
  "Intron" = "#4682B4",        # blue
  "Exon" = "grey70",           # lighter grey
  "UTR" = "black",
  "Intergenic" = "grey90",     # light grey
  "Other" = "darkred"
)

# Base pie chart
# Create pie chart (true pie)
p <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 5, xmin = 0, fill = Category)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  xlim(0, 7) +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  theme_void() +
  theme(legend.position = "none")

# Define tick start and end
label_df <- df %>%
  mutate(
    x_start = 5.05,
    x_end = 5.5,
    y = label_pos,
    label_x = 7  # label goes a bit further than tick end
  )

# Add labels + custom tick lines
final_plot <- p +
  geom_segment(
    data = label_df,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    color = "gray30",
    linewidth = 1
  ) +
  geom_text(
    data = label_df,
    aes(x = label_x, y = y, label = label),
    size = 6.5,
    hjust = 0
  )

svg(paste0(name, "_genomic_distribution.svg"), width = 6.5, height = 6.5, bg = "white")
print(final_plot)
dev.off()
}

# ---- GO ENRICHMENT ANALYSIS ----

for (name in names(peak_annot_list)) {
  cat("Running GO for", name, "\n")
  
  peak_anno_df <- as.data.frame(peak_annot_list[[name]])
  refseq_ids <- na.omit(unique(peak_anno_df$geneId))

  # Convert gene symbols to Entrez
  entrez_ids <- bitr(refseq_ids,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Gg.eg.db)

  if (nrow(entrez_ids) < 10) {
    cat("Too few genes for GO enrichment:", name, "\n")
    next
  }

  # Run GO enrichment (Biological Process)
  ego <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Gg.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)

  # Save table
  write.table(as.data.frame(ego),
              file = paste0(name, "_GO_enrichment.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Plot
  pdf(paste0(name, "_GO_barplot_new.pdf"), width = 8, height = 6)
  p <- barplot(ego, showCategory = 20, title = paste0("GO BP - ", name))
  print(p)
  dev.off()
}

# ---- GO ENRICHMENT ANALYSIS (ONLY PROMOTERS) ----

for (name in names(peak_annot_list)) {
  cat("Running GO for promoters in", name, "\n")
  
  peak_anno_df <- as.data.frame(peak_annot_list[[name]])

  # Keep only promoter annotations (up to 3kb from TSS)
  promoter_df <- peak_anno_df[grepl("Promoter \\(<=1kb\\)|Promoter \\(1-2kb\\)|Promoter \\(2-3kb\\)", peak_anno_df$annotation), ]
  
  # Extract gene IDs
  refseq_ids <- na.omit(unique(promoter_df$geneId))

  # Convert gene IDs to Entrez (check if SYMBOL or REFSEQ format)
  entrez_ids <- bitr(refseq_ids,
                     fromType = "SYMBOL",  # Change if your geneId is not SYMBOL
                     toType = "ENTREZID",
                     OrgDb = org.Gg.eg.db)

  if (nrow(entrez_ids) < 10) {
    cat("Too few promoter-associated genes for GO enrichment:", name, "\n")
    next
  }
  
  all_genes <- keys(org.Gg.eg.db, keytype = "ENTREZID")
  
  # Run GO enrichment (Biological Process)
  ego <- enrichGO(gene = entrez_ids$ENTREZID,
                  universe = all_genes,
                  OrgDb = org.Gg.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)

# Convert to data frame and order
ego_df <- as.data.frame(ego)
ego_df <- head(ego_df, 10)  # top 10 by adjusted p-value
ego_df <- ego_df[order(-ego_df$Count), ] 
ego_df <- ego_df %>%
  mutate(Description_wrapped = str_wrap(str_to_sentence(Description), width = 30))  # adjust width as needed
ego_df$Description_wrapped <- factor(ego_df$Description_wrapped, levels = rev(ego_df$Description_wrapped))  # reverse for horizontal bars

p <- ggplot(ego_df, aes(x = Description_wrapped, y = Count)) +
  geom_bar(stat = "identity", fill = "navy") +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, lineheight = 0.95, margin = margin(r = 10)),
    axis.text.x = element_text(size = 12)
  ) +
  labs(title = NULL, x = NULL, y = "Gene Count")
  
# Plot
pdf(paste0(name, "_GO_barplot_promoters.pdf"), width = 5, height = 4)
print(p)
dev.off()
}

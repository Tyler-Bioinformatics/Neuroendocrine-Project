# Neuroendocrine RNA-seq Metastasis Project
# Main R script version of the report workflow.

library(DESeq2)
library(tidyverse)
library(apeglm)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

counts <- read_tsv(
  "data/GSE98894_raw_counts_GRCh38.p13_NCBI.tsv",
  show_col_types = FALSE
) %>%
  column_to_rownames("GeneID")

meta_raw <- readLines("data/GSE98894_series_matrix.txt")

gsm_line <- meta_raw[grep("^!Sample_geo_accession", meta_raw)]
gsm_ids <- strsplit(gsm_line, "\t")[[1]][-1] |> gsub('"', "", x = _)

char_lines <- meta_raw[grep("!Sample_characteristics", meta_raw)]

site_line <- char_lines[grep("pancreas|intestin|ileum|rect", char_lines, ignore.case = TRUE)][1]
status_line <- char_lines[grep("primary|metast", char_lines, ignore.case = TRUE)][1]

site_vals <- strsplit(site_line, "\t")[[1]][-1] |> gsub('"', "", x = _) |> tolower()
status_vals <- strsplit(status_line, "\t")[[1]][-1] |> gsub('"', "", x = _) |> tolower()

meta_clean <- tibble(
  sample = gsm_ids,
  site_raw = site_vals,
  status_raw = status_vals
) %>%
  mutate(
    site = case_when(
      str_detect(site_raw, "pancreas") ~ "PanNET",
      str_detect(site_raw, "intestin|ileum") ~ "SI-NET",
      str_detect(site_raw, "rect") ~ "RE-NET",
      TRUE ~ NA_character_
    ),
    status = case_when(
      str_detect(status_raw, "primary") ~ "Primary",
      str_detect(status_raw, "metast") ~ "Metastasis",
      TRUE ~ NA_character_
    )
  ) %>%
  select(sample, site, status) %>%
  filter(site %in% c("PanNET", "SI-NET"))

counts <- counts[, meta_clean$sample]
meta_clean <- meta_clean %>% arrange(match(sample, colnames(counts)))
stopifnot(all(meta_clean$sample == colnames(counts)))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_clean,
  design = ~ site + status
)

dds <- dds[rowSums(counts(dds)) > 10, ]

vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

pca_data <- plotPCA(vsd, intgroup = c("site", "status"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_global <- ggplot(pca_data, aes(PC1, PC2, color = site, shape = status)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic() +
  ggtitle("PCA of GEP-NET Samples (VST)")

ggsave("figures/global_pca.png", p_global, width = 7, height = 5, dpi = 300)

dds_pan <- dds[, dds$site == "PanNET"]
dds_pan$site <- droplevels(dds_pan$site)
dds_pan$status <- relevel(dds_pan$status, ref = "Primary")
design(dds_pan) <- ~ status
dds_pan <- DESeq(dds_pan)

res_pan <- results(dds_pan, contrast = c("status", "Metastasis", "Primary"))
res_pan <- lfcShrink(dds_pan, coef = "status_Metastasis_vs_Primary", res = res_pan, type = "apeglm")

res_pan_df <- as.data.frame(res_pan) %>%
  rownames_to_column("GeneID") %>%
  arrange(padj)

write_csv(res_pan_df, "results/PanNET_DE_Metastasis_vs_Primary.csv")

dds_sinet <- dds[, dds$site == "SI-NET"]
dds_sinet$site <- droplevels(dds_sinet$site)
dds_sinet$status <- relevel(dds_sinet$status, ref = "Primary")
design(dds_sinet) <- ~ status
dds_sinet <- DESeq(dds_sinet)

res_sinet <- results(dds_sinet)

res_sinet_df <- as.data.frame(res_sinet) %>%
  rownames_to_column("GeneID") %>%
  arrange(padj)

write_csv(res_sinet_df, "results/SINET_DE_Metastasis_vs_Primary.csv")

padj_cutoff <- 0.05

res_pan_clean <- as.data.frame(res_pan)
res_sinet_clean <- as.data.frame(res_sinet)

res_pan_clean <- res_pan_clean[!is.na(res_pan_clean$padj), ]
res_sinet_clean <- res_sinet_clean[!is.na(res_sinet_clean$padj), ]

pan_sig <- res_pan_clean[res_pan_clean$padj < padj_cutoff, ]
sinet_sig <- res_sinet_clean[res_sinet_clean$padj < padj_cutoff, ]

shared_genes <- intersect(rownames(pan_sig), rownames(sinet_sig))

shared_up_met <- shared_genes[
  pan_sig[shared_genes, "log2FoldChange"] > 0 &
    sinet_sig[shared_genes, "log2FoldChange"] > 0
]

ego_shared <- enrichGO(
  gene = as.character(shared_up_met),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

write_csv(as.data.frame(ego_shared), "results/shared_upregulated_GO_enrichment.csv")

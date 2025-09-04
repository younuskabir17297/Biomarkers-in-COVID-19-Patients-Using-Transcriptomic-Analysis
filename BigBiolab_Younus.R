# ==== STEP 0 (REVISED): Setup work folder, install & load packages ====

# 0a) Create / set working directory
workdir <- "GSE270045_workdir"
if (!dir.exists(workdir)) dir.create(workdir)
setwd(workdir)

# 0b) Helper installer (CRAN + Bioconductor)
pkg_install <- function(pkgs_cran = character(), pkgs_bioc = character()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  need_cran <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
  need_bioc <- pkgs_bioc[!sapply(pkgs_bioc, requireNamespace, quietly = TRUE)]
  
  if (length(need_cran)) install.packages(need_cran, dependencies = TRUE)
  if (length(need_bioc)) BiocManager::install(need_bioc, ask = FALSE, update = TRUE)
}

# 0c) Install everything (explicitly include reactome.db BEFORE ReactomePA)
pkg_install(
  pkgs_cran = c("ggplot2","pheatmap","ggrepel","data.table","igraph"),
  pkgs_bioc = c("GEOquery","DESeq2","limma","sva","org.Hs.eg.db",
                "clusterProfiler","reactome.db","ReactomePA","STRINGdb","apeglm")
)

# 0d) Load libraries (load reactome.db before ReactomePA to avoid the error)
suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(limma)
  library(sva)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(data.table)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(reactome.db)   # <- ensure this is loaded first
  library(ReactomePA)
  library(STRINGdb)
  library(igraph)
})
cat("\nAll packages loaded successfully.\n")

#1
# ==== STEP 1: Download GEO supplementary files ====

geo_id <- "GSE270045"
dest_base <- file.path(getwd(), "geo_supp")
if (!dir.exists(dest_base)) dir.create(dest_base)

# Download (creates a folder GSE270045 under dest_base)
getGEOSuppFiles(geo_id, baseDir = dest_base)

# List all downloaded files (recursive)
supp_files <- list.files(dest_base, recursive = TRUE, full.names = TRUE)
length(supp_files); supp_files

# ==== STEP 2 (Windows-safe): Locate & load the count matrix (no zcat) ====

# You should already have 'supp_files' from STEP 1.
# If your R session was restarted, re-list:
# dest_base <- file.path(getwd(), "geo_supp")
# supp_files <- list.files(dest_base, recursive = TRUE, full.names = TRUE)

# 2a) Filter candidate files: prefer TSV/TXT/CSV (optionally .gz)
cand <- supp_files[grepl("\\.(tsv|txt|csv)(\\.gz)?$", supp_files, ignore.case = TRUE)]
# Remove obviously non-data files
cand <- cand[!grepl("readme|md5|supplementary.*xlsx|\\.pdf$", cand, ignore.case = TRUE)]

# Prefer files whose names hint "count" or "matrix"
pref <- cand[grepl("count|matrix|expr|expression", basename(cand), ignore.case = TRUE)]
if (length(pref) > 0) cand <- pref

# If multiple remain, pick the largest file (often the counts)
if (length(cand) == 0) stop("No suitable .tsv/.txt/.csv(.gz) supplementary file found.")
sizes <- file.info(cand)$size
counts_file <- cand[which.max(sizes)]
message("Using counts file: ", counts_file)

# 2b) Load the counts with base R (works on Windows; handles .gz)
if (grepl("\\.gz$", counts_file, ignore.case = TRUE)) {
  con <- gzfile(counts_file, open = "rt")
  counts_df <- read.delim(con, header = TRUE, check.names = FALSE)
  close(con)
} else {
  counts_df <- read.delim(counts_file, header = TRUE, check.names = FALSE)
}

# 2c) Put gene IDs into rownames (assume first column is gene id)
gene_col <- names(counts_df)[1]
rownames(counts_df) <- counts_df[[gene_col]]
counts_df[[gene_col]] <- NULL

# Keep only numeric columns (in case there are annotation columns)
is_num <- vapply(counts_df, is.numeric, logical(1))
if (!all(is_num)) {
  warning("Some non-numeric columns were found and removed: ",
          paste(names(counts_df)[!is_num], collapse = ", "))
  counts_df <- counts_df[, is_num, drop = FALSE]
}

# Coerce to integer matrix
counts <- as.matrix(counts_df)
mode(counts) <- "numeric"
# If counts are floating (shouldn't be), round:
if (any(counts %% 1 != 0, na.rm = TRUE)) counts <- round(counts)
storage.mode(counts) <- "integer"

# 2d) Quick sanity checks
cat("Matrix dimensions (genes x samples):", paste(dim(counts), collapse = " x "), "\n")
cat("First 5 sample columns:\n"); print(head(colnames(counts), 5))
stopifnot(ncol(counts) >= 30)   # expect ~36 samples
stopifnot(nrow(counts) >= 5000) # expect thousands of genes



# ==== STEP 3: Sample metadata (Group, Batch) ====

# Try to pull pheno from GEO (may take a moment)
gse <- try(getGEO(geo_id, GSEMatrix = TRUE), silent = TRUE)

# Initialize metadata with sample column names
sample_ids <- colnames(counts)
Group <- rep(NA_character_, length(sample_ids))
Batch <- rep(NA_character_, length(sample_ids))

if (!inherits(gse, "try-error")) {
  pheno <- pData(gse[[1]])
  # Map GEO GSM accessions to our columns if names match
  # Try to align by substring matching (robust to different naming)
  guess_map <- function(ids, geo_df){
    map <- rep(NA_character_, length(ids))
    for (i in seq_along(ids)){
      hit <- grep(ids[i], rownames(geo_df), value = TRUE)
      if (length(hit) == 1) map[i] <- hit
    }
    map
  }
  map_rows <- guess_map(sample_ids, pheno)
  
  # Attempt to extract Group from 'characteristics_ch1' fields
  if (!all(is.na(map_rows))) {
    grp_guess <- rep(NA_character_, length(sample_ids))
    bat_guess <- rep(NA_character_, length(sample_ids))
    for (i in seq_along(sample_ids)){
      r <- map_rows[i]
      if (!is.na(r)){
        chars <- as.character(unlist(pheno[r, grep("^characteristics", colnames(pheno), ignore.case = TRUE)]))
        # Search for keywords
        grp_tag <- grep("status|group|condition|disease|phenotype", chars, ignore.case = TRUE, value = TRUE)
        if (length(grp_tag)){
          # take first tag and extract value after ':'
          val <- sub(".*?:\\s*", "", grp_tag[1])
          grp_guess[i] <- val
        }
        bat_tag <- grep("batch|run|lane|flowcell|sequencing", chars, ignore.case = TRUE, value = TRUE)
        if (length(bat_tag)){
          valb <- sub(".*?:\\s*", "", bat_tag[1])
          bat_guess[i] <- valb
        }
      }
    }
    Group <- grp_guess
    Batch <- bat_guess
  }
}

# Fallbacks if parsing failed:
if (all(is.na(Group))) {
  message("Could not parse group from GEO; using fallback (EDIT if needed).")
  # EDIT HERE if your column order is different:
  Group <- c(rep("HealthyControl", 17), rep("LongCOVID", 19))
}
if (all(is.na(Batch))) {
  # If no batch info, just put all in one batch (we’ll handle surrogate batch later)
  Batch <- "Batch1"
}

# Clean up Group labels to factor
Group <- ifelse(grepl("long|covid|patient|case", Group, ignore.case = TRUE), "LongCOVID", Group)
Group <- ifelse(grepl("healthy|control|ctrl", Group, ignore.case = TRUE), "HealthyControl", Group)
Group <- factor(Group, levels = c("HealthyControl","LongCOVID"))
Batch <- factor(Batch)

coldata <- data.frame(Sample = sample_ids, Group = Group, Batch = Batch, row.names = sample_ids)
table(coldata$Group); table(coldata$Batch)


# ==== STEP 4: Filter low-expression genes ====

minCount   <- 10
minSamples <- ceiling(0.25 * ncol(counts))  # expressed in ≥25% samples
keep <- rowSums(counts >= minCount) >= minSamples
filtered_counts <- counts[keep, , drop = FALSE]
dim(filtered_counts)

# ==== STEP 5: DESeq2 dataset & model ====

# Use batch in the model if there is >1 level
if (length(levels(coldata$Batch)) > 1) {
  design_formula <- ~ Batch + Group
} else {
  design_formula <- ~ Group
}

dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData   = coldata,
                              design    = design_formula)
# Run DESeq2
dds <- DESeq(dds)

# ==== STEP 6: PCA before/after batch correction (visual only) ====
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# PCA before
pca1 <- prcomp(t(vsd_mat))
var_expl <- (pca1$sdev^2) / sum(pca1$sdev^2) * 100
df1 <- data.frame(PC1 = pca1$x[,1], PC2 = pca1$x[,2],
                  Group = coldata$Group, Batch = coldata$Batch, Sample = rownames(coldata))
ggplot(df1, aes(PC1, PC2, color = Group, shape = Batch, label = Sample)) +
  geom_point(size = 3) +
  labs(title = "PCA (before batch removal)",
       x = paste0("PC1 (", round(var_expl[1],1), "%)"),
       y = paste0("PC2 (", round(var_expl[2],1), "%)"))

# PCA after removing batch (or surrogate) effects
if (length(levels(coldata$Batch)) > 1) {
  vsd_rm <- removeBatchEffect(vsd_mat, batch = coldata$Batch, design = model.matrix(~ Group, coldata))
} else {
  # If no batch, estimate SVs and remove them for plotting
  mod  <- model.matrix(~ Group, coldata)
  mod0 <- model.matrix(~ 1, coldata)
  sv   <- sva::svaseq(as.matrix(filtered_counts), mod, mod0)
  # remove effect of SVs from vsd
  vsd_rm <- limma::removeBatchEffect(vsd_mat, covariates = sv$sv, design = mod)
}
pca2 <- prcomp(t(vsd_rm))
var_expl2 <- (pca2$sdev^2) / sum(pca2$sdev^2) * 100
df2 <- data.frame(PC1 = pca2$x[,1], PC2 = pca2$x[,2],
                  Group = coldata$Group)
ggplot(df2, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA (after batch/SV removal)",
       x = paste0("PC1 (", round(var_expl2[1],1), "%)"),
       y = paste0("PC2 (", round(var_expl2[2],1), "%)"))

# ==== STEP 7 (FIXED): DE results & significant DEGs ====

# 7a) Get the basic results first (this still uses contrast and is fine)
res <- results(dds, contrast = c("Group","LongCOVID","HealthyControl"), alpha = 0.05)

# 7b) Find the correct coefficient name for apeglm shrinkage
rn <- resultsNames(dds)
print(rn)  # so you can see what's available, e.g. "Group_LongCOVID_vs_HealthyControl"

# Try to auto-detect the coef for the LongCOVID vs HealthyControl comparison
coef_of_interest <- grep("Group.*LongCOVID.*vs.*HealthyControl", rn, value = TRUE)
if (length(coef_of_interest) != 1) {
  stop("Could not auto-detect coef. Choose one from resultsNames(dds), e.g.:
       coef_of_interest <- 'Group_LongCOVID_vs_HealthyControl'")
}

# 7c) LFC shrinkage with apeglm (requires coef=)
# (apeglm is installed in Step 0)
res_shrunk <- lfcShrink(dds, coef = coef_of_interest, type = "apeglm")

# If you ever get an apeglm-related error, you can FALL BACK to 'ashr' or 'normal':
# install.packages('ashr')  # once
# res_shrunk <- lfcShrink(dds, coef = coef_of_interest, type = "ashr")
# OR:
# res_shrunk <- lfcShrink(dds, contrast = c("Group","LongCOVID","HealthyControl"), type = "normal")

# 7d) Order by adjusted p-value and summarize
res_ord <- res_shrunk[order(res_shrunk$padj), ]
summary(res_ord)

# 7e) Significant DEGs: padj < 0.05 & |log2FC| > 1
sig <- subset(as.data.frame(res_ord), !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
nrow(sig)
head(sig)

# ==== STEP 8: Volcano plot ====

volc <- as.data.frame(res_shrunk)
volc$gene <- rownames(volc)
volc$negLog10Padj <- -log10(volc$padj)
volc$negLog10Padj[!is.finite(volc$negLog10Padj)] <- NA

volc$hit <- with(volc, !is.na(padj) & padj < 0.001 & abs(log2FoldChange) > 2)

ggplot(volc, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(alpha = !is.na(padj), color = hit)) +
  scale_color_manual(values = c("grey60","red")) +
  scale_alpha_manual(values = c(0.4, 0.9), guide = "none") +
  geom_hline(yintercept = -log10(0.001), linetype = 2) +
  geom_vline(xintercept = c(-2, 2), linetype = 2) +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 adjusted p-value") +
  theme_minimal()


# ==== STEP 9: Heatmap of top 10 up & top 10 down ====

sig_df <- sig
sig_df$gene <- rownames(sig_df)

top_up   <- head(sig_df[order(sig_df$log2FoldChange, decreasing = TRUE), ], 10)
top_down <- head(sig_df[order(sig_df$log2FoldChange, decreasing = FALSE), ], 10)
top_genes <- unique(c(top_up$gene, top_down$gene))

# Use VST values for heatmap
mat <- assay(vsd)[top_genes, , drop = FALSE]
# Row-scale
mat_scaled <- t(scale(t(mat)))

ann <- data.frame(Group = coldata$Group, Batch = coldata$Batch)
rownames(ann) <- rownames(coldata)

pheatmap(mat_scaled,
         annotation_col = ann,
         show_rownames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize_row = 8,
         main = "Top 10 Up & Top 10 Down (VST, row-scaled)")

# ==== STEP 10 (FIXED): Export results ====

# If sig_df isn't in your environment, recreate it safely from res_ord:
if (!exists("sig_df")) {
  sig_df <- subset(as.data.frame(res_ord),
                   !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
}

# Add gene column from rownames for clean CSVs
sig_df$gene <- rownames(sig_df)

# Pick only the columns that actually exist (apeglm may not have 'stat')
wanted_cols <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
have_cols   <- intersect(wanted_cols, colnames(sig_df))

# Full significant DEGs table
sig_out <- sig_df[, c("gene", have_cols), drop = FALSE]
write.csv(sig_out, "DEGs_LongCOVID_vs_Healthy.csv", row.names = FALSE)

# Top 15 up & down by log2FC (guard against NAs)
ord_up <- order(sig_df$log2FoldChange, decreasing = TRUE, na.last = NA)
ord_dn <- order(sig_df$log2FoldChange, decreasing = FALSE, na.last = NA)

top15_up   <- head(sig_df[ord_up, c("gene","log2FoldChange","padj")], 15)
top15_down <- head(sig_df[ord_dn, c("gene","log2FoldChange","padj")], 15)

top15_tbl <- rbind(
  cbind(Direction = "Up",   top15_up),
  cbind(Direction = "Down", top15_down)
)
write.csv(top15_tbl, "Top15_Up_and_Down_DEGs.csv", row.names = FALSE)

# Quick sanity prints
cat("Wrote:", nrow(sig_out), "significant DEGs to DEGs_LongCOVID_vs_Healthy.csv\n")
cat("Wrote top-15 up and down to Top15_Up_and_Down_DEGs.csv\n")
cat("Columns exported:", paste(colnames(sig_out), collapse = ", "), "\n")

# ==== 11A. Setup + STRING mapping + interactions ====
suppressPackageStartupMessages({
  library(STRINGdb); library(org.Hs.eg.db)
})

stopifnot(exists("res_shrunk"))
if (!exists("sig_df")) {
  sig_df <- subset(as.data.frame(res_shrunk),
                   !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
  sig_df$gene <- rownames(sig_df)
}
sym_pool <- if (exists("sym_pool")) unique(sym_pool) else unique(sig_df$gene)
if (length(sym_pool) < 2) stop("Not enough symbols; relax DEG filters in Step 7.")

options(timeout = 600, download.file.method = "libcurl")
cache_dir <- file.path(getwd(), "string_cache"); if (!dir.exists(cache_dir)) dir.create(cache_dir)

try_map <- function(score){
  db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = score, input_directory = cache_dir)
  mp <- try(db$map(data.frame(gene = sym_pool), "gene", removeUnmappedRows = TRUE), silent = TRUE)
  list(db = db, map = mp)
}

tmp <- try_map(400); if (inherits(tmp$map, "try-error") || is.null(tmp$map) || nrow(tmp$map) == 0) tmp <- try_map(200)
if (inherits(tmp$map, "try-error") || is.null(tmp$map) || nrow(tmp$map) == 0) {
  write.table(sym_pool, "DEG_symbols_for_Cytoscape.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  stop("STRING mapping failed; 'DEG_symbols_for_Cytoscape.txt' written (use Cytoscape stringApp).")
}

string_db <- tmp$db; mapped <- tmp$map; ids <- unique(mapped$STRING_id)
inter <- try(string_db$get_interactions(ids), silent = TRUE)
if (inherits(inter, "try-error") || is.null(inter)) {
  write.table(sym_pool, "DEG_symbols_for_Cytoscape.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  stop("STRING interactions download failed; 'DEG_symbols_for_Cytoscape.txt' written.")
}
cat("11A OK\n")

# ==== 11B. Normalize edges + graph + hubs ====
suppressPackageStartupMessages({ library(igraph) })

# choose the correct endpoint columns automatically
pick_edge_cols <- function(df){
  opts <- list(c("stringId_A","stringId_B"), c("from","to"),
               c("protein1","protein2"), c("preferredName_A","preferredName_B"))
  for (o in opts) if (all(o %in% colnames(df))) return(o)
  stop("Unknown edge columns in STRING table: ", paste(colnames(df), collapse=", "))
}
edge_cols <- pick_edge_cols(inter)

id2sym <- setNames(mapped$gene, mapped$STRING_id)
to_edges_sym <- function(inter, edge_cols, ids, id2sym){
  a <- inter[[edge_cols[1]]]; b <- inter[[edge_cols[2]]]
  if (all(edge_cols %in% c("preferredName_A","preferredName_B"))) {
    es <- data.frame(source = a, target = b)
  } else {
    keep <- (a %in% ids) & (b %in% ids)
    es <- data.frame(source = id2sym[a[keep]], target = id2sym[b[keep]])
    es <- es[!is.na(es$source) & !is.na(es$target), ]
  }
  es <- es[es$source != es$target, , drop = FALSE]
  unique(es)
}

edges_sym <- to_edges_sym(inter, edge_cols, ids, id2sym)
if (!nrow(edges_sym)) stop("No PPI edges among your mapped DEGs at this score.")

g <- graph_from_data_frame(edges_sym, directed = FALSE)
deg_vec <- sort(degree(g), decreasing = TRUE)
hubs <- names(deg_vec)[1:min(10, length(deg_vec))]
cat("Top hubs:", paste(hubs, collapse = ", "), "\n")
cat("11B OK\n")

# ==== 11C. Plot hubs + save files ====
suppressPackageStartupMessages({ library(ggplot2) })

# LFC lookup (works if res_shrunk rownames are SYMBOL; else try SYMBOL->ENSEMBL)
lfc_lookup <- setNames(res_shrunk$log2FoldChange, rownames(res_shrunk))
hub_lfc <- lfc_lookup[hubs]
if (any(is.na(hub_lfc))) {
  sym2ens <- suppressWarnings(try(mapIds(org.Hs.eg.db, keys = hubs, keytype = "SYMBOL", column = "ENSEMBL"), silent = TRUE))
  need <- is.na(hub_lfc); if (!inherits(sym2ens, "try-error")) hub_lfc[need] <- lfc_lookup[na.omit(sym2ens[hubs[need]])]
}

hub_df <- data.frame(Gene = hubs, log2FC = hub_lfc)
hub_df$Gene <- factor(hub_df$Gene, levels = hub_df$Gene[order(hub_df$log2FC, decreasing = TRUE)])

p <- ggplot(hub_df, aes(Gene, log2FC, fill = log2FC > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE"="red","FALSE"="blue"), guide = "none") +
  labs(title = "Hub genes (degree) — log2FC", x = "", y = "log2FC") +
  theme_minimal()
print(p)
ggsave("HubGenes_log2FC.png", width = 7, height = 5, dpi = 150)

# Cytoscape tables
nodes_out <- data.frame(id = unique(c(edges_sym$source, edges_sym$target)))
write.table(nodes_out, "PPI_nodes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(edges_sym, "PPI_edges.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nFiles written in:\n", getwd(), "\n")
print(list.files(pattern = "^(PPI_|HubGenes_)"))
cat("11C OK\n")

#  12A. Setup + map IDs to ENTREZ ----

# Recreate sig_df if needed (uses res_shrunk from Step 7)
if (!exists("sig_df")) {
  stopifnot(exists("res_shrunk"))
  sig_df <- subset(as.data.frame(res_shrunk),
                   !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
  sig_df$gene <- rownames(sig_df)
}
stopifnot(exists("filtered_counts"))
if (nrow(sig_df) == 0) stop("No DEGs at padj<0.05 & |log2FC|>1.")

# Function to detect ID types
detect_id_type <- function(ids) {
  ids <- unique(na.omit(ids))
  if (mean(grepl("^ENSG\\d+(\\.\\d+)?$", ids)) > 0.5) return("ENSEMBL")
  if (mean(grepl("^[0-9]+$", ids)) > 0.5) return("ENTREZID")
  "SYMBOL"
}
strip_version <- function(x) sub("\\.\\d+$", "", x)

# Function to map IDs to ENTREZ
map_to_entrez <- function(ids) {
  ids <- unique(na.omit(ids))
  if (!length(ids)) return(character(0))
  guess <- detect_id_type(ids)
  message("Detected ID type → ", guess)
  if (guess == "SYMBOL") {
    df <- suppressWarnings(tryCatch(bitr(ids, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                                    error = function(e) NULL))
    return(if (is.null(df)) character(0) else unique(df$ENTREZID))
  }
  if (guess == "ENSEMBL") {
    ids2 <- strip_version(ids)
    df <- suppressWarnings(tryCatch(bitr(ids2, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                                    error = function(e) NULL))
    if (!is.null(df) && nrow(df) > 0) return(unique(df$ENTREZID))
    ent <- suppressWarnings(tryCatch(mapIds(org.Hs.eg.db, keys = ids2, keytype = "ENSEMBL", column = "ENTREZID"),
                                     error = function(e) NULL))
    return(unique(na.omit(as.character(ent))))
  }
  unique(ids) # already ENTREZ
}

# Build lists
deg_ids <- sig_df$gene
bg_ids  <- rownames(filtered_counts)
deg_ent <- map_to_entrez(deg_ids)
bg_ent  <- map_to_entrez(bg_ids)

cat("DEGs mapped:", length(deg_ent), "\nBackground mapped:", length(bg_ent), "\n")
cat("12A OK\n")

# KEGG Enrichment
ekegg_res <- tryCatch(
  enrichKEGG(gene = deg_ent, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05),
  error = function(e) { message("KEGG error: ", conditionMessage(e)); NULL }
)

# KEGG results
if (!is.null(ekegg_res) && nrow(as.data.frame(ekegg_res)) > 0) {
  kegg_top10 <- head(as.data.frame(ekegg_res)[order(as.data.frame(ekegg_res)$p.adjust), ], 10)
  write.csv(kegg_top10, "Top10_KEGG_pathways.csv", row.names = FALSE)
  
  p1 <- barplot(ekegg_res, showCategory = 10, title = "KEGG: Top 10 (barplot)")
  ggsave("KEGG_top10_barplot.png", plot = p1, width = 8, height = 6, dpi = 150)
  
  p2 <- dotplot(ekegg_res, showCategory = 10, title = "KEGG: Top 10 (dotplot)")
  ggsave("KEGG_top10_dotplot.png", plot = p2, width = 8, height = 6, dpi = 150)
  
  cat("Saved KEGG CSV + plots.\n")
} else {
  cat("No significant KEGG pathways (or network blocked).\n")
}

cat("DEGs mapped to ENTREZ IDs:", length(deg_ent), "\n")

# Adjust p-value cutoff for KEGG
ekegg_res <- enrichKEGG(
  gene = deg_ent,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

cat("Significant DEGs:", nrow(sig), "\n")
head(deg_ent)

deg_ent <- c("7157", "1956")  # Example ENTREZ IDs for TP53 and EGFR
ekegg_res <- enrichKEGG(
  gene = deg_ent,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

cat("Background genes mapped:", length(bg_ent), "\n")
cat("Mapped DEGs to ENTREZ IDs:", length(deg_ent), "\n")

ekegg_res <- enrichKEGG(
  gene = deg_ent,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

deg_ent <- c("7157", "1956")  # Example ENTREZ IDs for TP53 and EGFR
ekegg_res <- enrichKEGG(
  gene = deg_ent,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

cat("DEGs mapped:", length(deg_ent), "\n")
cat("Background genes mapped:", length(bg_ent), "\n")

sig <- subset(as.data.frame(res_ord), !is.na(padj) & padj < 0.1 & abs(log2FoldChange) > 0.5)
summary(res_ord)

deg_ent <- c("7157", "1956")  # Example ENTREZ IDs for TP53 and EGFR
ekegg_res <- enrichKEGG(
  gene = deg_ent,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Adjust DEG filtering to capture more DEGs
sig <- subset(as.data.frame(res_ord), !is.na(padj) & padj < 0.1 & abs(log2FoldChange) > 0.5)

# Check the number of significant DEGs
nrow(sig)

deg_ids <- sig$gene  # Ensure this column contains gene symbols
deg_ent <- map_to_entrez(deg_ids)
cat("DEGs mapped:", length(deg_ent), "\n")

head(sig$gene)  # Check the first few gene symbols

deg_ids <- na.omit(sig$gene)
deg_ent <- map_to_entrez(deg_ids)
cat("DEGs mapped:", length(deg_ent), "\n")
head(sig$gene)

head(res_ord)
res_ord$gene <- rownames(res_ord)
head(res_ord$gene)

sig <- subset(as.data.frame(res_ord), !is.na(padj) & padj < 0.1 & abs(log2FoldChange) > 0.5)
head(sig)

deg_ent <- mapIds(org.Hs.eg.db, keys = sig$gene, keytype = "SYMBOL", column = "ENTREZID")
deg_ent <- na.omit(deg_ent)  # Remove NA values
cat("DEGs mapped:", length(deg_ent), "\n")

ekegg_res <- enrichKEGG(gene = deg_ent, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

ego_res <- enrichGO(gene = deg_ent, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = bg_ent, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

ere_res <- enrichPathway(gene = deg_ent, universe = bg_ent, pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

# KEGG Enrichment Barplot (Top 10 pathways)
if (!is.null(ekegg_res) && nrow(as.data.frame(ekegg_res)) > 0) {
  p1 <- barplot(ekegg_res, showCategory = 10, title = "KEGG: Top 10 Pathways (Barplot)")
  ggsave("KEGG_top10_barplot.png", plot = p1, width = 8, height = 6, dpi = 150)
  print(p1)  # Display the plot
}

# KEGG Enrichment Dotplot (Top 10 pathways)
if (!is.null(ekegg_res) && nrow(as.data.frame(ekegg_res)) > 0) {
  p2 <- dotplot(ekegg_res, showCategory = 10, title = "KEGG: Top 10 Pathways (Dotplot)")
  ggsave("KEGG_top10_dotplot.png", plot = p2, width = 8, height = 6, dpi = 150)
  print(p2)  # Display the plot
}

# Reactome Enrichment Barplot (Top 10 pathways)
if (!is.null(ere_res) && nrow(as.data.frame(ere_res)) > 0) {
  p1 <- barplot(ere_res, showCategory = 10, title = "Reactome: Top 10 Pathways (Barplot)")
  ggsave("Reactome_top10_barplot.png", plot = p1, width = 8, height = 6, dpi = 150)
  print(p1)  # Display the plot
}

# Reactome Enrichment Dotplot (Top 10 pathways)
if (!is.null(ere_res) && nrow(as.data.frame(ere_res)) > 0) {
  p2 <- dotplot(ere_res, showCategory = 10, title = "Reactome: Top 10 Pathways (Dotplot)")
  ggsave("Reactome_top10_dotplot.png", plot = p2, width = 8, height = 6, dpi = 150)
  print(p2)  # Display the plot
}


# Relax DEGs filtering (for more DEGs)
sig <- subset(as.data.frame(res_ord), !is.na(padj) & padj < 0.1 & abs(log2FoldChange) > 0.5)
deg_ids <- rownames(sig)  # Get DEG gene symbols after relaxing the filters

# Map DEGs to ENTREZ IDs
deg_ent <- mapIds(org.Hs.eg.db, keys = deg_ids, keytype = "SYMBOL", column = "ENTREZID")
deg_ent <- na.omit(deg_ent)  # Remove any NAs

cat("DEGs mapped:", length(deg_ent), "\n")

# Background genes (universe) from filtered counts
bg_ent <- mapIds(org.Hs.eg.db, keys = rownames(filtered_counts), keytype = "SYMBOL", column = "ENTREZID")
bg_ent <- na.omit(bg_ent)  # Remove NAs if present
cat("Number of background genes:", length(bg_ent), "\n")

# Perform GO enrichment (Molecular Function)
ego_res_mf <- enrichGO(gene = deg_ent, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                       ont = "MF", universe = bg_ent, pAdjustMethod = "BH", 
                       pvalueCutoff = 0.1, readable = TRUE)

# Perform GO enrichment (Cellular Component)
ego_res_cc <- enrichGO(gene = deg_ent, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                       ont = "CC", universe = bg_ent, pAdjustMethod = "BH", 
                       pvalueCutoff = 0.1, readable = TRUE)

# Check if any terms are found for Molecular Function
if (!is.null(ego_res_mf) && nrow(as.data.frame(ego_res_mf)) > 0) {
  print(ego_res_mf)
}

# Check if any terms are found for Cellular Component
if (!is.null(ego_res_cc) && nrow(as.data.frame(ego_res_cc)) > 0) {
  print(ego_res_cc)
}

# Barplot for the MF enrichment results
if (!is.null(ego_res_mf) && nrow(as.data.frame(ego_res_mf)) > 0) {
  p1 <- barplot(ego_res_mf, showCategory = 10, title = "GO: Molecular Function (Top Enriched Terms)")
  ggsave("GO_MolecularFunction_top_enriched_terms.png", plot = p1, width = 8, height = 6, dpi = 150)
  print(p1)  # Display the plot
}

# Dotplot for the MF enrichment results
if (!is.null(ego_res_mf) && nrow(as.data.frame(ego_res_mf)) > 0) {
  p2 <- dotplot(ego_res_mf, showCategory = 10, title = "GO: Molecular Function (Top Enriched Terms)")
  ggsave("GO_MolecularFunction_dotplot.png", plot = p2, width = 8, height = 6, dpi = 150)
  print(p2)  # Display the plot
}


# ==== Export top-10 pathways/terms from each enrichment ====

# Helper: write top-10 by adjusted p-value (falls back to raw p if ties)
save_top10 <- function(enrich_obj, out_stem){
  if (is.null(enrich_obj)) return(invisible(NULL))
  df <- try(as.data.frame(enrich_obj), silent = TRUE)
  if (inherits(df, "try-error") || is.null(df) || !nrow(df)) return(invisible(NULL))
  
  # sort by adjusted p-value then raw p-value
  ord <- order(df$p.adjust, df$pvalue, na.last = NA)
  top10 <- head(df[ord, , drop = FALSE], 10)
  
  # CSV
  write.csv(top10, paste0(out_stem, ".csv"), row.names = FALSE)
  
  # XLSX (optional)
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(top10, paste0(out_stem, ".xlsx"))
  }
  
  message("Wrote top 10 → ", out_stem, ".csv",
          if (requireNamespace("writexl", quietly = TRUE)) " and .xlsx" else "")
}

# Ensure GO:BP exists (you already computed MF/CC; BP is often required by briefs)
if (!exists("ego_res_bp")) {
  # Requires: deg_ent (ENTREZ IDs) and bg_ent (ENTREZ universe) defined earlier
  ego_res_bp <- tryCatch(
    enrichGO(gene = deg_ent, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "BP", universe = bg_ent, pAdjustMethod = "BH",
             pvalueCutoff = 0.05, readable = TRUE),
    error = function(e) NULL
  )
}

# Export top-10 tables for each enrichment you ran
save_top10(ekegg_res,  "Top10_KEGG_pathways")
save_top10(ere_res,    "Top10_Reactome_pathways")
save_top10(ego_res_bp, "Top10_GO_BP")
save_top10(ego_res_mf, "Top10_GO_MF")
save_top10(ego_res_cc, "Top10_GO_CC")





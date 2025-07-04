# ========================================================
# Microarray analysis using mouse data (Clariom S arrays)
# ========================================================
# This script performs loading and initial quality control
# of gene expression data from Affymetrix .CEL files.

# -----------------------------
# Load required libraries
# -----------------------------
library(oligo)                           # For reading and processing microarray data
library(clariomsmousetranscriptcluster.db)  # Annotation database for Clariom S Mouse arrays
library(AnnotationDbi)                  # Biological annotation tools
library(arrayQualityMetrics)           # Quality control reporting
library(FactoMineR)                    # Multivariate analysis (PCA, etc.)
library(factoextra)                    # Visualization of FactoMineR results
library(mixOmics)                      # Multivariate methods for omics data
library(limma)                         # Differential expression analysis
library(tidyverse)                     # Data manipulation and visualization
library(ggplot2)                       # Elegant data visualization
library(dplyr)                         # Data manipulation (tidyverse)
library(clusterProfiler)              # Functional enrichment (GO, KEGG, etc.)
library(org.Mm.eg.db)                 # Annotation for Mus musculus (mouse)
library(pathview)                     # KEGG pathway visualization
library(kableExtra)                   # Enhanced HTML tables for reporting
library(ggpubr)                       # Publication-ready ggplot2 plots with easy statistical annotations and multi‐plot layouts
library(ggrepel)                      # Improved text labels for ggplot2

# -------------------------------
# Load CEL files for each group
# -------------------------------
# Read .CEL files for KO and WT mice at 6 and 24 months
ko_mouses_6m <- list.files(path = "6m/CEL/ko", pattern = "[Cc][Ee][Ll]$", full.names = TRUE)
wt_mouses_6m <- list.files(path = "6m/CEL/wt", pattern = "[Cc][Ee][Ll]$", full.names = TRUE)
ko_mouses_24 <- list.files(path = "24m/CEL/ko", pattern = "[Cc][Ee][Ll]$", full.names = TRUE)
wt_mouses_24 <- list.files(path = "24m/CEL/wt", pattern = "[Cc][Ee][Ll]$", full.names = TRUE)

# Combine all .CEL file paths into a single vector
celFiles <- c(ko_mouses_6m, wt_mouses_6m, ko_mouses_24, wt_mouses_24)

# -------------------------------
# Read raw microarray data
# -------------------------------
microarray.raw.data <- read.celfiles(celFiles)

# -------------------------------
# Quality control of raw data
# -------------------------------
# Generate an HTML report with quality metrics and plots
arrayQualityMetrics(
  expressionset = microarray.raw.data,
  outdir = "QC_Report_Raw",       # Output directory
  force = TRUE,                   # Overwrite existing report if any
  do.logtransform = TRUE          # Apply log2 transformation if needed
)
# ===============================
# Data preprocessing
# ===============================

# ------------------------------
# Normalize the raw data using RMA
# ------------------------------
# RMA includes background correction, quantile normalization, and summarization
microarray.processed.data <- rma(microarray.raw.data)

# -------------------------------
# Boxplot of normalized intensities
# -------------------------------
# Asigna nombres personalizados a las muestras antes de hacer el boxplot
sampleID <- c("ko1_6m", "ko2_6m", "ko3_6m", "ko4_6m",
              "wt1_6m", "wt2_6m", "wt3_6m", "wt4_6m",
              "ko1_24m", "ko2_24m", "ko3_24m", "ko4_24m",
              "wt1_24m", "wt2_24m", "wt3_24m")

# Establece los nombres en el objeto ExpressionSet
colnames(microarray.processed.data) <- sampleID

# Ahora haz el boxplot con los nombres nuevos y más limpios
boxplot(microarray.processed.data,
        col = rainbow(length(sampleID)),
        las = 2,                      # Rotar etiquetas del eje x (vertical)
        ylab = "Fluorescence (R.U.)",
        main = "Normalized Expression (RMA)")

# --------------------------------------
# Extract normalized expression matrix
# --------------------------------------
expression.level <- exprs(microarray.processed.data)

# Inspect matrix
head(expression.level)
dim(expression.level)

# ------------------------------------------------
# Filter probes with low average expression (noise)
# ------------------------------------------------
mean_expression <- rowMeans(expression.level)
filtered_expression <- expression.level[mean_expression > 5, ]

# -----------------------------
# Rename columns with sample IDs
# -----------------------------
sampleID <- c("ko1_6m", "ko2_6m", "ko3_6m", "ko4_6m",
              "wt1_6m", "wt2_6m", "wt3_6m", "wt4_6m",
              "ko1_24m", "ko2_24m", "ko3_24m", "ko4_24m",
              "wt1_24m", "wt2_24m", "wt3_24m")

colnames(filtered_expression) <- sampleID

# Preview filtered matrix
head(filtered_expression)

# ==================================================
# Exploratory Analysis: Correlation between replicates
# ==================================================

# Compare KO replicates at 24 months, if you find it necessary you can do that to other conditions
x_min <- min(filtered_expression[,"ko1_24m"], na.rm = TRUE)
x_max <- max(filtered_expression[,"ko1_24m"], na.rm = TRUE)
y_min <- min(filtered_expression[,"ko2_24m"], na.rm = TRUE)
y_max <- max(filtered_expression[,"ko2_24m"], na.rm = TRUE)

# Scatter plot between two biological replicates
ggplot(filtered_expression, aes(x = ko1_24m, y = ko2_24m)) +
  geom_point(color = "grey", size = 0.5) +                      # puntos
  geom_smooth(method = "lm", se = FALSE, color = "red") +       # línea de regresión
  xlim(x_min, x_max) +                                          # límites eje X
  ylim(y_min, y_max) +                                          # límites eje Y
  stat_cor(                                                     # añade r de Pearson
    aes(label = ..r.label..),
    method  = "pearson",
    label.x = x_min,
    label.y = y_max
  ) +
  labs(
    x = "KO 24m - Réplica 1",
    y = "KO 24m - Réplica 2",
    title = "Correlación entre réplicas KO a 24 meses"
  ) +
  theme_minimal()

# ==============================
# Principal Component Analysis (PCA)
# ==============================

# Create a data frame for PCA: one row per sample
# Transpose expression matrix so that genes are columns and samples are rows
pca.gene.expression <- data.frame(colnames(filtered_expression), t(filtered_expression))
colnames(pca.gene.expression)[1] <- "Sample"  

# Run PCA using FactoMineR
# scale.unit = TRUE standardizes variables (genes) before analysis
# quali.sup = 1 tells PCA to treat the "Sample" column as a supplementary qualitative variable
res.pca <- PCA(pca.gene.expression, graph = FALSE, scale.unit = TRUE, quali.sup = 1)

# Hierarchical clustering on principal components (HCPC)
res.hcpc <- HCPC(res.pca, graph = FALSE)

# ------------------------
# Plot dendrogram of clusters
# ------------------------
fviz_dend(res.hcpc, 
          k = 4,                          # Number of clusters to cut
          cex = 0.75,                     # Font size
          palette = "jco",               # Color palette
          rect = TRUE,                   # Draw rectangles around clusters
          rect_fill = TRUE,              # Fill rectangles with color
          rect_border = "jco",           # Border color
          type = "rectangle",            # Dendrogram shape
          labels_track_height = 1400     # Adjust label alignment
)

# ==============================
# Partial Least Squares Discriminant Analysis (PLS-DA)
# ==============================

# Prepare the expression matrix for mixOmics:
# Transpose again: samples as rows, genes as columns
X <- t(filtered_expression)

# Create the class label vector (Y) as a factor
# Make sure it matches the sample order
Y <- factor(c(rep("KO_6m", 4),
              rep("WT_6m", 4),
              rep("KO_24m", 4),
              rep("WT_24m", 3)))  # Adjust if sample count changes

# Run PLS-DA with 3 components
plsda.res <- plsda(X, Y, ncomp = 3)

# Plot sample distribution along components 1 and 2
plotIndiv(plsda.res,
          comp = c(1, 2),           # Which components to display
          group = Y,                # Grouping variable
          legend = TRUE,            # Add legend
          ellipse = TRUE,           # Confidence ellipses
          title = "PLS-DA: Condition Separation")

# ==============================
# Identify important genes (loadings)
# ==============================

# For component 1: extract genes with highest absolute loadings. You can change comp = 2 or 3 to see other components.
genes_importantes_comp1 <- selectVar(plsda.res, comp = 1)$value
head(genes_importantes_comp1[order(-abs(genes_importantes_comp1[,1])), , drop = FALSE], 20) 


# ====================================
# Differential Expression: Linear Model Setup
# ====================================

# ------------------------------------------
# Design matrix for 4 experimental conditions
# ------------------------------------------
# We define the experimental groups as:
# Group 1: KO at 6 months
# Group 2: WT at 6 months
# Group 3: KO at 24 months
# Group 4: WT at 24 months
experimental.design <- model.matrix(~ -1 + factor(c(
  1, 1, 1, 1,   # KO_6m (4 samples)
  2, 2, 2, 2,   # WT_6m (4 samples)
  3, 3, 3, 3,   # KO_24m (4 samples)
  4, 4, 4       # WT_24m (3 samples)
)))

# Rename columns to reflect biological conditions
colnames(experimental.design) <- c("ko_6m", "wt_6m", "ko_24m", "wt_24m")

# Check the dimensions of expression matrix and design matrix
dim(filtered_expression)
dim(experimental.design)

# ------------------------------------------
# Fit the linear model to each gene (limma)
# ------------------------------------------
linear.fit <- lmFit(filtered_expression, experimental.design)

# ------------------------------------------
# Define contrasts between experimental groups
# ------------------------------------------
# Each contrast represents a specific biological comparison:
# - KO vs WT at 6 months
# - KO 24m vs KO 6m (aging effect in KO)
# - WT 24m vs WT 6m (aging effect in WT)
# - KO vs WT at 24 months
contrast.matrix <- makeContrasts(
  ko_6m - wt_6m,
  ko_24m - ko_6m,
  wt_24m - wt_6m,
  ko_24m - wt_24m,
  levels = c("ko_6m", "wt_6m", "ko_24m", "wt_24m")
)

# Check matrix dimensions for safety
dim(linear.fit$coefficients)
dim(contrast.matrix)

# Apply the contrasts to the fitted model
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)

# Use empirical Bayes shrinkage to improve statistical inference
contrast.results <- eBayes(contrast.linear.fit)

# ========================================================
# Create the background "universe" for enrichment analysis
# ========================================================

# We define the full gene universe as all probes that passed the expression filter
universe_ids <- rownames(filtered_expression)

# Optional: save IDs for external use or documentation
writeLines(universe_ids, "universe_ids.txt")

# ---------------------------------------------
# Map PROBEIDs to ENTREZIDs using annotation DB
# ---------------------------------------------
# Many probes (~10,000) may not have valid ENTREZID mappings
# despite attempts with other databases/identifiers
universe_mapped <- AnnotationDbi::select(
  clariomsmousetranscriptcluster.db,
  keys = universe_ids,        # Probe IDs from expression matrix
  columns = "ENTREZID",       # Column we want to retrieve
  keytype = "PROBEID"         # Type of input keys
)

# ---------------------------------------------
# Clean mapping: remove NAs and duplicates
# ---------------------------------------------
universe_clean <- universe_mapped %>%
  filter(!is.na(ENTREZID)) %>%              # Remove missing mappings
  distinct(PROBEID, .keep_all = TRUE)       # Keep only one entry per probe

# Extract unique ENTREZ IDs (gene universe)
universe_entrez <- unique(universe_clean$ENTREZID)

# ---------------------------------------------
# Summary of mapping success
# ---------------------------------------------
cat("Total filtered probes:", length(universe_ids), "\n")
cat("Mapped probes with ENTREZID:", nrow(universe_clean), "\n")
cat("Unique ENTREZIDs:", length(universe_entrez), "\n")


# ===================================================
# Get differential expression results for contrast 1
# (KO at 6 months vs WT at 6 months)
# ===================================================

ko_6m_vs_wt_6m <- topTable(contrast.results, number = Inf, coef = 1, sort.by = "none")

# Calculate -log10(p-value) for volcano plot
ko_6m_vs_wt_6m$negLogP <- -log10(ko_6m_vs_wt_6m$P.Value)

# Classify genes based on significance thresholds
ko_6m_vs_wt_6m$Significance <- "NS"  # Not Significant
ko_6m_vs_wt_6m$Significance[ko_6m_vs_wt_6m$logFC > 1.5 & ko_6m_vs_wt_6m$P.Value < 0.05] <- "Up"
ko_6m_vs_wt_6m$Significance[ko_6m_vs_wt_6m$logFC < -1.5 & ko_6m_vs_wt_6m$P.Value < 0.05] <- "Down"

################################## Volcano plot ##############################################

# Add PROBEID and ENTREZID to the differential expression table
de_df <- ko_6m_vs_wt_6m %>%
  tibble::rownames_to_column("PROBEID") %>%        # Turn row names into a PROBEID column
  inner_join(
    universe_clean %>% select(PROBEID, ENTREZID),
    by = "PROBEID"                                 # Keep only probes present in the universe
  )

# Select top 5 up- and top 5 down-regulated genes by logFC
top_up <- de_df %>%
  filter(Significance == "Up") %>%                 # Only up-regulated genes
  arrange(desc(logFC)) %>%                         # Sort by descending log fold-change
  slice_head(n = 5)                                # Take the top 5

top_down <- de_df %>%
  filter(Significance == "Down") %>%               # Only down-regulated genes
  arrange(logFC) %>%                               # Sort by ascending log fold-change
  slice_head(n = 5)                                # Take the top 5

top10 <- bind_rows(top_up, top_down)               # Combine the 10 selected genes

# Map ENTREZIDs to gene SYMBOLs for labeling
symbol_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = top10$ENTREZID,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
top10 <- left_join(top10, symbol_map, by = "ENTREZID")

# Create volcano plot highlighting and labeling the top 5 Up and top 5 Down genes
ggplot(de_df, aes(x = logFC, y = negLogP, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +                          # Plot all data points
  scale_color_manual(
    values = c("Up" = "red",                                    # Custom color for up-regulated
               "Down" = "blue",                                 # Custom color for down-regulated
               "NS" = "gray")                                   # Custom color for non-significant
  ) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +  # Add vertical lines at fold‐change cutoffs
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # Add horizontal line at p-value cutoff
  geom_text_repel(
    data          = top10,                                      # Only label the selected top genes
    aes(label     = SYMBOL),                                    # Use gene SYMBOL for labels
    show.legend   = FALSE,
    size          = 5,                                          # Increase text size
    fontface      = "bold",                                     # Bold font
    box.padding   = 0.8,                                        # Space around text boxes
    point.padding = 0.5,                                        # Minimum distance from points
    nudge_y       = 1,                                          # Nudge labels upward
    direction     = "both",                                     # Allow movement in x and y directions
    force         = 1.5,                                        # Repulsion force between labels
    max.overlaps  = Inf,                                        # Do not drop any labels
    segment.size  = 0.3,                                        # Thickness of connector lines
    segment.color = "black"                                     # Color of connector lines
  ) +
  scale_x_continuous(expand = expansion(add = c(0.3, 0.3))) +    # Add padding on x-axis so labels aren’t cut off
  scale_y_continuous(expand = expansion(add = c(0.8, 1))) +      # Add padding on y-axis so labels aren’t cut off
  coord_cartesian(clip = "off") +                                # Allow labels to extend beyond plotting area
  theme_minimal() +
  theme(
    plot.margin = unit(c(1, 3, 1, 1), "lines")                  # Extra margin to accommodate labels
  ) +
  labs(
    title = "Volcano Plot: KO 6m vs WT 6m",                     # Plot title
    x     = "log2 Fold Change",                                 # X-axis label
    y     = "-log10(p-value)",                                  # Y-axis label
    color = "Regulation"                                        # Legend title
  )

## Gene and pathway enrichment analysis for KO 6m vs WT 6m

# Full contrast table sorted by logFC
ko_6m_vs_wt_6m <- topTable(contrast.results, number = Inf, coef = 1, sort.by = "logFC")
fold.change.ko_6m_vs_wt_6m <- ko_6m_vs_wt_6m$logFC
genes.ids.ko_6m_vs_wt_6m <- rownames(ko_6m_vs_wt_6m)  # PROBEIDs

# Identify activated and repressed genes
activated_ids_ko_6m_wt_6m <- genes.ids.ko_6m_vs_wt_6m[fold.change.ko_6m_vs_wt_6m > 1.5]
repressed_ids_ko_6m_wt_6m <- genes.ids.ko_6m_vs_wt_6m[fold.change.ko_6m_vs_wt_6m < -1.5]

# Map PROBEID to ENTREZID for activated genes
mapped_activated <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                          keys = activated_ids_ko_6m_wt_6m,
                                          columns = "ENTREZID",
                                          keytype = "PROBEID")

entrez_ids_activated <- mapped_activated %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE) %>%
  pull(ENTREZID)

# GO enrichment: Biological Process
ego_activated <- enrichGO(
  gene          = entrez_ids_activated,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",  
  pAdjustMethod = "BH",
  universe      = universe_entrez,  
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)
barplot(ego_activated, showCategory = 10, title = "GO:BP - Activated KO 6m vs WT 6m")

# GO enrichment: Molecular Function
ego_activated_MF <- enrichGO(
  gene          = entrez_ids_activated,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",  
  pAdjustMethod = "BH",
  universe      = universe_entrez,  
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)
barplot(ego_activated_MF, showCategory = 10, title = "GO:MF - Activated KO 6m vs WT 6m")

# GO enrichment: Cellular Component
ego_activated_CC <- enrichGO(
  gene          = entrez_ids_activated,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",  
  pAdjustMethod = "BH",
  universe      = universe_entrez,  
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)
barplot(ego_activated_CC, showCategory = 10, title = "GO:CC - Activated KO 6m vs WT 6m")

# KEGG enrichment (activated genes)
ekegg_activated <- enrichKEGG(
  gene          = entrez_ids_activated,
  organism      = "mmu",  # Mus musculus
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)
barplot(ekegg_activated, showCategory = 15, title = "KEGG - Activated KO 6m vs WT 6m")

# Convert KEGG result to data frame
df.activated.enrich.kegg <- as.data.frame(ekegg_activated)
head(df.activated.enrich.kegg)

# Map all DE probe IDs to ENTREZID
mapped_all <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                    keys = genes.ids.ko_6m_vs_wt_6m,
                                    columns = "ENTREZID",
                                    keytype = "PROBEID")

# Assign ENTREZ IDs to logFC values
names(fold.change.ko_6m_vs_wt_6m) <- mapped_all$ENTREZID[match(rownames(ko_6m_vs_wt_6m), mapped_all$PROBEID)]
fold.change.ko_6m_vs_wt_6m <- fold.change.ko_6m_vs_wt_6m[!is.na(names(fold.change.ko_6m_vs_wt_6m))]

# Set non-activated genes to 0 for clean pathview
fold.change.ko_6m_vs_wt_6m_cleaned <- fold.change.ko_6m_vs_wt_6m
fold.change.ko_6m_vs_wt_6m_cleaned[!(names(fold.change.ko_6m_vs_wt_6m_cleaned) %in% entrez_ids_activated)] <- 0

# Visualize one KEGG pathway (adjust pathway ID as needed)
pathview(
  gene.data  = fold.change.ko_6m_vs_wt_6m_cleaned,
  pathway.id = "mmu05033",  
  species    = "mmu",
  gene.idtype = "ENTREZID",
  limit = list(gene = max(abs(fold.change.ko_6m_vs_wt_6m_cleaned)), cpd = 1)
)

# REPPRESSED GENES

# Map repressed PROBEIDs to ENTREZIDs
mapped_repressed <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                          keys = repressed_ids_ko_6m_wt_6m,
                                          columns = "ENTREZID",
                                          keytype = "PROBEID")

entrez_ids_repressed <- mapped_repressed %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE) %>%
  pull(ENTREZID)

# plots are not displayed when there are not enriched terms found
# GO enrichment: Biological Process
ego_repressed <- enrichGO(
  gene          = entrez_ids_repressed,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# GO enrichment: Molecular Function
ego_repressed_MF <- enrichGO(
  gene          = entrez_ids_repressed,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# GO enrichment: Cellular Component
ego_repressed_CC <- enrichGO(
  gene          = entrez_ids_repressed,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# KEGG enrichment (repressed genes)
ekegg_repressed <- enrichKEGG(
  gene          = entrez_ids_repressed,
  organism      = "mmu",  # Mus musculus
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)

# ===================================================
# Get differential expression results for contrast 4
# (KO at 24 months vs WT at 24 months)
# ===================================================

ko_24m_vs_wt_24m <- topTable(contrast.results, number = Inf, coef = 4, sort.by = "none")

# Calculate -log10(p-value) for volcano plot
ko_24m_vs_wt_24m$negLogP <- -log10(ko_24m_vs_wt_24m$P.Value)

# Classify genes by fold change and p-value
ko_24m_vs_wt_24m$Significance <- "NS"
ko_24m_vs_wt_24m$Significance[ko_24m_vs_wt_24m$logFC > 1.5 & ko_24m_vs_wt_24m$P.Value < 0.05] <- "Up"
ko_24m_vs_wt_24m$Significance[ko_24m_vs_wt_24m$logFC < -1.5 & ko_24m_vs_wt_24m$P.Value < 0.05] <- "Down"

######################################## Volcano plot####################################################
# Convert row names to a PROBEID column and left‐join to get ENTREZIDs
de24m_df <- ko_24m_vs_wt_24m %>%
  tibble::rownames_to_column("PROBEID") %>%            # row names → PROBEID
  left_join(
    universe_clean %>% select(PROBEID, ENTREZID),     # keep all DE probes
    by = "PROBEID"
  )

# Select top 5 up-regulated and top 5 down-regulated genes (with valid ENTREZID)
top_up24m <- de24m_df %>%
  filter(Significance == "Up", !is.na(ENTREZID)) %>%   # only Up with ENTREZID
  arrange(desc(logFC)) %>%                             # highest logFC first
  slice_head(n = 5)                                    # take top 5

top_down24m <- de24m_df %>%
  filter(Significance == "Down", !is.na(ENTREZID)) %>% # only Down with ENTREZID
  arrange(logFC) %>%                                   # most negative logFC first
  slice_head(n = 5)                                    # take top 5

top10_24m <- bind_rows(top_up24m, top_down24m)         # combine for labeling

# Map ENTREZID to gene SYMBOL for nicer labels
symbol_map24m <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = top10_24m$ENTREZID,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
top10_24m <- left_join(top10_24m, symbol_map24m, by = "ENTREZID")

# Build the volcano plot and add labels for the 10 genes
ggplot(de24m_df, aes(x = logFC, y = negLogP, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +                              # Plot all points
  scale_color_manual(values = c("Up"   = "red",                     # Custom color for up-regulated
                                "Down" = "blue",                    # Custom color for down-regulated
                                "NS"   = "gray")) +                # Custom color for non-significant
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +      # Add vertical lines at FC cutoffs
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +      # Add horizontal line at p-value cutoff
  geom_text_repel(
    data          = top10_24m,                                      # Only label the selected top genes
    aes(label     = SYMBOL),                                        # Use gene SYMBOL as label
    show.legend   = FALSE,
    size          = 5,                                              # Increase text size
    fontface      = "bold",                                         # Bold font for emphasis
    box.padding   = 0.8,                                            # Padding around text box
    point.padding = 0.5,                                            # Minimum distance from point
    nudge_y       = 1,                                              # Nudge labels upward
    direction     = "both",                                         # Allow movement in x and y directions
    force         = 1.5,                                            # Repulsion force between labels
    max.overlaps  = Inf,                                            # Do not drop any labels
    segment.size  = 0.3,                                            # Thickness of connector lines
    segment.color = "black"                                         # Color of connector lines
  ) +
  scale_x_continuous(expand = expansion(add = c(0.3, 0.3))) +        # Add padding on x-axis
  scale_y_continuous(expand = expansion(add = c(0.8, 1))) +         # Add padding on y-axis
  coord_cartesian(clip = "off") +                                   # Allow labels outside plot area
  theme_minimal() +
  theme(
    plot.margin = unit(c(1, 3, 1, 1), "lines")                     # Extra margin for labels
  ) +
  labs(
    title = "Volcano Plot: KO 24m vs WT 24m",# Plot title
    x     = "log2 Fold Change",                                     # X-axis label
    y     = "-log10(p-value)",                                      # Y-axis label
    color = "Regulation"                                            # Legend title
  )

## Functional enrichment: KO 24m vs WT 24m

# Sorted result table for fold change
ko_24m_vs_wt_24m <- topTable(contrast.results, number = 29129, coef = 4, sort.by = "logFC")
fold.change.ko_24m_vs_wt_24m <- ko_24m_vs_wt_24m$logFC
genes.ids.ko_24m_vs_wt_24m <- rownames(ko_24m_vs_wt_24m)

# Identify upregulated and downregulated genes
activated_ids_ko_24m_wt_24m <- genes.ids.ko_24m_vs_wt_24m[fold.change.ko_24m_vs_wt_24m > 1.5]
repressed_ids_ko_24m_wt_24m <- genes.ids.ko_24m_vs_wt_24m[fold.change.ko_24m_vs_wt_24m < -1.5]

# Map Transcript Cluster IDs → ENTREZIDs for upregulated genes
mapped_24m <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                    keys = activated_ids_ko_24m_wt_24m,
                                    columns = "ENTREZID",
                                    keytype = "PROBEID")

mapped_unique_24m <- mapped_24m %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE)

entrez_ids_activated_24m <- unique(mapped_unique_24m$ENTREZID)

# GO: Biological Process
ego_ko_24m_vs_wt_24m <- enrichGO(
  gene         = entrez_ids_activated_24m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)
barplot(ego_ko_24m_vs_wt_24m, showCategory = 10, title = "GO:BP - KO 24m vs WT 24m Activated")

# GO: Molecular Function
ego_ko_24m_vs_wt_24m_MF <- enrichGO(
  gene         = entrez_ids_activated_24m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)
barplot(ego_ko_24m_vs_wt_24m_MF, showCategory = 10, title = "GO:MF - KO 24m vs WT 24m Activated")

# GO: Cellular Component
ego_ko_24m_vs_wt_24m_CC <- enrichGO(
  gene         = entrez_ids_activated_24m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "CC",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)

# KEGG enrichment
ekegg_activated_24_24 <- enrichKEGG(
  gene          = entrez_ids_activated_24m,
  organism      = "mmu",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)
barplot(ekegg_activated_24_24, showCategory = 10, title = "KEGG - KO 24m vs WT 24m Activated")

df.activated_24_24.enrich.kegg <- as.data.frame(ekegg_activated_24_24)

# Assign ENTREZ IDs to all logFC values
mapped_all_24_24 <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                          keys = genes.ids.ko_24m_vs_wt_24m,
                                          columns = "ENTREZID",
                                          keytype = "PROBEID")

names(fold.change.ko_24m_vs_wt_24m) <- mapped_all$ENTREZID[match(rownames(ko_24m_vs_wt_24m), mapped_all$PROBEID)]
fold.change.ko_24m_vs_wt_24m <- fold.change.ko_24m_vs_wt_24m[!is.na(names(fold.change.ko_24m_vs_wt_24m))]

# Set non-upregulated genes to 0
fold.change.ko_24m_vs_wt_24m_cleaned <- fold.change.ko_24m_vs_wt_24m
fold.change.ko_24m_vs_wt_24m_cleaned[!(names(fold.change.ko_24m_vs_wt_24m_cleaned) %in% entrez_ids_activated_24m)] <- 0

# Visualize selected KEGG pathway
pathview(
  gene.data  = fold.change.ko_24m_vs_wt_24m_cleaned,
  pathway.id = "mmu04657",  
  species    = "mmu",
  gene.idtype = "ENTREZID",
  limit = list(gene = max(abs(fold.change.ko_24m_vs_wt_24m_cleaned)), cpd = 1)
)

# --- REPRESSED GENES ---

# Map downregulated probes
mapped_24m <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                    keys = repressed_ids_ko_24m_wt_24m,
                                    columns = "ENTREZID",
                                    keytype = "PROBEID")

mapped_unique_24m <- mapped_24m %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE)

entrez_ids_repressed_24m <- unique(mapped_unique_24m$ENTREZID)

# GO enrichment: BP, MF, CC
ego_ko_24m_vs_wt_24m_repressed <- enrichGO(
  gene         = entrez_ids_repressed_24m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)
barplot(ego_ko_24m_vs_wt_24m_repressed, showCategory = 10, title = "GO:BP - KO 24m vs WT 24m Repressed")

ego_ko_24m_vs_wt_24m_repressed_MF <- enrichGO(
  gene         = entrez_ids_repressed_24m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)
barplot(ego_ko_24m_vs_wt_24m_repressed_MF, showCategory = 10, title="GO:MF - KO 24m vs WT 24m Repressed")

ego_ko_24m_vs_wt_24m_repressed_CC <- enrichGO(
  gene         = entrez_ids_repressed_24m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "CC",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)

# KEGG enrichment for repressed genes in KO 24m vs WT 24m
ekegg_repressed_24_24 <- enrichKEGG(
  gene          = entrez_ids_repressed_24m,       # ENTREZ IDs of repressed genes
  organism      = "mmu",                          # Mouse organism code
  universe      = universe_entrez,                # Background universe
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)

# Barplot for enriched KEGG pathways
barplot(ekegg_repressed_24_24, showCategory = 10, title = "KEGG - KO 24m vs WT 24m Repressed")


# --------- Prepare data for Pathview visualization ----------

# Convert the enrichKEGG result into a data frame
df.repressed_24_24.enrich.kegg <- as.data.frame(ekegg_repressed_24_24)

# Display the first rows to inspect enriched pathways
head(df.repressed_24_24.enrich.kegg)

# Prepare the fold-change data for Pathview:
#    Make a copy of the original fold-change vector
fold.change.ko_24m_vs_wt_24m_cleaned <- fold.change.ko_24m_vs_wt_24m

# Set to zero all genes NOT in the list of repressed genes
#    names(...) returns the ENTREZIDs of each element in the vector
#    Any name not found in entrez_ids_repressed_24m is replaced with 0
fold.change.ko_24m_vs_wt_24m_cleaned[
  !(names(fold.change.ko_24m_vs_wt_24m_cleaned) %in% entrez_ids_repressed_24m)
] <- 0

# Call Pathview to visualize the "mmu04062" KEGG pathway
pathview(
  gene.data   = fold.change.ko_24m_vs_wt_24m_cleaned,             # log2FC data named by ENTREZID
  pathway.id  = "mmu04062",                                       # KEGG pathway ID to plot
  species     = "mmu",                                            # organism: Mus musculus
  gene.idtype = "ENTREZID",                                       # identifier type in gene.data
  limit       = list(
    gene = max(abs(fold.change.ko_24m_vs_wt_24m_cleaned)),        # max range for color scaling
    cpd  = 1                                                      # compound limit (not applicable here)
  )
)

# ===================================================
# Get differential expression results for contrast 2
# (KO at 24 months vs KO at 6 months)
# ===================================================

ko_24m_vs_ko_6m <- topTable(contrast.results, number = Inf, coef = 2, sort.by = "none")

# Calculate -log10(p-value)
ko_24m_vs_ko_6m$negLogP <- -log10(ko_24m_vs_ko_6m$P.Value)

# Classify genes based on significance thresholds
ko_24m_vs_ko_6m$Significance <- "NS"
ko_24m_vs_ko_6m$Significance[ko_24m_vs_ko_6m$logFC > 1.5 & ko_24m_vs_ko_6m$P.Value < 0.05] <- "Up"
ko_24m_vs_ko_6m$Significance[ko_24m_vs_ko_6m$logFC < -1.5 & ko_24m_vs_ko_6m$P.Value < 0.05] <- "Down"

# Volcano plot using ggplot2
# Annotate the DE results with PROBEID → ENTREZID (keeping all DE probes)
de_24v6_df <- ko_24m_vs_ko_6m %>%
  tibble::rownames_to_column("PROBEID") %>%           # convert rownames to PROBEID
  left_join(
    universe_clean %>% select(PROBEID, ENTREZID),    # add ENTREZID where available
    by = "PROBEID"
  )

# Pick the top 5 most up-regulated and top 5 most down-regulated genes
top_up_24v6 <- de_24v6_df %>%
  filter(Significance == "Up", !is.na(ENTREZID)) %>%  # only true Ups with ENTREZID
  arrange(desc(logFC)) %>%                            # highest logFC first
  slice_head(n = 5)                                   # take top 5

top_down_24v6 <- de_24v6_df %>%
  filter(Significance == "Down", !is.na(ENTREZID)) %>%# only true Downs with ENTREZID
  arrange(logFC) %>%                                  # most negative logFC first
  slice_head(n = 5)                                   # take top 5

top10_24v6 <- bind_rows(top_up_24v6, top_down_24v6)  # combine for labeling

# Map those ENTREZIDs to gene SYMBOLs for nicer labels
symbol_map_24v6 <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = top10_24v6$ENTREZID,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
top10_24v6 <- left_join(top10_24v6, symbol_map_24v6, by = "ENTREZID")

# Build the volcano plot, highlighting and labeling the top 10 genes
ggplot(de_24v6_df, aes(x = logFC, y = negLogP, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +                         # plot all points
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "NS" = "gray")    # custom colors
  ) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +  # FC cutoffs
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p-value cutoff
  geom_text_repel(
    data          = top10_24v6,                                 # label only top 10
    aes(label     = SYMBOL),                                    # use gene SYMBOL
    show.legend   = FALSE,                                      # exclude labels from legend
    size          = 5,                                          # larger text
    fontface      = "bold",                                     # bold font
    box.padding   = 0.8,                                        # padding around text
    point.padding = 0.5,                                        # min distance from points
    nudge_y       = 1,                                          # nudge labels upward
    direction     = "both",                                     # allow x & y adjustment
    force         = 1.5,                                        # stronger repulsion
    max.overlaps  = Inf,                                        # do not drop any labels
    segment.size  = 0.3,                                        # connector line width
    segment.color = "black"                                     # connector line color
  ) +
  scale_x_continuous(expand = expansion(add = c(0.3, 0.3))) +    # horizontal padding
  scale_y_continuous(expand = expansion(add = c(0.8, 1))) +      # vertical padding
  coord_cartesian(clip = "off") +                                # allow labels outside panel
  theme_minimal() +
  theme(
    plot.margin = unit(c(1, 3, 1, 1), "lines")                  # extra margin for labels
  ) +
  labs(
    title = "Volcano Plot: KO 24m vs KO 6m", 
    x     = "log2 Fold Change",
    y     = "-log10(p-value)",
    color = "Regulation"
  )
## Gene and pathway enrichment: KO 24m vs KO 6m

# Full results sorted by logFC
ko_24m_vs_ko_6m <- topTable(contrast.results, number = 29129, coef = 2, sort.by = "logFC")
fold.change.ko_24m_vs_ko_6m <- ko_24m_vs_ko_6m$logFC
genes.ids.ko_24m_vs_ko_6m <- rownames(ko_24m_vs_ko_6m)

# Identify activated and repressed genes
activated_ids_ko_24m_ko_6m <- genes.ids.ko_24m_vs_ko_6m[fold.change.ko_24m_vs_ko_6m > 1.5]
repressed_ids_ko_24m_ko_6m <- genes.ids.ko_24m_vs_ko_6m[fold.change.ko_24m_vs_ko_6m < -1.5]

# Map PROBEID (Transcript Cluster ID) → ENTREZID
mapped_ko_24m_ko_6m <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                             keys = activated_ids_ko_24m_ko_6m,
                                             columns = "ENTREZID",
                                             keytype = "PROBEID")

# Remove NAs and duplicates
mapped_unique_ko_24m_ko_6m <- mapped_ko_24m_ko_6m %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE)

# Get unique ENTREZ IDs
entrez_ids_activated_ko_24m_ko_6m <- unique(mapped_unique_ko_24m_ko_6m$ENTREZID)

# GO enrichment (Biological Process)
ego_ko_24m_vs_ko_6m <- enrichGO(
  gene         = entrez_ids_activated_ko_24m_ko_6m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)
barplot(ego_ko_24m_vs_ko_6m, showCategory = 10, title = "GO:BP - KO 24m vs KO 6m Activated")

# GO enrichment (Molecular Function)
ego_ko_24m_vs_ko_6m_MF <- enrichGO(
  gene         = entrez_ids_activated_ko_24m_ko_6m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)
barplot(ego_ko_24m_vs_ko_6m_MF, showCategory = 10, title = "GO:MF - KO 24m vs KO 6m Activated")

# GO enrichment (Cellular Component)
ego_ko_24m_vs_ko_6m_CC <- enrichGO(
  gene         = entrez_ids_activated_ko_24m_ko_6m,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "CC",
  pAdjustMethod= "BH",
  universe     = universe_entrez,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE
)
barplot(ego_ko_24m_vs_ko_6m_CC, showCategory = 10, title = "GO:CC - KO 24m vs KO 6m Activated")

# KEGG pathway enrichment
ekegg_activated_24_6_KO <- enrichKEGG(
  gene          = entrez_ids_activated_ko_24m_ko_6m,
  organism      = "mmu",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)
barplot(ekegg_activated_24_6_KO, showCategory = 10, title = "KEGG - KO 24m vs KO 6m Activated")

# Prepare data for pathview
mapped_all_24_6_KO <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                            keys = genes.ids.ko_24m_vs_ko_6m,
                                            columns = "ENTREZID",
                                            keytype = "PROBEID")
names(fold.change.ko_24m_vs_ko_6m) <- mapped_all$ENTREZID[match(rownames(ko_24m_vs_ko_6m), mapped_all$PROBEID)]
fold.change.ko_24m_vs_ko_6m <- fold.change.ko_24m_vs_ko_6m[!is.na(names(fold.change.ko_24m_vs_ko_6m))]

# Clean non-activated genes
fold.change.ko_24m_vs_ko_6m_cleaned <- fold.change.ko_24m_vs_ko_6m
fold.change.ko_24m_vs_ko_6m_cleaned[!(names(fold.change.ko_24m_vs_ko_6m_cleaned) %in% entrez_ids_activated_ko_24m_ko_6m)] <- 0

# Pathway visualization
pathview(
  gene.data  = fold.change.ko_24m_vs_ko_6m_cleaned,
  pathway.id = "mmu04613",  
  species    = "mmu",
  gene.idtype = "ENTREZID",
  limit = list(gene = max(abs(fold.change.ko_24m_vs_ko_6m_cleaned)), cpd = 1)
)

# --------------------
# REPPRESSED GENES
# --------------------

# Map repressed genes (Transcript Cluster IDs / PROBEID) to ENTREZ IDs
mapped_repressed_ko_24m_ko_6m <- AnnotationDbi::select(
  clariomsmousetranscriptcluster.db,
  keys = repressed_ids_ko_24m_ko_6m,
  columns = "ENTREZID",
  keytype = "PROBEID"
)

# Remove NAs and duplicate PROBEID entries
mapped_unique_repressed_ko_24m_ko_6m <- mapped_repressed_ko_24m_ko_6m %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE)

# Extract unique ENTREZ IDs for enrichment analysis
entrez_ids_repressed_ko_24m_ko_6m <- unique(mapped_unique_repressed_ko_24m_ko_6m$ENTREZID)


# ----------------------------
# GO Enrichment (BP, MF, CC)
# ----------------------------

# Biological Process (BP)
ego_repressed_ko_24m_ko_6m <- enrichGO(
  gene          = entrez_ids_repressed_ko_24m_ko_6m,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

barplot(ego_repressed_ko_24m_ko_6m,
        showCategory = 10,
        title = "GO:BP - KO 24m vs KO 6m Repressed")


# Molecular Function (MF)
ego_repressed_ko_24m_ko_6m_MF <- enrichGO(
  gene          = entrez_ids_repressed_ko_24m_ko_6m,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

barplot(ego_repressed_ko_24m_ko_6m_MF,
        showCategory = 10,
        title = "GO:MF - KO 24m vs KO 6m Repressed")


# Cellular Component (CC)
ego_repressed_ko_24m_ko_6m_CC <- enrichGO(
  gene          = entrez_ids_repressed_ko_24m_ko_6m,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

barplot(ego_repressed_ko_24m_ko_6m_CC,
        showCategory = 10,
        title = "GO:CC - KO 24m vs KO 6m Repressed")


# -----------------
# KEGG Enrichment
# -----------------

ekegg_repressed_ko_24_6 <- enrichKEGG(
  gene          = entrez_ids_repressed_ko_24m_ko_6m,
  organism      = "mmu",                   # Mus musculus
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)

barplot(ekegg_repressed_ko_24_6,
        showCategory = 10,
        title = "KEGG - KO 24m vs KO 6m Repressed")


# -------------------------
# KEGG Pathway Visualization (Repressed genes only)
# -------------------------

# Convert the enrichKEGG result for the KO 24m vs KO 6m repressed genes into a data frame
df.ekegg_repressed_ko_24_6 <- as.data.frame(ekegg_repressed_ko_24_6)

# Display the first few rows to inspect the enriched KEGG pathways
head(df.ekegg_repressed_ko_24_6)

# Prepare the fold-change data for Pathview:
#    Make a copy of the original fold-change vector for the 24m vs 6m KO comparison
fold.change.ko_24m_vs_ko_6m_cleaned <- fold.change.ko_24m_vs_ko_6m

# Set to zero all genes NOT in the list of repressed ENTREZ IDs
#    names(...) returns the ENTREZIDs for each element in the vector
#    Any name not found in entrez_ids_repressed_ko_24m_ko_6m is replaced with 0
fold.change.ko_24m_vs_ko_6m_cleaned[
  !(names(fold.change.ko_24m_vs_ko_6m_cleaned) %in% entrez_ids_repressed_ko_24m_ko_6m)
] <- 0

# Visualize the "mmu04721" KEGG pathway 
pathview(
  gene.data   = fold.change.ko_24m_vs_ko_6m_cleaned,           # log2FC data named by ENTREZID
  pathway.id  = "mmu04721",                                    # KEGG pathway ID to plot
  species     = "mmu",                                         # organism: Mus musculus
  gene.idtype = "ENTREZID",                                    # identifier type used in gene.data
  limit       = list(
    gene = max(abs(fold.change.ko_24m_vs_ko_6m_cleaned)),      # maximum range for color scaling
    cpd  = 1                                                   # compound limit (not used here)
  )
)

# Visualize another pathway, "mmu04724" 
pathview(
  gene.data   = fold.change.ko_24m_vs_ko_6m_cleaned,           # reuse the same filtered data
  pathway.id  = "mmu04724",                                    # another KEGG pathway ID
  species     = "mmu",                                         # Mus musculus
  gene.idtype = "ENTREZID",                                    # ENTREZID identifiers
  limit       = list(
    gene = max(abs(fold.change.ko_24m_vs_ko_6m_cleaned)),      # same color scale limit
    cpd  = 1                                                   # compound limit
  )
)


# ===================================================
# Get differential expression results for contrast 3
# (WT at 24 months vs WT at 6 months)
# ===================================================
## WT 24m vs WT 6m

# Full contrast table
wt_24m_vs_wt_6m <- topTable(contrast.results, number = Inf, coef = 3, sort.by = "none")

# Compute -log10(p-value) for volcano plot
wt_24m_vs_wt_6m$negLogP <- -log10(wt_24m_vs_wt_6m$P.Value)

# Classify genes by significance
wt_24m_vs_wt_6m$Significance <- "NS"  # Not Significant
wt_24m_vs_wt_6m$Significance[wt_24m_vs_wt_6m$logFC > 1.5 & wt_24m_vs_wt_6m$P.Value < 0.05] <- "Up"
wt_24m_vs_wt_6m$Significance[wt_24m_vs_wt_6m$logFC < -1.5 & wt_24m_vs_wt_6m$P.Value < 0.05] <- "Down"

# Volcano plot using ggplot2
# Convert row names to PROBEID and join to get ENTREZID (keeping all DE probes)
de_wt_df <- wt_24m_vs_wt_6m %>%
  tibble::rownames_to_column("PROBEID") %>%
  left_join(
    universe_clean %>% select(PROBEID, ENTREZID),
    by = "PROBEID"
  )

# Select top 5 up- and top 5 down-regulated genes (with valid ENTREZID)
top_up_wt <- de_wt_df %>%
  filter(Significance == "Up", !is.na(ENTREZID)) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 5)

top_down_wt <- de_wt_df %>%
  filter(Significance == "Down", !is.na(ENTREZID)) %>%
  arrange(logFC) %>%
  slice_head(n = 5)

top10_wt <- bind_rows(top_up_wt, top_down_wt)

# Map ENTREZIDs to SYMBOLs for labeling
symbol_map_wt <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = top10_wt$ENTREZID,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
top10_wt <- left_join(top10_wt, symbol_map_wt, by = "ENTREZID")

# Build the volcano plot with highlighted and labeled genes
ggplot(de_wt_df, aes(x = logFC, y = negLogP, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +                         # plot all points
  scale_color_manual(
    values = c("Up"   = "red",
               "Down" = "blue",
               "NS"   = "gray")
  ) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +  # fold-change cutoffs
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # p-value cutoff
  geom_text_repel(
    data          = top10_wt,                                   # only label top genes
    aes(label     = SYMBOL),
    show.legend   = FALSE,                                      # exclude text from legend
    size          = 5,                                          # larger text
    fontface      = "bold",                                     # bold font
    box.padding   = 0.8,                                        # padding around text
    point.padding = 0.5,                                        # min distance from points
    nudge_y       = 1,                                          # nudge labels upward
    direction     = "both",                                     # allow x & y adjustment
    force         = 1.5,                                        # stronger repulsion
    max.overlaps  = Inf,                                        # do not drop any
    segment.size  = 0.3,                                        # connector line width
    segment.color = "black"                                     # connector line color
  ) +
  scale_x_continuous(expand = expansion(add = c(0.3, 0.3))) +    # horizontal padding
  scale_y_continuous(expand = expansion(add = c(0.8, 1))) +      # vertical padding
  coord_cartesian(clip = "off") +                                # allow labels outside panel
  theme_minimal() +
  theme(
    plot.margin = unit(c(1, 3, 1, 1), "lines")                  # extra margin for labels
  ) +
  labs(
    title = "Volcano Plot: WT 24m vs WT 6m",
    x     = "log2 Fold Change",
    y     = "-log10(p-value)",
    color = "Regulation"
  )

## Gene and pathway enrichment - WT 24m vs WT 6m

# Get contrast result table sorted by logFC
wt_24m_vs_wt_6m <- topTable(contrast.results, number = 29129, coef = 3, sort.by = "logFC")
fold.change.wt_24m_vs_wt_6m <- wt_24m_vs_wt_6m$logFC
genes.ids.wt_24m_vs_wt_6m <- rownames(wt_24m_vs_wt_6m)

# Identify up- and downregulated genes
activated_ids <- genes.ids.wt_24m_vs_wt_6m[fold.change.wt_24m_vs_wt_6m > 1.5]
repressed_ids <- genes.ids.wt_24m_vs_wt_6m[fold.change.wt_24m_vs_wt_6m < -1.5]

# Map Probe IDs to ENTREZ IDs
mapped <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                keys = activated_ids,
                                columns = "ENTREZID",
                                keytype = "PROBEID")
mapped_unique <- mapped %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE)
entrez_ids_activated <- unique(mapped_unique$ENTREZID)

# GO enrichment for upregulated genes
ego_BP <- enrichGO(entrez_ids_activated, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",
                   pAdjustMethod = "BH", universe = universe_entrez,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
barplot(ego_BP, showCategory = 10, title = "GO:BP - WT 24m vs WT 6m Activated")

# Molecular Function
ego_MF <- enrichGO(entrez_ids_activated, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "MF",
                   pAdjustMethod = "BH", universe = universe_entrez,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
barplot(ego_MF, showCategory = 10, title = "GO:MF - WT 24m vs WT 6m Activated")

# Cellular Component
ego_CC <- enrichGO(entrez_ids_activated, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "CC",
                   pAdjustMethod = "BH", universe = universe_entrez,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
barplot(ego_CC, showCategory = 10, title = "GO:CC - WT 24m vs WT 6m Activated")

# KEGG pathway enrichment
ekegg_activated <- enrichKEGG(gene = entrez_ids_activated, organism = "mmu", universe = universe_entrez,
                              pvalueCutoff = 0.05, qvalueCutoff = 0.1)
barplot(ekegg_activated, showCategory = 10, title = "KEGG - WT 24m vs WT 6m Activated")

# Assign names and clean fold changes
mapped_all <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                    keys = genes.ids.wt_24m_vs_wt_6m,
                                    columns = "ENTREZID",
                                    keytype = "PROBEID")
names(fold.change.wt_24m_vs_wt_6m) <- mapped_all$ENTREZID[match(rownames(wt_24m_vs_wt_6m), mapped_all$PROBEID)]
fold.change.wt_24m_vs_wt_6m <- fold.change.wt_24m_vs_wt_6m[!is.na(names(fold.change.wt_24m_vs_wt_6m))]

# Zero out non-activated genes for Pathview
fold.change.cleaned <- fold.change.wt_24m_vs_wt_6m
fold.change.cleaned[!(names(fold.change.cleaned) %in% entrez_ids_activated)] <- 0

# Visualize specific KEGG pathway
pathview(
  gene.data  = fold.change.cleaned,
  pathway.id = "mmu04657", 
  species    = "mmu",
  gene.idtype = "ENTREZID",
  limit = list(gene = max(abs(fold.change.cleaned)), cpd = 1)
)

# REPPRESSED GENES

# Map repressed Probe IDs → ENTREZID
mapped_repressed <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                          keys = repressed_ids,
                                          columns = "ENTREZID",
                                          keytype = "PROBEID")

# Clean NAs and duplicates
mapped_unique_repressed <- mapped_repressed %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(PROBEID, .keep_all = TRUE)

# Get unique ENTREZ IDs
entrez_ids_repressed <- unique(mapped_unique_repressed$ENTREZID)

# GO enrichment for repressed genes
ego_BP_r <- enrichGO(
  gene          = entrez_ids_repressed,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)
barplot(ego_BP_r, showCategory = 10, title = "GO:BP - WT 24m vs WT 6m Repressed")

ego_MF_r <- enrichGO(
  gene          = entrez_ids_repressed,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)
barplot(ego_MF_r, showCategory = 10, title = "GO:MF - WT 24m vs WT 6m Repressed")

ego_CC_r <- enrichGO(
  gene          = entrez_ids_repressed,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)
barplot(ego_CC_r, showCategory = 10, title = "GO:CC - WT 24m vs WT 6m Repressed")

# KEGG enrichment
ekegg_repressed <- enrichKEGG(
  gene          = entrez_ids_repressed,
  organism      = "mmu",
  universe      = universe_entrez,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)
barplot(ekegg_repressed, showCategory = 10, title = "KEGG - WT 24m vs WT 6m Repressed")

# Map all genes used in contrast (PROBEID → ENTREZID) to align fold changes
mapped_all <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                    keys = genes.ids.wt_24m_vs_wt_6m,
                                    columns = "ENTREZID",
                                    keytype = "PROBEID")



# Convert the enrichKEGG result for WT 24m vs WT 6m into a data frame
df.ekegg_repressed_24_6_WT <- as.data.frame(ekegg_repressed)

# Display the first few rows to inspect enriched pathways
head(df.ekegg_repressed_24_6_WT)

# Prepare the fold-change data for Pathview:
# Create a copy of the original WT 24m vs WT 6m fold-change vector
fold.change.wt_24m_vs_wt_6m_cleaned <- fold.change.wt_24m_vs_wt_6m

# Zero out all genes NOT present in the list of differentially expressed ENTREZ IDs
#    names(...) returns the ENTREZIDs for each element in the vector
#    Any name not found in entrez_ids_wt_24m_vs_wt_6m is set to 0
fold.change.wt_24m_vs_wt_6m_cleaned[
  !(names(fold.change.wt_24m_vs_wt_6m_cleaned) %in% entrez_ids_repressed)
] <- 0

# Visualize the KEGG pathway with ID "mmu04020" (e.g., Calcium signaling pathway)
pathview(
  gene.data   = fold.change.wt_24m_vs_wt_6m_cleaned,           # fold-change data named by ENTREZID
  pathway.id  = "mmu04020",                                    # KEGG pathway ID to visualize
  species     = "mmu",                                         # organism: Mus musculus
  gene.idtype = "ENTREZID",                                    # identifier type in gene.data
  limit       = list(
    gene = max(abs(fold.change.wt_24m_vs_wt_6m_cleaned)),      # color scaling range based on absolute max
    cpd  = 1                                                   # compound limit (not used here)
  )
)

## Translation of genes

# Query the annotation database to get the PROBEID and ENTREZID for a specific gene symbol (e.g., "Lgals3")
symbol_to_probes <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                          keys = "Lgals3",  # Replace with any gene symbol of interest
                                          columns = c("PROBEID", "ENTREZID"),
                                          keytype = "SYMBOL")
symbol_to_probes


## Gene Lists per Contrast

# CONTRAST 1: KO 6m vs WT 6m

# Get gene symbols for activated probes
names_ko_6m_wt_6m_activated <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                     keys = activated_ids_ko_6m_wt_6m,
                                                     columns = "SYMBOL",
                                                     keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)

# Get gene symbols for repressed probes
names_ko_6m_wt_6m_repressed <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                     keys = repressed_ids_ko_6m_wt_6m,
                                                     columns = "SYMBOL",
                                                     keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)


# CONTRAST 2: KO 24m vs KO 6m

names_ko_24m_ko_6m_activated <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                      keys = activated_ids_ko_24m_ko_6m,
                                                      columns = "SYMBOL",
                                                      keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)

names_ko_24m_ko_6m_repressed <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                      keys = repressed_ids_ko_24m_ko_6m,
                                                      columns = "SYMBOL",
                                                      keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)


# CONTRAST 3: WT 24m vs WT 6m

# Activated
names_wt_24m_wt_6m_activated <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                      keys = activated_ids_wt_24m_vs_wt_6m,
                                                      columns = "SYMBOL",
                                                      keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)

# Repressed
names_wt_24m_wt_6m_repressed <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                      keys = repressed_ids_wt_24m_vs_wt_6m,
                                                      columns = "SYMBOL",
                                                      keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)


# CONTRAST 4: KO 24m vs WT 24m

names_ko_24m_wt_24m_activated <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                       keys = activated_ids_ko_24m_wt_24m,
                                                       columns = "SYMBOL",
                                                       keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)

names_ko_24m_wt_24m_repressed <- AnnotationDbi::select(clariomsmousetranscriptcluster.db,
                                                       keys = repressed_ids_ko_24m_wt_24m,
                                                       columns = "SYMBOL",
                                                       keytype = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(PROBEID, .keep_all = TRUE)


## Merge and Annotate

# Combine all contrasts into one table and add labels for contrast and regulation status
complet_table <- bind_rows(
  names_ko_6m_wt_6m_activated %>% mutate(Contrast = "KO6m_vs_WT6m", Status = "Activated"),
  names_ko_6m_wt_6m_repressed %>% mutate(Contrast = "KO6m_vs_WT6m", Status = "Repressed"),
  
  names_ko_24m_ko_6m_activated %>% mutate(Contrast = "KO24m_vs_KO6m", Status = "Activated"),
  names_ko_24m_ko_6m_repressed %>% mutate(Contrast = "KO24m_vs_KO6m", Status = "Repressed"),
  
  names_wt_24m_wt_6m_activated %>% mutate(Contrast = "WT24m_vs_WT6m", Status = "Activated"),
  names_wt_24m_wt_6m_repressed %>% mutate(Contrast = "WT24m_vs_WT6m", Status = "Repressed"),
  
  names_ko_24m_wt_24m_activated %>% mutate(Contrast = "KO24m_vs_WT24m", Status = "Activated"),
  names_ko_24m_wt_24m_repressed %>% mutate(Contrast = "KO24m_vs_WT24m", Status = "Repressed")
) %>%
  select(Contrast, Status, PROBEID, SYMBOL)

# Sort the table by contrast and regulation state
complet_table <- complet_table %>%
  arrange(Contrast, Status, SYMBOL)


## Render HTML Table

# Format the table using kable and style it for HTML export
html_table <- kable(complet_table, format = "html", escape = FALSE,
                    caption = "Activated and repressed genes with their symbols by contrast") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "left", font_size = 12)

# Save the table as an HTML file
save_kable(html_table, file = "contrasts_genes.html")

# Display the styled table in the viewer
html_table

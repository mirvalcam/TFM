###########################################################################
##  FULL WORKFLOW – PRESENCE/ABSENCE PROTEOMICS ANALYSIS                ##
##  Comparisons:                                                        ##
##  * Aging:        WT_24M vs WT_6M                                     ##
##  * Genotype:     KO_24M vs WT_24M                                    ##
##  Includes: GO & KEGG enrichment · Venn Diagrams · Presence table     ##
###########################################################################

## --------------------------- 1. LOAD LIBRARIES ------------------------ ##
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(VennDiagram)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(kableExtra)
})

## --------------------------- 2. IMPORT DATA --------------------------- ##
file <- "072822_Resultados_Proteomica liver Cuali.xlsx"
sheets <- excel_sheets(file)

# Read each replicate (sheet); keep only rows with Score ≥ 2
data <- lapply(sheets[1:16], function(sheet)
  read_excel(file, sheet = sheet) |> filter(Score >= 2))

names(data) <- sheets[1:16]

# Extract UNIPROT accession IDs per sample
ids_per_sample <- lapply(data, \(df) df$Accession)

## --------------------------- 3. GROUP BY CONDITION --------------------- ##
groups <- list(
  WT_6M  = unique(unlist(ids_per_sample[1:4])),
  KO_6M  = unique(unlist(ids_per_sample[5:8])),
  WT_24M = unique(unlist(ids_per_sample[9:12])),
  KO_24M = unique(unlist(ids_per_sample[13:16]))
)

## --------------------------- 4. SET OPERATIONS ------------------------- ##
# Compare presence/absence between groups

# Aging-related: gained or lost proteins in WT_24M vs WT_6M
gained_with_age  <- setdiff(groups$WT_24M, groups$WT_6M)
lost_with_age    <- setdiff(groups$WT_6M, groups$WT_24M)

# Genotype-related (at 24m): KO-specific or WT-specific proteins
recovered_in_KO  <- setdiff(groups$KO_24M, groups$WT_24M)
exclusive_to_WT  <- setdiff(groups$WT_24M, groups$KO_24M)

## --------------------------- 5. MAP UNIPROT → ENTREZ ------------------- ##
# Safe wrapper to avoid NAs
bitr_safe <- function(vec) bitr(vec, "UNIPROT", "ENTREZID", OrgDb = org.Mm.eg.db) |>
  filter(!is.na(ENTREZID))

# Perform mappings
ids_gained     <- bitr_safe(gained_with_age)
ids_lost       <- bitr_safe(lost_with_age)
ids_recovered  <- bitr_safe(recovered_in_KO)
ids_exclusive  <- bitr_safe(exclusive_to_WT)

## --------------------------- 6. ENRICHMENT ANALYSIS -------------------- ##
# GO enrichment for BP, MF, CC
run_go <- function(ids, ont) enrichGO(ids$ENTREZID, OrgDb = org.Mm.eg.db,
                                      keyType = "ENTREZID", ont = ont,
                                      pvalueCutoff = 0.05)

# Run for all sets and categories
ego_list <- list(
  gan_BP  = run_go(ids_gained,     "BP"),
  gan_MF  = run_go(ids_gained,     "MF"),
  gan_CC  = run_go(ids_gained,     "CC"),
  per_BP  = run_go(ids_lost,       "BP"),
  per_MF  = run_go(ids_lost,       "MF"),
  per_CC  = run_go(ids_lost,       "CC"),
  rec_BP  = run_go(ids_recovered,  "BP"),
  rec_MF  = run_go(ids_recovered,  "MF"),
  rec_CC  = run_go(ids_recovered,  "CC"),
  exc_BP  = run_go(ids_exclusive,  "BP"),
  exc_MF  = run_go(ids_exclusive,  "MF"),
  exc_CC  = run_go(ids_exclusive,  "CC")
)

# Define universe as all detected proteins (across all replicates)
univ_ids <- bitr_safe(unique(unlist(groups)))$ENTREZID |> unique()

# KEGG enrichment function
run_kegg <- function(ids, title) {
  ek <- enrichKEGG(ids$ENTREZID, organism = "mmu", universe = univ_ids,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  barplot(ek, showCategory = 15) + ggtitle(title)
  ek
}

# KEGG analyses for each protein group
kegg_gan <- run_kegg(ids_gained,     "Proteins gained with age")
kegg_per <- run_kegg(ids_lost,       "Proteins lost with age")
kegg_rec <- run_kegg(ids_recovered,  "Proteins recovered in KO24")
kegg_exc <- run_kegg(ids_exclusive,  "Proteins exclusive to WT24")

## --------------------------- 7. DOT/BAR PLOTS -------------------------- ##
# Create barplots for all GO results
plot_go <- function(ego, ttl) barplot(ego, showCategory = 15) + ggtitle(ttl)

plot_go(ego_list$gan_BP, "Exclusively in WT24m vs WT6m – BP")
plot_go(ego_list$gan_MF, "Exclusively in WT24m vs WT6m – MF")
plot_go(ego_list$gan_CC, "Exclusively in WT24m vs WT6m – CC")
plot_go(ego_list$per_BP, "Exclusively in WT6m vs WT24m – BP")
plot_go(ego_list$per_MF, "Exclusively in WT6m vs WT24m – MF")
plot_go(ego_list$per_CC, "Exclusively in WT6m vs WT24m – CC")
plot_go(ego_list$rec_BP, "Exclusively in KO24m vs WT24m – BP")
plot_go(ego_list$rec_MF, "Exclusively in KO24m vs WT24m – MF")
plot_go(ego_list$rec_CC, "Exclusively in KO24m vs WT24m – CC")
plot_go(ego_list$exc_BP, "Exclusively in WT24m vs KO24m – BP")
plot_go(ego_list$exc_MF, "Exclusively in WT24m vs KO24m – Exclusive – MF")
plot_go(ego_list$exc_CC, "Exclusively in WT24m vs KO24m – Exclusive – CC")


## --------------------------- 8. VENN DIAGRAMS -------------------------- ##
# Venn: aging comparison (WT)
venn.diagram(
  x              = list(WT_6M = groups$WT_6M, WT_24M = groups$WT_24M),
  category.names = c("WT 6 m", "WT 24 m"),
  filename       = "venn_aging_WT.png",
  imagetype      = "png", height = 2000, width = 2000, resolution = 500,
  fill           = c("skyblue", "orange"), alpha = 0.5, cex = 1.5, cat.cex = 1.3
)

# Venn: genotype comparison (24m)
venn.diagram(
  x              = list(WT_24M = groups$WT_24M, KO_24M = groups$KO_24M),
  category.names = c("WT 24 m", "KO 24 m"),
  filename       = "venn_genotype_24m.png",
  imagetype      = "png", height = 2000, width = 2000, resolution = 500,
  fill           = c("lightgreen", "tomato"), alpha = 0.5, cex = 1.5, cat.cex = 1.3
)

## --------------------------- 9. PRESENCE/ABSENCE TABLE ----------------- ##

# 1) Logical matrix of protein presence per group
presence_df <- tibble(
  UNIPROT = unique(unlist(groups))
) %>% 
  mutate(
    WT_6M   = UNIPROT %in% groups$WT_6M,
    KO_6M   = UNIPROT %in% groups$KO_6M,
    WT_24M  = UNIPROT %in% groups$WT_24M,
    KO_24M  = UNIPROT %in% groups$KO_24M
  )

# 2) Classification columns for interpretation
presence_df <- presence_df %>% 
  mutate(
    Age_Status = case_when(
      WT_6M  & !WT_24M ~ "Lost_with_age",
      !WT_6M & WT_24M  ~ "Gained_with_age",
      TRUE             ~ "Unchanged"
    ),
    KO_Rescue = case_when(
      Age_Status == "Lost_with_age" & KO_24M  ~ "Recovered_by_KO",
      Age_Status == "Lost_with_age" & !KO_24M ~ "Not_recovered",
      TRUE                                   ~ NA_character_
    ),
    Genotype24 = case_when(
      KO_24M & !WT_24M ~ "Exclusive_KO24",
      WT_24M & !KO_24M ~ "Exclusive_WT24",
      TRUE             ~ "Shared_24M"
    )
  ) %>% 
  arrange(UNIPROT)

# 3) Save full table as CSV
write.csv(presence_df, "proteins_presence_absence.csv", row.names = FALSE)

# 4) Preview first 100 rows in HTML format
html_tbl <- presence_df %>% 
  slice_head(n = 100) %>% 
  kable(format = "html", escape = FALSE,
        caption = "Presence/absence table of proteins") %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 12)

save_kable(html_tbl, "proteins_presence_absence_preview.html")

# 5) Confirmation message
cat("Table generated: proteins_presence_absence.csv\n")


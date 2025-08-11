# Chargement des packages nécessaires
if (!require("Biostrings")) install.packages("BiocManager"); BiocManager::install("Biostrings")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("patchwork")) install.packages("patchwork")
if (!require("stringr")) install.packages("stringr")

library(Biostrings)
library(tidyverse)
library(patchwork)
library(stringr)

# Répertoire de sortie
output_dir <- "/beegfs/project/horizon/PDV_Segments/Data/Metrics/"

# Fonction : parser les gènes à partir des headers
parse_gene_fasta <- function(file, species) {
  fasta <- readAAStringSet(file)
  headers <- names(fasta)

  df <- tibble(
    header = headers,
    segment = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    order = NA_integer_,
    species = species
  )

  for (i in seq_along(headers)) {
    h <- headers[i]
    parts <- strsplit(h, "_")[[1]]

    # Définir start et end selon l'espèce
    if (species == "Microplitis_demolitor") {
      start <- as.integer(parts[3])
      end <- as.integer(parts[4])
      segment <- parts[6]
    } else {
      start <- as.integer(parts[2])
      end <- as.integer(parts[3])
      segment <- parts[5]
    }

    # Extraire le rang du gène dans le segment (X/Y)
    match <- str_match(h, "_(\\d+)/\\d+_")
    order <- if (!is.na(match[1, 2])) as.integer(match[1, 2]) else NA_integer_

    df$start[i] <- min(start, end)
    df$end[i] <- max(start, end)
    df$segment[i] <- segment
    df$order[i] <- order
  }

  return(df)
}

# Fonction : parser les tailles de segments
get_segment_sizes <- function(file, species) {
  fasta <- readDNAStringSet(file)
  seg_names <- names(fasta)
  seg_names <- sapply(strsplit(seg_names, " "), function(x) x[1])

  tibble(
    segment = seg_names,
    segment_length = width(fasta),
    species = species
  )
}

# Fonction : stats sur distances entre gènes successifs
compute_distances_stats <- function(df) {
  df %>%
    group_by(species, segment) %>%
    arrange(order) %>%
    filter(n() > 1) %>%
    mutate(prev_end = lag(end)) %>%
    filter(!is.na(prev_end)) %>%
    mutate(distance = abs(start - prev_end)) %>%
    summarise(
      Min_gene_distance = min(distance),
      Max_gene_distance = max(distance),
      Median_gene_distance = median(distance),
      Mean_gene_distance = mean(distance),
      .groups = "drop"
    )
}

# Chemins
base_gene <- "/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Clean_Data/"
base_seg  <- "/beegfs/project/horizon/PDV_Segments/Data/Segments/"

species_list <- c("Cotesia_congregata", "Microplitis_demolitor", "Hyposoter_didymator")
colors <- c("Cotesia_congregata" = "skyblue", "Microplitis_demolitor" = "blue", "Hyposoter_didymator" = "purple")

# Lecture des gènes
genes_all <- bind_rows(
  parse_gene_fasta(file.path(base_gene, "Cotesia_congregata.fa"), "Cotesia_congregata"),
  parse_gene_fasta(file.path(base_gene, "Microplitis_demolitor.fa"), "Microplitis_demolitor"),
  parse_gene_fasta(file.path(base_gene, "Hyposoter_didymator.fa"), "Hyposoter_didymator")
)

# Lecture des tailles de segments
segment_sizes <- bind_rows(
  get_segment_sizes(file.path(base_seg, "Cotesia_congregata_segments.fa"), "Cotesia_congregata"),
  get_segment_sizes(file.path(base_seg, "Microplitis_demolitor_segments.fa"), "Microplitis_demolitor"),
  get_segment_sizes(file.path(base_seg, "Hyposoter_didymator_segments.fa"), "Hyposoter_didymator")
)

# Comptage des gènes par segment
gene_counts <- genes_all %>%
  group_by(species, segment) %>%
  summarise(Gene_count = n(), .groups = "drop")

# Statistiques sur les distances
dist_stats <- compute_distances_stats(genes_all)

# Fusion finale
all_metrics <- segment_sizes %>%
  rename(Segment = segment, Segment_length = segment_length, Species = species) %>%
  left_join(gene_counts, by = c("Species" = "species", "Segment" = "segment")) %>%
  left_join(dist_stats, by = c("Species" = "species", "Segment" = "segment")) %>%
  arrange(Species, Segment)

# Sauvegarde
write_tsv(all_metrics, file.path(output_dir, "PDV_all_metrics.tsv"))



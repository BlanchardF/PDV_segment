# Chargement des packages nécessaires
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("patchwork")) install.packages("patchwork")

library(tidyverse)
library(patchwork)

# Répertoire contenant le fichier TSV et où les plots seront sauvegardés
output_dir <- "/beegfs/project/horizon/PDV_Segments/Data/Metrics/"

# Définir les couleurs par espèce
colors <- c(
  Cotesia_congregata = "skyblue",
  Microplitis_demolitor = "blue",
  Hyposoter_didymator = "purple"
)

# Lecture du fichier de métriques
metrics <- read_tsv(file.path(output_dir, "PDV_all_metrics.tsv"), col_types = cols())

# Mise en forme longue pour les plots
plot_data <- metrics %>%
  pivot_longer(cols = c(Segment_length, Gene_count, Mean_gene_distance),
               names_to = "Plot_type", values_to = "Value") %>%
  filter(!is.na(Value))

# Définition des quantiles d'intérêt pour chaque métrique
quantiles_def <- list(
  Segment_length = c(0.05, 0.10, 0.20, 0.80, 0.90, 0.95),
  Mean_gene_distance = c(0.80, 0.90, 0.95),
  Gene_count = c(0.05, 0.10, 0.20)
)

# Calcul des quantiles à plat pour les lignes verticales sur le plot
quantiles_long <- bind_rows(lapply(names(quantiles_def), function(metric) {
  qs <- quantiles_def[[metric]]
  metric_data <- plot_data %>% filter(Plot_type == metric)
  map_dfr(qs, function(q) {
    metric_data %>%
      group_by(Species) %>%
      summarise(Value = quantile(Value, probs = q, na.rm = TRUE), .groups = "drop") %>%
      mutate(Plot_type = metric, Quantile = paste0("Q", q * 100))
  })
}))

# Génération d’un tableau large pour sortie
quantiles_wide <- quantiles_long %>%
  mutate(Quant = str_remove(Quantile, "^Q")) %>%
  pivot_wider(
    id_cols = c(Plot_type, Species),
    names_from = Quant,
    values_from = Value,
    names_sort = TRUE
  ) %>%
  select(Plot_type, Species, `5`, `10`, `20`, `80`, `90`, `95`)  # ordre demandé

# Sauvegarde du tableau des quantiles
write_tsv(quantiles_wide, file.path(output_dir, "PDV_plot_quantiles.tsv"))

# Fonction de génération d’un plot de densité avec lignes de quantiles
make_density_plot <- function(metric_name, title, xlab) {
  df <- filter(plot_data, Plot_type == metric_name)
  lines_df <- filter(quantiles_long, Plot_type == metric_name)

  ggplot(df, aes(x = Value, color = Species)) +
    geom_density(size = 1) +
    scale_color_manual(values = colors) +
    geom_vline(data = lines_df,
               aes(xintercept = Value, color = Species),
               linetype = "dashed", alpha = 0.5) +
    labs(title = title, x = xlab, y = "Density") +
    scale_x_continuous(n.breaks = 10) +
    theme_minimal()
}

# Création des trois graphiques
p1 <- make_density_plot("Segment_length", "Segment length distribution", "Segment length (bp)")
p2 <- make_density_plot("Mean_gene_distance", "Mean gene distance distribution", "Mean gene distance (bp)")
p3 <- make_density_plot("Gene_count", "Gene count per segment", "Number of genes per segment")

# Assemblage vertical
final_plot <- (p1 / p2 / p3) +
  plot_annotation(title = "Distributions of PDV segment properties")

# Export
ggsave(file.path(output_dir, "PDV_segments_distributions.svg"),
       plot = final_plot, width = 10, height = 12)


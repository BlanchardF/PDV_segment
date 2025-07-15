#!/bin/Rscript

cat("chargement package ")
#package 

package_name <- "karyoploteR"
if (!require(package_name, character.only = TRUE)) {
  install.packages(package_name)
  library(package_name, character.only = TRUE)
} else {
  library(package_name, character.only = TRUE)
}



package_name <- "dplyr"
if (!require(package_name, character.only = TRUE)) {
  install.packages(package_name)
  library(package_name, character.only = TRUE)
} else {
  library(package_name, character.only = TRUE)
}



package_name <- "GenomicRanges"
if (!require(package_name, character.only = TRUE)) {
  install.packages(package_name)
  library(package_name, character.only = TRUE)
} else {
  library(package_name, character.only = TRUE)
}




package_name <- "igraph"
if (!require(package_name, character.only = TRUE)) {
  install.packages(package_name)
  library(package_name, character.only = TRUE)
} else {
  library(package_name, character.only = TRUE)
}




# Fonction pour lire un fichier ou créer un data.frame vide
safe_read_table <- function(file, colnames) {
  if (file.info(file)$size == 0) {
    # Crée un data.frame vide avec les colonnes attendues
    df <- as.data.frame(matrix(ncol = length(colnames), nrow = 0))
    colnames(df) <- colnames
    return(df)
  } else {
    df <- read.table(file, header = FALSE, sep = "\t", quote="\"", comment.char = "")
    colnames(df) <- colnames
    return(df)
  }
}


package_name <- "argparse"
if (!require(package_name, character.only = TRUE)) {
  install.packages(package_name)
  library(package_name, character.only = TRUE)
} else {
  library(package_name, character.only = TRUE)
}


cat("  package charger  ")




################################
#### Chargement des données ####
################################

cat("   parser  ")

parser <- ArgumentParser(description= 'recup membres clusters of interest')

parser$add_argument('--input_Taille_contigs', '-iTc')
parser$add_argument('--input_Result_mmseqs2', '-iRm')
parser$add_argument('--input_HIM', '-iH', help= 'stat dataframe')
parser$add_argument('--input_DRJ', '-iD', help= '%CG')




parser$add_argument('--output_tab', '-ot')
parser$add_argument('--output_Plot', '-o')
parser$add_argument('--output_Plot_jpg', '-op')



xargs<- parser$parse_args()




########## chargement tableaux ##########
#Taille des contigs 
Taille_contig <- read.table(xargs$input_Taille_contigs, quote="\"", comment.char="")
#Resultats mmesqs search/convertalis 
result_mmseqs2 <- read.delim(xargs$Result_mmseqs2, header=FALSE)
# DRJ
DRJ <- safe_read_table(xargs$input_HIM, c("seqnames", "start", "end", "orientation"))
# HIM
HIM <- safe_read_table(xargs$input_DRJ, c("seqnames", "start", "end"))





#Variable 

#pourcentage de similariter entre les vrai gènes et ceux trouver par mmseqs 
tcov=0.60

#taille minimum pour qu'un contigs soit garder dans l'analyse 
Taille_contig_min=100000


#nombre de genes minimum pour former 
nb_gene_min=2

#distance entre gene 
dist_gene=10000

#taille min d'un cluster sans DRJ et HIM
min_cluster_gene=0

#taille min d'un cluster 
min_cluster=50000

#taille des estimation de segment a partir des DRJ/HIM 
taille_segment=20000


#renommer les header
colnames(Taille_contig)<- c("contig","taille")
colnames(result_mmseqs2)<- c("query","qlen","tlen","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qaln","tcov","qcov")
colnames(HIM)<-c("E-value","score","bias","contig","start","end")
colnames(DRJ)<-c("E-value","score","bias","contig","start","end")




# Filtrage des contigs de grande taille
contigs_grand <- Taille_contig[Taille_contig$taille > Taille_contig_min, ]

kp <- plotKaryotype(genome = contigs_grand, chromosomes = "all")




################
#### mmseqs ####
################
# Filtrage du data.frame des résultats MMseqs2
result_mmseqs2_filt <- result_mmseqs2[result_mmseqs2$query %in% contigs_grand$contig, ]

#filtrage qcov > 60%
result_mmseqs2_filt <- result_mmseqs2_filt[result_mmseqs2_filt$tcov > qcov, ]

# S'assurer que start <= end
start_fixed <- pmin(result_mmseqs2_filt$qstart, result_mmseqs2_filt$qend)
end_fixed <- pmax(result_mmseqs2_filt$qstart, result_mmseqs2_filt$qend)

# Création de l'objet GRanges
mmseqs2_gr <- GRanges(
  seqnames = result_mmseqs2_filt$query,
  ranges = IRanges(start = start_fixed, end = end_fixed),
  score = result_mmseqs2_filt$pident  # Pour coloration selon l'identité
)


genome_gr <- GRanges(
  seqnames = contigs_grand$contig,
  ranges = IRanges(start = 1, end = contigs_grand$taille)
)





######################################################################
#### Version avec les deux couleur différente en fonction du sens ####
######################################################################




################
#### mmseqs ####
################

# S'assurer que start <= end
start_fixed <- pmin(result_mmseqs2_filt$qstart, result_mmseqs2_filt$qend)
end_fixed <- pmax(result_mmseqs2_filt$qstart, result_mmseqs2_filt$qend)


# Création de l'objet GRanges
mmseqs2_gr <- GRanges(
  seqnames = result_mmseqs2_filt$query,
  ranges = IRanges(start = start_fixed, end = end_fixed),
  score = result_mmseqs2_filt$pident  # Pour coloration selon l'identité
)


genome_gr <- GRanges(
  seqnames = contigs_grand$contig,
  ranges = IRanges(start = 1, end = contigs_grand$taille)
)








####################
#### DRJ et HIM ####
####################

### DRJ ###
# Séparer les cas normaux et inversés
DRJ_normal <- DRJ[DRJ$start <= DRJ$end, ]
DRJ_reversed <- DRJ[DRJ$start > DRJ$end, ]

# Correction des valeurs start et end
DRJ_normal$start_fixed <- DRJ_normal$start
DRJ_normal$end_fixed <- DRJ_normal$end

DRJ_reversed$start_fixed <- DRJ_reversed$end
DRJ_reversed$end_fixed <- DRJ_reversed$start



safe_GRanges <- function(df, seq_col = "seqnames", start_col = "start", end_col = "end", extra_cols = list()) {
  if (nrow(df) == 0 || any(is.na(df[[start_col]])) || any(is.na(df[[end_col]]))) {
    return(GRanges())
  }
  gr <- GRanges(
    seqnames = df[[seq_col]],
    ranges = IRanges(start = df[[start_col]], end = df[[end_col]])
  )
  # Ajouter colonnes supplémentaires (ex : score)
  for (colname in names(extra_cols)) {
    mcols(gr)[[colname]] <- df[[extra_cols[[colname]]]]
  }
  return(gr)
}

# Exemple d'application aux DRJ normaux et inversés
DRJ_normal_gr <- safe_GRanges(DRJ_normal, seq_col = "contig", start_col = "start_fixed", end_col = "end_fixed", extra_cols = list(score = "score"))
DRJ_reversed_gr <- safe_GRanges(DRJ_reversed, seq_col = "contig", start_col = "start_fixed", end_col = "end_fixed", extra_cols = list(score = "score"))



# Inverser start et end si start > end pour HIM
HIM$start <- pmin(HIM$start, HIM$end)
HIM$end <- pmax(HIM$start, HIM$end)

# Créer l'objet GRanges pour HIM
HIM_gr <- safe_GRanges(HIM, seq_col = "contig", start_col = "start", end_col = "end")



###############################
#### Test nouvelle version ####
###############################
# Convertir les objets GRanges en data.frame
df_regions <- as.data.frame(mmseqs2_gr)

# Convertir les objets GRanges en data.frame
df_HIM <- as.data.frame(HIM_gr)




df_DRJ_normal <- as.data.frame(DRJ_normal_gr)
df_DRJ_reversed <- as.data.frame(DRJ_reversed_gr)

# Créer un objet karyotype avec les chromosomes détectés
kp <- plotKaryotype(genome = contigs_grand, chromosomes = "all")


# Ajouter les régions rouges (regions mmseqs2)
if (nrow(df_regions) > 0) {
  for (i in 1:nrow(df_regions)) {
    kpPlotRegions(kp, 
                  data = toGRanges(df_regions[i, c("seqnames", "start", "end")]), 
                  col = "red",  
                  border = "red", 
                  r0 = 0.06, r1 = 0.25)
  }
}

# Ajouter les HIM en vert
if (nrow(df_HIM) > 0) {
  df_HIM$end[is.na(df_HIM$end)] <- df_HIM$start[is.na(df_HIM$end)]
  for (i in 1:nrow(df_HIM)) {
    kpPlotRegions(kp, 
                  data = toGRanges(df_HIM[i, c("seqnames", "start", "end")]), 
                  col = "green",  
                  border = "green", 
                  r0 = 0.27, r1 = 0.45)
  }
}

# DRJ normaux (bleu clair)
if (nrow(df_DRJ_normal) > 0) {
  df_DRJ_normal <- df_DRJ_normal[!is.na(df_DRJ_normal$start) & !is.na(df_DRJ_normal$end), ]
  for (i in 1:nrow(df_DRJ_normal)) {
    kpPlotRegions(kp, 
                  data = toGRanges(df_DRJ_normal[i, c("seqnames", "start", "end")]), 
                  col = "#ADD8E6",  
                  border = "#ADD8E6", 
                  r0 = 0.52, r1 = 0.7)
  }
}

# DRJ inversés (bleu foncé)
if (nrow(df_DRJ_reversed) > 0) {
  df_DRJ_reversed <- df_DRJ_reversed[!is.na(df_DRJ_reversed$start) & !is.na(df_DRJ_reversed$end), ]
  for (i in 1:nrow(df_DRJ_reversed)) {
    kpPlotRegions(kp, 
                  data = toGRanges(df_DRJ_reversed[i, c("seqnames", "start", "end")]), 
                  col = "#00008B",  
                  border = "#00008B", 
                  r0 = 0.52, r1 = 0.7)
  }
}


# Ajouter les numéros de base
kpAddBaseNumbers(kp, tick.dist = 1000000, minor.tick.dist = 500000,
                 cex = 0.4, col = "black", tick.len = 5, units = "Mb",
                 digits = 2, add.units = TRUE)










##########################
#### Test autocluster ####
##########################


library(dplyr)

#### Clustering par proximité des gènes
cluster_genes <- function(df, threshold = dist_gene, min_genes = nb_gene_min ) {
  if (nrow(df) == 0) return(list())
  df$start <- as.integer(as.character(df$start))
  df$end <- as.integer(as.character(df$end))
  df <- df[order(df$seqnames, df$start), ]
  clusters <- list()
  
  for (contig in unique(df$seqnames)) {
    contig_df <- df %>% filter(seqnames == contig)
    if (nrow(contig_df) == 0) next
    current_cluster <- contig_df[1, , drop = FALSE]
    
    for (i in 2:nrow(contig_df)) {
      last_end <- current_cluster$end[nrow(current_cluster)]
      this_start <- contig_df$start[i]
      
      if (!is.na(this_start) && !is.na(last_end) && this_start - last_end <= threshold) {
        current_cluster <- rbind(current_cluster, contig_df[i, ])
      } else {
        if (nrow(current_cluster) >= min_genes) {
          clusters[[length(clusters) + 1]] <- current_cluster
        }
        current_cluster <- contig_df[i, , drop = FALSE]
      }
    }
    if (nrow(current_cluster) >= min_genes) {
      clusters[[length(clusters) + 1]] <- current_cluster
    }
  }
  
  return(clusters)
}


#### Clustering centré sur DRJ ou HIM
cluster_from_center <- function(df_center, window = taille_segment) {
  if (nrow(df_center) == 0) return(list())
  
  df_center <- df_center %>%
    mutate(
      seg_start = pmax(start - window, 0, na.rm = TRUE),
      seg_end = end + window
    ) %>%
    filter(!is.na(seqnames) & !is.na(seg_start) & !is.na(seg_end)) %>%
    arrange(seqnames, seg_start)
  
  merged_clusters <- list()
  
  for (contig in unique(df_center$seqnames)) {
    contig_df <- df_center[df_center$seqnames == contig, ]
    if (nrow(contig_df) == 0) next
    
    current_cluster <- contig_df[1, , drop = FALSE]
    
    for (i in 2:nrow(contig_df)) {
      this_start <- contig_df$seg_start[i]
      current_end <- current_cluster$seg_end[nrow(current_cluster)]
      
      if (!is.na(this_start) && !is.na(current_end) && this_start <= current_end) {
        current_cluster$seg_end[nrow(current_cluster)] <- max(current_end, contig_df$seg_end[i], na.rm = TRUE)
        current_cluster$end[nrow(current_cluster)] <- max(current_cluster$end[nrow(current_cluster)], contig_df$seg_end[i], na.rm = TRUE)
      } else {
        current_cluster$start <- current_cluster$seg_start
        current_cluster$end <- current_cluster$seg_end
        merged_clusters[[length(merged_clusters) + 1]] <- current_cluster
        current_cluster <- contig_df[i, , drop = FALSE]
      }
    }
    
    current_cluster$start <- current_cluster$seg_start
    current_cluster$end <- current_cluster$seg_end
    merged_clusters[[length(merged_clusters) + 1]] <- current_cluster
  }
  
  return(merged_clusters)
}


#### Résumé simple (sans segment de référence)
print_summary_no_ref <- function(clusters, methode) {
  if (length(clusters) == 0) {
    cat("\n", methode, ":\n", sep = "")
    cat("📦 Nombre total de clusters : 0\n")
    cat("🧬 Nombre total de bases dans les clusters : 0\n")
    return(invisible(NULL))
  }
  
  total_clusters <- length(clusters)
  
  total_bases <- sum(sapply(clusters, function(cl) {
    if (is.list(cl) && !is.data.frame(cl)) {
      # Liste avec start/end
      return(as.numeric(cl$end) - as.numeric(cl$start) + 1)
    } else if (is.data.frame(cl)) {
      # Data.frame avec colonnes start/end
      return(sum(as.numeric(cl$end) - as.numeric(cl$start) + 1))
    } else if (inherits(cl, "GRanges")) {
      # GRanges
      return(sum(width(cl)))
    } else {
      warning("Format de cluster non reconnu")
      return(0)
    }
  }))
  
  cat("\n", methode, ":\n", sep = "")
  cat("📦 Nombre total de clusters :", total_clusters, "\n")
  cat("🧬 Nombre total de bases dans les clusters :", total_bases, "\n")
}



#### Fusion de clusters multi-méthodes
merge_clusters_across_methods <- function(cluster_lists, min_length = 1000) {
  all_clusters <- do.call(rbind, lapply(cluster_lists, function(cl_list) {
    do.call(rbind, lapply(cl_list, function(cl) {
      cl <- cl[!is.na(cl$start) & !is.na(cl$end), ]
      if (nrow(cl) == 0) return(NULL)
      data.frame(
        seqnames = unique(cl$seqnames),
        start = suppressWarnings(min(cl$start, cl$seg_start, na.rm = TRUE)),
        end = suppressWarnings(max(cl$end, cl$seg_end, na.rm = TRUE))
      )
    }))
  }))
  
  if (is.null(all_clusters) || nrow(all_clusters) == 0) return(list())
  
  all_clusters <- all_clusters %>%
    filter(!is.na(start) & !is.na(end)) %>%
    arrange(seqnames, start)
  
  merged_clusters <- list()
  for (contig in unique(all_clusters$seqnames)) {
    contig_df <- all_clusters[all_clusters$seqnames == contig, ]
    if (nrow(contig_df) == 0) next
    
    current <- contig_df[1, , drop = FALSE]
    for (i in 2:nrow(contig_df)) {
      start_i <- contig_df$start[i]
      end_current <- current$end[nrow(current)]
      
      if (!is.na(start_i) && !is.na(end_current) && start_i <= end_current) {
        current$end[nrow(current)] <- max(end_current, contig_df$end[i], na.rm = TRUE)
      } else {
        merged_clusters[[length(merged_clusters) + 1]] <- current
        current <- contig_df[i, , drop = FALSE]
      }
    }
    merged_clusters[[length(merged_clusters) + 1]] <- current
  }
  
  merged_clusters <- Filter(function(cl) {
    cl_len <- cl$end - cl$start + 1
    !is.na(cl_len) && cl_len >= min_length
  }, merged_clusters)
  
  return(merged_clusters)
}


#### === APPLICATION === ####

# DRJ fusionné si applicable
df_DRJ <- if (exists("df_DRJ_normal") && exists("df_DRJ_reversed")) {
  rbind(df_DRJ_normal, df_DRJ_reversed)
} else {
  data.frame()
}

# Clustering
clusters_gene <- if (exists("df_regions")) cluster_genes(df_regions) else list()
clusters_DRJ <- if (exists("df_DRJ") && nrow(df_DRJ) > 0) cluster_from_center(df_DRJ) else list()
clusters_HIM <- if (exists("df_HIM") && nrow(df_HIM) > 0) cluster_from_center(df_HIM) else list()

# Affichage
print_summary_no_ref(clusters_gene, "Méthode 1 - Gènes proches")
print_summary_no_ref(clusters_DRJ, "Méthode 2 - DRJ centrés")
print_summary_no_ref(clusters_HIM, "Méthode 3 - HIM centrés")

# Vérifier si DRJ et HIM sont vides
is_DRJ_empty <- length(clusters_DRJ) == 0
is_HIM_empty <- length(clusters_HIM) == 0

# Appliquer la logique de fusion ou non
if (is_DRJ_empty && is_HIM_empty) {
  # Filtrer les clusters gènes à ≥ 1000 pb
  clusters_merged <- Filter(function(cl) {
    width <- sum(cl$end - cl$start + 1)
    width >= min_cluster_gene
  }, clusters_gene)
  
  print_summary_no_ref(clusters_merged, "Méthode 4 - Fusion (clusters = gènes seuls ≥ 1000 bp)")
} else {
  clusters_merged <- merge_clusters_across_methods(
    list(clusters_gene, clusters_DRJ, clusters_HIM),
    min_length = min_cluster 
  )
  print_summary_no_ref(clusters_merged, "Méthode 4 - Fusion (clusters ≥ 50000 bp)")
}










#### Création d'un tableau récapitulatif des clusters fusionnés
cluster_coords <- do.call(rbind, lapply(clusters_merged, function(cl) {
  data.frame(
    seqnames = unique(as.character(cl$seqnames)),   # En supposant un seul contig par cluster
    start = min(cl$start, na.rm = TRUE),
    end = max(cl$end, na.rm = TRUE),
    width = max(cl$end, na.rm = TRUE) - min(cl$start, na.rm = TRUE) + 1
  )
}))

# Affichage du tableau
print(cluster_coords)



# Optionnel : sauvegarde dans un fichier .tsv
# write.table(cluster_coords, file = "clusters_fusionnes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cluster_gr <- toGRanges(cluster_coords)

kpPlotRegions(kp, 
              data = cluster_gr,
              col = "purple", 
              border = "purple", 
              r0 = 0.75, r1 = 0.93)  # Ligne supérieure du plot




# Fonction pour extraire éléments (positions + nombre) dans un cluster donné
extract_elements_in_cluster <- function(cluster, elements_df, label = "gene") {
  if (nrow(elements_df) == 0) return(data.frame())
  
  # Filtrer éléments qui se chevauchent avec le cluster (sur le même contig)
  overlap <- elements_df %>%
    filter(seqnames == unique(cluster$seqnames) &
             !(end < cluster$start | start > cluster$end))
  
  if (nrow(overlap) == 0) {
    return(data.frame(
      n = 0,
      positions = NA_character_
    ))
  }
  
  positions <- paste0("[", overlap$start, "-", overlap$end, "]", collapse = ", ")
  
  data.frame(
    n = nrow(overlap),
    positions = positions
  )
}

# Construire le tableau résumé
cluster_summary <- do.call(rbind, lapply(clusters_merged, function(cl) {
  # Extraire coordonnées cluster
  cluster_info <- data.frame(
    seqnames = unique(as.character(cl$seqnames)),
    start = min(cl$start, na.rm = TRUE),
    end = max(cl$end, na.rm = TRUE)
  )
  
  # Gènes dans cluster
  genes_info <- extract_elements_in_cluster(cluster_info, df_regions, "gene")
  colnames(genes_info) <- c("n_genes", "positions_genes")
  
  # DRJ dans cluster
  DRJ_info <- extract_elements_in_cluster(cluster_info, df_DRJ, "DRJ")
  colnames(DRJ_info) <- c("n_DRJ", "positions_DRJ")
  
  # HIM dans cluster
  HIM_info <- extract_elements_in_cluster(cluster_info, df_HIM, "HIM")
  colnames(HIM_info) <- c("n_HIM", "positions_HIM")
  
  cbind(cluster_info, genes_info, DRJ_info, HIM_info)
}))







# Fonction pour extraire toutes les positions de type [start-end] depuis une colonne
extract_positions <- function(pos_column) {
  # Retourne NA s'il n'y a rien
  if (is.na(pos_column) || pos_column == "<NA>") return(integer(0))
  
  # Trouver tous les couples [start-end]
  matches <- gregexpr("\\[\\d+-\\d+\\]", pos_column)[[1]]
  if (matches[1] == -1) return(integer(0))
  
  sapply(regmatches(pos_column, gregexpr("\\d+-\\d+", pos_column))[[1]], function(x) {
    as.integer(unlist(strsplit(x, "-")))
  })
}

# Appliquer le recalcul du start et end
cluster_summary$adjusted_start <- NA_integer_
cluster_summary$adjusted_end <- NA_integer_

for (i in seq_len(nrow(cluster_summary))) {
  # Extraire les positions des gènes, DRJ et HIM
  genes <- extract_positions(cluster_summary$positions_genes[i])
  drj <- extract_positions(cluster_summary$positions_DRJ[i])
  him <- extract_positions(cluster_summary$positions_HIM[i])
  
  # Fusionner toutes les positions
  all_positions <- c(genes, drj, him)
  
  if (length(all_positions) > 0) {
    starts <- all_positions[seq(1, length(all_positions), by=2)]
    ends <- all_positions[seq(2, length(all_positions), by=2)]
    
    cluster_summary$adjusted_start[i] <- min(starts)
    cluster_summary$adjusted_end[i] <- max(ends)
  } else {
    # Si aucune position, on garde les anciennes
    cluster_summary$adjusted_start[i] <- cluster_summary$start[i]
    cluster_summary$adjusted_end[i] <- cluster_summary$end[i]
  }
}



cluster_summary$start <- cluster_summary$adjusted_start
cluster_summary$end <- cluster_summary$adjusted_end
cluster_summary$adjusted_start <- NULL
cluster_summary$adjusted_end <- NULL



print(cluster_summary)


# Enregistrer le tableau résumé
write.table(cluster_summary, xargs$output_tab, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Enregistrer le graphique
ggsave(filename = paste0(xargs$output_Plot), plot = Lb, width = 12, height = 10.08, dpi = 2400, device = "svg")
ggsave(filename = paste0(xargs$output_Plot_jpg), plot = Lb, width = 12, height = 10.08, dpi = 600)


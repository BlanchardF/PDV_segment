library(karyoploteR)
library(dplyr)
library(GenomicRanges)
library(igraph)

########## Cotesia sur Cotesia ##########
#Taille des contigs 
Taille_contig <- read.table("~/data/test/Taille_contig.tsv", quote="\"", comment.char="")
#Resultats mmesqs search/convertalis 
result_mmseqs2 <- read.delim("~/data/test/_result_mmseqs2.m8", header=FALSE)
#Positions des segments 
segment_positions <- read.table("~/data/test/segment_positions.tsv", quote="\"", comment.char="")
#HIM
HIM <- read.table("~/data/test/Cotesia_congregata-HIM-Cotesia_congregata_clean.tsv", quote="\"", comment.char="")
#DRJ
DRJ<- read.table("~/data/test/Cotesia_congregata-DRJ-Cotesia_congregata_clean.tsv", quote="\"", comment.char="")



########## Microplitis sur Cotesia ##########
#Taille des contigs 
Taille_contig <- read.table("~/data/test/Taille_contig.tsv", quote="\"", comment.char="")
#Resultats mmesqs search/convertalis 
result_mmseqs2 <- read.delim("~/data/test/Microplitis_Cotesia/_result_mmseqs2_test.m8", header=FALSE)
#Positions des segments 
segment_positions <- read.table("~/data/test/segment_positions.tsv", quote="\"", comment.char="")
#HIM
HIM <- read.table("~/data/test/Microplitis_Cotesia/Microplitis_demolitor-HIM-Cotesia_congregata_clean.tsv", quote="\"", comment.char="")
#DRJ
DRJ<- read.table("~/data/test/Microplitis_Cotesia/Microplitis_demolitor-DRJ-Cotesia_congregata_clean.tsv", quote="\"", comment.char="") ###pas bon 




########## Cotesia sur Microplitis ##########
#Taille des contigs 
Taille_contig <- read.table("~/data/test/Cotesia_Microplitis/Taille_contig.tsv", quote="\"", comment.char="")
#Resultats mmesqs search/convertalis 
result_mmseqs2 <- read.delim("~/data/test/Cotesia_Microplitis/_result_mmseqs2.m8", header=FALSE)
#Positions des segments 
segment_positions <- read.table("~/data/test/Cotesia_Microplitis/segment_positions_Microplitis_demolitor.tsv", quote="\"", comment.char="")
#HIM
HIM <- read.table("~/data/test/Cotesia_Microplitis/Cotesia_congregata-HIM-Microplitis_demolitor.tsv", quote="\"", comment.char="")
#DRJ
DRJ<- read.table("~/data/test/Cotesia_Microplitis/Cotesia_congregata-DRJ-Microplitis_demolitor.tsv", quote="\"", comment.char="")





########## Hyposoter_didymator sur Microplitis ##########
#Taille des contigs 
Taille_contig <- read.table("~/data/test/Cotesia_Microplitis/Taille_contig.tsv", quote="\"", comment.char="")
#Resultats mmesqs search/convertalis 
result_mmseqs2 <- read.delim("~/data/test/Cotesia_Microplitis/_result_mmseqs2.m8", header=FALSE)
#Positions des segments 
segment_positions <- read.table("~/data/test/Cotesia_Microplitis/segment_positions_Microplitis_demolitor.tsv", quote="\"", comment.char="")
#HIM
#HIM <- read.table("~/data/test/Cotesia_Microplitis/Cotesia_congregata-HIM-Microplitis_demolitor.tsv", quote="\"", comment.char="")
#DRJ
#DRJ<- read.table("~/data/test/Cotesia_Microplitis/Cotesia_congregata-DRJ-Microplitis_demolitor.tsv", quote="\"", comment.char="")






########## Hyposoter_didymator sur Cotesia ##########
#Taille des contigs 
Taille_contig <- read.table("~/data/test/Taille_contig.tsv", quote="\"", comment.char="")
#Resultats mmesqs search/convertalis 
result_mmseqs2 <- read.delim("~/data/test/Hyposoter_Cotesia/_result_mmseqs2.m8", header=FALSE)
#Positions des segments 
segment_positions <- read.table("~/data/test/segment_positions.tsv", quote="\"", comment.char="")
#HIM
#HIM <- read.table("~/data/test/Microplitis_Cotesia/Microplitis_demolitor-HIM-Cotesia_congregata_clean.tsv", quote="\"", comment.char="")
#DRJ
#DRJ<- read.table("~/data/test/Microplitis_Cotesia/Microplitis_demolitor-DRJ-Cotesia_congregata_clean.tsv", quote="\"", comment.char="") ###pas bon 



#renommer les header
colnames(Taille_contig)<- c("contig","taille")
colnames(result_mmseqs2)<- c("query","qlen","tlen","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qaln","tcov","qcov")
colnames(segment_positions)<-c("query","target","tstart","tend")
colnames(HIM)<-c("E-value","score","bias","contig","start","end")
colnames(DRJ)<-c("E-value","score","bias","contig","start","end")

head(Taille_contig)
head(result_mmseqs2)
head(segment_positions)
head(HIM)
head(DRJ)


# Filtrage des contigs de grande taille
contigs_grand <- Taille_contig[Taille_contig$taille > 100000, ]

kp <- plotKaryotype(genome = contigs_grand, chromosomes = "all")




################
#### mmseqs ####
################
# Filtrage du data.frame des résultats MMseqs2
result_mmseqs2_filt <- result_mmseqs2[result_mmseqs2$query %in% contigs_grand$contig, ]

#filtrage qcov > 80%
result_mmseqs2_filt <- result_mmseqs2_filt[result_mmseqs2_filt$tcov > 0.60, ]

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









##########################
#### Segment Positions ###
##########################

# Filtrer les segments sur les grands contigs
segment_positions_filt <- segment_positions[segment_positions$target %in% contigs_grand$contig, ]

# Créer l'objet GRanges
segment_gr <- GRanges(
  seqnames = segment_positions_filt$target,
  ranges = IRanges(start = segment_positions_filt$qstart, end = segment_positions_filt$qend)
)

# Replot tout (DRJ, HIM, mmseqs + segments)
kp <- plotKaryotype(genome = genome_gr, chromosomes = "all")




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

# Créer les objets GRanges séparément
DRJ_normal_gr <- GRanges(
  seqnames = DRJ_normal$contig,
  ranges = IRanges(start = DRJ_normal$start_fixed, end = DRJ_normal$end_fixed),
  score = DRJ_normal$score
)

DRJ_reversed_gr <- GRanges(
  seqnames = DRJ_reversed$contig,
  ranges = IRanges(start = DRJ_reversed$start_fixed, end = DRJ_reversed$end_fixed),
  score = DRJ_reversed$score
)




# Inverser start et end si start > end pour HIM
HIM$start <- pmin(HIM$start, HIM$end)
HIM$end <- pmax(HIM$start, HIM$end)

# Créer l'objet GRanges pour HIM
HIM_gr <- GRanges(
  seqnames = HIM$contig,
  ranges = IRanges(start = HIM$start, end = HIM$end),
  score = HIM$score  # Vous pouvez colorer selon le score si vous le souhaitez
)


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


# Ajouter les régions rouges (mmseqs2_gr) - Plus proche du chromosome
for (i in 1:nrow(df_regions)) {
  kpPlotRegions(kp, 
                data=toGRanges(df_regions[i, c("seqnames", "start", "end")]), 
                col="red",  
                border="red", 
                r0=0.06, r1=0.25)  # Proche du chromosome
}

# Ajouter les régions vertes (HIM_gr) - Juste au-dessus des rouges
for (i in 1:nrow(df_HIM)) {
  if (is.na(df_HIM$end[i])) { 
    df_HIM$end[i] <- df_HIM$start[i]
  }
  kpPlotRegions(kp, 
                data=toGRanges(df_HIM[i, c("seqnames", "start", "end")]), 
                col="green",  
                border="green", 
                r0=0.27, r1=0.45)  # Décalé au-dessus des rouges
} 


# Ajouter les régions bleues clair (DRJ normal) et bleu foncé (DRJ inversé)
for (i in 1:nrow(df_DRJ_normal)) {
  kpPlotRegions(kp, 
                data=toGRanges(df_DRJ_normal[i, c("seqnames", "start", "end")]), 
                col="#ADD8E6",  # Bleu clair
                border="#ADD8E6", 
                r0=0.52, r1=0.7)  
}

for (i in 1:nrow(df_DRJ_reversed)) {
  kpPlotRegions(kp, 
                data=toGRanges(df_DRJ_reversed[i, c("seqnames", "start", "end")]), 
                col="#00008B",  # Bleu foncé
                border="#00008B", 
                r0=0.52, r1=0.7)  
}

# Ajouter les numéros de base
kpAddBaseNumbers(kp, tick.dist = 1000000, minor.tick.dist = 500000,
                 cex = 0.4, col = "black", tick.len = 5, units = "Mb",
                 digits = 2, add.units = TRUE)




#### zoom 
#kp <- plotKaryotype(genome = genome_gr, chromosomes = "contig_2034", zoom = toGRanges("contig_2034:2300000-3100000"))
#kp <- plotKaryotype(genome = genome_gr, chromosomes = "contig_2034", zoom = toGRanges("contig_2034:19200000-19800000"))
#kp <- plotKaryotype(genome = genome_gr, chromosomes = "contig_2030", zoom = toGRanges("contig_2030:3500000-4000000"))
#kp <- plotKaryotype(genome = genome_gr, chromosomes = "contig_2030", zoom = toGRanges("contig_2030:10200000-10800000"))
#kp <- plotKaryotype(genome = genome_gr, chromosomes = "NC_068546.1", zoom = toGRanges("NC_068546.1:4000000-4500000"))


##########################
#### Segment Positions ###
##########################

# Filtrer les segments sur les grands contigs
segment_positions_filt <- segment_positions[segment_positions$target %in% contigs_grand$contig, ]

# Créer l'objet GRanges
segment_gr <- GRanges(
  seqnames = segment_positions_filt$target,
  ranges = IRanges(start = segment_positions_filt$tstart, end = segment_positions_filt$tend)
)

kpPlotRegions(kp, 
              data = segment_gr, 
              col = "pink", 
              border = "pink", 
              cex = 1.5, 
              r0=0, r1=0.05)  


























##########################
#### Test autocluster ####
##########################
######################################
#### Méthodes de clustering génomique
######################################
df_DRJ <- rbind(df_DRJ_normal, df_DRJ_reversed)

library(dplyr)

# Méthode 1 - Clustering par proximité des gènes
cluster_genes <- function(df, threshold = 10000, min_genes = 2) {
  df$start <- as.integer(as.character(df$start))
  df$end <- as.integer(as.character(df$end))
  df <- df[order(df$seqnames, df$start), ]
  clusters <- list()
  
  for (contig in unique(df$seqnames)) {
    contig_df <- df %>% filter(seqnames == contig)
    if (nrow(contig_df) == 0) next
    
    contig_df$start <- as.integer(as.character(contig_df$start))
    contig_df$end <- as.integer(as.character(contig_df$end))
    current_cluster <- contig_df[1, , drop = FALSE]
    
    for (i in 2:nrow(contig_df)) {
      last_end <- current_cluster$end[nrow(current_cluster)]
      this_start <- contig_df$start[i]
      
      if (!is.na(this_start) && !is.na(last_end)) {
        if (this_start - last_end <= threshold) {
          current_cluster <- rbind(current_cluster, contig_df[i, ])
        } else {
          if (nrow(current_cluster) >= min_genes) {
            clusters[[length(clusters) + 1]] <- current_cluster
          }
          current_cluster <- contig_df[i, , drop = FALSE]
        }
      } else {
        warning("start ou end manquant à l'itération ", i)
      }
    }
    
    if (nrow(current_cluster) >= min_genes) {
      clusters[[length(clusters) + 1]] <- current_cluster
    }
  }
  
  return(clusters)
}




# Méthode 2 et 3 - DRJ ou HIM centrés
cluster_from_center <- function(df_center, window = 20000) {
  df_center <- df_center %>%
    mutate(
      seg_start = pmax(start - window, 0, na.rm = TRUE),
      seg_end = end + window,
      start = pmax(start - window, 0, na.rm = TRUE),  # Ajouté pour compatibilité
      end = end + window                              # Ajouté pour compatibilité
    ) %>%
    filter(!is.na(seqnames) & !is.na(seg_start) & !is.na(seg_end)) %>%
    arrange(seqnames, seg_start)
  
  merged_clusters <- list()
  
  for (contig in unique(df_center$seqnames)) {
    contig_df <- df_center[df_center$seqnames == contig, ]
    if (nrow(contig_df) == 0) next
    
    merged <- contig_df[1, , drop = FALSE]
    
    for (i in 2:nrow(contig_df)) {
      seg_start_i <- contig_df$seg_start[i]
      seg_end_prev <- merged$seg_end[nrow(merged)]
      
      if (!is.na(seg_start_i) && !is.na(seg_end_prev) && seg_start_i <= seg_end_prev) {
        merged$seg_end[nrow(merged)] <- max(seg_end_prev, contig_df$seg_end[i], na.rm = TRUE)
        merged$end[nrow(merged)] <- max(merged$end[nrow(merged)], contig_df$seg_end[i], na.rm = TRUE)  # Mise à jour aussi de end
      } else {
        merged_clusters[[length(merged_clusters) + 1]] <- merged
        merged <- contig_df[i, , drop = FALSE]
      }
    }
    merged_clusters[[length(merged_clusters) + 1]] <- merged
  }
  
  return(merged_clusters)
}


  
cluster_from_center <- function(df_center, window = 20000) {
  # Crée les colonnes de segment élargi
  df_center <- df_center %>%
    mutate(
      seg_start = pmax(start - window, 0, na.rm = TRUE),
      seg_end = end + window
    ) %>%
    filter(!is.na(seqnames) & !is.na(seg_start) & !is.na(seg_end)) %>%
    arrange(seqnames, seg_start)
  
  # Initialise la liste des clusters
  merged_clusters <- list()
  
  for (contig in unique(df_center$seqnames)) {
    contig_df <- df_center[df_center$seqnames == contig, ]
    if (nrow(contig_df) == 0) next
    
    # Initialise le cluster en cours avec la première ligne
    current_cluster <- contig_df[1, , drop = FALSE]
    
    for (i in 2:nrow(contig_df)) {
      this_start <- contig_df$seg_start[i]
      current_end <- current_cluster$seg_end[nrow(current_cluster)]
      
      # Teste si chevauchement
      if (!is.na(this_start) && !is.na(current_end) && this_start <= current_end) {
        # Fusionne : on met à jour seg_end et end
        current_cluster$seg_end[nrow(current_cluster)] <- max(current_end, contig_df$seg_end[i], na.rm = TRUE)
        current_cluster$end[nrow(current_cluster)] <- max(current_cluster$end[nrow(current_cluster)], contig_df$seg_end[i], na.rm = TRUE)
      } else {
        # On ferme le cluster actuel, on le formate
        current_cluster$start <- current_cluster$seg_start
        current_cluster$end <- current_cluster$seg_end
        merged_clusters[[length(merged_clusters) + 1]] <- current_cluster
        current_cluster <- contig_df[i, , drop = FALSE]
      }
    }
    
    # Ajouter le dernier cluster du contig
    current_cluster$start <- current_cluster$seg_start
    current_cluster$end <- current_cluster$seg_end
    merged_clusters[[length(merged_clusters) + 1]] <- current_cluster
  }
  
  return(merged_clusters)
}




evaluate_coverage <- function(clusters, segment_positions) {
  covered <- 0
  cluster_hits <- logical(length(clusters))
  
  total_segment_bases <- sum(segment_positions$tend - segment_positions$tstart + 1, na.rm = TRUE)
  covered_bases <- 0
  total_cluster_bases <- 0
  false_positive_bases <- 0
  
  for (i in 1:nrow(segment_positions)) {
    seg <- segment_positions[i, ]
    for (j in seq_along(clusters)) {
      cl <- clusters[[j]]
      if (is.null(cl$start) || is.null(cl$end)) next
      contig <- unique(cl$seqnames)
      cl_start <- suppressWarnings(min(cl$start, na.rm = TRUE))
      cl_end <- suppressWarnings(max(cl$end, na.rm = TRUE))
      
      if (is.finite(cl_start) && is.finite(cl_end) &&
          seg$target == contig &&
          seg$tend >= cl_start &&
          seg$tstart <= cl_end) {
        covered <- covered + 1
        cluster_hits[j] <- TRUE
        
        # Bases communes (intersection)
        overlap_start <- max(cl_start, seg$tstart)
        overlap_end <- min(cl_end, seg$tend)
        if (overlap_start <= overlap_end) {
          covered_bases <- covered_bases + (overlap_end - overlap_start + 1)
        }
        
        break
      }
    }
  }
  
  for (j in seq_along(clusters)) {
    cl <- clusters[[j]]
    if (is.null(cl$start) || is.null(cl$end)) next
    cl_start <- suppressWarnings(min(cl$start, na.rm = TRUE))
    cl_end <- suppressWarnings(max(cl$end, na.rm = TRUE))
    if (is.finite(cl_start) && is.finite(cl_end)) {
      cl_size <- cl_end - cl_start + 1
      total_cluster_bases <- total_cluster_bases + cl_size
      if (!cluster_hits[j]) {
        false_positive_bases <- false_positive_bases + cl_size
      }
    }
  }
  
  return(list(
    total_segments = nrow(segment_positions),
    covered_segments = covered,
    uncovered_segments = nrow(segment_positions) - covered,
    uncovered_clusters = sum(!cluster_hits),
    total_clusters = length(clusters),
    total_segment_bases = total_segment_bases,
    covered_segment_bases = covered_bases,
    total_cluster_bases = total_cluster_bases,
    false_positive_bases = false_positive_bases
  ))
}


print_results <- function(result, methode) {
  cat("\n", methode, ":\n", sep = "")
  cat("🔢 Nombre total de segments :", result$total_segments, "\n")
  cat("✅ Segments couverts :", result$covered_segments, "\n")
  cat("📭 Segments non couverts :", result$uncovered_segments, "\n")
  cat("📦 Nombre total de clusters :", result$total_clusters, "\n")
  cat("❌ Clusters sans segment (faux positifs) :", result$uncovered_clusters, "\n")
  cat("🧬 Nombre total de bases dans les segments :", result$total_segment_bases, "\n")
  cat("🟩 Bases de segments couvertes :", result$covered_segment_bases, "\n")
  cat("📏 Nombre total de bases dans les clusters :", result$total_cluster_bases, "\n")
  cat("🚫 Bases des clusters hors segments (faux positifs) :", result$false_positive_bases, "\n")
}




# Méthode 1
clusters_gene <- cluster_genes(df_regions)
res1 <- evaluate_coverage(clusters_gene, segment_positions)

#Méthode 2
clusters_DRJ <- cluster_from_center(df_DRJ)
res2 <- evaluate_coverage(clusters_DRJ, segment_positions)

# Méthode 3
clusters_HIM <- cluster_from_center(df_HIM)
res3 <- evaluate_coverage(clusters_HIM, segment_positions)





print_results(res1, "Méthode 1 - Gènes proches")
print_results(res2, "Méthode 2 - DRJ centrés")
print_results(res3, "Méthode 3 - HIM centrés")






















merge_clusters_across_methods <- function(cluster_lists) {
  # Combiner tous les clusters dans un seul data.frame (start/end calculés proprement)
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
  
  # Nettoyage
  all_clusters <- all_clusters %>%
    filter(!is.na(start) & !is.na(end)) %>%
    arrange(seqnames, start)
  
  # Fusion des clusters qui se chevauchent
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
  
  return(merged_clusters)
}




clusters_merged <- merge_clusters_across_methods(list(clusters_gene, clusters_DRJ, clusters_HIM))
res_merged <- evaluate_coverage(clusters_merged, segment_positions)
print_results(res_merged, "Méthode 4 - Fusion de toutes les méthodes")


















summarize_clusters <- function(clusters, segment_positions) {
  summary_table <- data.frame()
  
  for (i in seq_along(clusters)) {
    cl <- clusters[[i]]
    
    # Extraire les bornes du cluster
    seqname <- unique(cl$seqnames)
    cl_start <- suppressWarnings(min(cl$start, cl$seg_start, na.rm = TRUE))
    cl_end <- suppressWarnings(max(cl$end, cl$seg_end, na.rm = TRUE))
    
    # Vérifier si le cluster couvre au moins un segment
    overlaps <- segment_positions %>%
      filter(target == seqname & !(tend < cl_start | tstart > cl_end))
    
    covers <- nrow(overlaps) > 0
    
    # Ajouter la ligne au tableau
    summary_table <- rbind(summary_table, data.frame(
      cluster_id = i,
      seqnames = seqname,
      start = cl_start,
      end = cl_end,
      length = cl_end - cl_start + 1,
      covers_segment = covers
    ))
  }
  
  return(summary_table)
}


summary_clusters <- summarize_clusters(clusters_merged, segment_positions)
head(summary_clusters)










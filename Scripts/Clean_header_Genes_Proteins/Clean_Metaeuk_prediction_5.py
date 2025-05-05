#!/bin/bash

# Fichier d'entrée
input_file="/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_4.fa"
# Fichier de sortie
output_file="/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_5.fa"


# Traitement des en-têtes
awk 'BEGIN {OFS=FS=" "} /^>/ { # Si la ligne commence par ">"
  $1="";                     # Supprime la première colonne
  gsub(/ /, "_", $0);        # Remplace les espaces par des underscores
  gsub(/\|/, "_", $0);       # Remplace les | par des underscores
  sub(/^_/, "", $0);         # Supprime le premier underscore
  print ">" $0               # Affiche la ligne modifiée avec ">"
  next
} {print}' "$input_file" > "$output_file"  # Copie les autres lignes sans modification
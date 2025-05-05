#!/bin/bash

# Chemins d'entrée et de sortie codés en dur
input_file="/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_3.fa"
output_file="/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_4.fa"

# Vérifier si le fichier d'entrée existe
if [[ ! -f "$input_file" ]]; then
  echo "Erreur : le fichier d'entrée $input_file n'existe pas."
  exit 1
fi

# Traitement des en-têtes pour ajouter le numéro de segment à la fin et remplacer les underscores par des tirets
awk '
/^>/ {
  split($0, a, "|")
  gsub("_", "-", a[2])  # Remplacer les underscores par des tirets dans la deuxième partie
  print $0 "|" a[2]
  next
} 
{print}' "$input_file" > "$output_file"

# Confirmation
echo "Traitement terminé. Le fichier modifié est enregistré sous $output_file."

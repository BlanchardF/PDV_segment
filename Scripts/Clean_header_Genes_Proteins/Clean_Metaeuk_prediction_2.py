#!/bin/bash

# Définir les fichiers d'entrée et de sortie
input_file="../../Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_1.fa"
output_file="../../Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_2.fa"


# Traiter le fichier d'entrée et générer la sortie
while IFS= read -r line; do
    # Vérifier si la ligne commence par ">"
    if [[ $line == ">"* ]]; then
        # Extraire la partie du header et les 7e et 8e éléments
        header_part=$(echo $line | cut -d '|' -f 7,8)
        # Remplacer les | par _
        header_part=$(echo $header_part | tr '|' '_')
        # Créer la nouvelle ligne d'en-tête avec les modifications
        new_header="${line} ${header_part}"
        echo $new_header >> "$output_file"
    else
        # Si ce n'est pas un header, juste afficher la ligne
        echo $line >> "$output_file"
    fi
done < "$input_file"

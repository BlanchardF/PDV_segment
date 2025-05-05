#!/bin/bash


#SBATCH --partition=normal
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=65G
#SBATCH --output=/beegfs/project/horizon/PDV_Segments/Scripts/Genes_Prediction/std_output_Genes_Prediction.txt
#SBATCH --error=/beegfs/project/horizon/PDV_Segments/Scripts/Genes_Prediction/std_error_Genes_Prediction.txt


# Définir les chemins
bin_dir="/beegfs/project/horizon/bin/miniconda3/bin/"
db="/beegfs/project/horizon/db/UniRef/UniRef_insecta90_viruses90_bacteria50.fa"
preds_dir="/beegfs/project/horizon/PDV_Segments/Data/Segments/Genes_Proteins/"

# Paramètres spécifiques
species=Cotesia_icipe
genome=/beegfs/project/horizon/PDV_Segments/Data/Segments/Cotesia_icipe_segments.fa

min_length=33
compressed=1
memory_limit=65G

# Construire les chemins pour les fichiers d'entrée et de sortie
species_fa=${preds_dir}${species}/${species}.fa
preds_results=${preds_dir}${species}_predsResults
temp_folder=${preds_dir}${species}_tempFolder


# Exécuter la commande
mkdir {$temp_folder}


${bin_dir}metaeuk easy-predict \
    ${genome} \
    ${db} \
    ${preds_results} \
    ${temp_folder} \
    --min-length ${min_length} \
    --compressed ${compressed} \
    --split-memory-limit ${memory_limit}

#!/bin/bash

# Chemins
SRA_TOOLKIT="/beegfs/data/soft/sratoolkit.2.10.5-ubuntu64/bin"
ACCESSION="ERR3829581"
OUTDIR="/beegfs/project/horizon/PDV_Segments/Data/Reads"
SRA_FILE="$OUTDIR/$ACCESSION/$ACCESSION.sra"

# Étape 1 : Télécharger le .sra avec limite augmentée
#echo ">>> Téléchargement de $ACCESSION"
#"$SRA_TOOLKIT/prefetch" --max-size 50000000 --output-directory "$OUTDIR" "$ACCESSION"

# Vérification du téléchargement
#if [[ ! -f "$SRA_FILE" ]]; then
#  echo "!!! ERREUR : Fichier $SRA_FILE introuvable après prefetch."
#  exit 1
#fi

# Étape 2 : Extraction des fichiers FASTQ
echo ">>> Extraction des FASTQ depuis $SRA_FILE"
"$SRA_TOOLKIT/fasterq-dump" --split-files --threads 8 -O "$OUTDIR" "$SRA_FILE"

# Étape 3 (optionnelle) : suppression du .sra
# echo ">>> Suppression du fichier .sra"
# rm -f "$SRA_FILE"

echo ">>> Terminé avec succès."


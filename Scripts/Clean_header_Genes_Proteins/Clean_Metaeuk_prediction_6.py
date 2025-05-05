#!/usr/bin/env python3
from collections import defaultdict
import re

# Fichier d'entrée et de sortie
input_file = "/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_5.fa"
output_file = "/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_6.fa"

# Fonction pour extraire les informations du header
def parse_header(header):
    match = re.match(r"^(.*?_)(\d+)_(\d+)_(.*?)_([^_]+)$", header)
    if match:
        scaffold, start, end, protein_name, segment = match.groups()
        return {
            "scaffold": scaffold,
            "start": int(start),
            "end": int(end),
            "protein_name": protein_name,
            "segment": segment,
        }
    return None

# Lecture et traitement des séquences
segments = defaultdict(list)
with open(input_file, "r") as infile:
    current_header = None
    current_sequence = []

    for line in infile:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                # Ajoute la séquence précédente au dictionnaire
                segments[current_segment].append({
                    "header": current_header,
                    "sequence": "".join(current_sequence),
                    "start": current_start,
                })

            # Analyse le nouveau header
            header_info = parse_header(line[1:])
            if header_info:
                current_header = line[1:]
                current_segment = header_info["segment"]
                current_start = header_info["start"]
                current_sequence = []
        else:
            current_sequence.append(line)

    # Ajoute la dernière séquence
    if current_header:
        segments[current_segment].append({
            "header": current_header,
            "sequence": "".join(current_sequence),
            "start": current_start,
        })

# Traitement et écriture des séquences modifiées
with open(output_file, "w") as outfile:
    for segment, proteins in segments.items():
        # Trie les protéines par position start
        proteins.sort(key=lambda x: x["start"])

        # Ajoute le rang au header
        total_proteins = len(proteins)
        for idx, protein in enumerate(proteins, start=1):
            new_header = f"{protein['header']}_{idx}/{total_proteins}"
            outfile.write(f">{new_header}\n")
            outfile.write(f"{protein['sequence']}\n")


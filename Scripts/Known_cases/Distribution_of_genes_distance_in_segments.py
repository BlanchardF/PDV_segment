import re
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def parse_fasta_headers(fasta_file):
    """Extrait les headers d'un fichier FASTA et retourne une liste."""
    headers = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line.strip()[1:])  # Retirer le '>'
    return headers

def extract_positions(headers):
    """Extrait les positions start, end et le segment depuis les headers."""
    segment_positions = defaultdict(list)
    
    for header in headers:
        parts = header.split('_')
        if len(parts) >= 5:
            try:
                start, end = int(parts[1]), int(parts[2])
                segment = parts[4]  # 5ᵉ élément (index 4)
                segment_positions[segment].append((start, end))
            except ValueError:
                continue  # Ignore les entrées mal formatées
                
    return segment_positions

def compute_gene_distances(segment_positions):
    """Calcule les distances entre gènes successifs dans chaque segment."""
    distances = []

    for segment, positions in segment_positions.items():
        # Trier les gènes par leur position de début
        positions.sort()
        for i in range(1, len(positions)):
            prev_end = positions[i - 1][1]
            curr_start = positions[i][0]
            if curr_start > prev_end:
                distance = curr_start - prev_end
                distances.append(distance)
            else:
                distances.append(0)  # gènes chevauchants ou collés

    return distances

def plot_distance_distribution(distances, output_file="distances_plot.png"):
    """Affiche et enregistre un histogramme de la distribution des distances."""
    plt.figure(figsize=(10, 6))
    plt.hist(distances, bins=50, color='skyblue', edgecolor='black')
    plt.title("Distribution des distances entre gènes")
    plt.xlabel("Distance (en bases)")
    plt.ylabel("Nombre d'occurrences")
    plt.grid(True)
    plt.tight_layout()
    
    # Sauvegarde
    plt.savefig(output_file, dpi=300)
    print(f"Plot enregistré sous : {output_file}")
    
    plt.show()

    # Affichage des stats
    mean_dist = np.mean(distances)
    median_dist = np.median(distances)
    print(f"Distance moyenne : {mean_dist:.2f} bases")
    print(f"Distance médiane : {median_dist:.2f} bases")

if __name__ == "__main__":
    fasta_file = "/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Clean_Data/Cotesia_congregata.fa"  # Remplace par le fichier voulu
    headers = parse_fasta_headers(fasta_file)
    segment_positions = extract_positions(headers)
    distances = compute_gene_distances(segment_positions)
    plot_distance_distribution(distances, output_file="distances_plot.png")


#!/usr/bin/env python3

import time
import os
import requests

# Fonction pour interroger l'API et obtenir des informations sur un identifiant
def get_info_from_api(identifier):
    url = f"https://rest.uniprot.org/uniprotkb/{identifier}.json"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Lève une exception en cas d'erreur HTTP
        data = response.json()  # Si la réponse est en format JSON
        return data
    except requests.exceptions.RequestException as e:
        print(f"Erreur lors de la récupération des données pour {identifier}: {e}")
        return None

# Remplacez le chemin du fichier FASTA
fasta_file = "../../Data/Genes_Proteins/Cotesia_icipe_predsResults.fa"
# Remplacez par le chemin du fichier de sortie
output_file = "../../Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_1.fa"

# Créer le répertoire de sortie s'il n'existe pas
output_dir = os.path.dirname(output_file)
os.makedirs(output_dir, exist_ok=True)

# Ouvrir et lire le fichier FASTA
with open(fasta_file, "r") as file, open(output_file, "w") as output:
    for line in file:
        if line.startswith(">"):  # Si c'est une ligne d'en-tête
            # Extraire l'identifiant sans le préfixe "UniRef90_" et sans underscore initial
            identifier = line.split("|")[0][9:].lstrip('_')  # Supprimer "UniRef90_" et l'underscore initial

            # Interroger l'API pour cet identifiant
            print(f"Récupération des informations pour l'identifiant {identifier}...")
            api_data = get_info_from_api(identifier)

            if api_data:
                # Extraire le premier identifiant de preuve (s'il existe)
                evidence_ids = api_data.get('organism', {}).get('evidences', [])
                first_evidence_id = evidence_ids[0]['id'] if evidence_ids else "Non disponible"

                # Ajouter l'identifiant de preuve au header
                new_header = f"{line.strip()} {first_evidence_id}"

                # Écrire la nouvelle ligne d'en-tête dans le fichier de sortie
                output.write(new_header + "\n")
            else:
                output.write(line)  # Si aucune donnée API, garder le header original
        else:
            output.write(line)  # Écrire les séquences sans modification

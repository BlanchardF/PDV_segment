#!/usr/bin/env python3

import requests
import xml.etree.ElementTree as ET

def get_protein_name_from_ncbi(accession):
    """
    Récupère le nom de la protéine pour un numéro d'accession NCBI.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "protein",
        "id": accession,
        "retmode": "xml"
    }
    
    response = requests.get(url, params=params)
    if response.status_code == 200:
        root = ET.fromstring(response.text)
        for feature in root.findall(".//GBFeature"):
            if feature.find(".//GBFeature_key").text == "Protein":
                for qualifier in feature.findall(".//GBFeature_quals/GBQualifier"):
                    if qualifier.find(".//GBQualifier_name").text == "product":
                        return qualifier.find(".//GBQualifier_value").text.replace(" ", "-")
    return "Nom-de-protéine-non-trouvé"

def update_fasta_with_protein_names(input_fasta, output_fasta):
    """
    Met à jour les en-têtes d'un fichier FASTA avec les noms de protéines correspondant aux numéros d'accession NCBI.
    """
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                header = line.strip()
                # Extraire le numéro d'accession (2ème élément après un espace)
                accession = header.split()[1]
                protein_name = get_protein_name_from_ncbi(accession)
                updated_header = f"{header} {protein_name}"
                outfile.write(updated_header + "\n")
            else:
                outfile.write(line)

# Chemins d'entrée et de sortie
input_fasta = "/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_2.fa"
output_fasta = "/beegfs/project/horizon/PDV_Segments/Data/Genes_Proteins/Tmp/Cotesia_icipe/Cotesia_icipe_3.fa"

# Exécution de la fonction principale
update_fasta_with_protein_names(input_fasta, output_fasta)

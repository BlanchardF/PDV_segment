#!/usr/bin/env python3

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract segment positions from MMseqs2 .m8 output with 100% coverage.")
    parser.add_argument("-in", "--input", required=True, help="Input .m8 file from mmseqs convertalis")
    parser.add_argument("-out", "--output", required=True, help="Output TSV file with Scaffold, Start, End, Segment_name")
    args = parser.parse_args()

    # Colonnes attendues dans le fichier convertalis
    cols = [
        "query", "qlen", "tlen", "target", "pident", "alnlen", "mismatch",
        "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qaln", "tcov"
    ]

    # Lecture du fichier .m8
    df = pd.read_csv(args.input, sep="\t", names=cols)

    # Filtrage sur tcov = 1.0 (100%)
    df = df[df["tcov"] == 1.0].copy()

    # Extraction et formatage des colonnes
    df["Start"] = df[["qstart", "qend"]].min(axis=1)
    df["End"] = df[["qstart", "qend"]].max(axis=1)

    df_out = df[["query", "Start", "End", "target"]].copy()
    df_out.columns = ["Scaffold", "Start", "End", "Segment_name"]

    # Sauvegarde du fichier final
    df_out.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()


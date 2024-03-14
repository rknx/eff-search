#!/usr/bin/env python3

# Anuj Sharma
# github.com/rknx
# v1, 03/01/2024

# Import required modules
import sys
from Bio import SeqIO

# Usage
if ( len(sys.argv) > 1 ) or ( sys.stdin.isatty() ):
    sys.exit("""gb2fasta.py: Extract all CDS from genbank file in fasta format.
    Usage: gb2fasta.py < input.gb [ > output.fasta ]
    Dependencies: Python3 and SeqIO library.""")

# Read genbank file from stdin
records = list(SeqIO.parse(sys.stdin, "genbank"))

# Initialize fasta dict
seqs = []

# Process genes
for record in records:
    for feature in record.features:
        if feature.type == "CDS":
            # --- Get sequence

            # # 1. For prokaryotes only, this will capture introns
            # seq = record[int(feature.location.start):int(feature.location.end)]
            # if feature.location.strand == -1: seq.seq = seq.seq.reverse_complement()

            # # 3. Iteration over exons for eukaryotes
            # seq = record[0:0]
            # for exon in feature.location.parts:
            #     seq += record[int(exon.start):int(exon.end)]
            # if feature.location.strand == -1: seq.seq = seq.seq.reverse_complement()
            
            # 2. Convenient function that works with both
            seq = feature.location.extract(record)


            # --- Get metadata

            # ID from locustag, use protein.id alternatively
            seq.id = feature.qualifiers['locus_tag'][0]
            if "ribosomal_slippage" in feature.qualifiers: seq.id += ":f"
            if "<" in str(feature.location.start) or ">" in str(feature.location.end): seq.id += ":t"
            
            # Protein functional annotation
            seq.description =  feature.qualifiers['product'][0]

            # --- Put everything together
            seqs.append(seq)

# Write to stdout
SeqIO.write(seqs, sys.stdout, 'fasta')
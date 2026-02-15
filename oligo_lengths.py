#!/usr/bin/env python3

from Bio import SeqIO
import glob

# pattern matching all tiles
pattern = "CVB3_full_proteome_DMS/A.Data/A2.synthetic_oligonucleotides_P2_P3/Tile_*_oligos.fa"

files = sorted(glob.glob(pattern))

lengths = []

for file in files:
    # read only first sequence
    first_record = next(SeqIO.parse(file, "fasta"))
    length = len(first_record.seq)

    lengths.append(length)

    print(f"{file}: {length} bp")

print("\nSummary:")
print(f"Number of tiles: {len(lengths)}")
print(f"Min length: {min(lengths)} bp")
print(f"Max length: {max(lengths)} bp")
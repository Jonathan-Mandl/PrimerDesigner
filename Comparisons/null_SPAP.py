import General.utils as GU
import time
import json
import tracemalloc
from pathlib import Path
import pandas as pd
import networkx as nx
from General.primer_graphs import create_primer_df, create_graph

from General.args import *
import sys
import csv


def export_null_paths_primers_py3(
    paths,
    template_5to3,
    out_csv_path,
    protein_name="SPAP",
    method="NullWeighted"
):
    """
    paths: list of paths, each path is ['s', (start,end,strand), ..., 'd']
    Writes a CSV of all primers from all null paths (Python 3 compatible).
    """

    with open(out_csv_path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "protein_name", "method", "null_path_id",
            "primer_order", "pair_id",
            "start", "end", "strand",
            "primer_seq_5to3", "length"
        ])

        for path_id, path in enumerate(paths):

            # strip source / sink
            core = path
            if core and core[0] == "s":
                core = core[1:]
            if core and core[-1] == "d":
                core = core[:-1]

            pair_id = -1

            for primer_order, node in enumerate(core):
                start, end, strand = node

                seq = primer_seq_from_template(
                    template_5to3, start, end, strand
                ).upper()

                if strand.lower() == "f":
                    pair_id += 1

                writer.writerow([
                    protein_name,
                    method,
                    path_id,
                    primer_order,
                    pair_id,
                    start,
                    end,
                    strand,
                    seq,
                    len(seq)
                ])

def primer_seq_from_template(template_5to3: str, start: int, end: int, strand: str) -> str:
    """
    start/end are 0-based, python slicing semantics [start:end)
    strand: 'f' or 'r'
    returns primer sequence in 5'->3' orientation as ordered.
    """
    subseq = template_5to3[start +len(GU.UPSTREAM_NT):end + len(GU.UPSTREAM_NT)]
    if strand.lower() == "f":
        return subseq
    elif strand.lower() == "r":
        return GU.revcomp(subseq)
    else:
        raise ValueError(f"Unknown strand: {strand}")


def main():
    protein_name = "SPAP"
     # Load config
    GU.init_config("configs/SPAP_experiment.json")

    mutreg_nt = GU.read_fasta("data/SPAP_reference.fa")

    sequence_nt = GU.UPSTREAM_NT + mutreg_nt + GU.DOWNSTREAM_NT

    args = get_args()

    t0 = time.time()

    # Build primer table 
    primer_df = create_primer_df(sequence_nt, args)

    # Build primer graph
    graph = create_graph(primer_df, len(mutreg_nt), args)

    sampled_paths = GU.sample_paths_dag_uniform(graph,  's', 'd', k=1000, max_tries=1000, seed=42)

    export_null_paths_primers_py3(
        paths=sampled_paths,
        template_5to3=sequence_nt,
        out_csv_path="results/null_paths_primers_SPAP.csv",
        protein_name=protein_name,
        method="NullWeighted"
    )
    print("✔ Null primer paths written")

if __name__ == "__main__":
    main()


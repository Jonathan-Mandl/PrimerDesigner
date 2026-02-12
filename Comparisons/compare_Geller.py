import sys
import time
import json
import tracemalloc
import random
from pathlib import Path
import pandas as pd
import networkx as nx

import General.utils as GU
from General.primer_graphs import create_primer_df, create_graph
from General.args import *

def export_primer_set(primer_df, nodes, protein_name, method):
    rows = []
    pair_id = -1

    for order, node in enumerate(nodes):
        r = primer_df.loc[node]

        start, end, strand = node
        seq = r["seq"].upper()
        length = len(seq)

        # Increment pair_id on forward primer
        if strand == "f":
            pair_id += 1

        rows.append({
            "protein_name": protein_name,
            "method": method,
            "primer_order": order,
            "pair_id": pair_id,
            "start": start,
            "end": end,
            "strand": strand,
            "primer_seq_5to3": seq,
            "length": length
        })

    return pd.DataFrame(rows)


def rc(seq: str) -> str:
    tbl = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(tbl)[::-1]


def find_unique(haystack: str, needle: str, prev_pos):
    """Return 0-based index of first occurrence at/after prev_pos; raise if not found."""
    p = haystack.find(needle, prev_pos)
    if p == -1:
        raise SystemExit(f"Primer {needle} not found!")
    return p


GELLER_PRIMERS = {
    1:  ("taccactactaggcaaagcatcact",           "ctcttggacctctactagacctgg"),
    2:  ("gcactacccaatttcgtttgaagga",           "ccgaatgcatttccaagctgttc"),
    3:  ("gaacagggagtgaaggactatgtg",            "cagccatagggattccgtaatattg"),
    4:  ("gtggctcaaacagaaggtgtca",              "gttcctggtcactttgggatgg"),
    5:  ("cacaatcgagcagagcgcg",                 "tttctcagcaagcgaccttcc"),
    6:  ("caagtcggtggcaacaaacttaatt",           "cggtgaggtgaacagaatgcc"),
    7:  ("catggctgccctagaagagaaa",              "aattgaccgggcaacactcatc"),
    8:  ("cccatgtcagtcaagacttgtgac",            "gtgcaacgctaattttgatctctct"),
    9:  ("caccgagatgtttagggagtacaat",           "ccactgacacaaatgtggtcaa"),
    10: ("ggctttcatttgcttacaggca",              "gcccacctgtcatagatgcc"),
    11: ("gaatatggcgagtttaccatgctg",            "gttgggaaacttgctggtgttaat"),
    12: ("ggttaatgaggcagtgctagca",              "caccttgctcatcattgaagtagtg"),
    13: ("cttctcagcagcactcctcaaa",              "gcgtagtggtccactgcttc"),
    14: ("acacgtggatgagtacatgctg",              "cttctctatggacctgagctcat"),
    15: ("cctgaacctaccaatggtgacttat",           "ccagacagggcttaagctagc"),
    16: ("agcatttgattactctgggtacgat",           "caccatatgcgatcatcctgaattg"),
    17: ("gtgtacaaagggattgacttggac",            "ttgattcgtgtatgtctttcatggg"),
    18: ("cttcctggtgcatcctgttatg",              "agcacagtagggttaagccaa"),
}


def main():

    GU.UPSTREAM_NT, GU.DOWNSTREAM_NT, GU.MAX_TM = GU.load_config("configs/geller_experiment.json")

    mutreg_nt = GU.read_fasta("data/geller_reference.fa")

    sequence_nt = GU.UPSTREAM_NT + mutreg_nt + GU.DOWNSTREAM_NT

    sequence_nt = sequence_nt.upper()

    # ---- Args (keep your sys.argv override logic) ----
    sys.argv = [
        sys.argv[0],
        "--file_path", "input_path",
        "--output", "output_path",
    ]
    args = get_args()

    args.oligo_lmin = 240
    args.oligo_lmax = 260

    t0 = time.time()

    # ---- Build primer table ----
    primer_df = create_primer_df(sequence_nt, args)

    competing_primers = []
    f_pos = 0

    for tile in sorted(GELLER_PRIMERS.keys()):
        fwd = GELLER_PRIMERS[tile][0].upper()
        rev = GELLER_PRIMERS[tile][1].upper()
        rc_rev = rc(rev)

        f_pos = find_unique(sequence_nt.upper(), fwd, f_pos)
        r_pos = find_unique(sequence_nt.upper(), rc_rev, f_pos)

        if f_pos is None:
            raise SystemExit(f"No sequence matched fwd primer {tile}!")
        else:
            competing_primers.append(
                (f_pos - len(GU.UPSTREAM_NT), f_pos + len(fwd) - len(GU.UPSTREAM_NT), "f")
            )

        if r_pos is None:
            raise SystemExit(f"No sequence matched rev primer {tile}!")
        else:
            competing_primers.append(
                (r_pos - len(GU.UPSTREAM_NT), r_pos + len(rev) - len(GU.UPSTREAM_NT), "r")
            )

    print("All primers were found in sequence!")

    geller_set = primer_df.loc[competing_primers].copy().reset_index()
    geller_efficiency = float(geller_set["efficiency"].sum())

    # ---- Build graph (time + peak mem) ----
    t_graph0 = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t_graph0
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    # ---- Longest path ----
    longest_path_t0 = time.time()

    full_path = GU.longest_path_dag(graph, "s", "d")
    primer_path_nodes = full_path[1:-1]  # strip s/d

    primer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
    primer_efficiency = float(primer_set["efficiency"].sum())

    longest_path_time = time.time() - longest_path_t0
    total_time = time.time() - t0
    # ---- Create results directory ----
    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- Summary CSV ----
    row = {
        "graph_nodes": len(graph.nodes),
        "graph_edges": len(graph.edges),
        "graph_time_sec": round(graph_time, 3),
        "graph_peak_mem_MB": round(graph_peak_mb, 1),
        "PD_single_efficiency": primer_efficiency,
        "Geller_efficiency": geller_efficiency,
        "PD_single_primers": len(primer_path_nodes),
        "QuickChange_primers": len(competing_primers),
        "shortest_path_time": longest_path_time,
        "total_time_sec": round(total_time, 3),
    }

    df = pd.DataFrame([row])
    summary_path = results_dir / "Geller_method_comparison.csv"
    df.to_csv(summary_path, index=False)

    # ---- Primer lists CSV ----
    pd_nodes = primer_path_nodes
    geller_nodes = competing_primers
    protein_name = "geller_reference"

    df_pd = export_primer_set(
        primer_df=primer_df,
        nodes=pd_nodes,
        protein_name=protein_name,
        method="PD_single",
    )

    df_geller = export_primer_set(
        primer_df=primer_df,
        nodes=geller_nodes,
        protein_name=protein_name,
        method="Geller",
    )

    df_all = pd.concat([df_pd, df_geller], ignore_index=True)

    primer_path = results_dir / "Geller_method_comparison_primers.csv"
    df_all.to_csv(primer_path, index=False)


if __name__ == "__main__":
    main()
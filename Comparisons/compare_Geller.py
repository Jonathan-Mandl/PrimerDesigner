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
import csv


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


def export_null_paths_primers_py3(
    paths,
    template_5to3,
    out_csv_path,
    protein_name="GELLER",
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
    args = get_args()

    args.oligo_lmin = 227
    args.oligo_lmax = 298

    t0 = time.time()

    # ---- Build primer table ----
    primer_df = create_primer_df(sequence_nt, args)

    competing_primers = []
    search_pos = 0  # pointer for searching forward primers

    for tile in sorted(GELLER_PRIMERS.keys()):
        fwd = GELLER_PRIMERS[tile][0].upper()
        rev = GELLER_PRIMERS[tile][1].upper()
        rc_rev = rc(rev)

        # find forward primer starting at search_pos
        f_hit = find_unique(sequence_nt, fwd, search_pos)

        # find reverse primer AFTER the forward primer
        r_hit = find_unique(sequence_nt, rc_rev, f_hit + len(fwd))

        # append correct positions (convert back to coords relative to mutreg start)
        competing_primers.append(
            (f_hit - len(GU.UPSTREAM_NT), f_hit + len(fwd) - len(GU.UPSTREAM_NT), "f")
        )
        competing_primers.append(
            (r_hit - len(GU.UPSTREAM_NT), r_hit + len(rev) - len(GU.UPSTREAM_NT), "r")
        )

        search_pos = r_hit + len(rev)

    print("All primers were found in sequence!")

    geller_set = primer_df.loc[competing_primers].copy().reset_index()
    geller_efficiency = float(geller_set["efficiency"].mean())

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
    primer_efficiency = float(primer_set["efficiency"].mean())

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
        "PD_avg_efficiency": primer_efficiency,
        "Geller_avg_efficiency": geller_efficiency,
        "PD_single_primers": len(primer_path_nodes),
        "QuickChange_primers": len(competing_primers),
        "longest_path_time": longest_path_time,
        "total_time_sec": round(total_time, 3),
    }

    df = pd.DataFrame([row])
    summary_path = results_dir / "Geller_method_comparison.csv"
    df.to_csv(summary_path, index=False)

    # ---- Primer lists CSV ----
    pd_nodes = primer_path_nodes
    geller_nodes = competing_primers
    protein_name = "geller"

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

    primer_path = results_dir / "Geller_method_primers.csv"
    df_all.to_csv(primer_path, index=False)

    # sample 1000 random paths for null distribution 
    sampled_paths = GU.sample_paths_dag_weighted(graph,  's', 'd', k=1000, max_tries=1000, seed=42)

    # save null paths
    export_null_paths_primers_py3(
        paths=sampled_paths,
        template_5to3=sequence_nt,
        out_csv_path="results/null_paths_primers_Geller.csv",
        protein_name=protein_name,
        method="NullWeighted"
    )
    print("✔ Null primer paths written")


if __name__ == "__main__":
    main()
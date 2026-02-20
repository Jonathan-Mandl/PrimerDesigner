import sys
import time
import tracemalloc
from pathlib import Path
import pandas as pd
from General.primer_graphs import create_primer_df, create_graph
import General.utils as GU
from General.args import *
import csv


def export_null_paths_primers_py3(
    paths,
    template_5to3,
    out_csv_path,
    protein_name="SPAP",
    method="NullWeighted",
    cfg = None
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
                    template_5to3, start, end, strand, cfg
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

def primer_seq_from_template(template_5to3: str, start: int, end: int, strand: str, cfg) -> str:
    """
    start/end are 0-based, python slicing semantics [start:end)
    strand: 'f' or 'r'
    returns primer sequence in 5'->3' orientation as ordered.
    """
    subseq = template_5to3[start +len(cfg.upstream):end + len(cfg.upstream)]
    if strand.lower() == "f":
        return subseq
    elif strand.lower() == "r":
        return GU.revcomp(subseq)
    else:
        raise ValueError(f"Unknown strand: {strand}")



# -----------------------------
# Export function (same as Geller/PrimalScheme)
# -----------------------------
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


def main():
    
    # Load config
    cfg = GU.load_config("configs/SPAP_experiment.json")

    # ---- Input sequence ----
    mutreg_nt = GU.read_fasta("data/SPAP_reference.fa")

    sequence_nt = cfg.upstream + mutreg_nt + cfg.downstream

    protein_name = "SPAP"

    args = get_args()

    t0 = time.time()

    # ---- Build primer table and graph ----
    primer_df = create_primer_df(sequence_nt, args, cfg)

    t_graph0 = time.time()
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_time = time.time() - t_graph0
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    longest_path_t0 = time.time()

    # find longest path in the graph (which corresponds to the primer set)
    full_path = GU.longest_path_dag(graph, "s", "d")
    primer_path_nodes = full_path[1:-1]

    # select the primers from the DataFrame that correspond to the longest path nodes
    primer_designer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
    primer_designer_efficiency = float(primer_designer_set["efficiency"].mean())

    longest_path_time = time.time() - longest_path_t0
    total_time = time.time() - t0

    # ---- Load QuickChange primers (from paper) ----
    quick_primers = [
        (31, 55, "f"), (211, 230, "r"), (183, 208, "f"), (357, 382, "r"),
        (322, 343, "f"), (499, 518, "r"), (471, 492, "f"), (644, 669, "r"),
        (612, 638, "f"), (784, 808, "r"), (752, 773, "f"), (925, 947, "r"),
        (892, 919, "f"), (1066, 1088, "r"), (988, 1009, "f"), (1156, 1180, "r"),
        (1132, 1154, "f"), (1303, 1324, "r"), (1261, 1281, "f"),
        (1432, 1459, "r"), (1406, 1427, "f"), (1563, 1586, "r"),
        (1529, 1550, "f"), (1639, 1662, "r"), (1613, 1637, "f"),
        (1730, 1751, "r"),
    ]

    for start, end, strand in quick_primers:
        quick_primers.append((start - len(cfg.upstream)-1, end - len(cfg.upstream), strand))

    # select only those primers from the graph that match the QuickChange primers
    quick_set = primer_df.loc[quick_primers].copy().reset_index()
    quick_efficiency = float(quick_set["efficiency"].mean())

    # ---- Load PrimalScheme primers (bed file) ----
    primal_scheme_primers = []
    with open("Comparisons/PrimalScheme_primer.bed") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[5]

            direction = "f" if strand == "+" else "r"

            primal_scheme_primers.append(
                (start - len(cfg.upstream), end - len(cfg.upstream), direction)
            )

    # select primers from primer_df that match the PrimalScheme's primers
    primal_scheme_set = primer_df.loc[primal_scheme_primers].copy().reset_index()
    primal_scheme_efficiency = float(primal_scheme_set["efficiency"].mean())

    # ---- Create results directory ----
    results_dir = Path("Results")
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- Summary CSV ----
    row = {
        "protein_name": protein_name,
        "graph_nodes": len(graph.nodes),
        "graph_edges": len(graph.edges),
        "graph_time_sec": round(graph_time, 3),
        "graph_peak_mem_MB": round(graph_peak_mb, 1),
        "PD_single_avg_efficiency": primer_designer_efficiency,
        "QuickChange_avg_efficiency": quick_efficiency,
        "PrimalScheme_avg_efficiency": primal_scheme_efficiency,
        "PD_single_primers": len(primer_path_nodes),
        "QuickChange_primers": len(quick_primers),
        "PrimalScheme_primers": len(primal_scheme_primers),
        "shortest_path_time": longest_path_time,
        "total_time_sec": round(total_time, 3),
    }

    pd.DataFrame([row]).to_csv(results_dir / "SpAP_comparison.csv", index=False)

    # ---- Primer lists CSV ----
    df_pd = export_primer_set(
        primer_df=primer_df,
        nodes=primer_path_nodes,
        protein_name=protein_name,
        method="PD_single",
    )

    df_quick = export_primer_set(
        primer_df=primer_df,
        nodes=quick_primers,
        protein_name=protein_name,
        method="QuickChange",
    )

    df_primal = export_primer_set(
        primer_df=primer_df,
        nodes=primal_scheme_primers,
        protein_name=protein_name,
        method="PrimalScheme",
    )

    pd.concat([df_pd, df_quick, df_primal], ignore_index=True).to_csv(
        results_dir / "SPAP_primers.csv",
        index=False,
    )

    sampled_paths = GU.sample_paths_dag_uniform(graph,  's', 'd', k=1000, max_tries=1000, seed=42)

    export_null_paths_primers_py3(
        paths=sampled_paths,
        template_5to3=sequence_nt,
        out_csv_path=results_dir / "null_paths_primers_SPAP.csv",
        protein_name=protein_name,
        method="NullWeighted", cfg=cfg
    )
    print("✔ Null primer paths written")

if __name__ == "__main__":
    main()
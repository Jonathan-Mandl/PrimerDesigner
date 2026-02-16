import sys
import time
import tracemalloc
from pathlib import Path
import pandas as pd
from General.primer_graphs import create_primer_df, create_graph
import General.utils as GU
from General.args import *

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
    cfg = GU.load_config("config.json")

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

    full_path = GU.longest_path_dag(graph, "s", "d")
    primer_path_nodes = full_path[1:-1]

    primer_designer_set = primer_df.loc[primer_path_nodes].copy().reset_index()
    primer_designer_efficiency = float(primer_designer_set["efficiency"].mean())

    longest_path_time = time.time() - longest_path_t0
    total_time = time.time() - t0

    # ---- QuickChange primers (from paper) ----
    primer_locations = [
        (31, 55, "f"), (211, 230, "r"), (183, 208, "f"), (357, 382, "r"),
        (322, 343, "f"), (499, 518, "r"), (471, 492, "f"), (644, 669, "r"),
        (612, 638, "f"), (784, 808, "r"), (752, 773, "f"), (925, 947, "r"),
        (892, 919, "f"), (1066, 1088, "r"), (988, 1009, "f"), (1156, 1180, "r"),
        (1132, 1154, "f"), (1303, 1324, "r"), (1261, 1281, "f"),
        (1432, 1459, "r"), (1406, 1427, "f"), (1563, 1586, "r"),
        (1529, 1550, "f"), (1639, 1662, "r"), (1613, 1637, "f"),
        (1730, 1751, "r"),
    ]

    quick_primers = []
    for start, end, strand in primer_locations:
        quick_primers.append((start - len(cfg.upstream)-1, end - len(cfg.upstream), strand))

    quick_set = primer_df.loc[quick_primers].copy().reset_index()
    quick_efficiency = float(quick_set["efficiency"].mean())

    # ---- Create results directory ----
    results_dir = Path("results")
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
        "PD_single_primers": len(primer_path_nodes),
        "QuickChange_primers": len(quick_primers),
        "shortest_path_time": longest_path_time,
        "total_time_sec": round(total_time, 3),
    }

    pd.DataFrame([row]).to_csv(results_dir / "QuickChange_comparison.csv", index=False)

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

    pd.concat([df_pd, df_quick], ignore_index=True).to_csv(
        results_dir / "QuickChange_comparison_primers.csv",
        index=False,
    )

if __name__ == "__main__":
    main()
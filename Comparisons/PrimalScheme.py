import sys
import time
import tracemalloc
from pathlib import Path
import pandas as pd
from General.primer_graphs import create_primer_df, create_graph
import General.utils as GU
from General.args import *


# -----------------------------
# Export function (same as Geller)
# -----------------------------
def export_primer_set(primer_df, nodes, protein_name, method):
    rows = []
    pair_id = -1

    for order, node in enumerate(nodes):
        r = primer_df.loc[node]

        start, end, strand = node

        seq = r["seq"].upper()
        length = len(seq)

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

    # ---- Config ----
    GU.UPSTREAM_NT, GU.DOWNSTREAM_NT, GU.MAX_TM = GU.load_config("configs/SPAP_experiment.json")

    # ---- Input sequence ----
    mutreg_nt = GU.read_fasta("data/SPAP_reference.fa")

    sequence_nt = GU.UPSTREAM_NT + mutreg_nt + GU.DOWNSTREAM_NT
    protein_name = "SPAP"

    # ---- Args ----
    sys.argv = [sys.argv[0], "--file_path", "input_path", "--output", "output_path"]
    args = get_args()

    t0 = time.time()

    # ---- Build primer table ----
    primer_df = create_primer_df(sequence_nt, args)

    # ---- Load PrimalScheme primers ----
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
                (start - len(GU.UPSTREAM_NT), end - len(GU.UPSTREAM_NT), direction)
            )

    primal_scheme_set = primer_df.loc[primal_scheme_primers].copy().reset_index()
    primal_scheme_efficiency = float(primal_scheme_set["efficiency"].sum())

    # ---- Graph + PD_single ----
    tracemalloc.start()
    graph = create_graph(primer_df, len(mutreg_nt), args)
    graph_peak_mb = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()

    full_path = GU.longest_path_dag(graph, "s", "d")
    PD_single_path = full_path[1:-1]

    primer_set = primer_df.loc[PD_single_path].copy().reset_index()
    PD_single_efficiency = float(primer_set["efficiency"].sum())

    total_time = time.time() - t0

    # ---- Create results directory ----
    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- Save summary ----
    summary = {
        "protein_name": protein_name,
        "PD_single_efficiency": PD_single_efficiency,
        "PrimalScheme_efficiency": primal_scheme_efficiency,
        "PD_single_primers": len(PD_single_path),
        "PrimalScheme_primers": len(primal_scheme_primers),
        "total_time_sec": round(total_time, 3),
        "graph_peak_mem_MB": round(graph_peak_mb, 1),
    }

    pd.DataFrame([summary]).to_csv(
        results_dir / "PrimalScheme_comparison.csv",
        index=False,
    )

    # ---- Save primer lists ----
    df_pd = export_primer_set(
        primer_df,
        PD_single_path,
        protein_name,
        "PD_single",
    )

    df_primal = export_primer_set(
        primer_df,
        primal_scheme_primers,
        protein_name,
        "PrimalScheme",
    )

    pd.concat([df_pd, df_primal], ignore_index=True).to_csv(
        results_dir / "PrimalScheme_comparison_primers.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
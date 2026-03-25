import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


# ── Step 1: Load ──────────────────────────────────────────────────────────────

def load_compounds(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    for col in ("compound_name", "smiles"):
        if col not in df.columns:
            sys.exit(f"ERROR: CSV must contain a '{col}' column.")

    if "pic50" not in df.columns and "ic50_nm" not in df.columns:
        sys.exit("ERROR: CSV must contain either 'pic50' or 'ic50_nm' column.")

    total = len(df)

    # Normalize activity
    if "pic50" not in df.columns:
        df["pic50"] = np.nan

    if "ic50_nm" in df.columns:
        ic50 = pd.to_numeric(df["ic50_nm"], errors="coerce")
        valid_ic50 = ic50 > 0
        computed = np.where(valid_ic50, 9 - np.log10(ic50.where(valid_ic50)), np.nan)
        df["pic50"] = df["pic50"].fillna(pd.Series(computed, index=df.index))

    df["pic50"] = pd.to_numeric(df["pic50"], errors="coerce")

    # Validate SMILES
    mols, smiles_valid = [], []
    for _, row in df.iterrows():
        smi = str(row["smiles"]).strip() if pd.notna(row["smiles"]) else ""
        mol = Chem.MolFromSmiles(smi) if smi else None
        mols.append(mol)
        smiles_valid.append(mol is not None)

    df["mol"] = mols
    df["_smiles_valid"] = smiles_valid
    n_bad_smiles = (~df["_smiles_valid"]).sum()
    n_bad_activity = df["_smiles_valid"] & df["pic50"].isna()
    n_bad_activity = n_bad_activity.sum()

    df = df[df["_smiles_valid"] & df["pic50"].notna()].copy()
    df = df[["compound_name", "smiles", "pic50", "mol"]].reset_index(drop=True)

    print(f"Loaded {total} rows, {len(df)} valid, "
          f"{n_bad_smiles} invalid/empty SMILES skipped, "
          f"{n_bad_activity} missing activity dropped.")
    return df


# ── Step 2: Tanimoto ──────────────────────────────────────────────────────────

def compute_tanimoto(df: pd.DataFrame, radius: int = 2, nbits: int = 2048) -> np.ndarray:
    fps = [
        AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        for mol in df["mol"]
    ]
    n = len(fps)
    matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        matrix[i, i] = 1.0
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            matrix[i, j] = sim
            matrix[j, i] = sim
    return matrix


# ── Step 3: SALI ──────────────────────────────────────────────────────────────

def compute_sali(tanimoto: np.ndarray, activities: np.ndarray) -> np.ndarray:
    n = len(activities)
    sali = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            t = tanimoto[i, j]
            if np.isclose(t, 1.0, atol=1e-12):
                sali[i, j] = np.nan
                sali[j, i] = np.nan
            else:
                val = abs(activities[i] - activities[j]) / (1.0 - t)
                sali[i, j] = val
                sali[j, i] = val
    return sali


# ── Step 4: Find cliffs ───────────────────────────────────────────────────────

def find_cliffs(
    sali: np.ndarray,
    labels: list[str],
    tanimoto: np.ndarray,
    activities: np.ndarray,
    threshold: float | None,
) -> pd.DataFrame:
    rows = []
    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            s = sali[i, j]
            if not np.isfinite(s):
                continue
            rows.append({
                "compound_i":  labels[i],
                "compound_j":  labels[j],
                "tanimoto":    round(tanimoto[i, j], 4),
                "pic50_i":     round(activities[i], 3),
                "pic50_j":     round(activities[j], 3),
                "delta_pic50": round(abs(activities[i] - activities[j]), 3),
                "sali":        round(s, 3),
            })

    df = pd.DataFrame(rows).sort_values("sali", ascending=False).reset_index(drop=True)

    if threshold is None:
        finite = df["sali"].dropna()
        if len(finite) > 1:
            threshold = finite.mean() + finite.std(ddof=1)
        elif len(finite) == 1:
            threshold = finite.iloc[0]
        else:
            threshold = np.inf

    df["is_cliff"] = df["sali"] > threshold
    return df, threshold


# ── Step 5: Heatmap ───────────────────────────────────────────────────────────

def plot_heatmap(sali: np.ndarray, labels: list[str], output_path: Path) -> None:
    n = len(labels)
    if n < 2:
        print("  Skipping heatmap: need at least 2 valid compounds.")
        return

    df_plot = pd.DataFrame(sali, index=labels, columns=labels)
    mask = np.isnan(sali)

    annot = n <= 15
    figsize = max(6, n * 1.1)
    fig, ax = plt.subplots(figsize=(figsize, figsize))

    # Gray background for NaN cells
    ax.set_facecolor("#cccccc")

    sns.heatmap(
        df_plot,
        mask=mask,
        ax=ax,
        cmap="YlOrRd",
        annot=annot,
        fmt=".1f" if annot else "",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "SALI score"},
    )
    ax.set_title("SALI Score Matrix — Activity Cliffs", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


# ── Step 6: CLI ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Activity cliff detection via SALI score.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True,
                        help="CSV with compound_name, smiles, and pic50 or ic50_nm")
    parser.add_argument("--radius", type=int, default=2, help="Morgan radius")
    parser.add_argument("--nbits", type=int, default=2048, help="Fingerprint bits")
    parser.add_argument("--cliff-threshold", type=float, default=None,
                        help="SALI threshold for cliff flag (default: mean + 1 SD)")
    parser.add_argument("--output-dir", default="output", help="Output directory")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load
    df = load_compounds(Path(args.input))
    n = len(df)
    if n == 0:
        sys.exit("ERROR: No valid compounds after filtering.")
    labels = df["compound_name"].tolist()
    activities = df["pic50"].to_numpy()

    # Tanimoto
    print(f"Computing Morgan fingerprints (radius={args.radius}, nBits={args.nbits})...")
    tanimoto = compute_tanimoto(df, radius=args.radius, nbits=args.nbits)

    # SALI
    print(f"Computing {n}x{n} SALI matrix...")
    sali = compute_sali(tanimoto, activities)

    # Save SALI matrix
    matrix_path = output_dir / "sali_matrix.csv"
    pd.DataFrame(sali, index=labels, columns=labels).to_csv(
        matrix_path, float_format="%.4f"
    )
    print(f"  Saved: {matrix_path}")

    # Find cliffs
    cliffs_df, threshold = find_cliffs(
        sali, labels, tanimoto, activities, args.cliff_threshold
    )
    cliffs_path = output_dir / "top_cliffs.csv"
    cliffs_df.to_csv(cliffs_path, index=False)
    print(f"  Saved: {cliffs_path}")

    # Heatmap
    plot_heatmap(sali, labels, output_dir / "sali_heatmap.png")

    # Summary
    n_cliffs = cliffs_df["is_cliff"].sum()
    print(f"\nSummary:")
    print(f"  Compounds : {n}")
    print(f"  Pairs     : {len(cliffs_df)}")
    print(f"  Threshold : {threshold:.3f} (SALI > this = cliff)")
    print(f"  Cliffs    : {n_cliffs}")

    print(f"\nAll pairs (ranked by SALI):")
    for _, row in cliffs_df.iterrows():
        cliff_flag = " *** CLIFF" if row["is_cliff"] else ""
        print(f"  {row['compound_i']:<20} vs {row['compound_j']:<20} "
              f"Tanimoto={row['tanimoto']:.3f}  "
              f"ΔSALI={row['sali']:.2f}  "
              f"ΔpIC50={row['delta_pic50']:.2f}{cliff_flag}")


if __name__ == "__main__":
    main()

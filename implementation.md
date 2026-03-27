# Phase 07 — Activity Cliff Detection (SALI Score)

**Version:** 1.1 (final as-built)
**Author:** Kerwyn Medrano
**Date:** 2026-03-25
**Track:** Track 1 — Cheminformatics Core
**Tier:** Micro (2–3 hrs)
**API Cost:** $0.00 — pure RDKit + pandas + seaborn + matplotlib
**Status:** ✅ Complete
**Repo:** `activity-cliffs/`

## As-Built Results

| Pair | Tanimoto | ΔpIC50 | SALI | Cliff? |
|---|---|---|---|---|
| adagrasib vs divarasib | 0.295 | 2.34 | **3.32** | ✅ YES |
| divarasib vs olomorasib | 0.246 | 2.33 | 3.09 | — |
| sotorasib vs divarasib | 0.336 | 1.80 | 2.71 | — |
| sotorasib vs olomorasib | 0.233 | 0.53 | 0.69 | — |
| sotorasib vs adagrasib | 0.153 | 0.54 | 0.64 | — |
| adagrasib vs olomorasib | 0.164 | 0.01 | 0.01 | — |

Auto threshold (mean + 1 SD): **3.196**. 1 cliff flagged from 6 pairs.

**Key insight:** Divarasib is 220× more potent than adagrasib (pIC50 11.54 vs 9.20) despite Tanimoto = 0.295. The CF₃-biaryl back-pocket substituent unique to divarasib is the primary structural driver of this cliff.

**Deviation from plan:** `seaborn` added to dependencies (heatmap NaN masking cleaner via `sns.heatmap` + `mask=` parameter than pure matplotlib). `scipy` added transitively.

---

---

## 1. Project Overview

### Goal

Identify activity cliffs in a compound set: pairs where a small structural
change causes a large potency jump. Uses the Structure-Activity Landscape Index
(SALI) score.

```
SALI(i,j) = |activity_i - activity_j| / (1 - Tanimoto(i,j))
```

High SALI = cliff. Low SALI = smooth SAR. Undefined when Tanimoto = 1.0.

```bash
python main.py --input data/kras_g12c_inhibitors.csv
```

Outputs:
- `output/sali_matrix.csv` — N×N SALI score matrix
- `output/sali_heatmap.png` — heatmap of SALI scores
- `output/top_cliffs.csv` — ranked cliff pairs

### What This Phase Teaches

| Concept | Detail |
|---|---|
| SALI formula | `|ΔpIC50| / (1 - Tanimoto)` — combines similarity + activity |
| Activity in pIC50 | Log-transform IC50 (nM) before computing differences |
| Cliff threshold | Typically SALI > mean + 1 SD flags a cliff |
| Upper-triangle iteration | Avoid double-counting pairs |
| Annotated heatmap | Overlay SALI values on similarity-ordered axes |

### Domain Context

Activity cliffs are medicinal chemistry's most instructive data points: they
reveal which structural features are load-bearing for potency. For KRAS G12C:
- The warhead (acrylamide vs fluorovinyl) is a known cliff driver
- Back-pocket substituents can cause 10–100x potency differences with 1-atom changes

---

## 2. Architecture

```
activity-cliffs/
├── main.py
├── requirements.txt
├── data/
│   └── kras_g12c_inhibitors.csv  (compound_name, smiles, pic50)
└── output/
    ├── sali_matrix.csv
    ├── sali_heatmap.png
    └── top_cliffs.csv
```

---

## Key Concepts

- **Structure-Activity Landscape Index (SALI)** quantifies activity cliffs as `|delta_pIC50| / (1 - Tanimoto)` -- high values indicate pairs where a small structural change causes a large potency shift
- **pIC50 transformation** converts IC50 values to a log scale (`pIC50 = -log10(IC50_M)`) so that activity differences are additive and comparable across orders of magnitude
- **Tanimoto similarity** from Morgan fingerprints measures structural resemblance; SALI amplifies activity differences between structurally similar compounds (denominator near zero)
- **Automatic cliff thresholding** using mean + 1 standard deviation of SALI scores provides an objective cutoff without manual tuning, though domain experts may override it
- **Upper-triangle iteration** over the N x N SALI matrix avoids double-counting symmetric pairs and excludes the diagonal (self-comparisons)
- **Divergence guard** is required when Tanimoto = 1.0 (identical compounds), where SALI is undefined -- these cells are set to NaN rather than infinity

---

## 3. Input CSV — Added Column

```
compound_name, smiles, pic50
```

`pic50` = -log10(IC50 in molar). If raw IC50 in nM is provided, convert:
`pIC50 = 9 - log10(IC50_nM)`.

---

## 4. Module Specification

### `load_compounds(path)` → pd.DataFrame
- Require: compound_name, smiles, pic50
- Validate SMILES, drop invalid
- Coerce pic50 to float, drop NaN

### `compute_tanimoto(df, radius, nbits)` → np.ndarray
- Morgan fingerprints → pairwise Tanimoto (same as Phase 05)

### `compute_sali(tanimoto_matrix, activities)` → np.ndarray
- N×N matrix; cell (i,j) = |activity_i - activity_j| / (1 - tanimoto[i,j])
- Set diagonal = 0
- Set cell = NaN where tanimoto = 1.0 (identical compounds)

### `find_cliffs(sali_matrix, labels, activities, threshold)` → pd.DataFrame
- Extract upper triangle, sort by SALI descending
- Flag as cliff if SALI > threshold (default: mean + 1 SD)
- Return: compound_i, compound_j, tanimoto, delta_pic50, sali, is_cliff

### `plot_heatmap(sali_matrix, labels, output_path)`
- Seaborn heatmap (not clustermap — preserve compound order)
- Annotate with SALI values if N ≤ 15

### `main()`
- `--input` (required)
- `--radius` (default 2), `--nbits` (default 2048)
- `--cliff-threshold` (default: auto = mean + 1 SD)
- `--output-dir` (default: output)

---

## 5. Verification Checklist

```bash
python main.py --input data/kras_g12c_inhibitors.csv

# Expected:
# - N×N SALI matrix computed
# - top_cliffs.csv: pairs ranked by SALI
# - Highest SALI pair flagged as cliff
# - sali_heatmap.png saved
```

---

## 6. Risks / Assumptions / Next Step

**Risks:**
- pIC50 values from different assays are not directly comparable — note in README
- SALI diverges for near-identical compounds (Tanimoto → 1) — NaN guard required

**Next step:** Phase 08 — Matched molecular pair (MMP) analysis.

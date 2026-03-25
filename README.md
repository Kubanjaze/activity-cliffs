# activity-cliffs — Phase 07

Identify activity cliffs in a compound set using the Structure-Activity Landscape Index (SALI).

## What it does

```
SALI(i,j) = |pIC50_i - pIC50_j| / (1 - Tanimoto(i,j))
```

High SALI = small structural change, large potency jump = activity cliff.
Low SALI = smooth SAR = structurally similar compounds have similar potency.

Pairs are flagged as cliffs when SALI > mean + 1 SD (auto threshold, overridable).

## Inputs / Outputs

Input CSV required columns:
- `compound_name`
- `smiles`
- `pic50` (preferred) OR `ic50_nm` (auto-converted: `pic50 = 9 - log10(ic50_nm)`)

Outputs:
- `output/sali_matrix.csv` — N×N SALI matrix (NaN on diagonal + identical pairs)
- `output/top_cliffs.csv` — all pairs ranked by SALI, with `is_cliff` flag
- `output/sali_heatmap.png` — annotated heatmap

## Setup

### RDKit (recommended via conda)
```bash
conda create -n activity-cliffs python=3.11 rdkit -c conda-forge
conda activate activity-cliffs
pip install -r requirements.txt
```

### pip only
```bash
pip install rdkit
pip install -r requirements.txt
```

## Run

```bash
python main.py --input data/kras_g12c_inhibitors.csv
python main.py --input data/kras_g12c_inhibitors.csv --cliff-threshold 15.0
python main.py --help
```

## Notes

- **Assay comparability**: pIC50 values from different assays (biochemical vs cellular) are not directly comparable. Cliffs detected across assay types should be interpreted with caution.
- **SALI singularity**: when two compounds have identical Morgan fingerprints (Tanimoto = 1.0), SALI is undefined → stored as NaN.
- **Empty SMILES / missing activity**: rows are skipped without crashing.
- **O(N²) complexity**: fast for any practical library size.

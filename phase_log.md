# Phase 07 — Activity Cliff Detection (SALI Score)
## Phase Log

**Status:** ✅ Complete
**Started:** 2026-03-25

---

## Log

### 2026-03-25 — Phase started
- Implementation plan written
- Build complete. 4 compounds, 6 pairs, 1 cliff detected
- Cliff: adagrasib vs divarasib (SALI=3.32, ΔpIC50=2.34, Tanimoto=0.295)
- Key insight: divarasib's 2.34 log unit potency advantage over adagrasib with only 0.295 Tanimoto similarity — the trifluoromethyl biaryl back-pocket substituent is a primary driver
- Committed to git

### 2026-03-26 — Documentation update
- Added "Key Concepts" section to implementation.md (6 bullet points covering SALI formula, pIC50 transformation, Tanimoto similarity, automatic thresholding, upper-triangle iteration, divergence guard)

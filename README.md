# anuran-ccm-ecology

Data and R code for the analyses of phylogenetic structure and acoustic niche overlap of anuran communities in 18 landscapes of the Cantareira-Mantiqueira Corridor (Atlantic Forest, Brazil).

Manuscript: *Forest cover is related to phylogenetic diversity but not to acoustic niche overlap of anuran communities in the Atlantic Forest, Brazil* (under review).

## Structure

```
data/       input data
scripts/    analysis scripts (run in numeric order)
results/    tables produced by the scripts
figures/    Figure 2
```

## Data

| File | Content |
|---|---|
| `comm.txt` | abundance matrix, 18 landscapes x 40 species (indirect acoustic abundance) |
| `phy.tre` | consensus tree (Newick); pruned to the 40 recorded species by the scripts |
| `traits.txt` | mean of four call parameters per species (90% duration, 5% frequency, 95% frequency, maximum frequency) |
| `calls_raw.csv` | call-level measurements (800 calls) behind `traits.txt` |
| `landscape_variables.csv` | variables of each 2-km landscape (see below) |
| `distance_matrix_km.csv` | distances between landscapes (km) |

Columns of `landscape_variables.csv`: `Cobertura_pct` (forest cover, %), `Heterogeneidade_SHDI` (Shannon index of land-cover classes), `Isolamento_m` (mean edge-to-edge distance to the nearest forest fragment, m), `Densidade_borda` (forest edge density: perimeter of forest fragments / buffer area). The water-body columns are not used in the analyses.

In `calls_raw.csv` one *Ischnocnema guentheri* call has the same value in the three frequency columns; it was excluded from the frequency means in `traits.txt`.

## Scripts

| Script | Content |
|---|---|
| `01_ses.R` | SES of PD, MPD and MNTD (picante, `taxa.labels` null model, 999 randomisations) |
| `02_univariate_models.R` | univariate models, multiple-testing corrections, Figure 2 |
| `03_model_selection.R` | AICc selection among the 11 candidate models per index |
| `04_diagnostics.R` | Shapiro-Wilk, Breusch-Pagan, HC3 and Moran's I of the top models |
| `05_acoustic_overlap.R` | Horn-Morisita overlap per landscape and regressions for hypothesis 3 |
| `06_phylogenetic_signal.R` | Blomberg's K of the call parameters |
| `run_all.R` | runs the scripts above |

Run from the project root:

```
Rscript scripts/run_all.R
```

`01_ses.R` needs the `picante` package. The SES results are also stored in `results/`, so scripts 02 to 06 run without it (`run_all.R` skips step 1 when `picante` is not installed). Other packages: `ape`, `lmtest`, `sandwich`, `ggplot2`, `patchwork`.

## Notes

- Random permutations use `set.seed(123)`. p-values close to 0.05 (for example K of the 90% duration) can change slightly with the seed.
- In `05_acoustic_overlap.R`, each call parameter is rescaled to 0-1 across species and a constant (`eps = 0.001`) is added before computing the Horn-Morisita index.
- SES.PD is computed with `include.root = FALSE`.

## Citation

Ferreira LF et al. (in prep.).

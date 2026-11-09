# HiFiTaxa_analysis

Benchmark harness, scoring scripts and figure code for the HiFiTaxa manuscript
(Garritano, Modolon and Warren, *HiFiTaxa: a reproducible analysis pipeline for full
length 16S rRNA gene and ITS amplicon sequencing from PacBio reads*). The pipeline
itself lives at https://github.com/alegarritano/HiFiTaxa; this repository holds
everything that was used to benchmark it and to make the manuscript's tables and
figures.

## Layout

| folder | contents |
|---|---|
| `benchmark/` | the job scripts that ran the four benchmark runs (`pbs/`) and the holdout splitters; see [benchmark/README.md](benchmark/README.md) |
| `analysis/` | the R scripts that score the runs and make the manuscript's tables and figures, and the identity job; see [analysis/README.md](analysis/README.md) |
| `data/` | the species list of the fungal mock that the scripts read as the truth (`fungi_103_species.txt`) and, in `identity/`, the identity of every holdout test sequence to its closest training sequence (output of `analysis/identity.pbs`, 13 October 2026) |
| `results/tables/` | the CSV tables the scripts produce from the four runs, the source of Supplementary Tables S1 to S12 |
| `manuscript/` | the eight figure files as submitted (`figures/`: Figure 1, Figure 2 and Supplementary Figures S1 to S6) and the supplementary tables workbook; Figure 1 is a hand drawn schematic, not produced by a script |

## Reproducing the tables and figures

1. Clone this repository. Install R 4.5 with ggplot2, dplyr, tidyr, stringr, readr,
   patchwork and scales.
2. Obtain the four pipeline run directories, which are not part of this repository
   (`results_atcc_16s`, `results_fungi_its`, `benchmark_holdout_16s_seed42`,
   `benchmark_holdout_its_seed42`; produced with the job scripts in `benchmark/pbs/`
   as described in `analysis/README.md`, or available from the corresponding author),
   and point `HIFITAXA_DATA` at the folder that holds them.
3. Run the scripts in `analysis/` in the order given in `analysis/README.md`, with
   `IDENT_DIR=data/identity` for script 16. Outputs go to `HIFITAXA_OUT`
   (default `analysis/regenerated/`); the tables should match `results/tables/` and the
   figures, once their fonts are embedded, the files in `manuscript/figures/`.

## Licence and citation

MIT (see `LICENSE`). Please cite the HiFiTaxa manuscript and this repository
(`CITATION.cff`).

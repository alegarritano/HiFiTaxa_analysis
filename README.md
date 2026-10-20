# HiFiTaxa_analysis

Benchmark harness, scoring scripts and figure code for the HiFiTaxa manuscript
(Garritano, Modolon and Warren, *HiFiTaxa: a reproducible analysis pipeline for full
length 16S rRNA gene and ITS amplicon sequencing from PacBio reads*). The pipeline
itself lives at https://github.com/alegarritano/HiFiTaxa; this repository holds
everything that was used to benchmark it and to make the manuscript's tables and
figures. It was split out of the pipeline repository on 20 October 2026 and the files
carry their history from there.

## Layout

| folder | contents |
|---|---|
| `benchmark/` | the Gadi job scripts that ran the four benchmark runs (`pbs/`), the holdout splitters, the reference depletion and primer scan tools and the mock exclusion lists; see [benchmark/README.md](benchmark/README.md) |
| `analysis/` | the R scripts `00_common.R` to `17_mock_per_library.R` that score the runs and draw the figures, and the Gadi identity jobs; see [analysis/README.md](analysis/README.md) |
| `data/identity/` | identity of every holdout test sequence to its closest training sequence (output of `analysis/identity_resume.pbs`, 13 October 2026) |
| `results/tables/` | the CSV tables the scripts produce from the four runs, the source of Supplementary Tables S1 to S13 |
| `results/figures/` | the figures the scripts produce, among them Figure 2 (`fig2_main`) and the supplementary figures |
| `results/logs/` | the console logs of the scoring runs behind the manuscript |
| `manuscript/` | the figure files and the supplementary tables workbook as submitted (Figure 1 is a hand drawn schematic, not produced by a script) |

## Reproducing the tables and figures

1. Clone this repository. Install R 4.5 with ggplot2, dplyr, tidyr, stringr, readr,
   patchwork and scales (ggalluvial, cowplot, ggplotify and ggbreak for some figures).
2. Obtain the four pipeline run directories, which are not part of this repository
   (`results_atcc_16s`, `results_fungi_its`, `benchmark_holdout_16s_seed42`,
   `benchmark_holdout_its_seed42`; produced with the job scripts in `benchmark/pbs/`
   as described in `analysis/README.md`, or available from the corresponding author),
   and point `HIFITAXA_DATA` at the folder that holds them.
3. Run the scripts in `analysis/` in the order given in `analysis/README.md`, with
   `IDENT_DIR=data/identity` for script 16. Outputs go to `HIFITAXA_OUT`
   (default `analysis/regenerated/`) and should match `results/`.

## The runs behind the manuscript

| run | data | pipeline version | date |
|---|---|---|---|
| 16S rRNA gene mock | eight ATCC MSA-1003 libraries of PacBio's Kinnex 16S Revio data set (F01 to F08 with R37 to R44) | HiFiTaxa v1.2.3 (commit 9d6dcf3) | 5 to 7 October 2026 |
| GTDB r232 holdout | 10% of GTDB SSU r232 (95,731 sequences) against the remaining 90% | development checkout between v1.1.0 and v1.2.0 | 3 to 5 July 2026 |
| UNITE holdout | 10% of UNITE 19.02.2025 (10,213 sequences) against the remaining 90% | development checkout (with commit 49aceef) | 17 July 2026 |
| fungal ITS mock | the 103 species mock of Tedersoo et al. 2026 (ENA PRJEB108994) | development checkout (with commit 49aceef) | 17 July 2026 |

The classifier code paths are unchanged between the July checkout and v1.2.3 apart
from commit 49aceef, which passes the BLCA identity threshold to BLAST.

## Licence and citation

MIT (see `LICENSE`). Please cite the HiFiTaxa manuscript and this repository
(`CITATION.cff`).

# DADA2 Shiny

*[Version française](README.fr.md)*

Point-and-click Shiny app for DADA2 amplicon sequencing analysis, from raw Illumina reads to
ASV tables, taxonomy and Krona plots. Designed for several markers (16S, 18S, rbcL, COI, ITS)
and for eDNA metabarcoding.

> The user interface is currently in French.

## Features

- Step-by-step DADA2 workflow with every parameter editable, or a **Run all** button that
  chains the steps once settings are validated.
- Primer removal with **cutadapt** (primer mixes and IUPAC codes supported), primer presence
  check in every orientation, primers typed in or loaded from a FASTA/CSV file.
- **Overlap calculator** that warns before filtering when `truncLen` values make read merging
  impossible.
- **Multiple sequencing runs** (error models learned run by run), single- or paired-end,
  error model for binned quality scores (NovaSeq/NextSeq).
- Contaminant removal with **decontam** and abundance filters (tag-jumps, minimum reads).
- Taxonomy with `assignTaxonomy` or **IDTAXA**; six reference databases downloaded
  automatically.
- **Relative abundance** plots (stacked bars or heatmap, at ASV level or any taxonomic rank)
  and interactive **Krona** charts, built on phyloseq.
- **Reproducibility**: each step's parameters are recorded, the session is saved after every
  step and can be resumed, and a standalone R script reproducing the analysis can be exported.
- Long steps run in the background: the interface stays responsive and the log updates live.

## Installation

```r
install.packages(c("shiny", "DT", "ggplot2", "future", "promises", "BiocManager"))
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "decontam", "DECIPHER"))
```

`phyloseq`, `Biostrings`, `decontam` and `DECIPHER` are optional: the matching features are
disabled when they are missing.

cutadapt (recommended for primer removal): `pip install cutadapt` or
`conda install -c bioconda cutadapt`. On Windows, use conda or WSL and give the full path to
the executable in the *Amorces* (primers) tab.

## Running the app

```r
shiny::runApp("path/to/dada2shiny")   # or open app.R in RStudio > Run App
```

Keep `app.R` next to the `pipeline/` folder. To debug (everything in the main R process, no
background computation): `options(dada2shiny.async = FALSE)` before `runApp()`.

## Repository layout

```
dada2shiny/
├── app.R                  Shiny interface
└── pipeline/
    ├── pipeline.R         DADA2 steps, independent of Shiny (usable in scripts)
    ├── viz.R              relative abundance and Krona functions (usable in scripts)
    └── presets.R          marker presets and reference databases
```

## Workflow

1. **Data**: R1/R2 files paired *by sample name*; multiple runs (one subfolder per run);
   single-end; sample metadata (CSV); resume a saved session.
2. **Quality**: quality profiles.
3. **Primers**: presets (16S V4, 16S V3-V4, 18S V9, diatom rbcL, COI Leray, ITS2), primer
   check, removal with cutadapt or by fixed length.
4. **Filtering**: `filterAndTrim`, with expected-overlap check and warnings.
5. **Error model**: learned per run; optional model for binned quality scores.
6. **Denoising and merging**: per run, then `mergeSequenceTables`; warning when < 50 % of
   reads merge.
7. **Length filter and chimeras**.
8. **Contaminants and filters** (optional): decontam (prevalence/frequency), tag-jump
   threshold, minimum reads per ASV and per sample.
9. **Taxonomy**: `assignTaxonomy` (DADA2 FASTA) or IDTAXA (.RData training set).

Then: read tracking through the pipeline, exports (ASV table, FASTA, phyloseq object,
reproducible R script, session) and exploration.

### Run all

The **Tout lancer** (run all) button, available from every tab, chains the steps with the
current settings. It starts from the first step not yet run or whose parameters changed,
stops at the step you choose, stops on the first error (completed steps stay saved), and skips
taxonomy when no reference is given. On a new dataset, run the first steps one by one
(quality profiles, error model, merge rate) before relying on it.

### Primers and presets

- **Primer file**: FASTA, where the direction is read from the name (`>rbcL_708F_1`, `>R3_2`,
  `>515F-Y`, `>806RB`, `>COI_fwd`…), or CSV/TSV with `sequence` and `direction` (F/R) columns.
  Use the CSV format when names do not carry the direction (e.g. `fITS7`, `ITS4`).
- **Personal presets**: save all settings from the primer to taxonomy tabs under a name and
  reuse them across projects. They are stored in R's user configuration folder
  (`tools::R_user_dir("dada2shiny", "config")`), outside the app, so they survive updates.

### Exploration

- **Relative abundance**: stacked bars or heatmap, per ASV (optionally labelled with their
  genus) or aggregated at any rank; low-abundance taxa grouped as "Autres" (others),
  unassigned ASVs as "Non assigné". Samples can be averaged by a metadata variable, split into
  facets, and ordered by name, Bray-Curtis similarity (UPGMA) or dominant taxon. PDF, PNG and
  CSV export.
- **Krona**: interactive chart generated directly in R (KronaTools not required), one dataset
  per sample, per group or pooled, with or without ASVs as the last level. Downloadable as
  HTML. The Krona engine is loaded from marbl.github.io, so an internet connection is needed to
  display it.
- Alpha diversity, Bray-Curtis NMDS, and the phyloseq object as RDS.

## Reference databases

Downloaded and cached automatically: SILVA 138.2 (16S), RDP trainset 19 (16S),
Greengenes2 2024.09 (16S), PR2 5.0.0 (18S), SILVA 138.2 Fungi 18S and 26S.

To download yourself and use as a local file (links in the taxonomy tab): Diat.barcode (rbcL),
UNITE (ITS), MIDORI2 (COI and other mitochondrial genes). There is no official BOLD reference
in DADA2 format.

## Outputs

After each step the output folder contains `dada2_session.rds` (to resume the project) and
`dada2_log.txt` (the log). The exported R script uses the parameters *actually used* at each
step and runs without the app.

## Limitations

- Presets are starting points for 2x250/2x300 Illumina runs: adjust them to your quality
  profiles. Check primer sequences against the original publications, especially the rbcL
  primer mix (Vasselon et al. 2017).
- Rank names must match the reference file used (e.g. Diat.barcode).
- Not designed for long reads (Nanopore/PacBio).

## Citation

If you use this tool, please also cite the methods it relies on: DADA2 (Callahan et al. 2016,
*Nature Methods*), cutadapt (Martin 2011, *EMBnet.journal*), decontam (Davis et al. 2018,
*Microbiome*), phyloseq (McMurdie & Holmes 2013, *PLoS ONE*), Krona (Ondov et al. 2011,
*BMC Bioinformatics*), IDTAXA (Murali et al. 2018, *Microbiome*) if used, and the reference
database used.

## License

MIT, see [LICENSE](LICENSE).

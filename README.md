# DADA2 Shiny

*[Version française](README.fr.md)*

Point-and-click Shiny app for DADA2 amplicon sequencing analysis, from raw Illumina reads to
ASV tables, taxonomy and Krona plots. Designed for several markers (16S, 18S, rbcL, COI, ITS)
and for eDNA metabarcoding.

> The user interface is currently in French.

The repository contains two apps that share the same R modules (`pipeline/`):
- **`app.R`** — Illumina amplicons with DADA2 (described below);
- **`nanopore/app.R`** — Oxford Nanopore amplicons, from FASTQ to a sequence table, entirely in
  R so it also runs on Windows (see [Nanopore app](#nanopore-app)).

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
- **PCR replicates**: detection of failed PCRs (composition too far from the other replicates of
  the same sample, method inspired by metabaR) and merging of replicates (sum or mean, with an
  optional "present in at least k replicates" rule).
- **Taxonomic filter** (like mothur's `remove.lineage`): remove taxa (chloroplasts,
  mitochondria…), keep only a target group (e.g. Bacillariophyta), drop unassigned ASVs.
- Taxonomy with `assignTaxonomy` or **IDTAXA**; six reference databases downloaded
  automatically.
- **Remote BLAST** of selected ASVs on NCBI servers, with results (best hit, identity, coverage,
  ties between species) shown next to the DADA2 assignment.
- **Mock community evaluation** (in the spirit of mothur's `seq.error`): expected sequences
  recovered, error variants, unexpected ASVs, residual error rate.
- **Reference database builder**: turns your own FASTA (+ taxonomy in headers or in a table)
  into DADA2 files (`assignTaxonomy`, `addSpecies`) and an IDTAXA classifier, with cleaning and
  consistency checks.
- **Rarefaction curves** and diversity averaged over repeated rarefactions (alpha diversity and
  Bray-Curtis NMDS).
- **Relative abundance** plots (stacked bars or heatmap, at ASV level or any taxonomic rank)
  and interactive **Krona** charts, built on phyloseq.
- **Reproducibility**: each step's parameters are recorded, the session is saved after every
  step and can be resumed, and a standalone R script reproducing the analysis can be exported.
- Long steps run in the background: the interface stays responsive and the log updates live.

## Installation

```r
install.packages(c("shiny", "DT", "ggplot2", "future", "promises", "curl", "BiocManager"))
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "pwalign", "decontam", "DECIPHER", "ShortRead"))
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
├── app.R                  Illumina / DADA2 app
├── nanopore/app.R         Oxford Nanopore app
└── pipeline/
    ├── pipeline.R         DADA2 steps, independent of Shiny (usable in scripts)
    ├── viz.R              relative abundance and Krona functions (usable in scripts)
    ├── blast.R            remote NCBI BLAST (usable in scripts)
    ├── curation.R         replicates, taxonomic filter, rarefaction, mock, reference builder
    ├── nanopore.R         nanopore QC, primer orientation, clustering, consensus and polishing
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
   threshold, minimum reads per ASV and per sample, PCR replicates (failed PCR detection and
   merging).
9. **Taxonomy**: `assignTaxonomy` (DADA2 FASTA) or IDTAXA (.RData training set).
10. **Taxonomic filter** (optional): remove or keep taxa, drop ASVs unassigned at a given rank.

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

### PCR replicates

Add a metadata column giving the biological sample of each PCR (e.g. `sample_id`) and select it
in the contaminants step. For each sample, the composition of each PCR is compared with the
centroid of its replicates (Bray-Curtis); a PCR is flagged as failed when this dissimilarity
exceeds the usual variation between replicates. The threshold is set automatically where the
within-sample and between-sample dissimilarity distributions cross (as in metabaR), or by hand.
The most distant PCR is removed first and centroids are recomputed, so a single failed PCR does
not drag down its good replicates; with only two replicates, both are flagged since the faulty
one cannot be identified. Replicates are then merged. ASV identifiers stay the same throughout
(ASV number in the chimera-free table).

### Mock community

Choose the mock sample(s) and load the expected sequences (FASTA, one per species; they can be
longer than the amplicon and in either orientation). Each ASV is aligned to the closest expected
sequence and classified as exact, error variant (up to N differences) or unexpected.

### Reference database builder

*Outil : base de référence* tab. Input: a FASTA, with taxonomy in the headers (with or without a
leading ID, separators `;`, `|` or `,`) or in a separate table (one column per rank, or a single
taxonomy column). Cleaning: uppercase, U→T, gaps removed, invalid characters and short
sequences dropped, prefixes such as `k__` removed, empty or "unknown" ranks truncated. Exact
duplicates are removed; identical sequences with different taxonomies are reduced to their
common ancestor (or dropped). Names attached to several parents are reported. Outputs (in
`refs/` of the output folder): an `assignTaxonomy` file, an `addSpecies` file if a Species rank
exists, and optionally an IDTAXA classifier (`LearnTaxa`). One click sends them to the taxonomy
tab. Note: DADA2's `assignTaxonomy` has no separate training step (it learns from the FASTA at
each run); IDTAXA is trained once and tends to make fewer over-classifications.

### BLAST

The **BLAST** tab sends selected ASVs to NCBI BLAST (Common URL API) in a single search, up to
100 sequences. Settings: database (`core_nt` or `nt`), megablast or blastn, number of hits,
optional Entrez restriction (e.g. `Bacillariophyta[Organism]`), exclusion of uncultured /
environmental sequences, and your e-mail address (requested by NCBI). The search runs in the
background (from 30 s to a few minutes). Results show, for each ASV, the best hit with its
identity and coverage, the species tied with the same score, and the DADA2 assignment for
comparison; all hits link to their GenBank record. Results are saved in the session, added to
the ASV table export, and downloadable as CSV. NCBI usage rules are enforced (at least 10 s
between requests, at most one status check per minute).

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
- **Rarefaction**: rarefaction curves (exact expected richness) with the chosen depth, and
  diversity averaged over repeated rarefactions (Observed, Shannon, Simpson with standard
  deviations; mean Bray-Curtis for the NMDS).
- Alpha diversity, Bray-Curtis NMDS, and the phyloseq object as RDS.

## Nanopore app

```r
shiny::runApp("path/to/dada2shiny/nanopore")
```

Input: basecalled FASTQ files (MinKNOW / Dorado), one per barcode or one sub-folder per barcode
(`fastq_pass/barcode01/…`). Steps:

1. **Quality**: read length and mean quality per sample (N50, histograms).
2. **Filter and primers**: length / mean-quality filter; primers (IUPAC, primer mixes, presets)
   are searched at both ends with an error tolerance, reads are reoriented and trimmed.
3. **Consensus**, per sample: reads are clustered (`DECIPHER::Clusterize`), groups that hide two
   close variants are split (positions where two alleles co-vary across reads), each group is
   aligned (`AlignSeqs`) and turned into a polished consensus:
   - quality-weighted votes, balanced between strands (nanopore errors are often strand-specific);
   - homopolymer length correction (most likely run length estimated from all reads);
   - positions covered by less than half of the reads (primer or adapter leftovers) are ignored.
   Two modes: *metabarcoding* (all groups with at least N reads, to keep rare taxa) and
   *barcoding* (majority groups only). Optional, experimental: **Dorado polish** (`dorado aligner`
   + `dorado polish`) when dorado is installed (requires Rsamtools).
4. **Table**: near-identical consensus sequences are merged across samples; chimeras are removed
   with DADA2 (pooled method, tolerating one difference since a rare chimera's consensus rests on
   few reads).
5. **Taxonomy** and taxonomic filter, **mock community** evaluation, relative abundances, Krona,
   exports (table, FASTA, phyloseq, reproducible R script).

No Python, Linux or WSL is needed (unlike NGSpeciesID / Medaka). On simulated R10-like data
(≈2 % errors, mostly homopolymer indels; 5 species including one rare and two 2 % apart;
reads in both orientations, noise and chimeras), every true sequence was recovered with zero
differences across repeated random subsamplings, and all spurious consensus sequences were
chimeras removed at the table step. This is a check on simulated data: validate on your own mock
community before relying on it.

A simulated test dataset (with the expected sequences) is in `nanopore/exemple/`.

Dependencies: `DECIPHER`, `ShortRead`, `Biostrings`, `dada2` (Bioconductor); optional
`Rsamtools` (Dorado polish) and `phyloseq`.

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
- The DADA2 app is not designed for long reads: use the nanopore app for Oxford Nanopore data.

## Citation

If you use this tool, please also cite the methods it relies on: DADA2 (Callahan et al. 2016,
*Nature Methods*), metabaR for the replicate approach (Zinger et al. 2021, *Methods in Ecology
and Evolution*), cutadapt (Martin 2011, *EMBnet.journal*), decontam (Davis et al. 2018,
*Microbiome*), phyloseq (McMurdie & Holmes 2013, *PLoS ONE*), Krona (Ondov et al. 2011,
*BMC Bioinformatics*), IDTAXA (Murali et al. 2018, *Microbiome*) if used, and the reference
database used.

## License

MIT, see [LICENSE](LICENSE).

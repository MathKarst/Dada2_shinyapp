# =============================================================================
#  presets.R — préréglages par marqueur et bases de référence
#  Les valeurs sont des points de départ pour du 2x250 / 2x300 Illumina :
#  à ajuster selon les profils de qualité (le calcul de chevauchement de
#  l'onglet Filtrage signale les combinaisons impossibles).
# =============================================================================

MARKERS <- list(
  custom = list(label = "Personnalisé (pas de préréglage)"),

  "16S_V4" = list(
    label = "16S V4 — 515F-Y / 806RB (procaryotes)",
    primersF = "GTGYCAGCMGCCGCGGTAA",
    primersR = "GGACTACNVGGGTWTCTAAT",
    amplen = 253, truncLen = c(220, 160), maxEE = c(2, 2), minLen = 50,
    ref = "silva", levels = "Kingdom,Phylum,Class,Order,Family,Genus"),

  "16S_V3V4" = list(
    label = "16S V3-V4 — 341F / 805R (procaryotes, 2x300)",
    primersF = "CCTACGGGNGGCWGCAG",
    primersR = "GACTACHVGGGTATCTAATCC",
    amplen = 426, truncLen = c(270, 210), maxEE = c(2, 2), minLen = 50,
    ref = "silva", levels = "Kingdom,Phylum,Class,Order,Family,Genus"),

  "18S_V9" = list(
    label = "18S V9 — 1391F / EukBr (eucaryotes, longueur variable)",
    primersF = "GTACACACCGCCCGTC",
    primersR = "TGATCCTTCTGCAGGTTCACCTAC",
    amplen = NA, truncLen = c(0, 0), maxEE = c(2, 2), minLen = 50,
    ref = "pr2", levels = "Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species"),

  # Mix d'amorces de Vasselon et al. (2017), Ecol. Indic. — à vérifier contre la publication
  "rbcL_diatoms" = list(
    label = "rbcL diatomées 312 pb — Diat_rbcL_708F / R3 (Vasselon et al. 2017)",
    primersF = "AGGTGAAGTAAAAGGTTCWTACTTAAA, AGGTGAAGTTAAAGGTTCWTAYTTAAA, AGGTGAAACTAAAGGTTCWTACTTAAA",
    primersR = "CCTTCTAATTTACCWACWACTG, CCTTCTAATTTACCWACAACAG",
    amplen = 263, truncLen = c(200, 170), maxEE = c(2, 2), minLen = 50,
    # rangs à adapter au fichier Diat.barcode utilisé
    ref = "path", levels = "Kingdom,Phylum,Class,Order,Family,Genus,Species"),

  "COI_Leray" = list(
    label = "COI — mlCOIintF / jgHCO2198 (Leray et al. 2013), 313 pb",
    primersF = "GGWACWGGWTGAACWGTWTAYCCYCC",
    primersR = "TANACYTCNGGRTGNCCRAARAAYCA",
    amplen = 313, truncLen = c(220, 200), maxEE = c(2, 2), minLen = 50,
    ref = "path", levels = "Kingdom,Phylum,Class,Order,Family,Genus,Species"),

  "ITS2" = list(
    label = "ITS2 — fITS7 / ITS4 (champignons, longueur variable)",
    primersF = "GTGARTCATCGAATCTTTG",
    primersR = "TCCTCCGCTTATTGATATGC",
    amplen = NA, truncLen = c(0, 0), maxEE = c(2, 2), minLen = 50,
    ref = "path", levels = "Kingdom,Phylum,Class,Order,Family,Genus,Species")
)
MARKER_CHOICES <- setNames(names(MARKERS), vapply(MARKERS, `[[`, "", "label"))

REFS <- list(
  silva = list(
    label   = "SILVA 138.2 — 16S procaryotes (jusqu'au genre ; espèces via addSpecies)",
    url     = "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz",
    species = "https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz",
    levels  = "Kingdom,Phylum,Class,Order,Family,Genus"),
  rdp19 = list(
    label   = "RDP trainset 19 — 16S procaryotes (jusqu'au genre)",
    url     = "https://zenodo.org/records/14168771/files/rdp_19_toGenus_trainset.fa.gz",
    species = NULL,   # pas de fichier assignSpecies pour RDP 19
    levels  = "Kingdom,Phylum,Class,Order,Family,Genus"),
  gg2 = list(
    label   = "Greengenes2 2024.09 — 16S procaryotes (jusqu'au genre)",
    url     = "https://zenodo.org/records/14169078/files/gg2_2024_09_toGenus_trainset.fa.gz",
    species = NULL,
    levels  = "Kingdom,Phylum,Class,Order,Family,Genus"),
  pr2 = list(
    label   = "PR2 5.0.0 — 18S eucaryotes (jusqu'à l'espèce)",
    url     = "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_dada2.fasta.gz",
    species = NULL,
    levels  = "Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species"),
  silva_fungi18S = list(
    label   = "SILVA 138.2 Fungi — 18S champignons (jusqu'au genre)",
    url     = "https://zenodo.org/records/15044434/files/SILVA_SSUfungi_nr99_v138_2_toGenus_trainset.fasta",
    species = NULL,
    levels  = "Kingdom,Phylum,Class,Order,Family,Genus"),
  silva_fungi26S = list(
    label   = "SILVA 138.2 Fungi — 26S/LSU champignons (jusqu'à l'espèce)",
    url     = "https://zenodo.org/records/15044501/files/SILVA_LSUfungi_nr99_v138_2_toSpecies_trainset.fasta",
    species = NULL,
    levels  = "Kingdom,Phylum,Class,Order,Family,Genus,Species")
)

# Bases sans adresse de téléchargement stable : à télécharger soi-même, puis
# « Fichier local ». Affiché dans l'onglet Taxonomie.
REF_LINKS <- list(
  list(name = "Diat.barcode (rbcL 263 pb, diatomées), format DADA2",
       url = "https://doi.org/10.57745/SE6GJH"),
  list(name = "UNITE (ITS, champignons et eucaryotes), General FASTA release",
       url = "https://unite.ut.ee/repository.php"),
  list(name = "MIDORI2 (COI et autres gènes mitochondriaux), format DADA2",
       url = "https://www.reference-midori.info/"),
  list(name = "Liste complète des bases formatées pour DADA2",
       url = "https://benjjneb.github.io/dada2/training.html")
)
REF_CHOICES <- setNames(names(REFS), vapply(REFS, `[[`, "", "label"))

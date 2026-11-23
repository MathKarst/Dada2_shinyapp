# =============================================================================
#  DADA2 Shiny — interface pas-à-pas pour DADA2 (Illumina, single- ou paired-end)
#
#  Lancer depuis R :  shiny::runApp("chemin/vers/dada2shiny")
#
#  Dépendances
#    CRAN         : shiny, DT, ggplot2, future, promises
#    Bioconductor : dada2
#    Optionnelles : phyloseq + Biostrings (exploration), decontam (contaminants),
#                   DECIPHER (IDTAXA) ; cutadapt (outil externe, retrait des amorces)
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(ggplot2)
})

PIPELINE_FILE <- normalizePath(file.path("pipeline", "pipeline.R"))
source(PIPELINE_FILE)
source(file.path("pipeline", "presets.R"))
source(file.path("pipeline", "viz.R"))
source(file.path("pipeline", "blast.R"))
source(file.path("pipeline", "curation.R"))

HAS <- function(p) requireNamespace(p, quietly = TRUE)
has_phyloseq <- HAS("phyloseq")

# Calcul en arrière-plan : l'interface reste utilisable pendant les étapes longues.
# Mettre options(dada2shiny.async = FALSE) avant runApp() pour tout exécuter
# dans le processus principal (utile pour déboguer).
ASYNC <- isTRUE(getOption("dada2shiny.async", TRUE)) && HAS("future") && HAS("promises")
if (ASYNC) {
  future::plan(future::multisession, workers = 2)
  options(future.globals.maxSize = Inf)
}
options(shiny.maxRequestSize = 50 * 1024^3)
theme_set(theme_bw())
IS_WINDOWS <- .Platform$OS.type == "windows"

# -----------------------------------------------------------------------------
# Étapes, entrées associées et construction des paramètres
# -----------------------------------------------------------------------------
STEPS <- c("data", "primers", "filter", "errors", "denoise", "chimera", "post", "taxonomy", "taxfilter")
STEP_LABELS <- c(data = "Données", primers = "Amorces", filter = "Filtrage",
                 errors = "Erreurs", denoise = "Débruitage", chimera = "Chimères",
                 post = "Contaminants", taxonomy = "Taxonomie", taxfilter = "Filtre taxo")
STEP_INPUTS <- list(
  data     = c("src", "fq_dir", "paired", "patF", "patR", "sep", "sidx", "multirun",
               "outdir", "threads"),
  primers  = c("primer_method", "primersF", "primersR", "cut_minlen", "discard_untrimmed",
               "cutadapt_path"),
  filter   = c("truncF", "truncR", "trimLF", "trimLR", "maxEEF", "maxEER", "truncQ",
               "minLen", "rmphix"),
  errors   = c("nbases", "randomize", "errfun"),
  denoise  = c("pool", "minOverlap", "maxMismatch", "trimOverhang", "justConcat"),
  chimera  = c("lenMin", "lenMax", "chim_method", "minFold"),
  post     = c("decontam", "neg_col", "neg_values", "conc_col", "threshold", "drop_controls",
               "rel_thresh", "min_reads_asv", "min_reads_sample",
               "repl_col", "repl_detect", "repl_auto", "repl_thresh", "repl_aggregate", "repl_min"),
  taxonomy = c("tax_method", "ref_src", "ref_preset", "ref_path", "taxLevels", "minBoot",
               "tryRC", "do_species", "sp_src", "sp_path"),
  taxfilter = c("tf_remove", "tf_keep", "tf_unassigned")
)
# Listes déroulantes dont les choix dépendent des données (métadonnées, rangs)
DYNAMIC_SELECTS <- c("neg_col", "conc_col", "repl_col", "tf_unassigned")
# Champs de l'état produits par chaque étape (effacés si l'étape est relancée)
STEP_FIELDS <- list(
  data = c("samples", "readlen"), primers = "cut_stats", filter = "out", errors = "errs",
  denoise = c("seqtab", "counts"), chimera = c("seqtab_len", "nochim"),
  post = c("final", "final_pcr", "decontam", "repl", "blast"), taxonomy = c("taxa", "boot"),
  taxfilter = c("final_tax", "tf_report", "ps"))

BUILD <- list(
  data = function(v) list(outdir = normalizePath(v$outdir, mustWork = FALSE),
                          threads = as.integer(v$threads), paired = v$paired),
  primers = function(v) list(method = v$primer_method, primersF = v$primersF,
                             primersR = v$primersR, cut_minlen = v$cut_minlen,
                             discard_untrimmed = v$discard_untrimmed, cutadapt = v$cutadapt_path),
  filter = function(v) list(truncLen = c(v$truncF, v$truncR), trimLeft = c(v$trimLF, v$trimLR),
                            maxEE = c(v$maxEEF, v$maxEER), truncQ = v$truncQ,
                            minLen = v$minLen, rm.phix = v$rmphix),
  errors = function(v) list(nbases = v$nbases, randomize = v$randomize, errfun = v$errfun),
  denoise = function(v) list(pool = v$pool, minOverlap = v$minOverlap, maxMismatch = v$maxMismatch,
                             trimOverhang = v$trimOverhang, justConcatenate = v$justConcat),
  chimera = function(v) list(lenMin = v$lenMin, lenMax = v$lenMax, method = v$chim_method,
                             minFold = v$minFold),
  post = function(v) list(decontam = v$decontam, neg_col = v$neg_col, neg_values = v$neg_values,
                          conc_col = v$conc_col, threshold = v$threshold,
                          drop_controls = v$drop_controls, rel_thresh = v$rel_thresh,
                          min_reads_asv = v$min_reads_asv, min_reads_sample = v$min_reads_sample,
                          repl_col = v$repl_col %||% "", repl_detect = v$repl_detect,
                          repl_auto = v$repl_auto, repl_thresh = v$repl_thresh,
                          repl_aggregate = v$repl_aggregate, repl_min = v$repl_min),
  taxfilter = function(v) list(remove = v$tf_remove, keep = v$tf_keep,
                               unassigned_rank = v$tf_unassigned %||% "")
)

# Préréglages personnels : sauvegardés hors du dossier de l'appli, donc conservés
# quand on met l'appli à jour
PRESET_FILE <- file.path(tools::R_user_dir("dada2shiny", "config"), "user_presets.rds")
PRESET_INPUTS <- c("primersF", "primersR", "primer_method", "cut_minlen", "discard_untrimmed",
                   "cutadapt_path", "amplen", STEP_INPUTS$filter, STEP_INPUTS$errors,
                   STEP_INPUTS$denoise, STEP_INPUTS$chimera,
                   "decontam", "neg_values", "threshold", "drop_controls", "rel_thresh",
                   "min_reads_asv", "min_reads_sample", STEP_INPUTS$taxonomy,
                   STEP_INPUTS$taxfilter)
load_user_presets <- function() {
  if (!file.exists(PRESET_FILE)) return(list())
  tryCatch(readRDS(PRESET_FILE), error = function(e) list())
}
save_user_presets <- function(up) {
  dir.create(dirname(PRESET_FILE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(up, PRESET_FILE)
}

primer_len <- function(x) {
  p <- tryCatch(parse_primers(x), error = function(e) character(0))
  if (length(p)) max(nchar(p)) else 0
}

# =============================================================================
# UI
# =============================================================================
num <- function(id, label, value, min = NA, max = NA, step = NA)
  numericInput(id, label, value = value, min = min, max = max, step = step)
run_btn <- function(id, label) actionButton(id, label, class = "btn-primary", icon = icon("play"))
help <- function(...) helpText(...)

ui <- navbarPage(
  title = "DADA2 Shiny", id = "nav",
  header = tagList(
    tags$style(HTML("
      .run-all-bar {margin: -12px 0 4px 0;}
      .run-all-bar > * {margin-right: 6px;}
      .run-all-bar .form-group {margin-bottom: 0;}
      .run-all-bar .checkbox {margin: 0;}
      .step-badges {margin: 0 0 12px 0;}
      .step-badges .badge {font-size: 12px; margin-right: 4px; padding: 5px 8px;}
      .b-done {background:#2e7d32;} .b-stale {background:#ef6c00;}
      .b-run {background:#1565c0;} .b-todo {background:#9e9e9e;}
      #log {max-height: 78vh; overflow-y: auto; font-size: 12px; white-space: pre-wrap;}
      .ovl-ok {color:#2e7d32;} .ovl-warn {color:#ef6c00;} .ovl-bad {color:#c62828; font-weight:bold;}
    ")),
    div(class = "container-fluid",
      div(class = "run-all-bar",
        actionButton("run_all", "Tout lancer", class = "btn-success btn-sm", icon = icon("forward")),
        span("jusqu'à"),
        div(style = "display:inline-block; width:170px; vertical-align:middle;",
            selectInput("run_until", NULL, setNames(STEPS[-1], STEP_LABELS[-1]),
                        selected = "taxfilter", width = "100%")),
        div(style = "display:inline-block; vertical-align:middle;",
            checkboxInput("run_from_start", "Tout recalculer depuis les amorces", FALSE)),
        actionButton("stop_all", "Arrêter après l'étape en cours", class = "btn-default btn-sm",
                     icon = icon("stop"))),
      uiOutput("status_bar"))
  ),

  # ------------------------------------------------------------------ Données
  tabPanel("1. Données",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Reprendre un projet"),
        fluidRow(
          column(7, fileInput("load_rds", NULL, accept = ".rds",
                              buttonLabel = "Session .rds…", placeholder = "dada2_session.rds")),
          column(5, actionButton("load_outdir", "Depuis le dossier de sortie", icon = icon("folder-open")))),
        hr(),
        h4("Nouveau projet"),
        radioButtons("src", "Source des FASTQ",
                     c("Dossier local" = "dir", "Upload" = "upload"), inline = TRUE),
        conditionalPanel("input.src == 'dir'",
          textInput("fq_dir", "Dossier contenant les FASTQ", value = getwd())),
        conditionalPanel("input.src == 'upload'",
          fileInput("fq_up", "Fichiers FASTQ", multiple = TRUE, accept = c(".fastq", ".fq", ".gz")),
          help("Préfère le dossier local pour de gros jeux de données.")),
        checkboxInput("multirun", "Plusieurs runs de séquençage (un sous-dossier par run)", FALSE),
        checkboxInput("paired", "Reads appariés (R1 + R2)", TRUE),
        fluidRow(
          column(6, textInput("patF", "Motif forward (regex)", "_R1_001")),
          column(6, conditionalPanel("input.paired", textInput("patR", "Motif reverse (regex)", "_R2_001")))),
        fluidRow(
          column(6, textInput("sep", "Séparateur du nom", "_")),
          column(6, num("sidx", "Position du nom d'échantillon", 1, min = 1))),
        textInput("outdir", "Dossier de sortie", file.path(getwd(), "dada2_out")),
        num("threads", "Threads", max(1, parallel::detectCores() - 1), min = 1),
        run_btn("scan", "Détecter les échantillons"),
        hr(),
        h4("Métadonnées (optionnel)"),
        fileInput("meta", "CSV : 1re colonne = noms d'échantillons", accept = c(".csv", ".tsv", ".txt")),
        help("Utilisées pour decontam (témoins négatifs, concentrations) et l'exploration.")
      ),
      mainPanel(
        uiOutput("data_summary"),
        DTOutput("files_tbl"),
        uiOutput("meta_summary")
      )
    )
  ),

  # ------------------------------------------------------------------ Qualité
  tabPanel("2. Qualité",
    fluidRow(
      column(5, selectizeInput("qsamples", "Échantillons (max 6)", choices = NULL, multiple = TRUE,
                               options = list(maxItems = 6))),
      column(2, num("qn", "Reads lus / fichier", 50000, min = 1000, step = 10000)),
      column(2, br(), checkboxInput("qagg", "Profil agrégé", FALSE)),
      column(3, br(), actionButton("qplot", "Tracer", class = "btn-primary"))),
    fluidRow(
      column(6, h4("Forward (R1)"), plotOutput("qF", height = 520)),
      column(6, h4("Reverse (R2)"), plotOutput("qR", height = 520)))
  ),

  # ------------------------------------------------------------------ Amorces
  tabPanel("3. Amorces",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Préréglage marqueur"),
        selectInput("marker", NULL, MARKER_CHOICES),
        actionButton("apply_marker", "Appliquer le préréglage", icon = icon("wand-magic-sparkles")),
        conditionalPanel("input.marker && input.marker.startsWith('user:')",
          actionButton("delete_preset", "Supprimer ce préréglage", icon = icon("trash"),
                       class = "btn-sm btn-link text-danger")),
        help("Préréglage intégré : amorces, longueur d'amplicon, troncature et référence. ",
             "Préréglage personnel : tous les réglages enregistrés (onglets Amorces à Taxonomie)."),
        fluidRow(
          column(7, textInput("preset_name", NULL, placeholder = "Nom, ex. rbcL rivières 2x250")),
          column(5, actionButton("save_preset", "Enregistrer mes réglages", icon = icon("floppy-disk"),
                                 class = "btn-sm"))),
        hr(),
        textAreaInput("primersF", "Amorce(s) forward (IUPAC, séparées par des virgules)", "", rows = 2),
        textAreaInput("primersR", "Amorce(s) reverse", "", rows = 2),
        fileInput("primer_file", "… ou charger un fichier d'amorces",
                  accept = c(".fasta", ".fa", ".fas", ".csv", ".tsv", ".txt"),
                  buttonLabel = "FASTA / CSV…"),
        help("FASTA : le sens est lu dans le nom (>rbcL_708F_1, >R3_2, >COI_fwd…). ",
             "CSV/TSV : colonnes « sequence » et « direction » (F/R)."),
        tableOutput("primer_file_tbl"),
        actionButton("check_primers", "Vérifier la présence des amorces", icon = icon("magnifying-glass")),
        hr(),
        radioButtons("primer_method", "Retrait des amorces",
          c("cutadapt (recommandé)" = "cutadapt",
            "Longueur fixe (trimLeft à l'étape Filtrage)" = "trimLeft",
            "Aucun (amorces déjà retirées)" = "none")),
        conditionalPanel("input.primer_method == 'cutadapt'",
          textInput("cutadapt_path", "Exécutable cutadapt", "cutadapt"),
          fluidRow(column(6, num("cut_minlen", "Longueur min. après coupe", 20, min = 1)),
                   column(6, br(), checkboxInput("discard_untrimmed", "Jeter les reads sans amorce", TRUE)))),
        conditionalPanel("input.primer_method == 'trimLeft'",
          help("trimLeft sera rempli avec la longueur des amorces. Ne convient pas aux amorces ",
               "précédées d'espaceurs, ni si les reads dépassent l'amplicon.")),
        run_btn("run_primers", "Valider / lancer cutadapt")
      ),
      mainPanel(
        h4("Présence des amorces (premiers reads des 3 premiers échantillons)"),
        help("Attendu avant retrait : F (sens) ≈ 100 % en R1, R (sens) ≈ 100 % en R2. ",
             "Les rev. comp. apparaissent quand les reads traversent l'amplicon. ",
             "Après cutadapt, tout doit être proche de 0 %."),
        DTOutput("primer_tbl"),
        uiOutput("cut_summary"),
        DTOutput("cut_tbl")
      )
    )
  ),

  # ------------------------------------------------------------------ Filtrage
  tabPanel("4. Filtrage",
    sidebarLayout(
      sidebarPanel(width = 4,
        num("amplen", "Longueur attendue de l'amplicon, sans amorces (pb, vide = variable)", NA, min = 1),
        h5("Troncature (0 = aucune)"),
        fluidRow(column(6, num("truncF", "truncLen F", 240, min = 0)),
                 column(6, conditionalPanel("input.paired", num("truncR", "truncLen R", 160, min = 0)))),
        h5("Coupe en 5'"),
        fluidRow(column(6, num("trimLF", "trimLeft F", 0, min = 0)),
                 column(6, conditionalPanel("input.paired", num("trimLR", "trimLeft R", 0, min = 0)))),
        h5("Erreurs attendues maximales"),
        fluidRow(column(6, num("maxEEF", "maxEE F", 2, min = 0, step = 0.5)),
                 column(6, conditionalPanel("input.paired", num("maxEER", "maxEE R", 2, min = 0, step = 0.5)))),
        fluidRow(column(6, num("truncQ", "truncQ", 2, min = 0)),
                 column(6, num("minLen", "minLen", 50, min = 1))),
        checkboxInput("rmphix", "Retirer PhiX", TRUE),
        help("maxN est fixé à 0 (DADA2 n'accepte pas les N)."),
        run_btn("run_filter", "Lancer le filtrage")
      ),
      mainPanel(
        h4("Contrôle de la longueur et du chevauchement"),
        uiOutput("overlap_ui"),
        hr(),
        uiOutput("filter_summary"),
        DTOutput("filter_tbl")
      )
    )
  ),

  # ------------------------------------------------------------------ Erreurs
  tabPanel("5. Modèle d'erreur",
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("errfun", "Fonction d'estimation",
          c("Standard (loessErrfun) — MiSeq, qualités complètes" = "loess",
            "Qualités binnées — NovaSeq / NextSeq / MiSeq récents" = "binned")),
        help("Si les points noirs des graphiques s'écartent fortement de la ligne rouge ",
             "avec seulement 4-6 valeurs de qualité, choisis « qualités binnées »."),
        num("nbases", "nbases", 1e8, min = 1e6, step = 1e7),
        checkboxInput("randomize", "Échantillons tirés au hasard", FALSE),
        run_btn("run_errors", "Apprendre les erreurs"),
        hr(),
        selectInput("err_run", "Run affiché", choices = NULL)
      ),
      mainPanel(fluidRow(
        column(6, h4("Forward"), plotOutput("errF_plot", height = 560)),
        column(6, h4("Reverse"), plotOutput("errR_plot", height = 560))))
    )
  ),

  # ------------------------------------------------------------------ Débruitage
  tabPanel("6. Débruitage & fusion",
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("pool", "Pooling (à l'intérieur de chaque run)",
                    c("Pseudo-pooling (recommandé eDNA)" = "pseudo",
                      "Par échantillon" = "FALSE", "Pooling complet (lent)" = "TRUE")),
        conditionalPanel("input.paired",
          hr(), h5("Fusion des paires"),
          num("minOverlap", "minOverlap", 12, min = 1),
          num("maxMismatch", "maxMismatch", 0, min = 0),
          checkboxInput("trimOverhang", "trimOverhang", TRUE),
          checkboxInput("justConcat", "justConcatenate (pas de chevauchement)", FALSE)),
        run_btn("run_denoise", "Débruiter + fusionner")
      ),
      mainPanel(
        uiOutput("denoise_summary"),
        h4("Distribution des longueurs des ASV"),
        plotOutput("len_plot", height = 330),
        DTOutput("denoise_tbl")
      )
    )
  ),

  # ------------------------------------------------------------------ Chimères
  tabPanel("7. Longueur & chimères",
    sidebarLayout(
      sidebarPanel(width = 3,
        h5("Gamme de longueur des ASV conservées"),
        fluidRow(column(6, num("lenMin", "Min (pb)", 1, min = 1)),
                 column(6, num("lenMax", "Max (pb)", 1000, min = 1))),
        uiOutput("len_hint"),
        hr(),
        selectInput("chim_method", "Méthode", c("consensus", "pooled", "per-sample")),
        num("minFold", "minFoldParentOverAbundance", 1, min = 1, step = 0.5),
        run_btn("run_chimera", "Filtrer + retirer les chimères")
      ),
      mainPanel(uiOutput("chim_summary"))
    )
  ),

  # ------------------------------------------------------------------ Contaminants
  tabPanel("8. Contaminants & filtres",
    sidebarLayout(
      sidebarPanel(width = 4,
        h5("decontam (nécessite des métadonnées)"),
        selectInput("decontam", "Méthode",
          c("Aucune" = "none", "Prévalence (témoins négatifs)" = "prevalence",
            "Fréquence (concentration ADN)" = "frequency", "Les deux (either)" = "either",
            "Combinée" = "combined")),
        conditionalPanel("['prevalence','either','combined'].includes(input.decontam)",
          selectInput("neg_col", "Colonne indiquant les témoins négatifs", choices = NULL),
          textInput("neg_values", "Valeurs = témoin négatif (virgules)", "TRUE, neg, negative, blank, control"),
          checkboxInput("drop_controls", "Retirer les témoins de la table finale", TRUE)),
        conditionalPanel("['frequency','either','combined'].includes(input.decontam)",
          selectInput("conc_col", "Colonne de concentration ADN", choices = NULL)),
        conditionalPanel("input.decontam != 'none'",
          num("threshold", "Seuil (threshold)", 0.1, min = 0, max = 1, step = 0.05)),
        if (!HAS("decontam")) div(class = "text-danger", "Package decontam non installé."),
        hr(),
        h5("Filtres d'abondance"),
        num("rel_thresh", "Mettre à 0 les occurrences < X % du total de l'ASV (tag-jumps ; 0 = off)",
            0, min = 0, step = 0.05),
        num("min_reads_asv", "Nombre min. de reads par ASV", 2, min = 1),
        num("min_reads_sample", "Nombre min. de reads par échantillon (0 = off)", 0, min = 0),
        hr(),
        h5("Réplicats PCR (optionnel)"),
        selectInput("repl_col", "Colonne indiquant l'échantillon biologique de chaque PCR",
                    c("(aucune : pas de réplicats)" = "")),
        conditionalPanel("input.repl_col != ''",
          checkboxInput("repl_detect", "Repérer et retirer les PCR aberrantes", TRUE),
          conditionalPanel("input.repl_detect",
            checkboxInput("repl_auto", "Seuil de dissimilarité automatique", TRUE),
            conditionalPanel("!input.repl_auto",
              num("repl_thresh", "Seuil de dissimilarité (Bray-Curtis au centroïde)", 0.5,
                  min = 0, max = 1, step = 0.05))),
          selectInput("repl_aggregate", "Fusionner les réplicats",
            c("Somme des reads" = "sum", "Moyenne (arrondie)" = "mean",
              "Ne pas fusionner (garder les PCR)" = "none")),
          conditionalPanel("input.repl_aggregate != 'none'",
            num("repl_min", "Garder une ASV dans un échantillon si présente dans au moins … réplicat(s)",
                1, min = 1, step = 1))),
        help("Méthode inspirée de metabaR (Zinger et al. 2021) : une PCR est aberrante si sa composition ",
             "s'écarte de celle des autres réplicats du même échantillon plus que la variation habituelle ",
             "entre réplicats."),
        run_btn("run_post", "Appliquer")
      ),
      mainPanel(
        uiOutput("post_summary"),
        tabsetPanel(
          tabPanel("decontam", br(), DTOutput("decontam_tbl")),
          tabPanel("Réplicats PCR", br(),
                   plotOutput("repl_plot", height = 330),
                   DTOutput("repl_tbl")))
      )
    )
  ),

  # ------------------------------------------------------------------ Taxonomie
  tabPanel("9. Taxonomie",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons("tax_method", "Classificateur",
          c("assignTaxonomy (RDP naive Bayes, FASTA DADA2)" = "assignTaxonomy",
            "IDTAXA (DECIPHER, trainingSet .RData)" = "idtaxa")),
        if (!HAS("DECIPHER")) help("DECIPHER non installé : IDTAXA indisponible."),
        radioButtons("ref_src", "Base de référence",
          c("Base prédéfinie (téléchargée)" = "preset", "Fichier local" = "path", "Upload" = "upload")),
        conditionalPanel("input.ref_src == 'preset'", selectInput("ref_preset", NULL, REF_CHOICES)),
        conditionalPanel("input.ref_src == 'path'", textInput("ref_path", "Chemin du fichier de référence", "")),
        conditionalPanel("input.ref_src == 'upload'", fileInput("ref_up", "Fichier de référence")),
        tags$details(
          tags$summary(class = "text-muted", "Autres bases (à télécharger puis « Fichier local »)"),
          tags$ul(lapply(REF_LINKS, function(r)
            tags$li(tags$a(href = r$url, target = "_blank", r$name)))),
          help("COI : pas de base BOLD officielle au format DADA2 (MIDORI2 ou base maison). ",
               "Adapte les rangs au fichier utilisé.")),
        textInput("taxLevels", "Rangs (virgules)", REFS$silva$levels),
        fluidRow(column(6, num("minBoot", "minBoot / seuil de confiance", 50, min = 0, max = 100)),
                 column(6, br(), checkboxInput("tryRC", "Tester le rev. comp.", FALSE))),
        hr(),
        checkboxInput("do_species", "addSpecies (correspondance exacte)", FALSE),
        conditionalPanel("input.do_species",
          radioButtons("sp_src", NULL, c("Fichier espèces de la base prédéfinie" = "preset",
                                         "Chemin local" = "path", "Upload" = "upload")),
          conditionalPanel("input.sp_src == 'path'", textInput("sp_path", "Chemin", "")),
          conditionalPanel("input.sp_src == 'upload'", fileInput("sp_up", "Fichier espèces"))),
        run_btn("run_taxonomy", "Assigner la taxonomie")
      ),
      mainPanel(uiOutput("tax_summary"), DTOutput("tax_tbl"))
    )
  ),

  # ------------------------------------------------------------------ Filtre taxonomique
  tabPanel("10. Filtre taxonomique",
    sidebarLayout(
      sidebarPanel(width = 4,
        textInput("tf_remove", "Retirer les ASV assignées à (virgules)", "",
                  placeholder = "ex. Chloroplast, Mitochondria"),
        textInput("tf_keep", "Ne garder que les ASV assignées à (virgules, vide = tout)", "",
                  placeholder = "ex. Bacillariophyta"),
        selectInput("tf_unassigned", "Retirer les ASV non assignées au rang",
                    c("(ne pas filtrer)" = "")),
        help("Les noms sont comparés exactement (sans tenir compte des majuscules) à tous les rangs. ",
             "Exemples : 16S → retirer « Chloroplast, Mitochondria, Eukaryota » ; ",
             "rbcL → garder « Bacillariophyta »."),
        run_btn("run_taxfilter", "Appliquer le filtre")
      ),
      mainPanel(uiOutput("tf_summary"), tableOutput("tf_report"))
    )
  ),

  # ------------------------------------------------------------------ BLAST
  tabPanel("BLAST",
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("BLAST NCBI des ASV sélectionnées"),
        uiOutput("blast_sel_info"),
        fluidRow(
          column(6, actionButton("blast_select_visible", "Sélectionner les lignes filtrées", class = "btn-sm")),
          column(6, actionButton("blast_clear", "Vider la sélection", class = "btn-sm"))),
        hr(),
        selectInput("blast_db", "Base NCBI",
          c("core_nt — nucléotides (base par défaut du NCBI)" = "core_nt",
            "nt — nucléotides complète (plus lente)" = "nt")),
        radioButtons("blast_algo", "Algorithme",
          c("megablast — séquences très proches, rapide" = "megablast",
            "blastn — plus sensible, plus lent" = "blastn")),
        num("blast_hits", "Nombre de hits par ASV", 10, min = 1, max = 100),
        textInput("blast_entrez", "Restreindre à (requête Entrez, optionnel)", "",
                  placeholder = "ex. Bacillariophyta[Organism]"),
        checkboxInput("blast_uncult", "Exclure les séquences « uncultured » / environnementales", TRUE),
        textInput("blast_email", "Ton adresse e-mail (demandée par le NCBI)", "",
                  placeholder = "prenom.nom@labo.fr"),
        actionButton("run_blast", "BLASTer la sélection", class = "btn-primary", icon = icon("magnifying-glass")),
        help("Les séquences partent toutes en une seule recherche (", BLAST_MAX_SEQS, " au maximum). ",
             "Comptez de 30 s à quelques minutes selon la charge du NCBI. ",
             "Pour de gros volumes, le NCBI demande de lancer les recherches le soir ou le week-end.")
      ),
      mainPanel(width = 9,
        uiOutput("blast_status"),
        tabsetPanel(id = "blast_tabs",
          tabPanel("ASV à sélectionner", br(),
                   help("Clique sur les lignes pour les sélectionner ; les filtres sous les en-têtes ",
                        "permettent par ex. de ne garder que les ASV non assignées au genre."),
                   DTOutput("blast_asv_tbl")),
          tabPanel("Synthèse par ASV", br(),
                   downloadButton("dl_blast_summary", "Synthèse (CSV)", class = "btn-sm"),
                   downloadButton("dl_blast_hits", "Tous les hits (CSV)", class = "btn-sm"),
                   br(), br(), DTOutput("blast_summary_tbl")),
          tabPanel("Tous les hits", br(), DTOutput("blast_hits_tbl"))))
    )
  ),

  # ------------------------------------------------------------------ Suivi
  tabPanel("Suivi des reads",
    plotOutput("track_plot", height = 420),
    DTOutput("track_tbl")
  ),

  # ------------------------------------------------------------------ Mock
  tabPanel("Communauté mock",
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Évaluer une communauté mock"),
        radioButtons("mock_table", "Table utilisée",
          c("Sans chimères (avant filtres)" = "nochim", "Table finale" = "final")),
        selectizeInput("mock_samples", "Échantillon(s) mock", choices = NULL, multiple = TRUE),
        fileInput("mock_fasta", "Séquences attendues (FASTA)", accept = c(".fa", ".fasta", ".fas", ".gz")),
        help("Une séquence par espèce de la communauté ; elles peuvent être plus longues que ",
             "l'amplicon (gène complet) et dans n'importe quel sens."),
        num("mock_maxdiff", "Différences max. pour une variante d'erreur", 3, min = 0, max = 20),
        actionButton("run_mock", "Évaluer", class = "btn-primary", icon = icon("vial"))
      ),
      mainPanel(width = 9,
        uiOutput("mock_status"),
        fluidRow(column(5, h4("Bilan"), tableOutput("mock_summary")),
                 column(7, h4("Par séquence attendue"), tableOutput("mock_perref"))),
        plotOutput("mock_plot", height = 380),
        h4("Par ASV"),
        downloadButton("dl_mock", "Tableau par ASV (CSV)", class = "btn-sm"), br(), br(),
        DTOutput("mock_asv"))
    )
  ),

  # ------------------------------------------------------------------ Outil : base de référence
  tabPanel("Outil : base de référence",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Préparer une base de référence perso"),
        fileInput("rf_fasta", "Séquences (FASTA, .gz accepté)", accept = c(".fa", ".fasta", ".fas", ".fna", ".gz", ".txt")),
        radioButtons("rf_source", "Où est la taxonomie ?",
          c("Dans les en-têtes du FASTA" = "header", "Dans un tableau à part (CSV/TSV)" = "table")),
        conditionalPanel("input.rf_source == 'header'",
          checkboxInput("rf_header_id", "L'en-tête commence par un identifiant (>ID taxo…)", TRUE)),
        conditionalPanel("input.rf_source == 'table'",
          fileInput("rf_table", "Tableau de taxonomie", accept = c(".csv", ".tsv", ".txt")),
          fluidRow(
            column(6, textInput("rf_table_id", "Colonne identifiant", "id")),
            column(6, textInput("rf_table_tax", "Colonne taxonomie (vide = une colonne par rang)", "")))),
        selectInput("rf_sep", "Séparateur des rangs", c("; (point-virgule)" = ";", "| (barre)" = "|", ", (virgule)" = ",")),
        textInput("rf_levels", "Rangs, du plus large au plus fin", "Kingdom,Phylum,Class,Order,Family,Genus,Species"),
        checkboxInput("rf_strip", "Retirer les préfixes k__, p__, d__… (UNITE, Greengenes)", TRUE),
        fluidRow(
          column(6, num("rf_minlen", "Longueur min. (pb)", 50, min = 1)),
          column(6, selectInput("rf_conflicts", "Même séquence, taxonomies différentes",
                                c("Ancêtre commun" = "lca", "Retirer" = "drop")))),
        checkboxInput("rf_species_train", "Garder le rang Species dans le fichier assignTaxonomy", FALSE),
        checkboxInput("rf_idtaxa", "Entraîner aussi un classificateur IDTAXA (DECIPHER)", FALSE),
        if (!HAS("DECIPHER")) help("DECIPHER non installé : option IDTAXA indisponible."),
        textInput("rf_name", "Nom de la base", "ma_base"),
        actionButton("run_ref", "Préparer la base", class = "btn-primary", icon = icon("database"))
      ),
      mainPanel(width = 8,
        uiOutput("ref_status"),
        uiOutput("ref_files"),
        fluidRow(column(5, h4("Nettoyage"), tableOutput("ref_report")),
                 column(7, h4("Rangs"), tableOutput("ref_ranks"))),
        h4("Noms rattachés à plusieurs parents"),
        help("À vérifier : ce sont souvent des homonymes ou des erreurs de la base d'origine."),
        DTOutput("ref_incons"))
    )
  ),

  # ------------------------------------------------------------------ Export
  tabPanel("Export",
    fluidRow(column(8,
      h4("Résultats"),
      p(downloadButton("dl_asv", "Table ASV + taxonomie (CSV)")),
      p(downloadButton("dl_fasta", "Séquences des ASV (FASTA)")),
      p(downloadButton("dl_track", "Suivi des reads (CSV)")),
      p(downloadButton("dl_boot", "Bootstraps / confiances taxonomiques (CSV)")),
      p(downloadButton("dl_decontam", "Résultats decontam (CSV)")),
      p(downloadButton("dl_ps", "Objet phyloseq (RDS)")),
      h4("Reproductibilité"),
      p(downloadButton("dl_script", "Script R reproductible (.R)")),
      p(downloadButton("dl_session", "Session complète (RDS, rechargeable dans l'appli)")),
      help("La session est aussi enregistrée automatiquement dans le dossier de sortie ",
           "(dada2_session.rds) après chaque étape, avec le journal (dada2_log.txt).")
    ))
  ),

  # ------------------------------------------------------------------ Exploration
  tabPanel("Exploration",
    if (!has_phyloseq) {
      div(class = "alert alert-warning", "Installe phyloseq pour activer cet onglet.")
    } else sidebarLayout(
      sidebarPanel(width = 2,
        textInput("exclude", "Échantillons à exclure (virgules, ex. Mock)", ""),
        actionButton("build_ps", "Mettre à jour", class = "btn-primary", icon = icon("rotate")),
        help("L'objet phyloseq est construit automatiquement à l'ouverture de l'onglet ",
             "(table finale + taxonomie + métadonnées)."),
        hr(),
        h5("Diversité alpha et NMDS"),
        selectInput("xvar", "Variable en X", choices = NULL),
        selectInput("colvar", "Couleur", choices = NULL)
      ),
      mainPanel(width = 10,
        tabsetPanel(id = "explo_tabs",
          tabPanel("Abondances relatives",
            br(),
            fluidRow(
              column(3, selectInput("ab_level", "Niveau", c("ASV (sans taxonomie)" = "ASV"))),
              column(3, conditionalPanel("input.ab_level == 'ASV'",
                selectInput("ab_label", "Annoter les ASV avec le rang", c("(rien)" = "")))),
              column(2, num("ab_top", "Nombre de taxons", 15, min = 1, max = 60)),
              column(2, radioButtons("ab_type", "Représentation",
                                     c("Barres empilées" = "bar", "Heatmap" = "heatmap"))),
              column(2, br(), downloadButton("dl_ab_pdf", "PDF", class = "btn-sm"),
                     downloadButton("dl_ab_png", "PNG", class = "btn-sm"),
                     downloadButton("dl_ab_csv", "CSV", class = "btn-sm"))),
            fluidRow(
              column(3, selectInput("ab_group", "Moyenner les échantillons par",
                                    c("(non : une barre par échantillon)" = ""))),
              column(3, selectInput("ab_facet", "Séparer en facettes selon", c("(aucune)" = ""))),
              column(3, selectInput("ab_order", "Ordre des échantillons",
                                    c("Nom" = "name", "Similarité (Bray-Curtis, UPGMA)" = "cluster",
                                      "Abondance du taxon dominant" = "dominant")))),
            help("Les taxons les moins abondants (moyenne sur les échantillons) sont regroupés dans « Autres » ; ",
                 "les ASV non assignées au rang choisi sont dans « Non assigné »."),
            plotOutput("ab_plot", height = "auto")),
          tabPanel("Krona",
            br(),
            fluidRow(
              column(3, radioButtons("kr_mode", "Jeux de données",
                c("Un par échantillon" = "sample", "Un par groupe (moyenne des %)" = "group",
                  "Tous les échantillons cumulés" = "all"))),
              column(3, conditionalPanel("input.kr_mode == 'group'",
                selectInput("kr_group", "Variable de regroupement", choices = NULL))),
              column(3, checkboxInput("kr_asv", "Afficher les ASV au dernier niveau", TRUE)),
              column(3, downloadButton("dl_krona", "Télécharger le Krona (HTML)"))),
            help("Chaque jeu de données se choisit dans la liste en haut à gauche du graphique. ",
                 "Le moteur Krona est chargé depuis marbl.github.io : une connexion internet est nécessaire ",
                 "pour l'afficher, y compris pour le fichier téléchargé."),
            uiOutput("krona_frame")),
          tabPanel("Raréfaction", br(),
            fluidRow(
              column(3, num("rare_depth", "Profondeur de raréfaction (reads)", NA, min = 1)),
              column(3, num("rare_iter", "Nombre de raréfactions répétées", 20, min = 1, max = 500)),
              column(3, br(), actionButton("run_rare", "Calculer la diversité raréfiée", class = "btn-primary"))),
            help("Courbes : richesse attendue selon le nombre de reads (formule exacte). Une courbe qui ",
                 "plafonne signifie que la profondeur suffit. La ligne verticale marque la profondeur choisie ; ",
                 "les échantillons en dessous sont exclus des analyses raréfiées."),
            plotOutput("rare_plot", height = 420),
            uiOutput("rare_info")),
          tabPanel("Diversité alpha", br(),
            radioButtons("alpha_method", NULL, inline = TRUE,
              c("Sans raréfaction (phyloseq)" = "raw", "Moyenne sur raréfactions répétées" = "rare")),
            plotOutput("rich_plot", height = 450),
            conditionalPanel("input.alpha_method == 'rare'",
              downloadButton("dl_alpha_rare", "Diversité raréfiée (CSV)", class = "btn-sm"))),
          tabPanel("NMDS Bray-Curtis", br(),
            radioButtons("nmds_method", NULL, inline = TRUE,
              c("Proportions" = "raw", "Bray-Curtis moyen sur raréfactions répétées" = "rare")),
            plotOutput("nmds_plot", height = 450)),
          tabPanel("Objet phyloseq", br(), verbatimTextOutput("ps_print"),
                   downloadButton("dl_ps2", "Télécharger (RDS)"))))
    )
  ),

  # ------------------------------------------------------------------ Journal
  tabPanel("Journal", verbatimTextOutput("log"))
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  state <- reactiveVal(list(params = list(), inputs = list()))
  busy <- reactiveVal(NULL)
  job_start <- reactiveVal(NULL)

  notify <- function(msg, type = "message", duration = 6)
    showNotification(msg, type = type, duration = duration)

  # Valeurs actuelles des entrées d'une étape
  snap <- function(step) {
    v <- lapply(setNames(STEP_INPUTS[[step]], STEP_INPUTS[[step]]), function(i) isolate(input[[i]]))
    v[!vapply(v, is.null, TRUE)]
  }
  done <- function(step, st = state()) !is.null(st$params[[step]])
  stale <- function(step, st = state()) {
    if (!done(step, st)) return(FALSE)
    a <- st$inputs[[step]]; b <- lapply(names(a), function(i) input[[i]])
    !identical(lapply(unname(a), as.character), lapply(b, as.character))
  }
  outdir <- function() state()$params$data$outdir
  threads <- function() state()$params$data$threads %||% 1L
  paired <- function() is_paired(state()$samples)

  # Efface l'étape donnée et les suivantes
  reset_from <- function(st, step) {
    for (s in STEPS[match(step, STEPS):length(STEPS)]) {
      for (f in STEP_FIELDS[[s]]) st[f] <- list(NULL)
      st$params[s] <- list(NULL); st$inputs[s] <- list(NULL)
    }
    if (match(step, STEPS) <= match("primers", STEPS) && !is.null(st$samples))
      st$samples$trimF <- st$samples$trimR <- NULL
    if (match(step, STEPS) <= match("filter", STEPS) && !is.null(st$samples))
      st$samples$filtF <- st$samples$filtR <- NULL
    st
  }
  save_session <- function(st = state()) {
    d <- st$params$data$outdir
    if (!is.null(d)) tryCatch(saveRDS(st, file.path(d, "dada2_session.rds")),
                              error = function(e) notify(paste("Session non sauvegardée :", conditionMessage(e)), "warning"))
  }
  log_file <- reactive({ d <- outdir(); if (is.null(d)) NULL else file.path(d, "dada2_log.txt") })
  log_line <- function(...) {
    lf <- isolate(log_file()); if (is.null(lf)) return()
    cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", file = lf, append = TRUE, sep = "")
  }

  # Vérifie qu'une étape peut être lancée
  can_run <- function(step, requires) {
    if (!is.null(busy())) {
      notify(paste0("Patiente : l'étape « ", STEP_LABELS[busy()], " » est en cours."), "warning")
      return(FALSE)
    }
    missing <- requires[!vapply(requires, done, TRUE)]
    if (length(missing)) {
      notify(paste0("Lance d'abord : ", paste(STEP_LABELS[missing], collapse = ", "), "."), "warning")
      return(FALSE)
    }
    v <- snap(step)
    bad <- names(v)[vapply(v, function(x) is.numeric(x) && anyNA(x), TRUE)]
    bad <- setdiff(bad, "amplen")
    if (length(bad)) {
      notify(paste("Valeur(s) numérique(s) manquante(s) :", paste(bad, collapse = ", ")), "error")
      return(FALSE)
    }
    TRUE
  }

  # Lance une étape (en arrière-plan si possible). `store(st, res)` range le résultat.
  # `then` : fonction appelée après une réussite (enchaînement « Tout lancer »)
  launch <- function(step, fun, args, params, store, then = NULL) {
    v <- snap(step)
    st <- reset_from(state(), step)
    state(st)
    busy(step); job_start(Sys.time())
    lf <- isolate(log_file())
    cat(sprintf("\n%s ===== %s =====\n", format(Sys.time(), "[%H:%M:%S]"), STEP_LABELS[step]),
        file = lf, append = TRUE)

    on_ok <- function(res) {
      st <- store(state(), res)
      st$params[[step]] <- params
      st$inputs[[step]] <- v
      state(st)
      el <- round(as.numeric(difftime(Sys.time(), job_start(), units = "mins")), 2)
      busy(NULL)
      log_line("✔ ", STEP_LABELS[step], " terminé (", el, " min)")
      save_session(st)
      notify(paste0(STEP_LABELS[step], " : terminé (", el, " min)"))
      if (!is.null(then)) then()
    }
    on_err <- function(e) {
      busy(NULL)
      msg <- conditionMessage(e)
      log_line("✖ ", STEP_LABELS[step], " : ", msg)
      showNotification(paste0(STEP_LABELS[step], " — ", msg), type = "error", duration = NULL)
      stop_chain(paste0("arrêt sur erreur à l'étape « ", STEP_LABELS[step], " »"))
    }
    if (ASYNC) {
      p <- promises::future_promise(run_job(fun, args, lf, PIPELINE_FILE), seed = TRUE)
      promises::then(p, onFulfilled = on_ok, onRejected = on_err)
    } else {
      res <- withProgress(message = STEP_LABELS[step], value = 0.5,
                          tryCatch(run_job(fun, args, lf, PIPELINE_FILE), error = function(e) e))
      if (inherits(res, "error")) on_err(res) else on_ok(res)
    }
    invisible(TRUE)
  }

  # =========================================================== Tout lancer
  RUN <- new.env()
  chain <- reactiveVal(NULL)     # NULL ou list(steps, i, until) pendant un enchaînement
  stop_chain <- function(why = NULL) {
    if (is.null(chain())) return(invisible())
    chain(NULL)
    if (!is.null(why)) {
      log_line("■ Tout lancer : ", why)
      notify(paste("Tout lancer :", why), "warning", NULL)
    }
  }
  # La taxonomie est sautée si aucune référence n'est indiquée
  tax_configured <- function() {
    switch(input$ref_src,
           preset = TRUE,
           path   = nzchar(input$ref_path %||% ""),
           upload = !is.null(input$ref_up))
  }
  chain_next <- function() {
    ch <- chain(); if (is.null(ch)) return()
    if (ch$i > length(ch$steps)) {
      chain(NULL)
      log_line("■ Tout lancer : terminé")
      return(notify("Tout lancer : toutes les étapes demandées sont terminées.", duration = NULL))
    }
    s <- ch$steps[ch$i]
    ch$i <- ch$i + 1; chain(ch)
    if (s == "taxfilter" && !done("taxonomy")) {
      log_line("Filtre taxonomique sauté : pas de taxonomie")
      return(chain_next())
    }
    if (s == "taxonomy" && !tax_configured()) {
      log_line("Taxonomie sautée : aucune référence indiquée")
      notify("Taxonomie sautée : aucune référence indiquée.", "warning")
      return(chain_next())
    }
    ok <- RUN[[s]](then = chain_next)
    if (!isTRUE(ok)) stop_chain(paste0("étape « ", STEP_LABELS[s], " » non lancée (voir le message)"))
  }
  observeEvent(input$run_all, {
    if (!is.null(busy()) || !is.null(chain())) return(notify("Une étape est déjà en cours.", "warning"))
    if (!done("data")) return(notify("Détecte d'abord les échantillons (onglet Données).", "warning"))
    if (stale("data")) return(notify("Les paramètres de l'onglet Données ont changé : relance la détection.", "warning"))
    pipe <- STEPS[-1]
    until <- match(input$run_until, pipe)
    # première étape à (re)faire : non faite ou paramètres modifiés
    todo <- which(!vapply(pipe, done, TRUE) | vapply(pipe, stale, TRUE))
    from <- if (isTRUE(input$run_from_start)) 1L else if (length(todo)) todo[1] else NA
    if (is.na(from) || from > until)
      return(notify("Rien à faire : ces étapes sont déjà à jour (coche « Tout recalculer » pour forcer).", "message"))
    steps <- pipe[from:until]
    log_line("■ Tout lancer : ", paste(STEP_LABELS[steps], collapse = " → "))
    chain(list(steps = steps, i = 1L))
    chain_next()
  })
  observeEvent(input$stop_all, {
    if (is.null(chain())) return()
    chain(NULL)
    log_line("■ Tout lancer : arrêt demandé (l'étape en cours se termine)")
    notify("Arrêt demandé : l'étape en cours va se terminer, les suivantes ne seront pas lancées.", "warning")
  })

  # ---- Barre de statut -------------------------------------------------------
  output$status_bar <- renderUI({
    st <- state(); b <- busy()
    if (!is.null(b)) invalidateLater(1000)
    badges <- lapply(STEPS, function(s) {
      if (identical(b, s)) {
        el <- as.integer(difftime(Sys.time(), job_start(), units = "secs"))
        span(class = "badge b-run", sprintf("⏳ %s %02d:%02d", STEP_LABELS[s], el %/% 60, el %% 60))
      } else if (stale(s, st)) span(class = "badge b-stale", title = "Paramètres modifiés depuis l'exécution",
                                     paste("⚠", STEP_LABELS[s]))
      else if (done(s, st)) span(class = "badge b-done", paste("✔", STEP_LABELS[s]))
      else span(class = "badge b-todo", STEP_LABELS[s])
    })
    ch <- chain()
    div(class = "step-badges", badges,
        if (!is.null(ch)) span(class = "text-primary",
                               sprintf(" Tout lancer : étape %d / %d", min(ch$i - 1, length(ch$steps)),
                                       length(ch$steps))),
        if (!ASYNC) span(class = "text-muted", " (mode synchrone : l'interface est bloquée pendant les calculs)"))
  })

  # ---- Journal en direct ---------------------------------------------------
  log_txt <- reactivePoll(1000, session,
    checkFunc = function() { lf <- log_file(); if (is.null(lf) || !file.exists(lf)) "" else file.info(lf)$size },
    valueFunc = function() {
      lf <- log_file()
      if (is.null(lf) || !file.exists(lf)) return("Aucun journal pour l'instant.")
      paste(utils::tail(readLines(lf, warn = FALSE), 800), collapse = "\n")
    })
  output$log <- renderText(log_txt())

  # =========================================================== 1. Données
  restore_inputs <- function(st) {
    for (s in names(st$inputs)) for (id in names(st$inputs[[s]])) {
      val <- st$inputs[[s]][[id]]
      if (id %in% DYNAMIC_SELECTS)
        updateSelectInput(session, id, choices = unique(c(if (id %in% c("repl_col", "tf_unassigned")) "",
                                                          colnames(st$meta %||% data.frame()),
                                                          colnames(st$taxa %||% matrix(nrow = 0, ncol = 0)), val)),
                          selected = val)
      else session$sendInputMessage(id, list(value = val))
    }
  }
  after_load <- function(st) {
    state(st)
    restore_inputs(st)
    if (!is.null(st$samples))
      updateSelectizeInput(session, "qsamples", choices = st$samples$sample,
                           selected = head(st$samples$sample, 2))
    if (!is.null(st$errs)) updateSelectInput(session, "err_run", choices = names(st$errs))
    update_meta_choices(st$meta)
    notify(sprintf("Session chargée : étapes faites = %s",
                   paste(STEP_LABELS[vapply(STEPS, done, TRUE, st = st)], collapse = ", ")))
  }
  load_session_file <- function(path) {
    if (!is.null(busy())) return(notify("Une étape est en cours.", "warning"))
    st <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.list(st) || is.null(st$params)) return(notify("Fichier de session invalide.", "error"))
    after_load(st)
  }
  observeEvent(input$load_rds, load_session_file(input$load_rds$datapath))
  observeEvent(input$load_outdir, {
    f <- file.path(input$outdir, "dada2_session.rds")
    if (!file.exists(f)) return(notify(paste("Pas de session dans", input$outdir), "error"))
    load_session_file(f)
  })

  observeEvent(input$scan, {
    if (!is.null(busy())) return(notify("Une étape est en cours.", "warning"))
    v <- snap("data")
    if (is.na(v$threads) || v$threads < 1) return(notify("Nombre de threads invalide.", "error"))
    p <- BUILD$data(v)
    dir.create(p$outdir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(p$outdir)) return(notify("Impossible de créer le dossier de sortie.", "error"))

    dir <- if (v$src == "upload") {
      if (is.null(input$fq_up)) return(notify("Aucun fichier uploadé.", "error"))
      d <- file.path(p$outdir, "raw_upload")
      dir.create(d, showWarnings = FALSE)
      file.copy(input$fq_up$datapath, file.path(d, input$fq_up$name), overwrite = TRUE)
      d
    } else v$fq_dir

    s <- tryCatch(detect_samples(dir, v$patF, if (isTRUE(v$paired)) v$patR else NULL,
                                 v$sep, v$sidx, isTRUE(v$multirun)),
                  error = function(e) { notify(conditionMessage(e), "error", NULL); NULL })
    if (is.null(s)) return()
    rl <- tryCatch(c(F = read_length_summary(head(s$fnF, 2))[["median"]],
                     R = if (is_paired(s)) read_length_summary(head(s$fnR, 2))[["median"]] else NA),
                   error = function(e) c(F = NA, R = NA))

    st <- reset_from(state(), "data")
    st$samples <- s; st$readlen <- rl
    st$params$data <- c(p, list(fq_dir = dir))
    st$inputs$data <- v
    if (!is.null(st$meta)) {
      st$params$data$meta_file <- file.path(p$outdir, "metadata.csv")
      write.csv(st$meta, st$params$data$meta_file)
    }
    state(st)
    cat(sprintf("%s ===== Nouveau projet : %d échantillons, %d run(s) =====\n",
                format(Sys.time(), "[%H:%M:%S]"), nrow(s), length(unique(s$run))),
        file = file.path(p$outdir, "dada2_log.txt"), append = TRUE)
    save_session(st)
    updateSelectizeInput(session, "qsamples", choices = s$sample, selected = head(s$sample, 2))
    notify(sprintf("%d échantillons détectés.", nrow(s)))
  })

  output$data_summary <- renderUI({
    st <- state()
    validate(need(!is.null(st$samples), "Détecte les échantillons ou reprends un projet."))
    s <- st$samples
    tags$p(sprintf("%d échantillons, %d run(s), %s. Longueur médiane des reads : %s pb (R1)%s.",
                   nrow(s), length(unique(s$run)),
                   if (is_paired(s)) "paired-end" else "single-end",
                   st$readlen[["F"]],
                   if (is_paired(s)) paste0(" / ", st$readlen[["R"]], " pb (R2)") else ""))
  })
  output$files_tbl <- renderDT({
    s <- state()$samples; req(s)
    df <- data.frame(Echantillon = s$sample, Run = s$run, Forward = basename(s$fnF),
                     Reverse = basename(s$fnR),
                     `Mo (F)` = round(file.size(s$fnF) / 1e6, 1), check.names = FALSE)
    datatable(df, options = list(pageLength = 15), rownames = FALSE)
  })

  # Métadonnées
  update_meta_choices <- function(meta) {
    cols <- colnames(meta %||% data.frame())
    updateSelectInput(session, "neg_col", choices = cols,
                      selected = isolate(input$neg_col) %||% cols[grep("neg|control|blank|type", cols, ignore.case = TRUE)[1]])
    updateSelectInput(session, "conc_col", choices = cols,
                      selected = isolate(input$conc_col) %||% cols[grep("conc|quant|dna", cols, ignore.case = TRUE)[1]])
    cur <- isolate(input$repl_col) %||% ""
    updateSelectInput(session, "repl_col", choices = c("(aucune : pas de réplicats)" = "", cols),
                      selected = if (cur %in% cols) cur else "")
  }
  observeEvent(input$meta, {
    f <- input$meta
    sep <- if (grepl("\\.(tsv|txt)$", f$name, ignore.case = TRUE)) "\t" else ","
    m <- tryCatch(read.csv(f$datapath, sep = sep, row.names = 1, check.names = FALSE),
                  error = function(e) { notify(paste("Lecture impossible :", conditionMessage(e)), "error"); NULL })
    if (is.null(m)) return()
    st <- state(); st$meta <- m; st$ps <- NULL   # l'objet phyloseq sera reconstruit
    if (!is.null(st$params$data)) {
      mf <- file.path(st$params$data$outdir, "metadata.csv")
      write.csv(m, mf); st$params$data$meta_file <- mf
    }
    state(st); save_session(st)
    update_meta_choices(m)
  })
  output$meta_summary <- renderUI({
    st <- state(); m <- st$meta
    if (is.null(m)) return(NULL)
    miss <- if (!is.null(st$samples)) setdiff(st$samples$sample, rownames(m)) else character(0)
    tagList(h4("Métadonnées"),
      p(sprintf("%d lignes, colonnes : %s", nrow(m), paste(colnames(m), collapse = ", "))),
      if (length(miss)) p(class = "text-danger",
                          paste("Échantillons sans métadonnées :", paste(miss, collapse = ", "))))
  })

  # =========================================================== 2. Qualité
  qdata <- eventReactive(input$qplot, {
    s <- state()$samples
    validate(need(!is.null(s), "Détecte d'abord les échantillons."),
             need(length(input$qsamples) > 0, "Choisis au moins un échantillon."))
    i <- match(input$qsamples, s$sample)
    list(F = s$fnF[i], R = s$fnR[i], agg = input$qagg, n = input$qn)
  })
  output$qF <- renderPlot({ q <- qdata(); plotQualityProfile(q$F, aggregate = q$agg, n = q$n) })
  output$qR <- renderPlot({
    q <- qdata(); validate(need(!anyNA(q$R), "Single-end."))
    plotQualityProfile(q$R, aggregate = q$agg, n = q$n)
  })

  # =========================================================== 3. Amorces
  # ---- Préréglages personnels (enregistrés dans le dossier de config de R) ----
  user_presets <- reactiveVal(load_user_presets())
  marker_choices <- function(up, selected = NULL) {
    ch <- list(`Préréglages intégrés` = MARKER_CHOICES)
    if (length(up)) ch$`Mes préréglages` <- setNames(paste0("user:", names(up)), names(up))
    updateSelectInput(session, "marker", choices = ch, selected = selected)
  }
  observeEvent(user_presets(), marker_choices(user_presets(), isolate(input$marker)), once = TRUE)

  observeEvent(input$save_preset, {
    nm <- trimws(input$preset_name %||% "")
    if (!nzchar(nm)) return(notify("Donne un nom au préréglage.", "warning"))
    vals <- lapply(setNames(PRESET_INPUTS, PRESET_INPUTS), function(i) input[[i]])
    vals <- vals[!vapply(vals, is.null, TRUE)]
    up <- user_presets(); up[[nm]] <- vals
    ok <- tryCatch({ save_user_presets(up); TRUE },
                   error = function(e) { notify(paste("Enregistrement impossible :", conditionMessage(e)), "error"); FALSE })
    if (!ok) return()
    user_presets(up)
    marker_choices(up, paste0("user:", nm))
    updateTextInput(session, "preset_name", value = "")
    notify(paste0("Préréglage « ", nm, " » enregistré (", PRESET_FILE, ")."))
  })
  observeEvent(input$delete_preset, {
    key <- input$marker
    if (!startsWith(key %||% "", "user:")) return()
    nm <- sub("^user:", "", key)
    up <- user_presets(); up[[nm]] <- NULL
    save_user_presets(up); user_presets(up)
    marker_choices(up, "custom")
    notify(paste0("Préréglage « ", nm, " » supprimé."))
  })

  observeEvent(input$apply_marker, {
    if (startsWith(input$marker, "user:")) {
      nm <- sub("^user:", "", input$marker)
      vals <- user_presets()[[nm]]
      if (is.null(vals)) return()
      for (id in names(vals)) session$sendInputMessage(id, list(value = vals[[id]]))
      return(notify(paste0("Préréglage « ", nm, " » appliqué.")))
    }
    m <- MARKERS[[input$marker]]
    if (is.null(m$primersF)) return()
    updateTextAreaInput(session, "primersF", value = m$primersF)
    updateTextAreaInput(session, "primersR", value = m$primersR)
    updateRadioButtons(session, "primer_method", selected = "cutadapt")
    updateNumericInput(session, "amplen", value = m$amplen)
    updateNumericInput(session, "truncF", value = m$truncLen[1])
    updateNumericInput(session, "truncR", value = m$truncLen[2])
    updateNumericInput(session, "trimLF", value = 0)
    updateNumericInput(session, "trimLR", value = 0)
    updateNumericInput(session, "maxEEF", value = m$maxEE[1])
    updateNumericInput(session, "maxEER", value = m$maxEE[2])
    updateNumericInput(session, "minLen", value = m$minLen)
    if (m$ref %in% names(REFS)) {
      updateRadioButtons(session, "ref_src", selected = "preset")
      updateSelectInput(session, "ref_preset", selected = m$ref)
    } else updateRadioButtons(session, "ref_src", selected = "path")
    updateTextInput(session, "taxLevels", value = m$levels)
    notify(paste("Préréglage appliqué :", m$label))
  })
  observeEvent(input$primer_method, {
    if (input$primer_method == "trimLeft") {
      updateNumericInput(session, "trimLF", value = primer_len(input$primersF))
      updateNumericInput(session, "trimLR", value = primer_len(input$primersR))
    } else {
      updateNumericInput(session, "trimLF", value = 0)
      updateNumericInput(session, "trimLR", value = 0)
    }
  }, ignoreInit = TRUE)

  primer_file_res <- reactiveVal(NULL)
  observeEvent(input$primer_file, {
    f <- input$primer_file
    res <- tryCatch(read_primer_file(f$datapath, f$name),
                    error = function(e) { notify(conditionMessage(e), "error", NULL); NULL })
    primer_file_res(res)
    if (is.null(res)) return()
    updateTextAreaInput(session, "primersF", value = paste(res$F, collapse = ", "))
    updateTextAreaInput(session, "primersR", value = paste(res$R, collapse = ", "))
    notify(sprintf("%d amorce(s) forward et %d reverse chargées.", length(res$F), length(res$R)))
  })
  output$primer_file_tbl <- renderTable({
    r <- primer_file_res(); req(r)
    r$table
  }, spacing = "xs")

  primer_res <- reactiveVal(NULL)
  observeEvent(input$check_primers, {
    st <- state()
    if (is.null(st$samples)) return(notify("Détecte d'abord les échantillons.", "warning"))
    res <- tryCatch({
      raw <- primer_check(st$samples, input$primersF, input$primersR)
      out <- cbind(Fichiers = "bruts", raw)
      if (!is.null(st$samples$trimF)) out <- rbind(out, cbind(Fichiers = "après cutadapt",
        primer_check(st$samples, input$primersF, input$primersR, which = "trimmed")))
      out
    }, error = function(e) { notify(conditionMessage(e), "error"); NULL })
    primer_res(res)
  })
  output$primer_tbl <- renderDT({
    req(primer_res())
    datatable(primer_res(), rownames = FALSE, options = list(dom = "t"))
  })

  RUN$primers <- function(then = NULL) {
    if (!can_run("primers", "data")) return()
    v <- snap("primers"); p <- BUILD$primers(v)
    if (p$method != "none" && !length(tryCatch(parse_primers(p$primersF), error = function(e) character(0))))
      return(notify("Renseigne les amorces (ou choisis « Aucun »).", "error"))
    if (p$method != "cutadapt") {
      st <- reset_from(state(), "primers")
      st$params$primers <- p; st$inputs$primers <- v
      state(st); save_session(st)
      log_line("Amorces : méthode « ", p$method, " » (pas de cutadapt)")
      notify("Méthode de retrait enregistrée.")
      if (!is.null(then)) then()
      return(TRUE)
    }
    st <- state()
    launch("primers", "run_cutadapt",
           list(samples = st$samples, outdir = outdir(), p = p, threads = threads()),
           params = p,
           then = then,
           store = function(st, res) { st$samples <- res$samples; st$cut_stats <- res$stats; st })
  }
  observeEvent(input$run_primers, RUN$primers())
  output$cut_summary <- renderUI({
    cs <- state()$cut_stats; req(cs)
    tagList(hr(), h4("cutadapt"),
      p(sprintf("%s reads en entrée → %s conservés (%.1f %%).",
                format(sum(cs$cut_in), big.mark = " "), format(sum(cs$cut_out), big.mark = " "),
                100 * sum(cs$cut_out) / max(1, sum(cs$cut_in)))),
      help("Relance « Vérifier la présence des amorces » pour contrôler les fichiers coupés."))
  })
  output$cut_tbl <- renderDT({
    cs <- state()$cut_stats; req(cs)
    cs$pct <- round(100 * cs$cut_out / pmax(cs$cut_in, 1), 1)
    datatable(cs, rownames = FALSE, options = list(pageLength = 10))
  })

  # =========================================================== 4. Filtrage
  output$overlap_ui <- renderUI({
    st <- state()
    validate(need(!is.null(st$readlen), "Détecte d'abord les échantillons."))
    pe <- is_paired(st$samples)
    method <- input$primer_method
    lpF <- primer_len(input$primersF); lpR <- primer_len(input$primersR)
    rl <- st$readlen
    # longueur des reads disponibles après retrait des amorces par cutadapt
    avF <- rl[["F"]] - if (method == "cutadapt") lpF else 0
    avR <- if (pe) rl[["R"]] - if (method == "cutadapt") lpR else 0 else NA
    tF <- input$truncF %||% 0; tR <- input$truncR %||% 0
    effF <- if (isTRUE(tF > 0)) tF else avF
    effR <- if (isTRUE(tR > 0)) tR else avR
    covF <- effF - (input$trimLF %||% 0)
    covR <- if (pe) effR - (input$trimLR %||% 0) else NA
    al <- input$amplen

    msg <- list()
    add <- function(cls, txt) msg[[length(msg) + 1]] <<- tags$li(class = cls, txt)
    add("", sprintf("Reads : %s pb (R1)%s%s.", rl[["F"]], if (pe) sprintf(" / %s pb (R2)", rl[["R"]]) else "",
                    if (method == "cutadapt") sprintf(" ; après retrait des amorces ≈ %s%s pb", avF,
                                                      if (pe) paste0(" / ", avR) else "") else ""))
    if (isTRUE(tF > avF)) add("ovl-bad", sprintf("truncLen F (%s) > longueur des reads (%s) : tous les reads seront éliminés.", tF, avF))
    if (pe && isTRUE(tR > avR)) add("ovl-bad", sprintf("truncLen R (%s) > longueur des reads (%s) : tous les reads seront éliminés.", tR, avR))
    if (method == "cutadapt" && (isTRUE(tF > 0) || isTRUE(tR > 0)))
      add("ovl-warn", "Avec cutadapt, les reads coupés n'ont pas tous la même longueur : un truncLen proche de la longueur maximale élimine les reads plus courts.")
    if (is.na(al)) {
      add("", "Longueur d'amplicon non renseignée (ou variable) : chevauchement non calculé. Pour un marqueur de longueur variable (ITS, V9), préfère truncLen = 0 avec cutadapt.")
    } else {
      if (pe && (is.na(covF) || is.na(covR))) {
        add("", "Renseigne truncLen et trimLeft pour calculer le chevauchement.")
      } else if (pe) {
        ov <- covF + covR - al
        cls <- if (ov < 12) "ovl-bad" else if (ov < 20) "ovl-warn" else "ovl-ok"
        add(cls, sprintf("Chevauchement attendu : %s + %s − %s = %s pb %s", covF, covR, al, ov,
                         if (ov < 12) "→ insuffisant, la fusion échouera (augmente truncLen ou utilise justConcatenate)."
                         else if (ov < 20) "→ juste (minOverlap = 12 par défaut)." else "→ OK."))
      } else if (isTRUE(covF < al)) add("ovl-warn", sprintf("Les reads (%s pb) ne couvrent pas tout l'amplicon (%s pb).", covF, al))
      if (method != "cutadapt" && isTRUE(rl[["F"]] > al + lpF + lpR))
        add("ovl-warn", sprintf("Les reads (%s pb) dépassent l'amplicon + amorces (%s pb) : ils contiennent l'amorce opposée et l'adaptateur. Utilise cutadapt.",
                                rl[["F"]], al + lpF + lpR))
    }
    tags$ul(msg)
  })

  RUN$filter <- function(then = NULL) {
    if (!can_run("filter", "data")) return()
    st <- state()
    if (!done("primers")) {   # pas d'étape amorces : on enregistre « aucun »
      st$params$primers <- list(method = "none"); st$inputs$primers <- snap("primers"); state(st)
    }
    v <- snap("filter"); p <- BUILD$filter(v)
    launch("filter", "run_filter",
           list(samples = st$samples, outdir = outdir(), p = p, threads = threads()),
           params = p,
           then = then,
           store = function(st, res) { st$samples <- res$samples; st$out <- res$out; st })
  }
  observeEvent(input$run_filter, RUN$filter())
  output$filter_summary <- renderUI({
    o <- state()$out
    validate(need(!is.null(o), "Règle les paramètres puis lance le filtrage."))
    tags$p(sprintf("%s → %s lectures (%.1f %%). %d / %d échantillons avec des lectures.",
                   format(sum(o[, 1]), big.mark = " "), format(sum(o[, 2]), big.mark = " "),
                   100 * sum(o[, 2]) / max(1, sum(o[, 1])), sum(o[, 2] > 0), nrow(o)))
  })
  output$filter_tbl <- renderDT({
    o <- state()$out; req(o)
    o <- as.data.frame(o); o$pct <- round(100 * o$reads.out / pmax(o$reads.in, 1), 1)
    datatable(o, options = list(pageLength = 15, order = list(list(3, "asc"))))
  })

  # =========================================================== 5. Erreurs
  RUN$errors <- function(then = NULL) {
    if (!can_run("errors", "filter")) return()
    v <- snap("errors"); p <- BUILD$errors(v)
    launch("errors", "run_learn_errors",
           list(samples = state()$samples, p = p, threads = threads()),
           params = p,
           then = then,
           store = function(st, res) { st$errs <- res; st })
  }
  observeEvent(input$run_errors, RUN$errors())
  observe({
    e <- state()$errs
    updateSelectInput(session, "err_run", choices = names(e) %||% character(0))
  })
  cur_err <- reactive({
    e <- state()$errs
    validate(need(!is.null(e), "Lance l'apprentissage des erreurs."))
    e[[input$err_run %||% names(e)[1]]] %||% e[[1]]
  })
  output$errF_plot <- renderPlot(plotErrors(cur_err()$F, nominalQ = TRUE))
  output$errR_plot <- renderPlot({
    e <- cur_err(); validate(need(!is.null(e$R), "Single-end."))
    plotErrors(e$R, nominalQ = TRUE)
  })

  # =========================================================== 6. Débruitage
  RUN$denoise <- function(then = NULL) {
    if (!can_run("denoise", "errors")) return()
    v <- snap("denoise"); p <- BUILD$denoise(v)
    st <- state()
    launch("denoise", "run_denoise",
           list(samples = st$samples, errs = st$errs, p = p, threads = threads()),
           params = p,
           then = then,
           store = function(st, res) {
             st$seqtab <- res$seqtab; st$counts <- res$counts
             lens <- nchar(colnames(res$seqtab))
             if (length(lens)) {
               updateNumericInput(session, "lenMin", value = min(lens))
               updateNumericInput(session, "lenMax", value = max(lens))
             }
             st
           })
  }
  observeEvent(input$run_denoise, RUN$denoise())
  output$denoise_summary <- renderUI({
    st <- state(); s <- st$seqtab
    validate(need(!is.null(s), "Lance le débruitage."))
    c <- st$counts
    tagList(
      p(sprintf("Table : %d échantillons × %d ASV, %s reads.", nrow(s), ncol(s),
                format(sum(s), big.mark = " "))),
      if (!all(is.na(c$merged))) {
        r <- 100 * sum(c$merged) / max(1, sum(c$denoisedF))
        p(class = if (r < 50) "ovl-bad" else if (r < 80) "ovl-warn" else "ovl-ok",
          sprintf("%.1f %% des reads débruités ont fusionné.", r))
      })
  })
  output$len_plot <- renderPlot({
    s <- state()$seqtab; req(s)
    df <- data.frame(len = nchar(colnames(s)), n_ASV = 1, reads = colSums(s))
    df <- aggregate(cbind(n_ASV, reads) ~ len, data = df, FUN = sum)
    ggplot(df, aes(len, reads)) + geom_col(fill = "steelblue") +
      geom_vline(xintercept = c(input$lenMin, input$lenMax), linetype = 2, colour = "red") +
      { if (!is.na(input$amplen)) geom_vline(xintercept = input$amplen, colour = "darkgreen") } +
      labs(x = "Longueur (pb)", y = "Nombre de reads",
           subtitle = "Rouge : bornes du filtre de longueur (onglet 7) ; vert : longueur attendue")
  })
  output$denoise_tbl <- renderDT({
    c <- state()$counts; req(c)
    if (!all(is.na(c$merged))) c$pct_merged <- round(100 * c$merged / pmax(c$denoisedF, 1), 1)
    datatable(c, rownames = FALSE, options = list(pageLength = 10))
  })

  # =========================================================== 7. Chimères
  output$len_hint <- renderUI({
    if (is.na(input$amplen)) return(NULL)
    help(sprintf("Longueur attendue : %s pb. Une fenêtre de ± 10-20 pb autour retire les produits non spécifiques.",
                 input$amplen))
  })
  RUN$chimera <- function(then = NULL) {
    if (!can_run("chimera", "denoise")) return()
    v <- snap("chimera"); p <- BUILD$chimera(v)
    launch("chimera", "run_chimera",
           list(seqtab = state()$seqtab, p = p, threads = threads()),
           params = p,
           then = then,
           store = function(st, res) { st$seqtab_len <- res$seqtab_len; st$nochim <- res$nochim; st })
  }
  observeEvent(input$run_chimera, RUN$chimera())
  output$chim_summary <- renderUI({
    st <- state()
    validate(need(!is.null(st$nochim), "Lance le filtre de longueur + chimères."))
    tags$p(sprintf("%d ASV après filtre de longueur → %d ASV non chimériques ; %.1f %% des reads conservés.",
                   ncol(st$seqtab_len), ncol(st$nochim), 100 * sum(st$nochim) / sum(st$seqtab_len)))
  })

  # =========================================================== 8. Contaminants
  RUN$post <- function(then = NULL) {
    if (!can_run("post", "chimera")) return()
    v <- snap("post"); p <- BUILD$post(v)
    st <- state()
    if (p$decontam != "none" && is.null(st$meta))
      return(notify("decontam nécessite des métadonnées (onglet Données).", "error"))
    launch("post", "run_postfilter",
           list(nochim = st$nochim, meta = st$meta, p = p),
           params = p,
           then = then,
           store = function(st, res) {
             st$final <- res$final; st$final_pcr <- res$final_pcr
             st$decontam <- res$decontam; st$repl <- res$repl; st
           })
  }
  observeEvent(input$run_post, RUN$post())
  output$post_summary <- renderUI({
    st <- state()
    validate(need(!is.null(st$final), "Optionnel : sans cette étape, la table sans chimères est utilisée telle quelle."))
    d <- st$decontam
    tagList(
      p(sprintf("Table finale : %d échantillons × %d ASV, %s reads (%.1f %% de la table sans chimères).",
                nrow(st$final), ncol(st$final), format(sum(st$final), big.mark = " "),
                100 * sum(st$final) / sum(st$nochim))),
      if (!is.null(d)) p(sprintf("decontam : %d contaminants identifiés.", sum(d$contaminant, na.rm = TRUE))),
      if (!is.null(st$repl)) p(sprintf("Réplicats : %s%s",
        if (length(st$repl$removed)) paste0(length(st$repl$removed), " PCR aberrante(s) retirée(s) (",
                                            paste(st$repl$removed, collapse = ", "), ")") else "aucune PCR retirée",
        if (isTRUE(st$repl$aggregated)) sprintf(" ; réplicats fusionnés en %d échantillons.", nrow(st$final)) else ".")))
  })
  output$repl_plot <- renderPlot({
    r <- state()$repl
    validate(need(!is.null(r$table), "Choisis une colonne d'échantillon biologique et active la détection des PCR aberrantes."))
    dw <- r$table$dissimilarite[!is.na(r$table$dissimilarite)]
    validate(need(length(dw) >= 2, "Pas assez de réplicats pour comparer."))
    df <- rbind(data.frame(d = dw, type = "Entre réplicats d'un même échantillon"),
                if (length(r$db)) data.frame(d = r$db, type = "Entre échantillons différents"))
    p <- ggplot(df, aes(d, fill = type)) + geom_density(alpha = 0.5, colour = NA) +
      scale_fill_manual(values = c("#4E79A7", "#E15759"), name = NULL) +
      labs(x = "Dissimilarité de Bray-Curtis (au centroïde des réplicats)", y = "Densité") +
      theme(legend.position = "top") + xlim(0, 1)
    if (!is.na(r$threshold %||% NA)) p <- p + geom_vline(xintercept = r$threshold, linetype = 2) +
      annotate("text", x = r$threshold, y = Inf, vjust = 1.5, hjust = -0.05,
               label = sprintf("seuil %s = %.2f", if (isTRUE(r$auto)) "auto" else "fixé", r$threshold))
    p
  })
  output$repl_tbl <- renderDT({
    r <- state()$repl; req(r$table)
    datatable(r$table[order(-r$table$aberrante, -r$table$dissimilarite), ], rownames = FALSE,
              options = list(pageLength = 15))
  })
  output$decontam_tbl <- renderDT({
    d <- state()$decontam; req(d)
    d <- d[order(-d$contaminant, d$p), ]
    d$sequence <- substr(d$sequence, 1, 40)
    datatable(d, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  final_tab <- reactive({ st <- state(); st$final_tax %||% st$final %||% st$nochim })
  # table avant filtre taxonomique (celle sur laquelle la taxonomie est assignée)
  pretax_tab <- reactive({ st <- state(); st$final %||% st$nochim })
  # métadonnées au niveau des échantillons de la table finale
  sample_meta <- reactive({
    st <- state(); m <- st$meta
    if (is.null(m)) return(NULL)
    if (isTRUE(st$repl$aggregated) && st$repl$col %in% colnames(m)) collapse_meta(m, st$repl$col) else m
  })

  # =========================================================== 9. Taxonomie
  observeEvent(input$ref_preset, {
    updateTextInput(session, "taxLevels", value = REFS[[input$ref_preset]]$levels)
  }, ignoreInit = TRUE)

  keep_upload <- function(fi) {   # copie dans outdir/refs pour garder une trace
    d <- file.path(outdir(), "refs"); dir.create(d, showWarnings = FALSE)
    dest <- file.path(d, fi$name); file.copy(fi$datapath, dest, overwrite = TRUE); dest
  }

  RUN$taxonomy <- function(then = NULL) {
    if (!can_run("taxonomy", "chimera")) return()
    v <- snap("taxonomy")
    if (v$tax_method == "idtaxa" && !HAS("DECIPHER")) return(notify("DECIPHER n'est pas installé.", "error"))
    if (v$tax_method == "idtaxa" && v$ref_src == "preset")
      return(notify("IDTAXA nécessite un trainingSet .RData (fichier local ou upload).", "error"))
    ref <- switch(v$ref_src,
                  preset = REFS[[v$ref_preset]]$url,
                  path   = v$ref_path,
                  upload = if (!is.null(input$ref_up)) keep_upload(input$ref_up))
    if (is.null(ref) || !nzchar(ref)) return(notify("Indique une base de référence.", "error"))
    species <- NULL
    if (isTRUE(v$do_species)) {
      species <- switch(v$sp_src,
                        preset = if (v$ref_src == "preset") REFS[[v$ref_preset]]$species,
                        path   = v$sp_path,
                        upload = if (!is.null(input$sp_up)) keep_upload(input$sp_up))
      if (is.null(species) || !nzchar(species))
        return(notify("Pas de fichier espèces disponible pour cette référence.", "error"))
    }
    p <- list(method = v$tax_method, ref = ref, levels = v$taxLevels, minBoot = v$minBoot,
              tryRC = v$tryRC, species = species)
    launch("taxonomy", "run_taxonomy",
           list(seqs = colnames(pretax_tab()), p = p, threads = threads()),
           params = p,
           then = then,
           store = function(st, res) { st$taxa <- res$taxa; st$boot <- res$boot; st })
  }
  observeEvent(input$run_taxonomy, RUN$taxonomy())

  # =========================================================== 10. Filtre taxonomique
  observe({
    tx <- state()$taxa
    cur <- isolate(input$tf_unassigned) %||% ""
    ranks <- colnames(tx %||% matrix(nrow = 0, ncol = 0))
    updateSelectInput(session, "tf_unassigned", choices = c("(ne pas filtrer)" = "", ranks),
                      selected = if (cur %in% ranks) cur else "")
  })
  RUN$taxfilter <- function(then = NULL) {
    if (!can_run("taxfilter", "taxonomy")) return(FALSE)
    v <- snap("taxfilter"); p <- BUILD$taxfilter(v)
    launch("taxfilter", "run_taxfilter",
           list(final = pretax_tab(), taxa = state()$taxa, p = p),
           params = p, then = then,
           store = function(st, res) { st$final_tax <- res$final_tax; st$tf_report <- res$report; st })
  }
  observeEvent(input$run_taxfilter, RUN$taxfilter())
  output$tf_summary <- renderUI({
    st <- state()
    validate(need(!is.null(st$taxa), "Assigne d'abord la taxonomie (onglet 9)."))
    validate(need(!is.null(st$final_tax), "Optionnel : sans filtre, toutes les ASV sont gardées. Laisse les champs vides et clique sur « Appliquer » pour valider l'étape."))
    pre <- pretax_tab(); ft <- st$final_tax
    p(sprintf("%d / %d ASV gardées, %s / %s reads (%.1f %%), %d échantillons.",
              ncol(ft), ncol(pre), format(sum(ft), big.mark = " "), format(sum(pre), big.mark = " "),
              100 * sum(ft) / max(1, sum(pre)), nrow(ft)))
  })
  output$tf_report <- renderTable({ r <- state()$tf_report; req(r, nrow(r) > 0); r })

  # identifiants stables : numéro de l'ASV dans la table sans chimères
  asv_ids <- reactive(paste0("ASV", match(colnames(final_tab()), colnames(state()$nochim))))

  output$tax_summary <- renderUI({
    tx <- state()$taxa
    validate(need(!is.null(tx), "Choisis une référence puis lance l'assignation."))
    pct <- round(100 * colMeans(!is.na(tx)), 1)
    tags$p(tags$b("% d'ASV assignées : "), paste(sprintf("%s %s %%", colnames(tx), pct), collapse = " · "))
  })
  output$tax_tbl <- renderDT({
    tx <- state()$taxa; req(tx)
    ft <- final_tab()
    df <- data.frame(ASV = asv_ids(), reads = colSums(ft), as.data.frame(tx[colnames(ft), , drop = FALSE]),
                     check.names = FALSE, row.names = NULL)
    datatable(df, filter = "top", rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
  })

  # =========================================================== BLAST
  blast_busy <- reactiveVal(NULL)          # heure de lancement pendant une recherche
  blast_n <- reactiveVal(0)

  blast_asv_df <- reactive({
    ft <- final_tab()
    validate(need(!is.null(ft), "Termine d'abord l'étape Chimères : la table d'ASV est nécessaire."))
    tx <- state()$taxa
    df <- data.frame(ASV = asv_ids(), reads = colSums(ft), longueur = nchar(colnames(ft)),
                     check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(tx)) df <- cbind(df, as.data.frame(tx[colnames(ft), , drop = FALSE], stringsAsFactors = FALSE))
    b <- state()$blast
    if (!is.null(b)) {
      sm <- blast_summary(b$hits, b$asv)
      df$BLAST <- ifelse(df$ASV %in% sm$ASV,
                         paste0(sm$meilleur_hit[match(df$ASV, sm$ASV)], " (",
                                sm$identite[match(df$ASV, sm$ASV)], " %)"), "")
      df$BLAST <- sub(" \\(NA %\\)$", "", df$BLAST)
    }
    df$sequence <- colnames(ft)
    rownames(df) <- NULL
    df
  })
  output$blast_asv_tbl <- renderDT({
    df <- blast_asv_df()
    df$sequence <- paste0(substr(df$sequence, 1, 25), "…")
    datatable(df, filter = "top", rownames = FALSE, selection = "multiple",
              options = list(pageLength = 25, scrollX = TRUE))
  }, server = TRUE)
  blast_proxy <- dataTableProxy("blast_asv_tbl")
  observeEvent(input$blast_select_visible, {
    rows <- input$blast_asv_tbl_rows_all
    if (length(rows) > BLAST_MAX_SEQS)
      notify(sprintf("%d lignes filtrées : seules les %d premières sont sélectionnées (limite par recherche).",
                     length(rows), BLAST_MAX_SEQS), "warning")
    selectRows(blast_proxy, head(rows, BLAST_MAX_SEQS))
  })
  observeEvent(input$blast_clear, selectRows(blast_proxy, NULL))
  blast_selected <- reactive({
    df <- blast_asv_df()
    df[input$blast_asv_tbl_rows_selected %||% integer(0), , drop = FALSE]
  })
  output$blast_sel_info <- renderUI({
    n <- nrow(blast_selected())
    p(tags$b(n), " ASV sélectionnée(s)",
      if (n > BLAST_MAX_SEQS) span(class = "text-danger", sprintf(" — maximum %d par recherche", BLAST_MAX_SEQS)))
  })

  observeEvent(input$run_blast, {
    if (!is.null(blast_busy())) return(notify("Une recherche BLAST est déjà en cours.", "warning"))
    sel <- blast_selected()
    if (!nrow(sel)) return(notify("Sélectionne au moins une ASV dans le tableau.", "warning"))
    if (nrow(sel) > BLAST_MAX_SEQS) return(notify(sprintf("Maximum %d ASV par recherche.", BLAST_MAX_SEQS), "error"))
    if (is.na(input$blast_hits) || input$blast_hits < 1) return(notify("Nombre de hits invalide.", "error"))
    if (!HAS("curl") && !nzchar(Sys.getenv("DADA2SHINY_BLAST_MOCK")))
      return(notify("Installe le package curl : install.packages('curl').", "error"))
    seqs <- setNames(sel$sequence, sel$ASV)
    p <- list(database = input$blast_db, megablast = input$blast_algo == "megablast",
              max_hits = as.integer(input$blast_hits), entrez = input$blast_entrez,
              exclude_uncultured = input$blast_uncult, email = trimws(input$blast_email),
              timeout_min = 30)
    blast_busy(Sys.time()); blast_n(length(seqs))
    lf <- isolate(log_file())
    cat(sprintf("\n%s ===== BLAST NCBI (%d ASV) =====\n", format(Sys.time(), "[%H:%M:%S]"), length(seqs)),
        file = lf, append = TRUE)
    updateTabsetPanel(session, "blast_tabs", selected = "Synthèse par ASV")

    on_ok <- function(res) {
      blast_busy(NULL)
      st <- state(); old <- st$blast
      res$hits <- res$hits[res$hits$ASV %in% res$asv, , drop = FALSE]
      # les nouveaux résultats remplacent ceux des mêmes ASV, les autres sont conservés
      if (!is.null(old)) {
        res$hits <- rbind(old$hits[!old$hits$ASV %in% res$asv, , drop = FALSE], res$hits)
        res$asv <- union(setdiff(old$asv, res$asv), res$asv)
        res$rids <- c(old$rids %||% old$rid, res$rid)
      } else res$rids <- res$rid
      res$params <- p
      st$blast <- res; state(st); save_session(st)
      log_line("✔ BLAST NCBI terminé (RID ", res$rid, ")")
      notify(sprintf("BLAST terminé : %d hits.", sum(res$hits$ASV %in% names(seqs))))
    }
    on_err <- function(e) {
      blast_busy(NULL)
      log_line("✖ BLAST NCBI : ", conditionMessage(e))
      showNotification(paste("BLAST —", conditionMessage(e)), type = "error", duration = NULL)
    }
    if (ASYNC) {
      pr <- promises::future_promise(run_job("run_ncbi_blast", list(seqs = seqs, p = p), lf, PIPELINE_FILE),
                                     seed = TRUE)
      promises::then(pr, onFulfilled = on_ok, onRejected = on_err)
    } else {
      res <- withProgress(message = "BLAST NCBI", value = 0.5,
                          tryCatch(run_job("run_ncbi_blast", list(seqs = seqs, p = p), lf, PIPELINE_FILE),
                                   error = function(e) e))
      if (inherits(res, "error")) on_err(res) else on_ok(res)
    }
    invisible(NULL)
  })

  output$blast_status <- renderUI({
    t0 <- blast_busy()
    if (!is.null(t0)) {
      invalidateLater(1000)
      el <- as.integer(difftime(Sys.time(), t0, units = "secs"))
      rid_line <- tail(grep("BLAST NCBI : RID", strsplit(log_txt(), "\n")[[1]], value = TRUE), 1)
      rid <- if (length(rid_line)) sub(".*RID ([A-Z0-9]+).*", "\\1", rid_line) else NULL
      return(div(class = "alert alert-info",
        sprintf("⏳ Recherche BLAST en cours pour %d ASV — %02d:%02d. ", blast_n(), el %/% 60, el %% 60),
        if (!is.null(rid)) tagList("RID ", tags$a(href = paste0(NCBI_BLAST_URL, "?CMD=Get&RID=", rid),
                                                  target = "_blank", rid), " (consultable sur le site du NCBI). "),
        "Tu peux continuer à utiliser l'appli."))
    }
    b <- state()$blast
    if (is.null(b)) return(NULL)
    div(class = "text-muted",
        sprintf("Résultats BLAST pour %d ASV (dernière recherche : %s, RID ", length(b$asv),
                format(b$time, "%d/%m %H:%M")),
        tags$a(href = paste0(NCBI_BLAST_URL, "?CMD=Get&RID=", b$rid), target = "_blank", b$rid),
        "). Le NCBI ne conserve les résultats en ligne que temporairement.")
  })

  blast_summary_df <- reactive({
    b <- state()$blast
    validate(need(!is.null(b), "Sélectionne des ASV dans le premier onglet puis lance « BLASTer la sélection »."))
    sm <- blast_summary(b$hits, b$asv)
    # comparaison avec l'assignation DADA2
    tx <- state()$taxa; ft <- final_tab()
    if (!is.null(tx) && !is.null(ft)) {
      seq_of <- setNames(colnames(ft), asv_ids())
      txs <- tx[seq_of[sm$ASV], , drop = FALSE]
      last <- apply(txs, 1, function(v) { v <- v[!is.na(v)]; if (length(v)) paste0(tail(v, 1), " (", names(v)[length(v)], ")") else "Non assigné" })
      sm <- cbind(sm[, 1, drop = FALSE], DADA2 = unname(last), sm[, -1])
    }
    sm
  })
  link_acc <- function(acc) ifelse(is.na(acc), "",
    sprintf('<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s" target="_blank">%s</a>', acc, acc))
  output$blast_summary_tbl <- renderDT({
    sm <- blast_summary_df(); sm$accession <- link_acc(sm$accession)
    datatable(sm, escape = FALSE, rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })
  output$blast_hits_tbl <- renderDT({
    b <- state()$blast
    validate(need(!is.null(b), "Pas encore de résultats BLAST."))
    h <- b$hits; h$accession <- link_acc(h$accession)
    datatable(h, escape = FALSE, rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })
  output$dl_blast_summary <- downloadHandler("blast_synthese.csv",
    function(file) write.csv(blast_summary_df(), file, row.names = FALSE))
  output$dl_blast_hits <- downloadHandler("blast_hits.csv", function(file) {
    b <- state()$blast; req(b); write.csv(b$hits, file, row.names = FALSE)
  })

  # =========================================================== Outils en arrière-plan
  tool_busy <- reactiveVal(NULL)
  tool_start <- reactiveVal(NULL)
  run_tool <- function(label, fun, args, on_ok) {
    if (!is.null(tool_busy())) return(notify(paste0("« ", tool_busy(), " » est déjà en cours."), "warning"))
    tool_busy(label); tool_start(Sys.time())
    lf <- isolate(log_file()) %||% file.path(tempdir(), "dada2_log.txt")
    cat(sprintf("\n%s ===== %s =====\n", format(Sys.time(), "[%H:%M:%S]"), label), file = lf, append = TRUE)
    ok <- function(res) { tool_busy(NULL); on_ok(res); notify(paste(label, ": terminé")) }
    ko <- function(e) {
      tool_busy(NULL)
      cat(format(Sys.time(), "[%H:%M:%S] "), "✖ ", label, " : ", conditionMessage(e), "\n", file = lf, append = TRUE, sep = "")
      showNotification(paste(label, "—", conditionMessage(e)), type = "error", duration = NULL)
    }
    if (ASYNC) {
      pr <- promises::future_promise(run_job(fun, args, lf, PIPELINE_FILE), seed = TRUE)
      promises::then(pr, onFulfilled = ok, onRejected = ko)
    } else {
      res <- withProgress(message = label, value = 0.5,
                          tryCatch(run_job(fun, args, lf, PIPELINE_FILE), error = function(e) e))
      if (inherits(res, "error")) ko(res) else ok(res)
    }
    invisible(NULL)
  }
  tool_status <- function(label) renderUI({
    if (!identical(tool_busy(), label)) return(NULL)
    invalidateLater(1000)
    el <- as.integer(difftime(Sys.time(), tool_start(), units = "secs"))
    div(class = "alert alert-info", sprintf("⏳ %s en cours — %02d:%02d (détails dans le Journal)", label, el %/% 60, el %% 60))
  })

  # =========================================================== Communauté mock
  mock_tab <- reactive({ st <- state(); if (identical(input$mock_table, "final")) final_tab() else st$nochim })
  observe({
    t <- mock_tab(); cur <- isolate(input$mock_samples)
    sn <- rownames(t %||% matrix(nrow = 0, ncol = 0))
    updateSelectizeInput(session, "mock_samples", choices = sn,
                         selected = intersect(cur %||% sn[grepl("mock", sn, ignore.case = TRUE)], sn))
  })
  mock_res <- reactiveVal(NULL)
  observeEvent(input$run_mock, {
    t <- mock_tab()
    if (is.null(t)) return(notify("Termine d'abord l'étape Chimères.", "warning"))
    if (!length(input$mock_samples)) return(notify("Choisis le ou les échantillons mock.", "warning"))
    if (is.null(input$mock_fasta)) return(notify("Charge le FASTA des séquences attendues.", "warning"))
    refs <- tryCatch(read_fasta_simple(input$mock_fasta$datapath), error = function(e) NULL)
    if (!length(refs)) return(notify("FASTA illisible.", "error"))
    ids <- paste0("ASV", match(colnames(t), colnames(state()$nochim)))
    run_tool("Évaluation du mock", "mock_eval",
             list(counts = t, mock_samples = input$mock_samples, ref_seqs = refs, asv_ids = ids,
                  max_diff = input$mock_maxdiff %||% 3),
             function(res) mock_res(res))
  })
  output$mock_status <- tool_status("Évaluation du mock")
  output$mock_summary <- renderTable({ r <- mock_res(); validate(need(!is.null(r), "Choisis l'échantillon mock, charge les séquences attendues puis clique sur « Évaluer »."));  r$summary })
  output$mock_perref <- renderTable({ req(mock_res())$per_ref })
  output$mock_plot <- renderPlot({
    r <- mock_res(); req(r)
    a <- r$asv; a$ASV <- factor(a$ASV, levels = a$ASV)
    ggplot(a, aes(ASV, reads, fill = classe)) + geom_col() + scale_y_log10() +
      scale_fill_manual(values = setNames(c("#59A14F", "#F28E2B", "#E15759"),
                                          c("Exacte", sprintf("Variante d'erreur (1-%d diff.)", r$max_diff),
                                            "Inattendue (contaminant, chimère…)")), name = NULL) +
      labs(x = NULL, y = "Reads (échelle log)", title = "ASV du mock, de la plus à la moins abondante") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7), legend.position = "top")
  })
  output$mock_asv <- renderDT({
    r <- mock_res(); req(r)
    a <- r$asv; a$sequence <- paste0(substr(a$sequence, 1, 30), "…")
    datatable(a, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })
  output$dl_mock <- downloadHandler("mock_evaluation.csv", function(file) {
    r <- mock_res(); req(r); write.csv(r$asv, file, row.names = FALSE)
  })

  # =========================================================== Outil : base de référence
  ref_res <- reactiveVal(NULL)
  observeEvent(input$run_ref, {
    if (is.null(input$rf_fasta)) return(notify("Charge le FASTA des séquences.", "warning"))
    if (input$rf_source == "table" && is.null(input$rf_table)) return(notify("Charge le tableau de taxonomie.", "warning"))
    if (isTRUE(input$rf_idtaxa) && !HAS("DECIPHER")) return(notify("DECIPHER n'est pas installé.", "error"))
    keep_name <- function(fi) { d <- file.path(tempdir(), "ref_in"); dir.create(d, showWarnings = FALSE)
                                f <- file.path(d, fi$name); file.copy(fi$datapath, f, overwrite = TRUE); f }
    od <- file.path(normalizePath(input$outdir, mustWork = FALSE), "refs")
    p <- list(tax_source = input$rf_source, header_id = input$rf_header_id, sep = input$rf_sep,
              table = if (!is.null(input$rf_table)) keep_name(input$rf_table),
              table_id = input$rf_table_id, table_tax = input$rf_table_tax, levels = input$rf_levels,
              strip_prefix = input$rf_strip, min_len = input$rf_minlen, conflicts = input$rf_conflicts,
              species_in_train = input$rf_species_train, idtaxa = input$rf_idtaxa, name = input$rf_name)
    run_tool("Préparation de la base de référence", "format_reference",
             list(fasta = keep_name(input$rf_fasta), p = p, outdir = od),
             function(res) ref_res(res))
  })
  output$ref_status <- tool_status("Préparation de la base de référence")
  output$ref_report <- renderTable({ r <- ref_res(); validate(need(!is.null(r), "Charge tes séquences, indique où est la taxonomie, puis « Préparer la base »."));  r$report })
  output$ref_ranks <- renderTable({ req(ref_res())$ranks })
  output$ref_incons <- renderDT({ r <- ref_res(); req(r); datatable(r$inconsistencies, rownames = FALSE, options = list(pageLength = 10)) })
  output$ref_files <- renderUI({
    r <- ref_res(); req(r)
    f <- r$files
    tagList(
      h4("Fichiers produits"),
      tags$ul(lapply(names(f), function(k) tags$li(tags$code(f[[k]])))),
      actionButton("use_ref", "Utiliser cette base dans l'onglet Taxonomie", icon = icon("arrow-right"),
                   class = "btn-success"),
      help("Les fichiers sont dans le sous-dossier « refs » du dossier de sortie."), hr())
  })
  observeEvent(input$use_ref, {
    r <- ref_res(); req(r)
    use_idtaxa <- !is.null(r$files$idtaxa)
    updateRadioButtons(session, "tax_method", selected = if (use_idtaxa) "idtaxa" else "assignTaxonomy")
    updateRadioButtons(session, "ref_src", selected = "path")
    updateTextInput(session, "ref_path", value = if (use_idtaxa) r$files$idtaxa else r$files$train)
    updateTextInput(session, "taxLevels", value = paste(r$levels, collapse = ","))
    if (!is.null(r$files$species)) {
      updateCheckboxInput(session, "do_species", value = TRUE)
      updateRadioButtons(session, "sp_src", selected = "path")
      updateTextInput(session, "sp_path", value = r$files$species)
    }
    updateNavbarPage(session, "nav", selected = "9. Taxonomie")
    notify("Base sélectionnée dans l'onglet Taxonomie.")
  })

  # =========================================================== Suivi
  track <- reactive(make_track(state()))
  output$track_tbl <- renderDT({
    tr <- track(); validate(need(!is.null(tr), "Rien à suivre pour l'instant."))
    datatable(tr, options = list(pageLength = 25, order = list(list(ncol(tr), "asc"))))
  })
  output$track_plot <- renderPlot({
    tr <- track(); req(tr, ncol(tr) > 2)
    steps <- setdiff(colnames(tr), "pct_conserve")
    long <- do.call(rbind, lapply(steps, function(s)
      data.frame(sample = rownames(tr), step = s, pct = 100 * tr[[s]] / pmax(tr[[1]], 1))))
    long$step <- factor(long$step, levels = steps)
    ggplot(long, aes(step, pct, group = sample)) +
      geom_line(alpha = 0.4) + geom_point(size = 1) +
      labs(x = NULL, y = "% des reads d'entrée", title = "Reads conservés à chaque étape (une ligne par échantillon)") +
      ylim(0, 100)
  })

  # =========================================================== Export
  asv_table <- reactive({
    ft <- final_tab(); req(ft)
    tx <- state()$taxa
    df <- data.frame(ASV = asv_ids(), sequence = colnames(ft), length = nchar(colnames(ft)),
                     check.names = FALSE)
    if (!is.null(tx)) df <- cbind(df, as.data.frame(tx[colnames(ft), , drop = FALSE]))
    b <- state()$blast
    if (!is.null(b)) {   # meilleur hit BLAST pour les ASV blastées
      sm <- blast_summary(b$hits, b$asv)
      i <- match(df$ASV, sm$ASV)
      df$BLAST_meilleur_hit <- sm$meilleur_hit[i]
      df$BLAST_identite <- sm$identite[i]
      df$BLAST_couverture <- sm$couverture[i]
      df$BLAST_accession <- sm$accession[i]
    }
    df <- cbind(df, as.data.frame(t(ft), check.names = FALSE))
    rownames(df) <- NULL
    df
  })
  output$dl_asv <- downloadHandler("asv_table.csv", function(file) write.csv(asv_table(), file, row.names = FALSE))
  output$dl_fasta <- downloadHandler("asv_sequences.fasta", function(file) {
    ft <- final_tab(); req(ft); writeLines(paste0(">", asv_ids(), "\n", colnames(ft)), file)
  })
  output$dl_track <- downloadHandler("track_reads.csv", function(file) { req(track()); write.csv(track(), file) })
  output$dl_boot <- downloadHandler("taxonomy_confidence.csv", function(file) {
    b <- state()$boot; req(b)
    write.csv(data.frame(ASV = asv_ids(), b[colnames(final_tab()), , drop = FALSE], row.names = NULL), file, row.names = FALSE)
  })
  output$dl_decontam <- downloadHandler("decontam.csv", function(file) {
    d <- state()$decontam; req(d); write.csv(d, file, row.names = FALSE)
  })
  output$dl_ps <- downloadHandler("phyloseq.rds", function(file) { req(state()$ps); saveRDS(state()$ps, file) })
  output$dl_script <- downloadHandler("dada2_pipeline.R", function(file) {
    req(state()$samples); writeLines(generate_script(state()), file)
  })
  output$dl_session <- downloadHandler("dada2_session.rds", function(file) saveRDS(state(), file))

  # =========================================================== Exploration
  if (has_phyloseq) {
    build_ps <- function(quiet = FALSE) {
      st <- state(); ft <- final_tab()
      if (is.null(ft)) {
        if (!quiet) notify("Termine d'abord l'étape Chimères.", "warning")
        return(invisible(FALSE))
      }
      ids <- setNames(asv_ids(), colnames(ft))   # mêmes identifiants que dans l'export
      ps <- tryCatch({
        sm <- sample_meta()
        samdf <- if (!is.null(sm)) {
          common <- intersect(rownames(ft), rownames(sm))
          if (!length(common)) stop("Aucun nom commun entre métadonnées et table ASV.")
          ft <- ft[common, , drop = FALSE]
          sm[common, , drop = FALSE]
        } else data.frame(Echantillon = rownames(ft), row.names = rownames(ft))
        parts <- list(phyloseq::otu_table(ft, taxa_are_rows = FALSE), phyloseq::sample_data(samdf))
        if (!is.null(st$taxa)) parts <- c(parts, list(phyloseq::tax_table(st$taxa[colnames(ft), , drop = FALSE])))
        ps <- do.call(phyloseq::phyloseq, parts)
        excl <- parse_values(input$exclude)
        if (length(excl)) ps <- phyloseq::prune_samples(!phyloseq::sample_names(ps) %in% excl, ps)
        ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) > 0, ps)
        ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
        if (HAS("Biostrings")) {
          dna <- Biostrings::DNAStringSet(phyloseq::taxa_names(ps))
          names(dna) <- phyloseq::taxa_names(ps)
          ps <- phyloseq::merge_phyloseq(ps, dna)
        }
        phyloseq::taxa_names(ps) <- unname(ids[phyloseq::taxa_names(ps)])
        ps
      }, error = function(e) { notify(conditionMessage(e), "error"); NULL })
      if (is.null(ps)) return(invisible(FALSE))
      st$ps <- ps; state(st)

      vars <- colnames(phyloseq::sample_data(ps))
      keep_sel <- function(cur, choices, default) if (!is.null(cur) && cur %in% choices) cur else default
      updateSelectInput(session, "xvar", choices = vars, selected = keep_sel(input$xvar, vars, vars[1]))
      updateSelectInput(session, "colvar", choices = c("(aucune)" = "", vars), selected = keep_sel(input$colvar, vars, ""))
      updateSelectInput(session, "ab_group", choices = c("(non : une barre par échantillon)" = "", vars),
                        selected = keep_sel(input$ab_group, vars, ""))
      updateSelectInput(session, "ab_facet", choices = c("(aucune)" = "", vars),
                        selected = keep_sel(input$ab_facet, vars, ""))
      updateSelectInput(session, "kr_group", choices = vars, selected = keep_sel(input$kr_group, vars, vars[1]))
      ranks <- phyloseq::rank_names(ps) %||% character(0)
      def_rank <- if ("Family" %in% ranks) "Family" else if (length(ranks)) ranks[min(5, length(ranks))] else "ASV"
      updateSelectInput(session, "ab_level", choices = c("ASV (sans taxonomie)" = "ASV", ranks),
                        selected = keep_sel(input$ab_level, c("ASV", ranks), def_rank))
      def_lab <- if ("Genus" %in% ranks) "Genus" else ""
      updateSelectInput(session, "ab_label", choices = c("(rien)" = "", ranks),
                        selected = keep_sel(input$ab_label, ranks, def_lab))
      if (!quiet) notify(sprintf("phyloseq : %d échantillons, %d ASV.", phyloseq::nsamples(ps), phyloseq::ntaxa(ps)))
      invisible(TRUE)
    }
    observeEvent(input$build_ps, build_ps())
    # Construction automatique à l'ouverture de l'onglet (ou après une étape qui a invalidé l'objet)
    observeEvent(list(input$nav, state()$taxa, final_tab()), {
      if (identical(input$nav, "Exploration") && is.null(state()$ps) && !is.null(final_tab()))
        build_ps(quiet = TRUE)
    }, ignoreNULL = FALSE)

    colv <- reactive(if (nzchar(input$colvar %||% "")) input$colvar else NULL)
    ps <- reactive({
      p <- state()$ps
      validate(need(!is.null(p), "Pas encore de table à explorer : termine au moins l'étape Chimères."))
      p
    })
    ps_parts <- reactive({
      p <- ps()
      counts <- as(phyloseq::otu_table(p), "matrix")
      if (phyloseq::taxa_are_rows(p)) counts <- t(counts)
      tt <- phyloseq::tax_table(p, errorIfNULL = FALSE)
      list(counts = counts, tax = if (is.null(tt)) NULL else as(tt, "matrix"),
           ids = phyloseq::taxa_names(p), meta = as(phyloseq::sample_data(p), "data.frame"))
    })

    output$ps_print <- renderPrint(ps())
    output$dl_ps2 <- downloadHandler("phyloseq.rds", function(file) saveRDS(ps(), file))
    # ---- Raréfaction ----------------------------------------------------------
    observe({   # profondeur par défaut : le plus petit échantillon
      d <- tryCatch(ps_parts(), error = function(e) NULL); req(d)
      if (is.na(isolate(input$rare_depth) %||% NA))
        updateNumericInput(session, "rare_depth", value = min(rowSums(d$counts)))
    })
    rare_curves <- reactive(rarefaction_curves(ps_parts()$counts))
    output$rare_plot <- renderPlot({
      rc <- rare_curves(); d <- ps_parts()
      cv <- colv()
      if (!is.null(cv)) rc$groupe <- as.character(d$meta[rc$sample, cv])
      ends <- rc[!duplicated(rc$sample, fromLast = TRUE), ]
      p <- ggplot(rc, aes(depth, richness, group = sample)) +
        { if (is.null(cv)) geom_line(colour = "#4E79A7", alpha = 0.7) else geom_line(aes(colour = groupe), alpha = 0.8) } +
        geom_text(data = ends, aes(label = sample), hjust = -0.1, size = 2.8, colour = "grey30") +
        labs(x = "Nombre de reads", y = "Nombre d'ASV attendu", colour = cv) +
        scale_x_continuous(expand = expansion(mult = c(0.01, 0.12)))
      if (!is.na(input$rare_depth %||% NA)) p <- p + geom_vline(xintercept = input$rare_depth, linetype = 2)
      p
    })
    output$rare_info <- renderUI({
      d <- ps_parts(); dep <- input$rare_depth
      req(!is.na(dep %||% NA))
      below <- rownames(d$counts)[rowSums(d$counts) < dep]
      tagList(
        p(sprintf("Profondeur %s reads : %d / %d échantillons gardés.", format(dep, big.mark = " "),
                  nrow(d$counts) - length(below), nrow(d$counts))),
        if (length(below)) p(class = "text-warning", "Exclus (moins de reads) : ", paste(below, collapse = ", ")))
    })
    rare <- eventReactive(input$run_rare, {
      d <- ps_parts()
      validate(need(!is.na(input$rare_depth %||% NA), "Indique une profondeur de raréfaction."))
      res <- withProgress(message = "Raréfactions répétées", value = 0.5,
                          tryCatch(rarefied_diversity(d$counts, input$rare_depth, input$rare_iter %||% 20),
                                   error = function(e) e))
      validate(need(!inherits(res, "error"), if (inherits(res, "error")) conditionMessage(res)))
      res
    })
    rare_or_msg <- function() {
      validate(need(input$run_rare > 0, "Lance d'abord « Calculer la diversité raréfiée » dans l'onglet Raréfaction."))
      rare()
    }
    alpha_rare_df <- reactive({
      r <- rare_or_msg(); m <- ps_parts()$meta
      a <- r$alpha
      long <- do.call(rbind, lapply(c("Observed", "Shannon", "Simpson"), function(k)
        data.frame(sample = a$sample, metric = k, value = a[[k]], sd = a[[paste0(k, "_sd")]])))
      cbind(long, m[long$sample, , drop = FALSE])
    })

    output$rich_plot <- renderPlot({
      req(input$xvar)
      if (identical(input$alpha_method, "rare")) {
        df <- alpha_rare_df(); r <- rare()
        ggplot(df, aes(.data[[input$xvar]], value, colour = if (!is.null(colv())) .data[[colv()]])) +
          geom_pointrange(aes(ymin = value - sd, ymax = value + sd), position = position_jitter(width = 0.15, seed = 1)) +
          facet_wrap(~metric, scales = "free_y") +
          labs(y = "Moyenne ± écart-type", colour = colv(),
               title = sprintf("Diversité alpha moyenne sur %d raréfactions à %s reads", r$n_iter,
                               format(r$depth, big.mark = " "))) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      } else
        phyloseq::plot_richness(ps(), x = input$xvar, measures = c("Observed", "Shannon", "Simpson"), color = colv())
    })
    output$dl_alpha_rare <- downloadHandler("diversite_rarefiee.csv", function(file)
      write.csv(rare()$alpha, file, row.names = FALSE))

    output$nmds_plot <- renderPlot({
      validate(need(phyloseq::nsamples(ps()) >= 3, "Au moins 3 échantillons nécessaires."))
      if (identical(input$nmds_method, "rare")) {
        r <- rare_or_msg()
        validate(need(attr(r$bray, "Size") >= 3, "Au moins 3 échantillons à la profondeur choisie."))
        ord <- vegan::metaMDS(r$bray, trace = 0)
        sc <- as.data.frame(vegan::scores(ord, display = "sites"))
        names(sc)[1:2] <- c("NMDS1", "NMDS2")
        m <- ps_parts()$meta; cv <- colv()
        if (!is.null(cv)) sc$groupe <- as.character(m[rownames(sc), cv])
        ggplot(sc, aes(NMDS1, NMDS2)) +
          { if (is.null(cv)) geom_point(size = 3, colour = "#4E79A7") else geom_point(aes(colour = groupe), size = 3) } +
          labs(colour = cv, title = sprintf("NMDS sur Bray-Curtis moyen (%d raréfactions à %s reads) — stress = %.3f",
                                            r$n_iter, format(r$depth, big.mark = " "), ord$stress))
      } else {
        prop <- phyloseq::transform_sample_counts(ps(), function(x) x / sum(x))
        ord <- phyloseq::ordinate(prop, method = "NMDS", distance = "bray", trace = 0)
        phyloseq::plot_ordination(prop, ord, color = colv(),
                                  title = sprintf("NMDS Bray-Curtis (stress = %.3f)", ord$stress))
      }
    })

    # ---- Abondances relatives ------------------------------------------------
    ab <- reactive({
      d <- ps_parts()
      lvl <- input$ab_level %||% "ASV"
      validate(need(lvl == "ASV" || !is.null(d$tax), "Assigne d'abord la taxonomie (onglet 9), ou choisis le niveau ASV."))
      grp_var <- input$ab_group %||% ""
      groups <- if (nzchar(grp_var)) d$meta[[grp_var]] else NULL
      top <- max(1, input$ab_top %||% 15, na.rm = TRUE)
      label <- if (nzchar(input$ab_label %||% "")) input$ab_label else NULL
      mat <- abund_matrix(d$counts, d$tax, d$ids, lvl, label, groups, top)
      full <- abund_matrix(d$counts, d$tax, d$ids, lvl, label, groups, .Machine$integer.max)
      facet_var <- input$ab_facet %||% ""
      facet <- if (nzchar(facet_var) && is.null(groups))
        data.frame(sample = rownames(d$meta), facet = as.character(d$meta[[facet_var]])) else NULL
      list(mat = mat, order = sample_order(full, input$ab_order %||% "name"), facet = facet,
           legend = if (lvl == "ASV") "ASV" else lvl,
           title = sprintf("Abondance relative — %s (%s premiers)", if (lvl == "ASV") "ASV" else lvl,
                           min(top, nrow(full))))
    })
    ab_plot_obj <- reactive({
      a <- ab()
      plot_abund(a$mat, a$order, input$ab_type, a$facet, a$title, a$legend)
    })
    ab_height <- function() {
      a <- tryCatch(ab(), error = function(e) NULL)
      if (is.null(a)) return(400)
      if (identical(input$ab_type, "heatmap")) max(350, 24 * nrow(a$mat) + 160) else max(500, 18 * nrow(a$mat) + 120)
    }
    output$ab_plot <- renderPlot(ab_plot_obj(), height = ab_height)
    ab_width <- function() max(7, 0.28 * ncol(ab()$mat) + 5)
    output$dl_ab_pdf <- downloadHandler("abondances_relatives.pdf", function(file)
      ggsave(file, ab_plot_obj(), width = ab_width(), height = ab_height() / 80, limitsize = FALSE))
    output$dl_ab_png <- downloadHandler("abondances_relatives.png", function(file)
      ggsave(file, ab_plot_obj(), width = ab_width(), height = ab_height() / 80, dpi = 300, limitsize = FALSE))
    output$dl_ab_csv <- downloadHandler("abondances_relatives.csv", function(file) {
      m <- ab()$mat
      write.csv(data.frame(Taxon = rownames(m), round(m, 4), check.names = FALSE), file, row.names = FALSE)
    })

    # ---- Krona ---------------------------------------------------------------
    krona_dir <- file.path(tempdir(), paste0("krona_", session$token))
    dir.create(krona_dir, showWarnings = FALSE)
    addResourcePath(basename(krona_dir), krona_dir)
    session$onSessionEnded(function() {
      try(removeResourcePath(basename(krona_dir)), silent = TRUE)
      unlink(krona_dir, recursive = TRUE)
    })
    krona_lines <- reactive({
      d <- ps_parts()
      mode <- input$kr_mode %||% "sample"
      groups <- if (mode == "group") {
        validate(need(nzchar(input$kr_group %||% ""), "Choisis une variable de regroupement."))
        d$meta[[input$kr_group]]
      }
      ds <- krona_datasets(d$counts, mode, groups)
      krona_html(ds, d$tax, d$ids, add_asv = isTRUE(input$kr_asv) || is.null(d$tax),
                 magnitude_label = if (mode == "group") "% moyen" else "Reads")
    })
    output$krona_frame <- renderUI({
      h <- krona_lines()
      f <- file.path(krona_dir, "krona.html")
      writeLines(enc2utf8(h), f, useBytes = TRUE)
      tagList(
        if (is.null(ps_parts()$tax)) p(class = "text-warning",
          "Pas de taxonomie : le Krona ne montre que les ASV. Assigne la taxonomie (onglet 9) pour l'arbre complet."),
        tags$iframe(src = sprintf("%s/krona.html?v=%s", basename(krona_dir), as.numeric(Sys.time())),
                    style = "width:100%; height:780px; border:1px solid #ddd;"))
    })
    output$dl_krona <- downloadHandler("krona.html", function(file)
      writeLines(enc2utf8(krona_lines()), file, useBytes = TRUE))
  }
}

shinyApp(ui, server)

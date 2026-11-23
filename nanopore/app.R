# =============================================================================
#  Nanopore Shiny — amplicons Oxford Nanopore : du FASTQ à la table de séquences
#  Entièrement en R (fonctionne sous Windows) ; Dorado polish en option.
#
#  Lancer : shiny::runApp("chemin/vers/dada2shiny/nanopore")
#  Dépendances : shiny, DT, ggplot2, future, promises ; Bioconductor : dada2,
#                DECIPHER, ShortRead, Biostrings ; optionnel : phyloseq, Rsamtools
# =============================================================================

suppressPackageStartupMessages({
  library(shiny); library(DT); library(ggplot2)
})

REPO <- normalizePath("..")
PIPELINE_FILE <- file.path(REPO, "pipeline", "pipeline.R")
for (f in c("pipeline.R", "curation.R", "nanopore.R", "viz.R", "presets.R"))
  source(file.path(REPO, "pipeline", f))

HAS <- function(p) requireNamespace(p, quietly = TRUE)
ASYNC <- isTRUE(getOption("dada2shiny.async", TRUE)) && HAS("future") && HAS("promises")
if (ASYNC) { future::plan(future::multisession, workers = 2); options(future.globals.maxSize = Inf) }
options(shiny.maxRequestSize = 50 * 1024^3)
theme_set(theme_bw())

DORADO <- Sys.which("dorado")
STEPS <- c("data", "prep", "consensus", "table", "taxonomy", "taxfilter")
LABELS <- c(data = "Données", prep = "Filtre & amorces", consensus = "Consensus", table = "Table",
            taxonomy = "Taxonomie", taxfilter = "Filtre taxo")
FIELDS <- list(data = c("samples", "qc"), prep = "prep", consensus = c("cons", "clusters"),
               table = c("tab", "cons_tab", "chimeras"), taxonomy = c("taxa", "boot"),
               taxfilter = c("tab_tax", "tf_report"))

num <- function(id, label, value, min = NA, max = NA, step = NA)
  numericInput(id, label, value = value, min = min, max = max, step = step)
help <- function(...) helpText(...)
run_btn <- function(id, label) actionButton(id, label, class = "btn-primary", icon = icon("play"))

# =============================================================================
# UI
# =============================================================================
ui <- navbarPage(
  "Nanopore Shiny", id = "nav",
  header = tagList(
    tags$style(HTML("
      .badges .badge {font-size:12px; margin-right:4px; padding:5px 8px;}
      .b-done{background:#2e7d32;} .b-run{background:#1565c0;} .b-todo{background:#9e9e9e;}
      #log {max-height:78vh; overflow-y:auto; font-size:12px; white-space:pre-wrap;}")),
    div(class = "container-fluid",
        div(style = "margin:-12px 0 4px 0;",
            actionButton("run_all", "Tout lancer", class = "btn-success btn-sm", icon = icon("forward")),
            span(class = "text-muted", " enchaîne les étapes non faites avec les réglages actuels")),
        uiOutput("status"))),

  # ------------------------------------------------------------------ Données
  tabPanel("1. Données",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Reprendre un projet"),
        actionButton("load_outdir", "Recharger la session du dossier de sortie", icon = icon("folder-open")),
        hr(), h4("Nouveau projet"),
        textInput("fq_dir", "Dossier des FASTQ (sortie MinKNOW / Dorado)", getwd()),
        checkboxInput("recursive", "Un sous-dossier par code-barres (fastq_pass/barcode01/…)", TRUE),
        help("Sinon : un fichier FASTQ (.fastq ou .fastq.gz) par échantillon, nommé comme l'échantillon."),
        textInput("outdir", "Dossier de sortie", file.path(getwd(), "nanopore_out")),
        num("threads", "Threads", max(1, parallel::detectCores() - 1), min = 1),
        run_btn("scan", "Détecter les échantillons"),
        hr(),
        fileInput("meta", "Métadonnées (CSV, 1re colonne = échantillon ; optionnel)", accept = c(".csv", ".tsv", ".txt"))
      ),
      mainPanel(uiOutput("data_info"), DTOutput("samples_tbl"))
    )),

  # ------------------------------------------------------------------ QC
  tabPanel("2. Qualité",
    fluidRow(column(3, num("qc_n", "Reads lus par échantillon", 5000, min = 500, step = 1000)),
             column(3, br(), actionButton("run_qc", "Calculer", class = "btn-primary"))),
    fluidRow(column(6, plotOutput("qc_len", height = 420)), column(6, plotOutput("qc_q", height = 420))),
    DTOutput("qc_tbl")),

  # ------------------------------------------------------------------ Filtre & amorces
  tabPanel("3. Filtre & amorces",
    sidebarLayout(
      sidebarPanel(width = 4,
        selectInput("marker", "Préréglage d'amorces", MARKER_CHOICES),
        actionButton("apply_marker", "Appliquer", icon = icon("wand-magic-sparkles"), class = "btn-sm"),
        textAreaInput("primersF", "Amorce(s) forward (IUPAC, virgules)", "", rows = 2),
        textAreaInput("primersR", "Amorce(s) reverse", "", rows = 2),
        help("Laisser vide si les amorces sont déjà retirées et les reads orientés."),
        fluidRow(column(6, num("primer_err", "Erreurs tolérées dans l'amorce", 0.2, min = 0, max = 0.4, step = 0.05)),
                 column(6, num("window", "Fenêtre de recherche aux extrémités (pb)", 150, min = 30))),
        checkboxInput("require_both", "Exiger les deux amorces (recommandé)", TRUE),
        hr(), h5("Filtre des reads (avant découpe)"),
        fluidRow(column(4, num("min_len", "Long. min", 200, min = 1)),
                 column(4, num("max_len", "Long. max", 2000, min = 1)),
                 column(4, num("min_q", "Q moyenne min", 10, min = 0))),
        run_btn("run_prep", "Filtrer, orienter et découper")
      ),
      mainPanel(uiOutput("prep_info"), DTOutput("prep_tbl"), plotOutput("prep_len", height = 300))
    )),

  # ------------------------------------------------------------------ Consensus
  tabPanel("4. Consensus",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons("mode", "Type de données",
          c("Metabarcoding (mélanges : tous les groupes assez fournis)" = "metabarcoding",
            "Barcoding (un spécimen : groupes majoritaires)" = "barcoding")),
        fluidRow(column(6, num("min_reads", "Reads min. par consensus", 10, min = 3)),
                 column(6, conditionalPanel("input.mode == 'barcoding'",
                          num("min_frac", "Part min. des reads", 0.1, min = 0, max = 1, step = 0.05)))),
        num("cluster_cutoff", "Distance de regroupement des reads", 0.08, min = 0.01, max = 0.3, step = 0.01),
        help("Environ 2 à 3 fois le taux d'erreur des reads : 0,08 pour du hac/sup R10.4.1, plus pour des reads plus bruités."),
        checkboxInput("split", "Séparer les variantes proches à l'intérieur d'un groupe", TRUE),
        conditionalPanel("input.split",
          fluidRow(column(6, num("split_freq", "Fréquence min. de la 2e variante", 0.2, min = 0.05, max = 0.5, step = 0.05)),
                   column(6, num("split_sites", "Positions min. qui co-varient", 3, min = 1)))),
        hr(), h5("Polissage"),
        checkboxInput("hp_correct", "Corriger la longueur des homopolymères", TRUE),
        radioButtons("polisher", NULL,
          c("Polisseur R (votes pondérés par la qualité, équilibrés entre brins)" = "R",
            "Dorado polish (expérimental, nécessite dorado)" = "dorado")),
        conditionalPanel("input.polisher == 'dorado'",
          textInput("dorado_path", "Exécutable dorado", if (nzchar(DORADO)) DORADO else "dorado"),
          textInput("dorado_model", "Modèle de polissage (vide = automatique)", "")),
        hr(),
        fluidRow(column(6, num("max_reads", "Reads max. par échantillon", 5000, min = 100, step = 500)),
                 column(6, num("max_msa", "Reads max. par alignement", 150, min = 20, step = 10))),
        run_btn("run_consensus", "Construire les consensus")
      ),
      mainPanel(uiOutput("cons_info"), plotOutput("clus_plot", height = 320), DTOutput("cons_tbl"))
    )),

  # ------------------------------------------------------------------ Table
  tabPanel("5. Table",
    sidebarLayout(
      sidebarPanel(width = 4,
        num("merge_cutoff", "Fusionner entre échantillons les consensus distants de moins de", 0.01, min = 0, max = 0.05, step = 0.005),
        help("0 = ne fusionner que les consensus identiques ; 0,01 ≈ 4 différences sur 400 pb."),
        checkboxInput("remove_chim", "Retirer les chimères (DADA2)", TRUE),
        conditionalPanel("input.remove_chim",
          selectInput("chim_method", "Méthode", c("Toutes les séquences ensemble (pooled)" = "pooled",
                                                  "Échantillon par échantillon (consensus)" = "consensus")),
          checkboxInput("one_off", "Tolérer une différence avec les parents (conseillé en nanopore)", TRUE)),
        run_btn("run_table", "Construire la table")
      ),
      mainPanel(uiOutput("table_info"), plotOutput("table_plot", height = 300), DTOutput("table_tbl"))
    )),

  # ------------------------------------------------------------------ Taxonomie
  tabPanel("6. Taxonomie",
    sidebarLayout(
      sidebarPanel(width = 4,
        radioButtons("tax_method", "Classificateur", c("assignTaxonomy" = "assignTaxonomy", "IDTAXA (.RData)" = "idtaxa")),
        radioButtons("ref_src", "Base de référence", c("Prédéfinie" = "preset", "Fichier local" = "path")),
        conditionalPanel("input.ref_src == 'preset'", selectInput("ref_preset", NULL, REF_CHOICES)),
        conditionalPanel("input.ref_src == 'path'", textInput("ref_path", "Chemin", "")),
        textInput("taxLevels", "Rangs", REFS$silva$levels),
        num("minBoot", "minBoot / confiance", 50, min = 0, max = 100),
        checkboxInput("tryRC", "Tester le rev. comp.", TRUE),
        help("Les bases perso se préparent avec l'outil « base de référence » de l'appli DADA2."),
        run_btn("run_taxonomy", "Assigner"),
        hr(), h5("Filtre taxonomique (optionnel)"),
        textInput("tf_remove", "Retirer (virgules)", ""),
        textInput("tf_keep", "Ne garder que (virgules)", ""),
        selectInput("tf_unassigned", "Retirer les non assignés au rang", c("(ne pas filtrer)" = "")),
        run_btn("run_taxfilter", "Appliquer le filtre")
      ),
      mainPanel(uiOutput("tax_info"), DTOutput("tax_tbl"))
    )),

  # ------------------------------------------------------------------ Mock
  tabPanel("Communauté mock",
    sidebarLayout(
      sidebarPanel(width = 3,
        selectizeInput("mock_samples", "Échantillon(s) mock", choices = NULL, multiple = TRUE),
        fileInput("mock_fasta", "Séquences attendues (FASTA)"),
        num("mock_maxdiff", "Différences max. pour une variante d'erreur", 3, min = 0),
        actionButton("run_mock", "Évaluer", class = "btn-primary")),
      mainPanel(fluidRow(column(5, tableOutput("mock_summary")), column(7, tableOutput("mock_perref"))),
                DTOutput("mock_asv")))),

  # ------------------------------------------------------------------ Exploration
  tabPanel("Exploration",
    tabsetPanel(
      tabPanel("Abondances relatives", br(),
        fluidRow(
          column(3, selectInput("ab_level", "Niveau", c("Séquence (sans taxonomie)" = "ASV"))),
          column(2, num("ab_top", "Nombre de taxons", 15, min = 1)),
          column(2, radioButtons("ab_type", NULL, c("Barres" = "bar", "Heatmap" = "heatmap"))),
          column(3, selectInput("ab_group", "Moyenner par", c("(non)" = ""))),
          column(2, br(), downloadButton("dl_ab", "PDF", class = "btn-sm"))),
        plotOutput("ab_plot", height = "auto")),
      tabPanel("Krona", br(),
        fluidRow(column(4, radioButtons("kr_mode", NULL, c("Un jeu par échantillon" = "sample", "Tout cumulé" = "all"), inline = TRUE)),
                 column(4, downloadButton("dl_krona", "Télécharger (HTML)", class = "btn-sm"))),
        help("Le moteur Krona est chargé depuis internet."),
        uiOutput("krona_frame")))),

  # ------------------------------------------------------------------ Export
  tabPanel("Export",
    h4("Résultats"),
    p(downloadButton("dl_table", "Table + taxonomie (CSV)")),
    p(downloadButton("dl_fasta", "Séquences (FASTA)")),
    p(downloadButton("dl_cons", "Détail des consensus par échantillon (CSV)")),
    p(downloadButton("dl_ps", "Objet phyloseq (RDS)")),
    h4("Reproductibilité"),
    p(downloadButton("dl_script", "Script R reproductible")),
    p(downloadButton("dl_session", "Session (RDS)")),
    help("La session et le journal sont aussi enregistrés dans le dossier de sortie après chaque étape.")),

  tabPanel("Journal", verbatimTextOutput("log"))
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {
  state <- reactiveVal(list(params = list()))
  busy <- reactiveVal(NULL); t_start <- reactiveVal(NULL); chain <- reactiveVal(NULL)
  notify <- function(m, type = "message", duration = 6) showNotification(m, type = type, duration = duration)
  done <- function(s, st = state()) !is.null(st$params[[s]])
  outdir <- function() state()$params$data$outdir
  threads <- function() state()$params$data$threads %||% 1L
  log_file <- reactive({ d <- outdir(); if (is.null(d)) NULL else file.path(d, "nanopore_log.txt") })
  log_line <- function(...) { lf <- isolate(log_file()); if (!is.null(lf))
    cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", file = lf, append = TRUE, sep = "") }
  reset_from <- function(st, step) {
    for (s in STEPS[match(step, STEPS):length(STEPS)]) {
      for (f in FIELDS[[s]]) st[f] <- list(NULL)
      st$params[s] <- list(NULL)
    }
    st
  }
  save_session <- function(st = state()) {
    d <- st$params$data$outdir
    if (!is.null(d)) try(saveRDS(st, file.path(d, "nanopore_session.rds")), silent = TRUE)
  }

  launch <- function(step, fun, args, params, store, then = NULL) {
    if (!is.null(busy())) { notify("Une étape est déjà en cours.", "warning"); return(FALSE) }
    st <- reset_from(state(), step); state(st)
    busy(step); t_start(Sys.time())
    lf <- isolate(log_file())
    cat(sprintf("\n%s ===== %s =====\n", format(Sys.time(), "[%H:%M:%S]"), LABELS[step]), file = lf, append = TRUE)
    ok <- function(res) {
      st <- store(state(), res); st$params[[step]] <- params; state(st)
      busy(NULL); save_session(st)
      log_line("✔ ", LABELS[step], " terminé")
      notify(paste(LABELS[step], ": terminé"))
      if (!is.null(then)) then()
    }
    ko <- function(e) {
      busy(NULL); chain(NULL)
      log_line("✖ ", LABELS[step], " : ", conditionMessage(e))
      showNotification(paste(LABELS[step], "—", conditionMessage(e)), type = "error", duration = NULL)
    }
    if (ASYNC) {
      promises::then(promises::future_promise(run_job(fun, args, lf, PIPELINE_FILE), seed = TRUE),
                     onFulfilled = ok, onRejected = ko)
    } else {
      res <- withProgress(message = LABELS[step], value = 0.5,
                          tryCatch(run_job(fun, args, lf, PIPELINE_FILE), error = function(e) e))
      if (inherits(res, "error")) ko(res) else ok(res)
    }
    TRUE
  }
  need_done <- function(prev) {
    if (!done(prev)) { notify(paste0("Lance d'abord l'étape « ", LABELS[prev], " »."), "warning"); return(FALSE) }
    TRUE
  }

  output$status <- renderUI({
    st <- state(); b <- busy()
    if (!is.null(b)) invalidateLater(1000)
    div(class = "badges", lapply(STEPS, function(s) {
      if (identical(b, s)) {
        el <- as.integer(difftime(Sys.time(), t_start(), units = "secs"))
        span(class = "badge b-run", sprintf("⏳ %s %02d:%02d", LABELS[s], el %/% 60, el %% 60))
      } else span(class = paste("badge", if (done(s, st)) "b-done" else "b-todo"),
                  paste(if (done(s, st)) "✔" else "", LABELS[s]))
    }))
  })
  log_txt <- reactivePoll(1000, session,
    function() { lf <- log_file(); if (is.null(lf) || !file.exists(lf)) "" else file.info(lf)$size },
    function() { lf <- log_file(); if (is.null(lf) || !file.exists(lf)) "Aucun journal." else
      paste(utils::tail(readLines(lf, warn = FALSE), 800), collapse = "\n") })
  output$log <- renderText(log_txt())

  # ------------------------------------------------------------------ Données
  observeEvent(input$scan, {
    if (!is.null(busy())) return(notify("Une étape est en cours.", "warning"))
    s <- tryCatch(nano_detect(input$fq_dir, recursive = input$recursive),
                  error = function(e) { notify(conditionMessage(e), "error", NULL); NULL })
    if (is.null(s)) return()
    od <- normalizePath(input$outdir, mustWork = FALSE)
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    st <- reset_from(state(), "data")
    st$samples <- s
    st$params$data <- list(outdir = od, threads = as.integer(input$threads), fq_dir = input$fq_dir)
    state(st); save_session(st)
    cat(sprintf("%s ===== Nouveau projet : %d échantillons =====\n", format(Sys.time(), "[%H:%M:%S]"), nrow(s)),
        file = file.path(od, "nanopore_log.txt"), append = TRUE)
    updateSelectizeInput(session, "mock_samples", choices = s$sample)
    notify(sprintf("%d échantillons détectés.", nrow(s)))
  })
  observeEvent(input$load_outdir, {
    f <- file.path(input$outdir, "nanopore_session.rds")
    if (!file.exists(f)) return(notify("Pas de session dans ce dossier.", "error"))
    st <- readRDS(f); state(st)
    if (!is.null(st$samples)) updateSelectizeInput(session, "mock_samples", choices = st$samples$sample)
    notify("Session rechargée.")
  })
  observeEvent(input$meta, {
    m <- tryCatch(read.csv(input$meta$datapath, row.names = 1, check.names = FALSE,
                           sep = if (grepl("\\.(tsv|txt)$", input$meta$name)) "\t" else ","),
                  error = function(e) NULL)
    if (is.null(m)) return(notify("Métadonnées illisibles.", "error"))
    st <- state(); st$meta <- m; state(st); save_session(st)
    updateSelectInput(session, "ab_group", choices = c("(non)" = "", colnames(m)))
  })
  output$data_info <- renderUI({
    s <- state()$samples
    validate(need(!is.null(s), "Indique le dossier des FASTQ puis « Détecter les échantillons »."))
    p(sprintf("%d échantillons, %d fichiers FASTQ.", nrow(s), sum(s$n_files)))
  })
  output$samples_tbl <- renderDT({
    s <- state()$samples; req(s)
    datatable(data.frame(Echantillon = s$sample, Fichiers = s$n_files,
                         `Mo` = round(vapply(s$file, function(x) sum(file.size(.files_of(x))), 0) / 1e6, 1),
                         check.names = FALSE), rownames = FALSE, options = list(pageLength = 25))
  })

  # ------------------------------------------------------------------ QC
  observeEvent(input$run_qc, {
    s <- state()$samples
    if (is.null(s)) return(notify("Détecte d'abord les échantillons.", "warning"))
    qc <- withProgress(message = "Qualité des reads", value = 0.5, nano_qc(s, input$qc_n))
    st <- state(); st$qc <- qc; state(st)
  })
  output$qc_len <- renderPlot({
    qc <- state()$qc; validate(need(!is.null(qc), "Clique sur « Calculer »."))
    ggplot(qc, aes(length)) + geom_histogram(bins = 80, fill = "#4E79A7") +
      facet_wrap(~sample, scales = "free_y") + labs(x = "Longueur (pb)", y = "Reads", title = "Longueur des reads") +
      geom_vline(xintercept = c(input$min_len, input$max_len), linetype = 2, colour = "red")
  })
  output$qc_q <- renderPlot({
    qc <- state()$qc; req(qc)
    ggplot(qc, aes(meanq)) + geom_histogram(bins = 60, fill = "#F28E2B") + facet_wrap(~sample, scales = "free_y") +
      geom_vline(xintercept = input$min_q, linetype = 2, colour = "red") +
      labs(x = "Qualité moyenne (Phred)", y = "Reads", title = "Qualité moyenne des reads")
  })
  output$qc_tbl <- renderDT({
    qc <- state()$qc; req(qc)
    n50 <- function(x) { x <- sort(x, decreasing = TRUE); x[which(cumsum(x) >= sum(x) / 2)[1]] }
    d <- do.call(rbind, lapply(split(qc, qc$sample), function(x) data.frame(
      Echantillon = x$sample[1], Reads_lus = nrow(x), Longueur_mediane = stats::median(x$length),
      N50 = n50(x$length), Q_mediane = stats::median(x$meanq),
      `% Q>=min` = round(100 * mean(x$meanq >= input$min_q), 1), check.names = FALSE)))
    datatable(d, rownames = FALSE, options = list(dom = "t", pageLength = 100))
  })

  # ------------------------------------------------------------------ Filtre & amorces
  observeEvent(input$apply_marker, {
    m <- MARKERS[[input$marker]]; if (is.null(m$primersF)) return()
    updateTextAreaInput(session, "primersF", value = m$primersF)
    updateTextAreaInput(session, "primersR", value = m$primersR)
    if (m$ref %in% names(REFS)) { updateRadioButtons(session, "ref_src", selected = "preset")
                                  updateSelectInput(session, "ref_preset", selected = m$ref) }
    else updateRadioButtons(session, "ref_src", selected = "path")
    updateTextInput(session, "taxLevels", value = m$levels)
    if (!is.na(m$amplen %||% NA)) {
      updateNumericInput(session, "min_len", value = round(m$amplen * 0.8))
      updateNumericInput(session, "max_len", value = round(m$amplen * 1.6 + 150))
    }
    notify(paste("Préréglage :", m$label))
  })
  RUN <- new.env()
  RUN$prep <- function(then = NULL) {
    if (!need_done("data")) return(FALSE)
    p <- list(min_len = input$min_len, max_len = input$max_len, min_q = input$min_q,
              primersF = input$primersF, primersR = input$primersR, primer_max_err = input$primer_err,
              window = input$window, require_both = input$require_both)
    if (nzchar(trimws(p$primersF)) != nzchar(trimws(p$primersR)))
      { notify("Indique les deux amorces, ou aucune.", "error"); return(FALSE) }
    launch("prep", "nano_prepare_all", list(samples = state()$samples, outdir = outdir(), p = p), p,
           function(st, res) { st$prep <- res; st }, then)
  }
  observeEvent(input$run_prep, RUN$prep())
  output$prep_info <- renderUI({
    pr <- state()$prep; validate(need(!is.null(pr), "Règle les filtres et les amorces, puis lance l'étape."))
    p(sprintf("%s reads → %s après filtre → %s orientés et découpés (%.1f %%).",
              format(sum(pr$input), big.mark = " "), format(sum(pr$filtered), big.mark = " "),
              format(sum(pr$oriented), big.mark = " "), 100 * sum(pr$oriented) / max(1, sum(pr$input))))
  })
  output$prep_tbl <- renderDT({
    pr <- state()$prep; req(pr)
    pr$pct <- round(100 * pr$oriented / pmax(pr$input, 1), 1)
    datatable(pr[, c("sample", "input", "filtered", "oriented", "pct")], rownames = FALSE, options = list(pageLength = 25))
  })
  output$prep_len <- renderPlot({
    pr <- state()$prep; req(pr)
    L <- do.call(rbind, lapply(seq_len(nrow(pr)), function(i) {
      if (!file.exists(pr$prep_file[i])) return(NULL)
      s <- ShortRead::FastqSampler(pr$prep_file[i], 2000); on.exit(close(s))
      data.frame(sample = pr$sample[i], length = Biostrings::width(ShortRead::yield(s)))
    }))
    req(L)
    ggplot(L, aes(length, fill = sample)) + geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
      labs(x = "Longueur après découpe des amorces (pb)", y = "Reads", fill = NULL)
  })

  # ------------------------------------------------------------------ Consensus
  RUN$consensus <- function(then = NULL) {
    if (!need_done("prep")) return(FALSE)
    if (input$polisher == "dorado" && !nzchar(Sys.which(input$dorado_path)) && !file.exists(input$dorado_path))
      { notify("dorado introuvable : choisis le polisseur R ou indique le chemin de dorado.", "error"); return(FALSE) }
    if (input$polisher == "dorado" && !HAS("Rsamtools"))
      { notify("Dorado polish nécessite le package Rsamtools.", "error"); return(FALSE) }
    p <- list(mode = input$mode, cluster_cutoff = input$cluster_cutoff, min_reads = input$min_reads,
              min_frac = input$min_frac, max_reads = input$max_reads, max_msa = input$max_msa,
              split = input$split, split_freq = input$split_freq, split_min_sites = input$split_sites,
              hp_correct = input$hp_correct, polisher = input$polisher, dorado = input$dorado_path,
              dorado_model = input$dorado_model, threads = threads(), seed = 1)
    launch("consensus", "nano_consensus_all", list(prep = state()$prep, p = p), p,
           function(st, res) { st$cons <- res$consensus; st$clusters <- res$clusters; st }, then)
  }
  observeEvent(input$run_consensus, RUN$consensus())
  output$cons_info <- renderUI({
    c <- state()$cons; validate(need(!is.null(state()$params$consensus), "Lance la construction des consensus."))
    if (is.null(c)) return(p(class = "text-danger", "Aucun consensus : baisse le nombre minimal de reads ou vérifie l'étape précédente."))
    p(sprintf("%d consensus sur %d échantillons ; %s reads représentés.", nrow(c), length(unique(c$sample)),
              format(sum(c$reads), big.mark = " ")))
  })
  output$clus_plot <- renderPlot({
    cl <- state()$clusters; req(cl)
    cl <- do.call(rbind, lapply(split(cl, cl$sample), function(x) { x$rang <- seq_len(nrow(x)); x }))
    ggplot(cl[cl$rang <= 30, ], aes(rang, size, fill = kept)) + geom_col() + facet_wrap(~sample, scales = "free_y") +
      scale_fill_manual(values = c(`TRUE` = "#59A14F", `FALSE` = "#BAB0AC"), labels = c(`TRUE` = "gardé", `FALSE` = "écarté"), name = NULL) +
      labs(x = "Groupes de reads (du plus au moins fourni)", y = "Reads", title = "Taille des groupes (30 premiers)")
  })
  output$cons_tbl <- renderDT({
    c <- state()$cons; req(c)
    c$sequence <- paste0(substr(c$sequence, 1, 30), "…")
    datatable(c, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })

  # ------------------------------------------------------------------ Table
  RUN$table <- function(then = NULL) {
    if (!need_done("consensus")) return(FALSE)
    if (is.null(state()$cons)) { notify("Aucun consensus.", "error"); return(FALSE) }
    p <- list(merge_cutoff = input$merge_cutoff, remove_chimeras = input$remove_chim,
              chim_method = input$chim_method, allow_one_off = input$one_off)
    launch("table", "nano_build_table",
           list(cons = state()$cons, merge_cutoff = p$merge_cutoff, remove_chimeras = p$remove_chimeras,
                threads = threads(), chim_method = p$chim_method, allow_one_off = p$allow_one_off), p,
           function(st, res) { st$tab <- res$table; st$cons_tab <- res$consensus; st$chimeras <- res$chimeras; st }, then)
  }
  observeEvent(input$run_table, RUN$table())
  final_tab <- reactive({ st <- state(); st$tab_tax %||% st$tab })
  seq_ids <- reactive(paste0("SEQ", match(colnames(final_tab()), colnames(state()$tab))))
  output$table_info <- renderUI({
    t <- state()$tab; validate(need(!is.null(t), "Construis la table."))
    p(sprintf("%d échantillons × %d séquences, %s reads ; %d chimère(s) retirée(s).", nrow(t), ncol(t),
              format(sum(t), big.mark = " "), length(state()$chimeras)))
  })
  output$table_plot <- renderPlot({
    t <- state()$tab; req(t)
    df <- data.frame(length = nchar(colnames(t)), reads = colSums(t))
    ggplot(df, aes(length, reads)) + geom_col(fill = "#4E79A7") + labs(x = "Longueur (pb)", y = "Reads", title = "Longueur des séquences")
  })
  output$table_tbl <- renderDT({
    t <- final_tab(); req(t)
    datatable(data.frame(ID = seq_ids(), longueur = nchar(colnames(t)), t(t), check.names = FALSE, row.names = NULL),
              rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })
  observe({ t <- final_tab(); if (!is.null(t)) updateSelectizeInput(session, "mock_samples", choices = rownames(t),
                                                                    selected = isolate(input$mock_samples)) })

  # ------------------------------------------------------------------ Taxonomie
  observeEvent(input$ref_preset, updateTextInput(session, "taxLevels", value = REFS[[input$ref_preset]]$levels), ignoreInit = TRUE)
  tax_configured <- function() input$ref_src == "preset" || nzchar(input$ref_path)
  RUN$taxonomy <- function(then = NULL) {
    if (!need_done("table")) return(FALSE)
    if (!tax_configured()) { notify("Indique une base de référence.", "warning"); return(FALSE) }
    p <- list(method = input$tax_method,
              ref = if (input$ref_src == "preset") REFS[[input$ref_preset]]$url else input$ref_path,
              levels = input$taxLevels, minBoot = input$minBoot, tryRC = input$tryRC)
    launch("taxonomy", "run_taxonomy", list(seqs = colnames(state()$tab), p = p, threads = threads()), p,
           function(st, res) { st$taxa <- res$taxa; st$boot <- res$boot; st }, then)
  }
  observeEvent(input$run_taxonomy, RUN$taxonomy())
  observe({
    tx <- state()$taxa; r <- colnames(tx %||% matrix(nrow = 0, ncol = 0))
    updateSelectInput(session, "tf_unassigned", choices = c("(ne pas filtrer)" = "", r))
    updateSelectInput(session, "ab_level", choices = c("Séquence (sans taxonomie)" = "ASV", r),
                      selected = if ("Family" %in% r) "Family" else "ASV")
  })
  RUN$taxfilter <- function(then = NULL) {
    if (!need_done("taxonomy")) return(FALSE)
    p <- list(remove = input$tf_remove, keep = input$tf_keep, unassigned_rank = input$tf_unassigned %||% "")
    launch("taxfilter", "run_taxfilter", list(final = state()$tab, taxa = state()$taxa, p = p), p,
           function(st, res) { st$tab_tax <- res$final_tax; st$tf_report <- res$report; st }, then)
  }
  observeEvent(input$run_taxfilter, RUN$taxfilter())
  output$tax_info <- renderUI({
    tx <- state()$taxa; validate(need(!is.null(tx), "Choisis une base et lance l'assignation."))
    pct <- round(100 * colMeans(!is.na(tx[colnames(final_tab()), , drop = FALSE])), 1)
    tagList(p(tags$b("% assignées : "), paste(colnames(tx), pct, "%", collapse = " · ")),
            if (!is.null(state()$tab_tax)) p(sprintf("Filtre taxonomique : %d / %d séquences gardées.",
                                                       ncol(state()$tab_tax), ncol(state()$tab))))
  })
  output$tax_tbl <- renderDT({
    tx <- state()$taxa; req(tx); t <- final_tab()
    datatable(data.frame(ID = seq_ids(), reads = colSums(t), tx[colnames(t), , drop = FALSE], check.names = FALSE, row.names = NULL),
              rownames = FALSE, filter = "top", options = list(pageLength = 20, scrollX = TRUE))
  })

  # ------------------------------------------------------------------ Tout lancer
  chain_next <- function() {
    ch <- chain(); if (is.null(ch)) return()
    todo <- ch[!vapply(ch, done, TRUE)]
    if (!length(todo)) { chain(NULL); return(notify("Tout lancer : terminé.", duration = NULL)) }
    s <- todo[1]
    if (s == "taxonomy" && !tax_configured()) { chain(NULL); return(notify("Taxonomie sautée : pas de référence.", "warning")) }
    if (s == "taxfilter") { chain(NULL); return(notify("Tout lancer : terminé (le filtre taxonomique est à lancer à la main si besoin).", duration = NULL)) }
    if (!isTRUE(RUN[[s]](then = chain_next))) chain(NULL)
  }
  observeEvent(input$run_all, {
    if (!done("data")) return(notify("Détecte d'abord les échantillons.", "warning"))
    chain(c("prep", "consensus", "table", "taxonomy", "taxfilter")); chain_next()
  })

  # ------------------------------------------------------------------ Mock
  mock_res <- reactiveVal(NULL)
  observeEvent(input$run_mock, {
    t <- state()$tab
    if (is.null(t) || !length(input$mock_samples) || is.null(input$mock_fasta))
      return(notify("Il faut la table, le ou les échantillons mock et le FASTA attendu.", "warning"))
    res <- tryCatch(withProgress(message = "Évaluation du mock", value = 0.5,
                                 mock_eval(t, input$mock_samples, read_fasta_simple(input$mock_fasta$datapath),
                                           paste0("SEQ", seq_len(ncol(t))), input$mock_maxdiff)),
                    error = function(e) { notify(conditionMessage(e), "error"); NULL })
    mock_res(res)
  })
  output$mock_summary <- renderTable({ r <- mock_res(); validate(need(!is.null(r), "Charge les séquences attendues puis « Évaluer »."));  r$summary })
  output$mock_perref <- renderTable(req(mock_res())$per_ref)
  output$mock_asv <- renderDT({ a <- req(mock_res())$asv; a$sequence <- paste0(substr(a$sequence, 1, 30), "…")
    datatable(a, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE)) })

  # ------------------------------------------------------------------ Exploration
  ab <- reactive({
    t <- final_tab(); validate(need(!is.null(t), "Construis d'abord la table."))
    tx <- state()$taxa; if (!is.null(tx)) tx <- tx[colnames(t), , drop = FALSE]
    lvl <- input$ab_level %||% "ASV"
    validate(need(lvl == "ASV" || !is.null(tx), "Assigne d'abord la taxonomie."))
    g <- input$ab_group %||% ""
    groups <- if (nzchar(g) && !is.null(state()$meta)) state()$meta[rownames(t), g] else NULL
    mat <- abund_matrix(t, tx, seq_ids(), lvl, if (lvl == "ASV" && !is.null(tx) && "Genus" %in% colnames(tx)) "Genus", groups, input$ab_top)
    list(mat = mat, title = paste("Abondance relative —", if (lvl == "ASV") "séquences" else lvl))
  })
  ab_plot <- reactive(plot_abund(ab()$mat, colnames(ab()$mat), input$ab_type, NULL, ab()$title))
  output$ab_plot <- renderPlot(ab_plot(), height = function() if (identical(input$ab_type, "heatmap")) 450 else 520)
  output$dl_ab <- downloadHandler("abondances.pdf", function(f) ggsave(f, ab_plot(), width = 10, height = 7))
  krona_dir <- file.path(tempdir(), paste0("kronanano_", session$token)); dir.create(krona_dir, showWarnings = FALSE)
  addResourcePath(basename(krona_dir), krona_dir)
  krona_lines <- reactive({
    t <- final_tab(); validate(need(!is.null(t), "Construis d'abord la table."))
    tx <- state()$taxa; if (!is.null(tx)) tx <- tx[colnames(t), , drop = FALSE]
    krona_html(krona_datasets(t, input$kr_mode), tx, seq_ids(), add_asv = TRUE)
  })
  output$krona_frame <- renderUI({
    writeLines(enc2utf8(krona_lines()), file.path(krona_dir, "krona.html"), useBytes = TRUE)
    tags$iframe(src = sprintf("%s/krona.html?v=%s", basename(krona_dir), as.numeric(Sys.time())),
                style = "width:100%; height:760px; border:1px solid #ddd;")
  })
  output$dl_krona <- downloadHandler("krona.html", function(f) writeLines(enc2utf8(krona_lines()), f, useBytes = TRUE))

  # ------------------------------------------------------------------ Export
  out_table <- reactive({
    t <- final_tab(); req(t); tx <- state()$taxa
    df <- data.frame(ID = seq_ids(), sequence = colnames(t), longueur = nchar(colnames(t)), check.names = FALSE)
    if (!is.null(tx)) df <- cbind(df, tx[colnames(t), , drop = FALSE])
    cbind(df, t(t), row.names = NULL)
  })
  output$dl_table <- downloadHandler("table_nanopore.csv", function(f) write.csv(out_table(), f, row.names = FALSE))
  output$dl_fasta <- downloadHandler("sequences_nanopore.fasta", function(f) {
    t <- final_tab(); req(t); writeLines(paste0(">", seq_ids(), "\n", colnames(t)), f) })
  output$dl_cons <- downloadHandler("consensus_par_echantillon.csv", function(f) {
    c <- state()$cons_tab %||% state()$cons; req(c); write.csv(c, f, row.names = FALSE) })
  output$dl_ps <- downloadHandler("phyloseq_nanopore.rds", function(f) {
    req(HAS("phyloseq")); t <- final_tab(); req(t)
    parts <- list(phyloseq::otu_table(t, taxa_are_rows = FALSE))
    m <- state()$meta
    parts <- c(parts, list(phyloseq::sample_data(if (!is.null(m)) m[rownames(t), , drop = FALSE]
                                                 else data.frame(Echantillon = rownames(t), row.names = rownames(t)))))
    if (!is.null(state()$taxa)) parts <- c(parts, list(phyloseq::tax_table(state()$taxa[colnames(t), , drop = FALSE])))
    ps <- do.call(phyloseq::phyloseq, parts)
    phyloseq::taxa_names(ps) <- seq_ids()
    saveRDS(ps, f)
  })
  output$dl_script <- downloadHandler("pipeline_nanopore.R", function(f) writeLines(nano_script(state(), REPO), f))
  output$dl_session <- downloadHandler("nanopore_session.rds", function(f) saveRDS(state(), f))
}

shinyApp(ui, server)

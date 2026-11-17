# =============================================================================
#  pipeline.R — fonctions du pipeline DADA2, indépendantes de Shiny
#
#  Chaque étape est une fonction pure : elle prend des données + une liste de
#  paramètres et renvoie une liste de résultats. L'avancement est signalé par
#  message(). On peut donc les utiliser hors de l'appli (script, cluster…).
# =============================================================================

suppressPackageStartupMessages(library(dada2))

`%||%` <- function(a, b) if (is.null(a)) b else a

is_paired <- function(samples) !all(is.na(samples$fnR))

# dada()/mergePairs() renvoient un objet nu (pas une liste) pour un seul échantillon
as_named_list <- function(x, nm) {
  if (inherits(x, "dada") || is.data.frame(x)) setNames(list(x), nm) else x
}

getN <- function(x) sum(getUniques(x))

# -----------------------------------------------------------------------------
# Exécution d'une étape avec journal en direct dans un fichier
# (utilisé par l'appli, en processus séparé ou non)
# -----------------------------------------------------------------------------
run_job <- function(fun, args, log_file, pipeline_file) {
  env <- new.env()
  sys.source(pipeline_file, envir = env)
  blast_file <- file.path(dirname(pipeline_file), "blast.R")
  if (file.exists(blast_file)) sys.source(blast_file, envir = env)
  f <- get(fun, envir = env)

  con <- file(log_file, open = "a")
  stamp <- function() format(Sys.time(), "[%H:%M:%S] ")
  sink(con, type = "output")
  on.exit({ sink(type = "output"); close(con) }, add = TRUE)

  withCallingHandlers(
    tryCatch(do.call(f, args), error = function(e) {
      cat(stamp(), "ERREUR : ", conditionMessage(e), "\n", sep = "")
      stop(e)
    }),
    message = function(m) {
      cat(conditionMessage(m)); flush(con)
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      cat(stamp(), "Avertissement : ", conditionMessage(w), "\n", sep = ""); flush(con)
      invokeRestart("muffleWarning")
    })
}

# -----------------------------------------------------------------------------
# 1. Détection des échantillons
# -----------------------------------------------------------------------------
#' @param multirun  si TRUE, chaque sous-dossier de `dir` est un run de séquençage
#' @param patR      motif reverse ; NULL ou "" = single-end
detect_samples <- function(dir, patF, patR = NULL, sep = "_", idx = 1, multirun = FALSE) {
  if (!dir.exists(dir)) stop("Dossier introuvable : ", dir)
  paired <- !is.null(patR) && nzchar(patR)

  name_of <- function(f) vapply(strsplit(basename(f), sep, fixed = TRUE),
                                function(p) if (length(p) >= idx) p[idx] else NA_character_, "")
  run_of <- function(f) if (multirun) basename(dirname(f)) else rep("run1", length(f))

  fnF <- sort(list.files(dir, pattern = patF, full.names = TRUE, recursive = multirun))
  if (!length(fnF)) stop("Aucun fichier ne correspond au motif forward « ", patF, " ».")
  df <- data.frame(sample = name_of(fnF), run = run_of(fnF), fnF = fnF, stringsAsFactors = FALSE)
  if (anyNA(df$sample))
    stop("Impossible d'extraire le nom d'échantillon de : ",
         paste(basename(fnF[is.na(df$sample)]), collapse = ", "))
  keyF <- paste(df$run, df$sample, sep = "/")
  if (anyDuplicated(keyF))
    stop("Noms d'échantillons dupliqués (ajuste séparateur/position) : ",
         paste(unique(keyF[duplicated(keyF)]), collapse = ", "))

  if (paired) {
    fnR <- sort(list.files(dir, pattern = patR, full.names = TRUE, recursive = multirun))
    if (!length(fnR)) stop("Aucun fichier ne correspond au motif reverse « ", patR, " ».")
    if (length(intersect(fnF, fnR))) stop("Les motifs forward et reverse sélectionnent les mêmes fichiers.")
    keyR <- paste(run_of(fnR), name_of(fnR), sep = "/")
    noR <- setdiff(keyF, keyR); noF <- setdiff(keyR, keyF)
    if (length(noR) || length(noF))
      stop("Paires incomplètes.",
           if (length(noR)) paste0(" Sans reverse : ", paste(noR, collapse = ", "), "."),
           if (length(noF)) paste0(" Sans forward : ", paste(noF, collapse = ", "), "."))
    df$fnR <- fnR[match(keyF, keyR)]
  } else df$fnR <- NA_character_

  if (anyDuplicated(df$sample)) {
    message("Mêmes noms d'échantillons dans plusieurs runs : noms préfixés par le run.")
    df$sample <- paste(df$run, df$sample, sep = "__")
  }
  df
}

# Longueur des reads sur les premières lectures de quelques fichiers
read_length_summary <- function(files, n = 2000) {
  files <- files[!is.na(files)]
  w <- unlist(lapply(files, function(f) {
    fs <- ShortRead::FastqStreamer(f, n = n)
    on.exit(close(fs))
    Biostrings::width(ShortRead::sread(ShortRead::yield(fs)))
  }))
  if (!length(w)) return(c(median = NA, max = NA))
  c(median = stats::median(w), max = max(w))
}

# -----------------------------------------------------------------------------
# 2. Amorces
# -----------------------------------------------------------------------------
parse_primers <- function(x) {
  p <- toupper(trimws(unlist(strsplit(x %||% "", "[,;[:space:]]+"))))
  p <- gsub("I", "N", p)                      # inosine -> N
  p <- p[nzchar(p)]
  bad <- grepl("[^ACGTRYSWKMBDHVN]", p)
  if (any(bad)) stop("Caractères non IUPAC dans les amorces : ", paste(p[bad], collapse = ", "))
  p
}
rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))

#' Devine le sens d'une amorce d'après son nom : « F », « R » ou NA.
#' Reconnaît fwd/forward/rev/reverse, _F/_R, 708F, R3, 515F-Y…
primer_direction <- function(name) {
  tok <- tolower(unlist(strsplit(name, "[^A-Za-z0-9]+")))
  tok <- tok[nzchar(tok)]
  has_num <- grepl("[0-9]", tok)
  isF <- any(tok %in% c("f", "fw", "fwd", "forward") | (has_num & grepl("f[0-9]*$", tok)) |
               grepl("^f[0-9]+$", tok) | grepl("^[0-9]+f[a-z]$", tok))     # ex. 515FB
  isR <- any(tok %in% c("r", "rv", "rev", "reverse") | (has_num & grepl("r[0-9]*$", tok)) |
               grepl("^r[0-9]+$", tok) | grepl("^[0-9]+r[a-z]$", tok))     # ex. 806RB
  if (isF && !isR) "F" else if (isR && !isF) "R" else NA_character_
}

#' Lit un fichier d'amorces.
#'  - FASTA : le sens est déduit du nom (ex. >rbcL_708F_1, >R3_2, >COI_fwd)
#'  - CSV/TSV : une colonne de séquences (nom contenant « seq ») et, de préférence,
#'    une colonne « direction » (F/R, forward/reverse) ; sinon le sens est
#'    déduit d'une colonne de noms.
#' @return list(F = séquences forward, R = séquences reverse, table = détail)
read_primer_file <- function(path, filename = basename(path)) {
  first <- readLines(path, n = 1, warn = FALSE)
  if (length(first) && grepl("^>", first)) {
    l <- trimws(readLines(path, warn = FALSE)); l <- l[nzchar(l)]
    h <- grep("^>", l)
    nm <- trimws(sub("^>", "", l[h]))
    ends <- c(h[-1] - 1, length(l))
    sq <- vapply(seq_along(h), function(i)
      if (ends[i] > h[i]) paste(l[(h[i] + 1):ends[i]], collapse = "") else "", "")
    dir <- vapply(nm, primer_direction, "", USE.NAMES = FALSE)
  } else {
    seps <- c("\t", ";", ",")
    counts <- vapply(seps, function(s) lengths(regmatches(first, gregexpr(s, first, fixed = TRUE))), 0)
    sep <- seps[which.max(counts)]
    df <- utils::read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE,
                            quote = "\"", comment.char = "", check.names = FALSE)
    cn <- tolower(names(df))
    sc <- which(grepl("seq", cn))[1]
    if (is.na(sc)) stop("Colonne de séquences introuvable : son nom doit contenir « seq ».")
    dc <- which(cn %in% c("direction", "dir", "sens", "orientation", "strand", "fr"))[1]
    nc <- setdiff(which(grepl("name|nom|id|primer|amorce", cn)), sc)[1]
    sq <- as.character(df[[sc]])
    nm <- if (!is.na(nc)) as.character(df[[nc]]) else paste0("amorce_", seq_along(sq))
    dir <- if (!is.na(dc)) {
      d <- toupper(substr(trimws(as.character(df[[dc]])), 1, 1))
      ifelse(d %in% c("F", "R"), d, NA_character_)
    } else vapply(nm, primer_direction, "", USE.NAMES = FALSE)
  }
  sq <- toupper(gsub("[[:space:]]", "", sq))
  if (!length(sq) || any(!nzchar(sq))) stop("Séquence vide ou fichier sans amorce.")
  if (anyNA(dir))
    stop("Sens (forward/reverse) non reconnu pour : ", paste(nm[is.na(dir)], collapse = ", "),
         ". Ajoute _F / _R au nom (FASTA) ou une colonne « direction » contenant F ou R (CSV).")
  parse_primers(paste(sq, collapse = ","))   # vérifie les caractères IUPAC
  list(F = sq[dir == "F"], R = sq[dir == "R"],
       table = data.frame(Nom = nm, Sens = dir, Sequence = sq, Longueur = nchar(sq)))
}

#' Proportion de reads contenant chaque amorce, dans chaque orientation
primer_check <- function(samples, primersF, primersR, n_samples = 3, n_reads = 2000,
                         which = c("raw", "trimmed")) {
  which <- match.arg(which)
  pF <- parse_primers(primersF); pR <- parse_primers(primersR)
  colF <- if (which == "trimmed") "trimF" else "fnF"
  colR <- if (which == "trimmed") "trimR" else "fnR"
  s <- utils::head(samples, n_samples)
  read_first <- function(f) {
    if (is.na(f) || !file.exists(f)) return(NULL)
    fs <- ShortRead::FastqStreamer(f, n = n_reads); on.exit(close(fs))
    ShortRead::sread(ShortRead::yield(fs))
  }
  readsF <- do.call(c, lapply(s[[colF]], read_first))
  readsR <- if (is_paired(samples)) do.call(c, lapply(s[[colR]], read_first)) else NULL
  pct <- function(primers, reads) {
    if (is.null(reads) || !length(primers)) return(NA_real_)
    hit <- Reduce(`|`, lapply(primers, function(p)
      Biostrings::vcountPattern(p, reads, fixed = FALSE) > 0))
    round(100 * mean(hit), 1)
  }
  orients <- list(
    "Amorce F (sens)" = pF, "Amorce F (rev. comp.)" = rc(pF),
    "Amorce R (sens)" = pR, "Amorce R (rev. comp.)" = rc(pR))
  orients <- orients[lengths(orients) > 0]
  data.frame(
    Recherche = names(orients),
    `% reads R1` = vapply(orients, pct, 0, reads = readsF),
    `% reads R2` = vapply(orients, pct, 0, reads = readsR),
    check.names = FALSE, row.names = NULL)
}

#' Retrait des amorces avec cutadapt (doit être installé)
run_cutadapt <- function(samples, outdir, p, threads = 1) {
  exe <- p$cutadapt %||% "cutadapt"
  ver <- tryCatch(suppressWarnings(system2(exe, "--version", stdout = TRUE, stderr = TRUE)),
                  error = function(e) NULL)
  if (is.null(ver) || !is.null(attr(ver, "status")))
    stop("cutadapt introuvable (« ", exe, " »). Installe-le (pip/conda) ou indique son chemin.")
  message("cutadapt version ", ver[1])

  pF <- parse_primers(p$primersF); pR <- parse_primers(p$primersR)
  if (!length(pF) && !length(pR)) stop("Aucune amorce renseignée.")
  paired <- is_paired(samples)
  flag <- function(f, x) if (length(x)) c(rbind(f, x)) else character(0)
  flags <- c(flag("-g", pF), flag("-a", rc(pR)))
  if (paired) flags <- c(flags, flag("-G", pR), flag("-A", rc(pF)))
  common <- c(flags, "-n", 2, "-m", p$cut_minlen %||% 20, "-j", threads,
              if (isTRUE(p$discard_untrimmed)) "--discard-untrimmed")

  d <- file.path(outdir, "cutadapt")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  samples$trimF <- file.path(d, paste0(samples$sample, "_R1_trim.fastq.gz"))
  samples$trimR <- if (paired) file.path(d, paste0(samples$sample, "_R2_trim.fastq.gz")) else NA_character_

  num <- function(lines, pattern) {
    l <- grep(pattern, lines, value = TRUE)[1]
    if (is.na(l)) return(NA_real_)
    as.numeric(gsub("[^0-9]", "", sub("^[^:]*:\\s*([0-9,]+).*$", "\\1", l)))
  }
  stats <- data.frame(sample = samples$sample, cut_in = NA_real_, cut_out = NA_real_)
  for (i in seq_len(nrow(samples))) {
    io <- if (paired)
      c("-o", shQuote(samples$trimF[i]), "-p", shQuote(samples$trimR[i]),
        shQuote(samples$fnF[i]), shQuote(samples$fnR[i]))
    else c("-o", shQuote(samples$trimF[i]), shQuote(samples$fnF[i]))
    res <- suppressWarnings(system2(exe, c(common, io), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(res, "status")) && attr(res, "status") != 0)
      stop("cutadapt a échoué pour ", samples$sample[i], " :\n", paste(utils::tail(res, 10), collapse = "\n"))
    stats$cut_in[i]  <- num(res, "^Total (read pairs|reads) processed")
    stats$cut_out[i] <- num(res, "^(Pairs|Reads) written")
    message(sprintf("%s : %s -> %s reads", samples$sample[i],
                    format(stats$cut_in[i], big.mark = " "), format(stats$cut_out[i], big.mark = " ")))
  }
  list(samples = samples, stats = stats)
}

# -----------------------------------------------------------------------------
# 3. Filtrage
# -----------------------------------------------------------------------------
run_filter <- function(samples, outdir, p, threads = 1) {
  paired <- is_paired(samples)
  inF <- if (!is.null(samples$trimF)) samples$trimF else samples$fnF
  inR <- if (!is.null(samples$trimR)) samples$trimR else samples$fnR

  d <- file.path(outdir, "filtered")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  filtF <- file.path(d, paste0(samples$sample, "_F_filt.fastq.gz"))
  filtR <- file.path(d, paste0(samples$sample, "_R_filt.fastq.gz"))
  unlink(c(filtF, filtR))
  mt <- if (.Platform$OS.type == "windows") FALSE else threads

  args <- list(truncQ = p$truncQ, minLen = p$minLen, maxN = 0, rm.phix = p$rm.phix,
               compress = TRUE, multithread = mt, verbose = TRUE)
  out <- if (paired) {
    do.call(filterAndTrim, c(list(inF, filtF, inR, filtR, truncLen = p$truncLen,
                                  trimLeft = p$trimLeft, maxEE = p$maxEE), args))
  } else {
    do.call(filterAndTrim, c(list(inF, filtF, truncLen = p$truncLen[1],
                                  trimLeft = p$trimLeft[1], maxEE = p$maxEE[1]), args))
  }
  out <- matrix(out, ncol = 2, dimnames = list(samples$sample, c("reads.in", "reads.out")))
  samples$filtF <- ifelse(file.exists(filtF), filtF, NA_character_)
  samples$filtR <- if (paired) ifelse(file.exists(filtR), filtR, NA_character_) else NA_character_
  empty <- is.na(samples$filtF) | (paired & is.na(samples$filtR))
  if (any(empty)) message("Échantillons sans lecture après filtrage (exclus de la suite) : ",
                          paste(samples$sample[empty], collapse = ", "))
  message(sprintf("Filtrage : %s / %s lectures conservées (%.1f %%)",
                  format(sum(out[, 2]), big.mark = " "), format(sum(out[, 1]), big.mark = " "),
                  100 * sum(out[, 2]) / max(1, sum(out[, 1]))))
  list(samples = samples, out = out)
}

# -----------------------------------------------------------------------------
# 4. Modèle d'erreur
# -----------------------------------------------------------------------------
#' Fonction d'estimation d'erreur pour qualités « binnées » (NovaSeq, NextSeq,
#' MiSeq récents) : loess pondéré + taux d'erreur forcés décroissants avec la
#' qualité. Adapté de la discussion benjjneb/dada2#1307.
binned_qual_errfun <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1) / tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
        mod.lo <- stats::loess(rlogp ~ q, df, weights = log10(tot + 1), span = 2)
        pred <- stats::predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred))); minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  est[est > 0.25] <- 0.25
  est[est < 1e-7] <- 1e-7
  # contrainte : le taux d'erreur ne peut pas augmenter avec la qualité
  est <- t(apply(est, 1, function(r) rev(cummax(rev(r)))))
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ],
               est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ],
               est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
               est[10:12, ], 1 - colSums(est[10:12, ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  err
}

usable <- function(samples) {
  ok <- !is.na(samples$filtF)
  if (is_paired(samples)) ok <- ok & !is.na(samples$filtR)
  ok
}

run_learn_errors <- function(samples, p, threads = 1) {
  fun <- if (identical(p$errfun, "binned")) binned_qual_errfun else dada2::loessErrfun
  ok <- usable(samples)
  runs <- unique(samples$run[ok])
  lapply(setNames(runs, runs), function(r) {
    s <- samples[ok & samples$run == r, ]
    message("== Run ", r, " (", nrow(s), " échantillons) — forward")
    e <- list(F = learnErrors(s$filtF, nbases = p$nbases, randomize = p$randomize,
                              errorEstimationFunction = fun, multithread = threads, verbose = TRUE))
    if (is_paired(samples)) {
      message("== Run ", r, " — reverse")
      e$R <- learnErrors(s$filtR, nbases = p$nbases, randomize = p$randomize,
                         errorEstimationFunction = fun, multithread = threads, verbose = TRUE)
    }
    e
  })
}

# -----------------------------------------------------------------------------
# 5. Débruitage + fusion (run par run, puis fusion des tables)
# -----------------------------------------------------------------------------
run_denoise <- function(samples, errs, p, threads = 1) {
  pool <- switch(as.character(p$pool), "TRUE" = TRUE, "pseudo" = "pseudo", FALSE)
  paired <- is_paired(samples)
  ok <- usable(samples)
  runs <- unique(samples$run[ok])
  counts <- list(); tables <- list()
  for (r in runs) {
    s <- samples[ok & samples$run == r, ]
    nm <- s$sample
    fF <- setNames(s$filtF, nm)
    message("== Run ", r, " : débruitage forward")
    dF <- as_named_list(dada(fF, err = errs[[r]]$F, pool = pool, multithread = threads), nm)
    cnt <- data.frame(sample = nm, denoisedF = vapply(dF, getN, 0), row.names = NULL)
    if (paired) {
      fR <- setNames(s$filtR, nm)
      message("== Run ", r, " : débruitage reverse")
      dR <- as_named_list(dada(fR, err = errs[[r]]$R, pool = pool, multithread = threads), nm)
      message("== Run ", r, " : fusion des paires")
      mg <- as_named_list(mergePairs(dF, fF, dR, fR, minOverlap = p$minOverlap,
                                     maxMismatch = p$maxMismatch, trimOverhang = p$trimOverhang,
                                     justConcatenate = p$justConcatenate, verbose = TRUE), nm)
      cnt$denoisedR <- vapply(dR, getN, 0)
      cnt$merged <- vapply(mg, getN, 0)
      low <- cnt$merged / pmax(cnt$denoisedF, 1) < 0.5
      if (any(low)) message("ATTENTION : moins de 50 % des reads fusionnés pour ",
                            paste(nm[low], collapse = ", "),
                            " — chevauchement insuffisant (truncLen trop courts ?)")
      tables[[r]] <- makeSequenceTable(mg)
    } else {
      cnt$denoisedR <- NA; cnt$merged <- NA
      tables[[r]] <- makeSequenceTable(dF)
    }
    counts[[r]] <- cnt
  }
  seqtab <- if (length(tables) > 1) {
    message("Fusion des tables de ", length(tables), " runs")
    mergeSequenceTables(tables = unname(tables))
  } else tables[[1]]
  message(sprintf("Table : %d échantillons x %d ASV", nrow(seqtab), ncol(seqtab)))
  list(seqtab = seqtab, counts = do.call(rbind, counts))
}

# -----------------------------------------------------------------------------
# 6. Longueur + chimères
# -----------------------------------------------------------------------------
run_chimera <- function(seqtab, p, threads = 1) {
  lens <- nchar(colnames(seqtab))
  st_len <- seqtab[, lens >= p$lenMin & lens <= p$lenMax, drop = FALSE]
  if (!ncol(st_len)) stop("Aucune ASV dans la gamme de longueurs ", p$lenMin, "-", p$lenMax, " pb.")
  message(ncol(seqtab) - ncol(st_len), " ASV hors gamme de longueur retirées (",
          round(100 * (1 - sum(st_len) / sum(seqtab)), 2), " % des reads)")
  nochim <- removeBimeraDenovo(st_len, method = p$method,
                               minFoldParentOverAbundance = p$minFold,
                               multithread = threads, verbose = TRUE)
  message(sprintf("%d ASV non chimériques ; %.1f %% des reads conservés",
                  ncol(nochim), 100 * sum(nochim) / sum(st_len)))
  list(seqtab_len = st_len, nochim = nochim)
}

# -----------------------------------------------------------------------------
# 7. Contaminants (decontam) + filtres d'abondance
# -----------------------------------------------------------------------------
run_postfilter <- function(nochim, meta, p) {
  st <- nochim
  cont_tbl <- NULL
  if (!identical(p$decontam, "none")) {
    if (!requireNamespace("decontam", quietly = TRUE))
      stop("Le package decontam n'est pas installé (BiocManager::install('decontam')).")
    if (is.null(meta)) stop("decontam nécessite des métadonnées (onglet Données).")
    m <- meta[match(rownames(st), rownames(meta)), , drop = FALSE]
    if (anyNA(rownames(m))) stop("Échantillons absents des métadonnées : ",
                                 paste(rownames(st)[is.na(rownames(m))], collapse = ", "))
    neg <- NULL; conc <- NULL
    if (p$decontam %in% c("prevalence", "either", "both", "combined")) {
      v <- m[[p$neg_col]]
      neg <- if (is.logical(v)) v else tolower(trimws(as.character(v))) %in%
        tolower(trimws(parse_values(p$neg_values)))
      neg[is.na(neg)] <- FALSE
      if (!any(neg)) stop("Aucun témoin négatif trouvé dans la colonne « ", p$neg_col, " ».")
      message(sum(neg), " témoins négatifs identifiés")
    }
    if (p$decontam %in% c("frequency", "either", "both", "combined")) {
      conc <- suppressWarnings(as.numeric(m[[p$conc_col]]))
      if (anyNA(conc) || any(conc <= 0)) stop("La colonne de concentration doit être numérique et > 0 pour tous les échantillons.")
    }
    res <- decontam::isContaminant(st, neg = neg, conc = conc, method = p$decontam,
                                   threshold = p$threshold)
    cont_tbl <- data.frame(sequence = rownames(res), res, row.names = NULL)
    message(sum(res$contaminant, na.rm = TRUE), " ASV identifiées comme contaminants (",
            round(100 * sum(st[, which(res$contaminant)]) / sum(st), 2), " % des reads)")
    st <- st[, !res$contaminant %in% TRUE, drop = FALSE]
    if (isTRUE(p$drop_controls) && !is.null(neg)) {
      st <- st[!neg, , drop = FALSE]
      message("Témoins négatifs retirés de la table")
    }
  }
  if (p$rel_thresh > 0) {
    thr <- colSums(st) * p$rel_thresh / 100
    low <- sweep(st, 2, thr, "<") & st > 0
    message(sum(low), " occurrences < ", p$rel_thresh, " % du total de leur ASV mises à 0 (",
            sum(st[low]), " reads)")
    st[low] <- 0
  }
  if (p$min_reads_asv > 1) {
    keep <- colSums(st) >= p$min_reads_asv
    message(sum(!keep), " ASV avec moins de ", p$min_reads_asv, " reads retirées")
    st <- st[, keep, drop = FALSE]
  }
  st <- st[, colSums(st) > 0, drop = FALSE]
  if (p$min_reads_sample > 0) {
    keep <- rowSums(st) >= p$min_reads_sample
    if (any(!keep)) message("Échantillons sous ", p$min_reads_sample, " reads retirés : ",
                            paste(rownames(st)[!keep], collapse = ", "))
    st <- st[keep, , drop = FALSE]
    st <- st[, colSums(st) > 0, drop = FALSE]
  }
  message(sprintf("Table finale : %d échantillons x %d ASV, %s reads",
                  nrow(st), ncol(st), format(sum(st), big.mark = " ")))
  list(final = st, decontam = cont_tbl)
}

parse_values <- function(x) {
  v <- trimws(unlist(strsplit(x %||% "", ",")))
  v[nzchar(v)]
}

# -----------------------------------------------------------------------------
# 8. Taxonomie
# -----------------------------------------------------------------------------
download_ref <- function(url, cache = tools::R_user_dir("dada2shiny", "cache")) {
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cache, basename(url))
  if (!file.exists(dest) || file.size(dest) == 0) {
    message("Téléchargement de ", basename(url), " ...")
    old <- options(timeout = max(3600, getOption("timeout"))); on.exit(options(old))
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  } else message("Référence en cache : ", dest)
  dest
}

run_taxonomy <- function(seqs, p, threads = 1) {
  levels <- parse_values(p$levels)
  if (grepl("^https?://", p$ref)) p$ref <- download_ref(p$ref)
  if (!is.null(p$species) && grepl("^https?://", p$species)) p$species <- download_ref(p$species)
  if (!file.exists(p$ref)) stop("Fichier de référence introuvable : ", p$ref)
  message("Référence : ", p$ref, "\nMéthode : ", p$method, " ; ", length(seqs), " ASV")

  if (p$method == "idtaxa") {
    if (!requireNamespace("DECIPHER", quietly = TRUE))
      stop("Le package DECIPHER n'est pas installé (BiocManager::install('DECIPHER')).")
    e <- new.env(); load(p$ref, envir = e)
    ts <- if (exists("trainingSet", envir = e)) e$trainingSet else get(ls(e)[1], envir = e)
    if (!inherits(ts, "Taxa")) stop("Le fichier .RData ne contient pas de trainingSet IDTAXA.")
    ids <- DECIPHER::IdTaxa(Biostrings::DNAStringSet(seqs), ts,
                            strand = if (isTRUE(p$tryRC)) "both" else "top",
                            threshold = p$minBoot, processors = threads, verbose = FALSE)
    fix <- function(x, what) {
      v <- x[[what]]
      if (length(x$taxon) && tolower(x$taxon[1]) %in% c("root", "rootrank")) v <- v[-1]
      length(v) <- length(levels)
      v
    }
    taxa <- t(vapply(ids, fix, character(length(levels)), what = "taxon"))
    taxa[!is.na(taxa) & startsWith(taxa, "unclassified_")] <- NA
    boot <- t(vapply(ids, fix, numeric(length(levels)), what = "confidence"))
    dimnames(taxa) <- dimnames(boot) <- list(seqs, levels)
  } else {
    tx <- assignTaxonomy(seqs, p$ref, minBoot = p$minBoot, tryRC = isTRUE(p$tryRC),
                         taxLevels = levels, outputBootstraps = TRUE, multithread = threads)
    taxa <- tx$tax; boot <- tx$boot
  }

  if (!is.null(p$species) && nzchar(p$species)) {
    if ("Species" %in% colnames(taxa)) {
      message("Rang « Species » déjà présent : addSpecies ignoré.")
    } else {
      if (!file.exists(p$species)) stop("Fichier espèces introuvable : ", p$species)
      taxa <- addSpecies(taxa, p$species, tryRC = isTRUE(p$tryRC))
    }
  }
  pct <- round(100 * colMeans(!is.na(taxa)), 1)
  message("% d'ASV assignées : ", paste(colnames(taxa), pct, collapse = " | "))
  list(taxa = taxa, boot = boot)
}

# -----------------------------------------------------------------------------
# Suivi des lectures
# -----------------------------------------------------------------------------
make_track <- function(st) {
  sn <- st$samples$sample
  pick <- function(v, nm) if (is.null(v)) NULL else unname(v[match(sn, nm)])
  cols <- list(
    input     = if (!is.null(st$out)) pick(st$out[, 1], rownames(st$out)),
    cutadapt  = if (!is.null(st$cut_stats)) pick(st$cut_stats$cut_out, st$cut_stats$sample),
    filtered  = if (!is.null(st$out)) pick(st$out[, 2], rownames(st$out)),
    denoisedF = if (!is.null(st$counts)) pick(st$counts$denoisedF, st$counts$sample),
    denoisedR = if (!is.null(st$counts) && is_paired(st$samples)) pick(st$counts$denoisedR, st$counts$sample),
    merged    = if (!is.null(st$counts) && is_paired(st$samples)) pick(st$counts$merged, st$counts$sample),
    length_ok = if (!is.null(st$seqtab_len)) pick(rowSums(st$seqtab_len), rownames(st$seqtab_len)),
    nonchim   = if (!is.null(st$nochim)) pick(rowSums(st$nochim), rownames(st$nochim)),
    final     = if (!is.null(st$final)) pick(rowSums(st$final), rownames(st$final)))
  if (!is.null(st$cut_stats)) cols$input <- pick(st$cut_stats$cut_in, st$cut_stats$sample)
  cols <- Filter(Negate(is.null), cols)
  if (!length(cols)) return(NULL)
  tr <- as.data.frame(cols, row.names = sn)
  tr[is.na(tr)] <- 0
  last <- tr[[ncol(tr)]]
  tr$pct_conserve <- round(100 * last / pmax(tr[[1]], 1), 1)
  tr
}

# -----------------------------------------------------------------------------
# Script R reproductible à partir des paramètres réellement utilisés
# -----------------------------------------------------------------------------
generate_script <- function(st) {
  P <- st$params
  dp <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  s <- st$samples
  paired <- is_paired(s)
  L <- c(
    "# =============================================================================",
    "#  Pipeline DADA2 — script généré par DADA2 Shiny",
    paste0("#  Date : ", format(Sys.time(), "%Y-%m-%d %H:%M"),
           " ; dada2 ", as.character(utils::packageVersion("dada2")),
           " ; ", R.version.string),
    "#  Les paramètres sont ceux réellement utilisés à chaque étape.",
    "# =============================================================================",
    "library(dada2)", "",
    paste0("outdir  <- ", dp(P$data$outdir)),
    paste0("threads <- ", P$data$threads),
    "dir.create(outdir, recursive = TRUE, showWarnings = FALSE)", "",
    "# ---- Échantillons détectés --------------------------------------------------",
    "samples <- data.frame(",
    paste0("  sample = ", dp(s$sample), ","),
    paste0("  run    = ", dp(s$run), ","),
    paste0("  fnF    = ", dp(s$fnF), ","),
    paste0("  fnR    = ", dp(s$fnR), ")"),
    "inF <- samples$fnF; inR <- samples$fnR", "")

  pr <- P$primers
  if (!is.null(pr) && identical(pr$method, "cutadapt")) {
    L <- c(L,
      "# ---- Retrait des amorces (cutadapt) -----------------------------------------",
      paste0("primersF <- ", dp(parse_primers(pr$primersF))),
      paste0("primersR <- ", dp(parse_primers(pr$primersR))),
      "rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))",
      "flags <- c(c(rbind(\"-g\", primersF)), c(rbind(\"-a\", rc(primersR)))",
      if (paired) "           , c(rbind(\"-G\", primersR)), c(rbind(\"-A\", rc(primersF)))",
      "           )",
      paste0("flags <- c(flags, \"-n\", 2, \"-m\", ", pr$cut_minlen, ", \"-j\", threads",
             if (isTRUE(pr$discard_untrimmed)) ", \"--discard-untrimmed\"", ")"),
      "cut_dir <- file.path(outdir, \"cutadapt\"); dir.create(cut_dir, showWarnings = FALSE)",
      "trimF <- file.path(cut_dir, paste0(samples$sample, \"_R1_trim.fastq.gz\"))",
      if (paired) "trimR <- file.path(cut_dir, paste0(samples$sample, \"_R2_trim.fastq.gz\"))",
      "for (i in seq_along(inF)) {",
      if (paired)
        paste0("  system2(", dp(pr$cutadapt), ", c(flags, \"-o\", shQuote(trimF[i]), \"-p\", shQuote(trimR[i]),",
               " shQuote(inF[i]), shQuote(inR[i])))")
      else paste0("  system2(", dp(pr$cutadapt), ", c(flags, \"-o\", shQuote(trimF[i]), shQuote(inF[i])))"),
      "}",
      "inF <- trimF", if (paired) "inR <- trimR", "")
  }

  f <- P$filter
  if (!is.null(f)) {
    L <- c(L,
      "# ---- Filtrage ---------------------------------------------------------------",
      "filtF <- file.path(outdir, \"filtered\", paste0(samples$sample, \"_F_filt.fastq.gz\"))",
      if (paired) "filtR <- file.path(outdir, \"filtered\", paste0(samples$sample, \"_R_filt.fastq.gz\"))",
      if (paired)
        paste0("out <- filterAndTrim(inF, filtF, inR, filtR, truncLen = ", dp(f$truncLen),
               ", trimLeft = ", dp(f$trimLeft), ", maxEE = ", dp(f$maxEE),
               ",\n                     truncQ = ", f$truncQ, ", minLen = ", f$minLen,
               ", maxN = 0, rm.phix = ", f$rm.phix,
               ", compress = TRUE,\n                     multithread = ",
               "if (.Platform$OS.type == \"windows\") FALSE else threads)")
      else
        paste0("out <- filterAndTrim(inF, filtF, truncLen = ", f$truncLen[1],
               ", trimLeft = ", f$trimLeft[1], ", maxEE = ", f$maxEE[1],
               ",\n                     truncQ = ", f$truncQ, ", minLen = ", f$minLen,
               ", maxN = 0, rm.phix = ", f$rm.phix,
               ", compress = TRUE,\n                     multithread = ",
               "if (.Platform$OS.type == \"windows\") FALSE else threads)"),
      "rownames(out) <- samples$sample",
      paste0("ok <- file.exists(filtF)", if (paired) " & file.exists(filtR)"), "")
  }

  e <- P$errors; d <- P$denoise
  if (!is.null(e)) {
    if (identical(e$errfun, "binned")) {
      fun_src <- deparse(binned_qual_errfun)
      L <- c(L, "# Fonction d'erreur pour qualités binnées (NovaSeq/NextSeq), d'après dada2#1307",
             paste("binned_qual_errfun <-", fun_src[1]), fun_src[-1], "")
    }
    errfun <- if (identical(e$errfun, "binned")) "binned_qual_errfun" else "loessErrfun"
    L <- c(L,
      "# ---- Modèle d'erreur, débruitage et fusion, run par run -------------------",
      "tables <- list()",
      "for (r in unique(samples$run[ok])) {",
      "  i <- ok & samples$run == r",
      "  fF <- setNames(filtF[i], samples$sample[i])",
      paste0("  errF <- learnErrors(fF, nbases = ", format(e$nbases, scientific = TRUE),
             ", randomize = ", e$randomize, ", errorEstimationFunction = ", errfun,
             ", multithread = threads)"),
      if (paired) "  fR <- setNames(filtR[i], samples$sample[i])",
      if (paired)
        paste0("  errR <- learnErrors(fR, nbases = ", format(e$nbases, scientific = TRUE),
               ", randomize = ", e$randomize, ", errorEstimationFunction = ", errfun,
               ", multithread = threads)"))
    if (!is.null(d)) {
      pool <- switch(as.character(d$pool), "TRUE" = "TRUE", "pseudo" = "\"pseudo\"", "FALSE")
      L <- c(L,
        paste0("  dadaF <- dada(fF, err = errF, pool = ", pool, ", multithread = threads)"),
        if (paired) paste0("  dadaR <- dada(fR, err = errR, pool = ", pool, ", multithread = threads)"),
        if (paired)
          paste0("  mergers <- mergePairs(dadaF, fF, dadaR, fR, minOverlap = ", d$minOverlap,
                 ", maxMismatch = ", d$maxMismatch, ",\n                        trimOverhang = ",
                 d$trimOverhang, ", justConcatenate = ", d$justConcatenate, ")"),
        paste0("  tables[[r]] <- makeSequenceTable(", if (paired) "mergers" else "dadaF", ")"))
    }
    L <- c(L, "}",
      "seqtab <- if (length(tables) > 1) mergeSequenceTables(tables = unname(tables)) else tables[[1]]", "")
  }

  ch <- P$chimera
  if (!is.null(ch)) {
    L <- c(L,
      "# ---- Longueur + chimères ---------------------------------------------------",
      paste0("lens <- nchar(colnames(seqtab))"),
      paste0("seqtab.len <- seqtab[, lens >= ", ch$lenMin, " & lens <= ", ch$lenMax, ", drop = FALSE]"),
      paste0("seqtab.nochim <- removeBimeraDenovo(seqtab.len, method = \"", ch$method,
             "\", minFoldParentOverAbundance = ", ch$minFold, ", multithread = threads, verbose = TRUE)"),
      "final <- seqtab.nochim", "")
  }

  pf <- P$post
  if (!is.null(pf)) {
    L <- c(L, "# ---- Contaminants et filtres d'abondance ----------------------------------")
    if (!identical(pf$decontam, "none")) {
      L <- c(L,
        paste0("meta <- read.csv(", dp(P$data$meta_file %||% "metadata.csv"), ", row.names = 1, check.names = FALSE)"),
        "meta <- meta[rownames(final), , drop = FALSE]",
        if (pf$decontam %in% c("prevalence", "either", "both", "combined"))
          paste0("neg <- tolower(as.character(meta[[", dp(pf$neg_col), "]])) %in% ",
                 dp(tolower(parse_values(pf$neg_values)))) else "neg <- NULL",
        if (pf$decontam %in% c("frequency", "either", "both", "combined"))
          paste0("conc <- as.numeric(meta[[", dp(pf$conc_col), "]])") else "conc <- NULL",
        paste0("contam <- decontam::isContaminant(final, neg = neg, conc = conc, method = \"",
               pf$decontam, "\", threshold = ", pf$threshold, ")"),
        "final <- final[, !contam$contaminant %in% TRUE, drop = FALSE]",
        if (isTRUE(pf$drop_controls) && !is.null(pf$neg_col)) "if (!is.null(neg)) final <- final[!neg, , drop = FALSE]")
    }
    if (pf$rel_thresh > 0) L <- c(L,
      paste0("final[sweep(final, 2, colSums(final) * ", pf$rel_thresh, " / 100, \"<\")] <- 0  # tag-jumps"))
    if (pf$min_reads_asv > 1) L <- c(L,
      paste0("final <- final[, colSums(final) >= ", pf$min_reads_asv, ", drop = FALSE]"))
    L <- c(L, "final <- final[, colSums(final) > 0, drop = FALSE]")
    if (pf$min_reads_sample > 0) L <- c(L,
      paste0("final <- final[rowSums(final) >= ", pf$min_reads_sample, ", , drop = FALSE]"),
      "final <- final[, colSums(final) > 0, drop = FALSE]")
    L <- c(L, "")
  }

  tx <- P$taxonomy
  if (!is.null(tx)) {
    L <- c(L, "# ---- Taxonomie ----------------------------------------------------------------")
    if (tx$method == "idtaxa") {
      L <- c(L,
        paste0("load(", dp(tx$ref), ")  # objet trainingSet"),
        paste0("ids <- DECIPHER::IdTaxa(Biostrings::DNAStringSet(colnames(final)), trainingSet, strand = \"",
               if (isTRUE(tx$tryRC)) "both" else "top", "\", threshold = ", tx$minBoot,
               ", processors = threads)"),
        paste0("ranks <- ", dp(parse_values(tx$levels))),
        "taxa <- t(sapply(ids, function(x) { t <- x$taxon[-1]; length(t) <- length(ranks); t }))",
        "taxa[!is.na(taxa) & startsWith(taxa, \"unclassified_\")] <- NA",
        "dimnames(taxa) <- list(colnames(final), ranks)")
    } else {
      L <- c(L,
        paste0("taxa <- assignTaxonomy(final, ", dp(tx$ref), ", minBoot = ", tx$minBoot,
               ", tryRC = ", isTRUE(tx$tryRC), ",\n                       taxLevels = ",
               dp(parse_values(tx$levels)), ", multithread = threads)"))
    }
    if (!is.null(tx$species) && nzchar(tx$species))
      L <- c(L, paste0("taxa <- addSpecies(taxa, ", dp(tx$species), ")"))
    L <- c(L, "")
  }

  L <- c(L,
    "# ---- Export ------------------------------------------------------------------",
    "asv_id <- paste0(\"ASV\", seq_len(ncol(final)))",
    "writeLines(paste0(\">\", asv_id, \"\\n\", colnames(final)), file.path(outdir, \"asv_sequences.fasta\"))",
    "tab <- data.frame(ASV = asv_id, sequence = colnames(final), t(final), check.names = FALSE)",
    if (!is.null(tx)) "tab <- cbind(tab[, 1:2], taxa[colnames(final), , drop = FALSE], tab[, -(1:2), drop = FALSE])",
    "write.csv(tab, file.path(outdir, \"asv_table.csv\"), row.names = FALSE)")
  paste(L[!is.na(L)], collapse = "\n")
}

# =============================================================================
#  curation.R — réplicats PCR, filtre taxonomique, raréfaction, communauté mock
#               et préparation de bases de référence. Indépendant de Shiny.
#
#  Convention : `counts` = matrice échantillons (ou PCR) x ASV, colonnes = séquences.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# Dissimilarité de Bray-Curtis entre les lignes de deux matrices de proportions
.bc_rows <- function(A, B) {
  A <- A / pmax(rowSums(A), 1e-12); B <- B / pmax(rowSums(B), 1e-12)
  outer(seq_len(nrow(A)), seq_len(nrow(B)),
        Vectorize(function(i, j) 1 - sum(pmin(A[i, ], B[j, ]))))
}

# -----------------------------------------------------------------------------
# 1. Réplicats PCR (inspiré de metabaR::pcrslayer, Zinger et al. 2021)
# -----------------------------------------------------------------------------
#' Seuil de dissimilarité : intersection des densités intra-échantillon (dw) et
#' inter-échantillons (db), au-delà de la médiane de dw. NA si non déterminable.
replicate_threshold <- function(dw, db) {
  dw <- dw[is.finite(dw)]; db <- db[is.finite(db)]
  if (length(dw) < 3 || length(db) < 3) return(NA_real_)
  dW <- stats::density(dw, from = 0, to = 1, n = 512)
  dB <- stats::density(db, from = 0, to = 1, n = 512)
  d <- dW$y - dB$y
  cross <- which(d[-length(d)] >= 0 & d[-1] < 0) + 1
  x <- dW$x[cross]
  x <- x[x > stats::median(dw)]
  if (length(x)) x[1] else NA_real_
}

#' Repère les PCR dont la composition s'écarte trop des autres réplicats du même
#' échantillon (distance au centroïde des réplicats > seuil). Itératif : après
#' retrait des PCR aberrantes, les centroïdes sont recalculés.
#' @param groups échantillon biologique de chaque PCR (même ordre que les lignes)
#' @param threshold seuil fixé ; NA = automatique
detect_failed_pcrs <- function(counts, groups, threshold = NA, max_iter = 10) {
  g <- as.character(groups); g[is.na(g) | !nzchar(g)] <- rownames(counts)[is.na(g) | !nzchar(g)]
  rel <- counts / pmax(rowSums(counts), 1)
  pcr <- rownames(counts)
  failed <- setNames(rep(FALSE, length(pcr)), pcr)
  iter <- setNames(rep(NA_integer_, length(pcr)), pcr)
  dw_first <- NULL; db_first <- NULL; thr <- threshold; auto <- is.na(threshold)

  for (it in seq_len(max_iter)) {
    ok <- !failed
    gg <- g[ok]; rr <- rel[ok, , drop = FALSE]
    n_rep <- table(gg)
    multi <- names(n_rep)[n_rep >= 2]
    if (!length(multi)) break
    cen <- rowsum(rr, gg, reorder = FALSE)
    cen <- cen / as.vector(n_rep[rownames(cen)])
    dw <- vapply(seq_len(nrow(rr)), function(i)
      if (gg[i] %in% multi) 1 - sum(pmin(rr[i, ] / max(sum(rr[i, ]), 1e-12),
                                         cen[gg[i], ] / max(sum(cen[gg[i], ]), 1e-12))) else NA_real_, 0)
    names(dw) <- rownames(rr)
    if (it == 1) {
      dw_first <- dw
      db_m <- .bc_rows(cen, cen)
      db_first <- db_m[upper.tri(db_m)]
      if (auto) thr <- replicate_threshold(dw[!is.na(dw)], db_first)
      if (is.na(thr)) break
    }
    over <- !is.na(dw) & dw > thr
    if (!any(over)) break
    # par échantillon, on ne retire que la PCR la plus éloignée, puis on recalcule
    # le centroïde ; avec 2 réplicats seulement, on ne peut pas savoir lequel est
    # mauvais : les deux sont signalés
    bad <- unlist(lapply(unique(gg[over]), function(s) {
      i <- which(gg == s)
      if (length(i) == 2) return(names(dw)[i])
      names(dw)[i][which(dw[i] == max(dw[i]))]
    }))
    failed[bad] <- TRUE; iter[bad] <- it
  }
  tab <- data.frame(PCR = pcr, echantillon = g, n_replicats = as.vector(table(g)[g]),
                    dissimilarite = round(unname((dw_first %||% setNames(rep(NA, length(pcr)), pcr))[pcr]), 3),
                    aberrante = unname(failed), iteration = unname(iter), stringsAsFactors = FALSE)
  if (is.na(thr)) message("Réplicats : seuil non déterminable automatiquement (trop peu de réplicats ",
                          "ou distributions qui ne se séparent pas) — aucune PCR retirée. Fixe un seuil à la main si besoin.")
  else message(sprintf("Réplicats : seuil de dissimilarité %s = %.3f ; %d PCR aberrante(s) : %s",
                       if (auto) "automatique" else "fixé", thr, sum(failed),
                       if (any(failed)) paste(pcr[failed], collapse = ", ") else "aucune"))
  list(table = tab, db = db_first %||% numeric(0), threshold = thr, auto = auto)
}

#' Fusionne les réplicats PCR par échantillon biologique
#' @param fun "sum" (somme des reads) ou "mean" (moyenne arrondie)
#' @param min_reps une ASV n'est gardée dans un échantillon que si elle est
#'   présente dans au moins min_reps réplicats (plafonné au nombre de réplicats)
aggregate_replicates <- function(counts, groups, fun = c("sum", "mean"), min_reps = 1) {
  fun <- match.arg(fun)
  g <- as.character(groups); g[is.na(g) | !nzchar(g)] <- rownames(counts)[is.na(g) | !nzchar(g)]
  agg <- rowsum(counts, g, reorder = FALSE)
  n <- as.vector(table(g)[rownames(agg)])
  if (fun == "mean") agg <- round(agg / n)
  if (min_reps > 1) {
    pres <- rowsum((counts > 0) * 1, g, reorder = FALSE)
    need <- pmin(min_reps, n)
    low <- sweep(pres, 1, need, "<")
    message(sprintf("Réplicats : %d occurrences présentes dans moins de %d réplicat(s) mises à 0 (%s reads)",
                    sum(low & agg > 0), min_reps, format(sum(agg[low]), big.mark = " ")))
    agg[low] <- 0
  }
  agg <- agg[, colSums(agg) > 0, drop = FALSE]
  message(sprintf("Réplicats fusionnés (%s) : %d PCR -> %d échantillons", fun, nrow(counts), nrow(agg)))
  storage.mode(agg) <- "integer"
  agg
}

#' Métadonnées au niveau de l'échantillon biologique : une ligne par échantillon,
#' colonnes gardées seulement si elles sont constantes entre réplicats
collapse_meta <- function(meta, col) {
  g <- as.character(meta[[col]]); g[is.na(g)] <- rownames(meta)[is.na(g)]
  keep <- vapply(meta, function(v) all(tapply(as.character(v), g, function(x) length(unique(x)) <= 1)), TRUE)
  out <- meta[!duplicated(g), keep, drop = FALSE]
  rownames(out) <- g[!duplicated(g)]
  out$n_replicats <- as.vector(table(g)[rownames(out)])
  out
}

# -----------------------------------------------------------------------------
# 2. Filtre taxonomique (équivalent de remove.lineage de mothur)
# -----------------------------------------------------------------------------
#' @param p list(remove = "Chloroplast, Mitochondria", keep = "Bacillariophyta",
#'               unassigned_rank = "" ou nom de rang)
#' Les noms sont comparés exactement (sans tenir compte de la casse) à tous les rangs.
run_taxfilter <- function(final, taxa, p) {
  vals <- function(x) tolower(trimws(unlist(strsplit(x %||% "", ","))))
  rm_v <- vals(p$remove); rm_v <- rm_v[nzchar(rm_v)]
  kp_v <- vals(p$keep); kp_v <- kp_v[nzchar(kp_v)]
  tx <- taxa[colnames(final), , drop = FALSE]
  low <- tolower(tx)
  hit <- function(v) if (!length(v)) rep(FALSE, nrow(tx)) else apply(low, 1, function(r) any(r %in% v))
  rm <- hit(rm_v)
  kp <- if (length(kp_v)) hit(kp_v) else rep(TRUE, nrow(tx))
  un <- if (nzchar(p$unassigned_rank %||% "") && p$unassigned_rank %in% colnames(tx))
    is.na(tx[, p$unassigned_rank]) else rep(FALSE, nrow(tx))
  reads <- colSums(final)
  report <- data.frame(
    Critere = c(if (length(rm_v)) paste("Retirés :", paste(rm_v, collapse = ", ")),
                if (length(kp_v)) paste("Hors groupe gardé :", paste(kp_v, collapse = ", ")),
                if (any(un)) paste("Non assignés au rang", p$unassigned_rank)),
    stringsAsFactors = FALSE)
  crit <- list(if (length(rm_v)) rm, if (length(kp_v)) !kp, if (any(un)) un)
  crit <- Filter(Negate(is.null), crit)
  if (length(crit)) {
    report$ASV <- vapply(crit, sum, 0)
    report$reads <- vapply(crit, function(x) sum(reads[x]), 0)
  }
  sel <- kp & !rm & !un
  out <- final[, sel, drop = FALSE]
  empty <- rowSums(out) == 0
  if (any(empty)) message("Échantillons vides après filtre taxonomique (retirés) : ",
                          paste(rownames(out)[empty], collapse = ", "))
  out <- out[!empty, , drop = FALSE]
  message(sprintf("Filtre taxonomique : %d / %d ASV gardées (%.1f %% des reads)",
                  ncol(out), ncol(final), 100 * sum(out) / max(1, sum(final))))
  list(final_tax = out, report = report)
}

# -----------------------------------------------------------------------------
# 3. Raréfaction
# -----------------------------------------------------------------------------
#' Courbes de raréfaction : richesse attendue (formule exacte d'Hurlbert)
rarefaction_curves <- function(counts, n_points = 40) {
  do.call(rbind, lapply(rownames(counts), function(s) {
    x <- counts[s, ]; x <- x[x > 0]; N <- sum(x)
    if (N == 0) return(NULL)
    depths <- unique(round(seq(1, N, length.out = n_points)))
    rich <- vapply(depths, function(n) sum(1 - exp(lchoose(N - x, n) - lchoose(N, n))), 0)
    data.frame(sample = s, depth = depths, richness = rich, stringsAsFactors = FALSE)
  }))
}

#' Diversité moyenne sur des raréfactions répétées (esprit de q2-boots)
#' @return list(alpha = data.frame(sample, Observed, Shannon, Simpson + écarts-types),
#'              bray = dist moyenne, dropped = échantillons sous la profondeur)
rarefied_diversity <- function(counts, depth, n_iter = 20, seed = 1) {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Le package vegan est nécessaire (installé avec phyloseq).")
  keep <- rowSums(counts) >= depth
  dropped <- rownames(counts)[!keep]
  cts <- counts[keep, , drop = FALSE]
  if (nrow(cts) < 2) stop("Moins de 2 échantillons atteignent la profondeur ", depth, ".")
  set.seed(seed)
  obs <- sha <- sim <- matrix(NA_real_, nrow(cts), n_iter, dimnames = list(rownames(cts), NULL))
  bray <- matrix(0, nrow(cts), nrow(cts))
  for (i in seq_len(n_iter)) {
    r <- vegan::rrarefy(cts, depth)
    obs[, i] <- rowSums(r > 0)
    sha[, i] <- vegan::diversity(r, "shannon")
    sim[, i] <- vegan::diversity(r, "simpson")
    bray <- bray + as.matrix(vegan::vegdist(r, "bray"))
  }
  bray <- bray / n_iter
  dimnames(bray) <- list(rownames(cts), rownames(cts))
  alpha <- data.frame(sample = rownames(cts),
                      Observed = rowMeans(obs), Observed_sd = apply(obs, 1, stats::sd),
                      Shannon = rowMeans(sha), Shannon_sd = apply(sha, 1, stats::sd),
                      Simpson = rowMeans(sim), Simpson_sd = apply(sim, 1, stats::sd),
                      row.names = NULL)
  list(alpha = alpha, bray = stats::as.dist(bray), dropped = dropped, depth = depth, n_iter = n_iter)
}

# -----------------------------------------------------------------------------
# 4. Communauté mock (esprit de seq.error de mothur)
# -----------------------------------------------------------------------------
read_fasta_simple <- function(path) {
  l <- readLines(path, warn = FALSE)                  # gère aussi les .gz
  l <- trimws(l); l <- l[nzchar(l)]
  h <- grep("^>", l)
  if (!length(h)) stop("Fichier FASTA sans en-tête (>).")
  ends <- c(h[-1] - 1, length(l))
  seqs <- vapply(seq_along(h), function(i)
    if (ends[i] > h[i]) paste(l[(h[i] + 1):ends[i]], collapse = "") else "", "")
  setNames(seqs, sub("^>\\s*", "", l[h]))
}

# Fonctions d'alignement par paires : déplacées de Biostrings vers pwalign
# (Bioconductor >= 3.19, retirées de Biostrings >= 2.77). On prend pwalign s'il
# est installé et exporte la fonction, sinon Biostrings (anciennes versions).
.pw <- function(f) {
  if (requireNamespace("pwalign", quietly = TRUE) && f %in% getNamespaceExports("pwalign"))
    getExportedValue("pwalign", f) else getExportedValue("Biostrings", f)
}
.aligner <- function() .pw("pairwiseAlignment")
.nindel_sum <- function(aln) {
  ni <- .pw("nindel")(aln)
  ni@insertion[, "WidthSum"] + ni@deletion[, "WidthSum"]
}
.nmismatch <- function(aln) .pw("nmismatch")(aln)

#' Nombre de différences (substitutions + indels) entre chaque ASV et chaque
#' séquence attendue ; l'ASV est alignée en entier, la référence localement
#' (les références peuvent être plus longues que l'amplicon), dans les deux sens.
mock_distances <- function(asv_seqs, ref_seqs) {
  pa <- .aligner()
  sm <- .pw("nucleotideSubstitutionMatrix")(match = 1, mismatch = -3, baseOnly = FALSE)
  refs <- Biostrings::DNAStringSet(ref_seqs)
  refs_rc <- Biostrings::reverseComplement(refs)
  D <- matrix(NA_real_, length(asv_seqs), length(ref_seqs),
              dimnames = list(names(asv_seqs), names(ref_seqs)))
  for (i in seq_along(asv_seqs)) {
    q <- Biostrings::DNAStringSet(rep(asv_seqs[[i]], length(refs)))
    best <- rep(Inf, length(refs))
    for (subj in list(refs, refs_rc)) {
      aln <- pa(q, subj, type = "global-local", substitutionMatrix = sm,
                gapOpening = 5, gapExtension = 2)
      best <- pmin(best, .nmismatch(aln) + .nindel_sum(aln))
    }
    D[i, ] <- best
  }
  D
}

#' Évalue une communauté mock
#' @param counts table échantillons x ASV ; @param mock_samples échantillon(s) mock
#' @param ref_seqs séquences attendues (vecteur nommé) ; @param max_diff seuil
#'   de différences pour considérer une ASV comme une variante d'erreur d'une
#'   séquence attendue plutôt qu'une séquence inattendue
mock_eval <- function(counts, mock_samples, ref_seqs, asv_ids = colnames(counts), max_diff = 3) {
  x <- colSums(counts[mock_samples, , drop = FALSE])
  ids <- setNames(asv_ids, colnames(counts))
  x <- x[x > 0]
  if (!length(x)) stop("Aucune ASV dans l'échantillon mock.")
  ref_seqs <- toupper(gsub("[^A-Za-z]", "", ref_seqs))
  ref_seqs <- gsub("U", "T", ref_seqs)
  names(ref_seqs) <- sub("\\s.*$", "", names(ref_seqs))
  message(sprintf("Mock : %d ASV comparées à %d séquences attendues…", length(x), length(ref_seqs)))
  D <- mock_distances(setNames(names(x), ids[names(x)]), ref_seqs)
  best_ref <- colnames(D)[apply(D, 1, which.min)]
  best_d <- apply(D, 1, min)
  cls <- ifelse(best_d == 0, "Exacte",
                ifelse(best_d <= max_diff, sprintf("Variante d'erreur (1-%d diff.)", max_diff),
                       "Inattendue (contaminant, chimère…)"))
  asv <- data.frame(ASV = ids[names(x)], reads = as.vector(x),
                    pct_reads = round(100 * as.vector(x) / sum(x), 2),
                    reference_proche = best_ref, differences = best_d,
                    identite = round(100 * (1 - best_d / nchar(names(x))), 2),
                    classe = cls, longueur = nchar(names(x)), sequence = names(x),
                    stringsAsFactors = FALSE, row.names = NULL)
  asv <- asv[order(-asv$reads), ]
  per_ref <- do.call(rbind, lapply(names(ref_seqs), function(r) {
    a <- asv[asv$reference_proche == r, , drop = FALSE]
    data.frame(reference = r,
               retrouvee = if (any(a$differences == 0)) "oui (exacte)"
                           else if (nrow(a) && min(a$differences) <= max_diff) "approchée" else "non",
               n_ASV_exactes = sum(a$differences == 0),
               n_variantes = sum(a$differences > 0 & a$differences <= max_diff),
               reads = sum(a$reads[a$differences <= max_diff]),
               stringsAsFactors = FALSE)
  }))
  ok <- asv$differences <= max_diff
  err_rate <- sum(asv$reads[ok] * asv$differences[ok]) / max(1, sum(asv$reads[ok] * asv$longueur[ok]))
  summary <- data.frame(
    Indicateur = c("Séquences attendues", "Retrouvées exactement", "Retrouvées (exactes ou approchées)",
                   "ASV dans le mock", "ASV exactes", "ASV variantes d'erreur", "ASV inattendues",
                   "% reads sur ASV exactes", "% reads sur variantes d'erreur", "% reads sur ASV inattendues",
                   "Taux d'erreur résiduel (par base)"),
    Valeur = as.character(c(length(ref_seqs), sum(per_ref$retrouvee == "oui (exacte)"), sum(per_ref$retrouvee != "non"),
               nrow(asv), sum(asv$differences == 0), sum(asv$differences > 0 & ok), sum(!ok),
               round(100 * sum(asv$reads[asv$differences == 0]) / sum(asv$reads), 2),
               round(100 * sum(asv$reads[asv$differences > 0 & ok]) / sum(asv$reads), 2),
               round(100 * sum(asv$reads[!ok]) / sum(asv$reads), 2),
               signif(err_rate, 3))),
    stringsAsFactors = FALSE)
  message(sprintf("Mock : %d/%d séquences retrouvées exactement, %d ASV inattendues",
                  sum(per_ref$retrouvee == "oui (exacte)"), length(ref_seqs), sum(!ok)))
  list(summary = summary, asv = asv, per_ref = per_ref, max_diff = max_diff)
}

# -----------------------------------------------------------------------------
# 5. Préparation d'une base de référence (formats DADA2 et IDTAXA)
# -----------------------------------------------------------------------------
.na_tax <- c("", "na", "nan", "none", "unknown", "unclassified", "unidentified", "unassigned",
             "incertae_sedis", "incertae sedis", "uncultured")

#' @param fasta chemin du FASTA (éventuellement .gz)
#' @param p liste de paramètres :
#'   tax_source "header" (taxonomie dans l'en-tête) ou "table" (fichier séparé)
#'   header_id  TRUE si l'en-tête commence par un identifiant suivi d'un espace
#'   sep        séparateur des rangs (";" par défaut ; "|" ou "," possibles)
#'   table      chemin du tableau (CSV/TSV) ; table_id = colonne identifiant ;
#'              table_tax = colonne unique de taxonomie (avec sep) OU vide pour
#'              utiliser les colonnes nommées comme les rangs
#'   levels     noms des rangs, du plus large au plus fin (ex. "Kingdom,...,Species")
#'   strip_prefix TRUE pour retirer les préfixes de type k__, p__, d__
#'   min_len    longueur minimale
#'   conflicts  "lca" (ancêtre commun) ou "drop" pour une même séquence à
#'              taxonomies différentes
#'   species_in_train TRUE pour garder le rang Species dans le fichier assignTaxonomy
#'   idtaxa     TRUE pour entraîner aussi un classificateur IDTAXA (DECIPHER)
#'   name       nom de base des fichiers produits
format_reference <- function(fasta, p, outdir) {
  levels <- trimws(unlist(strsplit(p$levels, ","))); levels <- levels[nzchar(levels)]
  if (length(levels) < 2) stop("Indique au moins deux rangs taxonomiques.")
  sep <- p$sep %||% ";"
  seqs <- read_fasta_simple(fasta)
  n_in <- length(seqs)
  message("Base de référence : ", n_in, " séquences lues")

  # --- taxonomie
  hdr <- names(seqs)
  if (identical(p$tax_source, "table")) {
    first <- readLines(p$table, n = 1, warn = FALSE)
    s <- if (grepl("\t", first)) "\t" else if (grepl(";", first) && !identical(sep, ";")) ";" else ","
    tb <- utils::read.table(p$table, header = TRUE, sep = s, quote = "\"", comment.char = "",
                            stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    idc <- p$table_id %||% names(tb)[1]
    if (!idc %in% names(tb)) stop("Colonne identifiant « ", idc, " » absente du tableau.")
    ids <- sub("\\s.*$", "", hdr)
    m <- match(ids, tb[[idc]])
    if (all(is.na(m))) stop("Aucun identifiant du FASTA ne correspond à la colonne « ", idc, " ».")
    if (anyNA(m)) message(sum(is.na(m)), " séquences sans taxonomie dans le tableau (retirées)")
    tax <- if (nzchar(p$table_tax %||% "")) {
      lapply(strsplit(tb[[p$table_tax]][m], sep, fixed = TRUE), trimws)
    } else {
      miss <- setdiff(levels, names(tb))
      if (length(miss)) stop("Colonnes de rangs absentes du tableau : ", paste(miss, collapse = ", "))
      lapply(seq_along(m), function(i) if (is.na(m[i])) NA_character_ else unlist(tb[m[i], levels]))
    }
    tax[is.na(m)] <- list(NA_character_)
  } else {
    ids <- if (isTRUE(p$header_id)) sub("\\s.*$", "", hdr) else paste0("ref", seq_along(hdr))
    txt <- if (isTRUE(p$header_id)) sub("^\\S+\\s*", "", hdr) else hdr
    tax <- lapply(strsplit(txt, sep, fixed = TRUE), trimws)
  }
  clean_tax <- function(v) {
    v <- as.character(v)
    if (isTRUE(p$strip_prefix)) v <- sub("^[a-z]__", "", v)
    v <- gsub("_", " ", v); v <- trimws(v)
    v[tolower(v) %in% .na_tax | is.na(v)] <- NA
    length(v) <- length(levels)
    k <- which(is.na(v))[1]                     # tronque au premier rang manquant
    if (!is.na(k)) v[k:length(v)] <- NA
    v
  }
  T <- do.call(rbind, lapply(tax, clean_tax))
  colnames(T) <- levels

  # --- séquences
  sq <- toupper(gsub("[-.[:space:]]", "", seqs)); sq <- gsub("U", "T", sq)
  bad <- grepl("[^ACGTRYSWKMBDHVN]", sq)
  short <- nchar(sq) < (p$min_len %||% 50)
  notax <- is.na(T[, 1])
  keep <- !bad & !short & !notax
  message(sprintf("Retirées : %d caractères invalides, %d trop courtes (< %d pb), %d sans taxonomie",
                  sum(bad), sum(short & !bad), p$min_len %||% 50, sum(notax & !bad & !short)))
  sq <- sq[keep]; T <- T[keep, , drop = FALSE]; ids <- ids[keep]

  # --- doublons et conflits
  key_tax <- apply(T, 1, paste, collapse = ";")
  dup_exact <- duplicated(paste(sq, key_tax))
  sq <- sq[!dup_exact]; T <- T[!dup_exact, , drop = FALSE]; ids <- ids[!dup_exact]
  n_conf <- 0
  dseq <- unique(sq[duplicated(sq)])
  if (length(dseq)) {
    n_conf <- length(dseq)
    drop <- rep(FALSE, length(sq))
    for (s in dseq) {
      i <- which(sq == s)
      if (identical(p$conflicts, "drop")) { drop[i] <- TRUE; next }
      lca <- T[i[1], ]
      for (j in i[-1]) { d <- which(T[j, ] != lca | is.na(T[j, ]) != is.na(lca))[1]
                         if (!is.na(d)) lca[d:length(lca)] <- NA }
      T[i[1], ] <- lca; drop[i[-1]] <- TRUE
    }
    sq <- sq[!drop]; T <- T[!drop, , drop = FALSE]; ids <- ids[!drop]
    message(sprintf("%d doublons exacts retirés ; %d séquences identiques à taxonomies différentes %s",
                    sum(dup_exact), n_conf, if (identical(p$conflicts, "drop")) "retirées" else "ramenées à leur ancêtre commun"))
  } else message(sum(dup_exact), " doublons exacts retirés")
  if (!length(sq)) stop("Aucune séquence valide après nettoyage.")

  # --- cohérence de la hiérarchie : un même nom sous plusieurs parents
  incons <- do.call(rbind, lapply(2:ncol(T), function(k) {
    d <- unique(data.frame(nom = T[, k], parent = T[, k - 1], stringsAsFactors = FALSE))
    d <- d[!is.na(d$nom), ]
    multi <- names(which(table(d$nom) > 1))
    if (!length(multi)) return(NULL)
    data.frame(rang = colnames(T)[k], nom = multi,
               parents = vapply(multi, function(x) paste(unique(d$parent[d$nom == x]), collapse = " | "), ""),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  if (!is.null(incons) && nrow(incons))
    message(nrow(incons), " noms rattachés à plusieurs parents (voir le rapport)")

  # --- fichiers
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  name <- gsub("[^A-Za-z0-9_.-]", "_", p$name %||% "ma_base")
  has_sp <- "Species" %in% levels
  train_cols <- if (has_sp && !isTRUE(p$species_in_train)) setdiff(levels, "Species") else levels
  tax_str <- apply(T[, train_cols, drop = FALSE], 1, function(v) paste0(paste(v[!is.na(v)], collapse = ";"), ";"))
  f_train <- file.path(outdir, paste0(name, "_assignTaxonomy.fa.gz"))
  con <- gzfile(f_train, "w"); writeLines(paste0(">", tax_str, "\n", sq), con); close(con)
  files <- list(train = f_train)
  message("Écrit : ", f_train)

  if (has_sp && "Genus" %in% levels) {
    sp <- T[, "Species"]; ge <- T[, "Genus"]
    ok <- !is.na(sp) & !is.na(ge)
    binom <- ifelse(grepl(" ", sp), sp, paste(ge, sp))
    binom <- vapply(strsplit(binom, " "), function(w) paste(utils::head(w, 2), collapse = " "), "")
    if (any(ok)) {
      f_sp <- file.path(outdir, paste0(name, "_assignSpecies.fa.gz"))
      con <- gzfile(f_sp, "w"); writeLines(paste0(">", ids[ok], " ", binom[ok], "\n", sq[ok]), con); close(con)
      files$species <- f_sp
      message("Écrit : ", f_sp, " (", sum(ok), " séquences avec genre + espèce)")
    }
  }

  problems <- NULL
  if (isTRUE(p$idtaxa)) {
    if (!requireNamespace("DECIPHER", quietly = TRUE)) stop("DECIPHER n'est pas installé : IDTAXA indisponible.")
    message("Entraînement IDTAXA (LearnTaxa)… cela peut prendre plusieurs minutes")
    tx_idt <- paste0("Root;", apply(T, 1, function(v) paste(v[!is.na(v)], collapse = ";")))
    trainingSet <- DECIPHER::LearnTaxa(train = Biostrings::DNAStringSet(sq), taxonomy = tx_idt, verbose = FALSE)
    f_idt <- file.path(outdir, paste0(name, "_IDTAXA.RData"))
    save(trainingSet, file = f_idt)
    files$idtaxa <- f_idt
    problems <- length(trainingSet$problemSequences$Index %||% integer(0))
    message("Écrit : ", f_idt, if (problems) paste0(" (", problems, " séquences signalées comme problématiques par LearnTaxa)"))
  }

  assigned <- colSums(!is.na(T))
  report <- data.frame(
    Etape = c("Séquences lues", "Caractères invalides", "Trop courtes", "Sans taxonomie",
              "Doublons exacts", "Conflits de taxonomie", "Séquences finales"),
    Nombre = c(n_in, sum(bad), sum(short & !bad), sum(notax & !bad & !short),
               sum(dup_exact), n_conf, length(sq)))
  list(files = files, report = report,
       ranks = data.frame(Rang = names(assigned), Sequences_assignees = as.vector(assigned),
                          Taxons_distincts = vapply(colnames(T), function(k) length(unique(stats::na.omit(T[, k]))), 0)),
       inconsistencies = incons %||% data.frame(rang = character(0), nom = character(0), parents = character(0)),
       levels = train_cols, idtaxa_problems = problems)
}

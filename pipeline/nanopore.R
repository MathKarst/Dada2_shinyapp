# =============================================================================
#  nanopore.R — amplicons nanopore : QC, orientation, regroupement, consensus
#               et polissage, entièrement en R (fonctionne sous Windows).
#
#  Principe (proche de NGSpeciesID) :
#    1. filtre longueur / qualité moyenne
#    2. recherche des amorces aux deux extrémités, remise dans le même sens, découpe
#    3. par échantillon : regroupement des reads (DECIPHER::Clusterize), séparation
#       des variantes proches (positions où deux allèles co-varient), alignement
#       multiple, consensus « poli » : votes pondérés par la qualité, équilibrés
#       entre brins, correction de la longueur des homopolymères
#    4. option : polissage par Dorado (dorado aligner + dorado polish) si installé
#    5. table échantillons x séquences consensus (fusion des consensus identiques
#       ou quasi identiques entre échantillons), chimères retirées avec DADA2
#
#  Dépendances : Biostrings, ShortRead, DECIPHER ; dada2 (chimères) ;
#                Rsamtools (seulement pour Dorado polish)
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# Qualités d'un ShortReadQ en liste d'entiers (Phred)
.qint <- function(r) methods::as(Biostrings::PhredQuality(Biostrings::quality(Biostrings::quality(r))), "IntegerList")
# Qualité moyenne « vraie » d'un read : moyenne des probabilités d'erreur
.mean_q <- function(r) {
  vapply(.qint(r), function(x) if (length(x)) -10 * log10(mean(10^(-x / 10))) else NA_real_, 0)
}

# -----------------------------------------------------------------------------
# 1. Détection des fichiers et QC
# -----------------------------------------------------------------------------
nano_detect <- function(dir, pattern = "\\.(fastq|fq)(\\.gz)?$", recursive = FALSE) {
  if (!dir.exists(dir)) stop("Dossier introuvable : ", dir)
  f <- sort(list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive))
  if (!length(f)) stop("Aucun FASTQ trouvé dans ", dir)
  # nom d'échantillon : nom du fichier sans extension, ou nom du dossier pour
  # l'organisation MinKNOW (fastq_pass/barcode01/*.fastq.gz)
  base <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(f))
  parent <- basename(dirname(f))
  multi <- duplicated(parent) | duplicated(parent, fromLast = TRUE)
  name <- if (recursive && all(multi | grepl("^barcode", parent))) parent else base
  df <- data.frame(sample = name, file = f, stringsAsFactors = FALSE)
  # plusieurs fichiers par code-barres (MinKNOW) : regroupés par échantillon
  df <- do.call(rbind, lapply(split(df, df$sample), function(d)
    data.frame(sample = d$sample[1], file = paste(d$file, collapse = "|"), n_files = nrow(d))))
  rownames(df) <- NULL
  df
}
.files_of <- function(x) strsplit(x, "|", fixed = TRUE)[[1]]

.read_fq <- function(files, n_max = Inf) {
  rs <- lapply(files, function(f) {
    if (is.finite(n_max)) {
      s <- ShortRead::FastqSampler(f, n = n_max); on.exit(close(s)); ShortRead::yield(s)
    } else ShortRead::readFastq(f)
  })
  if (length(rs) == 1) rs[[1]] else Reduce(ShortRead::append, rs)
}

#' Longueur et qualité moyenne des reads (échantillon de n_max reads par échantillon)
nano_qc <- function(samples, n_max = 5000) {
  do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
    r <- .read_fq(.files_of(samples$file[i]), n_max)
    if (!length(r)) return(NULL)
    data.frame(sample = samples$sample[i], length = Biostrings::width(r),
               meanq = round(.mean_q(r), 1))
  }))
}

# -----------------------------------------------------------------------------
# 2. Amorces : recherche, orientation, découpe
# -----------------------------------------------------------------------------
.rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
# Fonctions d'alignement par paires : déplacées de Biostrings vers pwalign
# (Bioconductor >= 3.19, retirées de Biostrings >= 2.77). On prend pwalign s'il
# est installé et exporte la fonction, sinon Biostrings (anciennes versions).
.pw <- function(f) {
  if (requireNamespace("pwalign", quietly = TRUE) && f %in% getNamespaceExports("pwalign"))
    getExportedValue("pwalign", f) else getExportedValue("Biostrings", f)
}
# matrice calculée à la première utilisation (pas au chargement du fichier)
.SM_cache <- new.env()
.SM <- function() {
  if (is.null(.SM_cache$m))
    .SM_cache$m <- .pw("nucleotideSubstitutionMatrix")(match = 1, mismatch = -1, baseOnly = FALSE)
  .SM_cache$m
}

# Meilleure position d'une liste d'amorces dans des fenêtres de séquence.
# Score maximal = nombre de bases non ambiguës ; une erreur coûte environ 2 points.
.find_primer <- function(windows, primers, max_err) {
  best <- list(ok = rep(FALSE, length(windows)), start = rep(NA_integer_, length(windows)),
               end = rep(NA_integer_, length(windows)), score = rep(-Inf, length(windows)))
  pa <- .pw("pairwiseAlignment")
  for (pr in primers) {
    aln <- pa(windows, Biostrings::DNAString(pr), type = "local-global", substitutionMatrix = .SM(),
              gapOpening = 1, gapExtension = 1)
    sc <- .pw("score")(aln)
    maxsc <- sum(strsplit(pr, "")[[1]] %in% c("A", "C", "G", "T"))
    ok <- sc >= maxsc - 2 * ceiling(max_err * nchar(pr))
    better <- ok & sc > best$score
    pat <- .pw("pattern")(aln)
    best$ok[better] <- TRUE
    best$start[better] <- Biostrings::start(pat)[better]
    best$end[better] <- Biostrings::end(pat)[better]
    best$score[better] <- sc[better]
  }
  best
}

#' Oriente et découpe des reads (ShortReadQ) entre les amorces.
#' @return list(reads = ShortReadQ orientés et découpés, strand = "+"/"-")
nano_orient_trim <- function(reads, primersF, primersR, max_err = 0.2, window = 150,
                             require_both = TRUE) {
  pF <- primersF; pR <- primersR
  L <- Biostrings::width(reads)
  w <- pmin(window, L)
  s5 <- Biostrings::subseq(ShortRead::sread(reads), 1, w)
  s3 <- Biostrings::subseq(ShortRead::sread(reads), L - w + 1, L)
  F5 <- .find_primer(s5, pF, max_err); R3 <- .find_primer(s3, .rc(pR), max_err)   # sens direct
  R5 <- .find_primer(s5, pR, max_err); F3 <- .find_primer(s3, .rc(pF), max_err)   # sens inverse
  nf <- F5$ok + R3$ok; nr <- R5$ok + F3$ok
  sf <- pmax(F5$score, 0) + pmax(R3$score, 0); sr <- pmax(R5$score, 0) + pmax(F3$score, 0)
  fwd <- nf > nr | (nf == nr & sf >= sr)
  n_ok <- ifelse(fwd, nf, nr)
  keep <- if (require_both) n_ok == 2 else n_ok >= 1
  st <- ifelse(fwd, ifelse(F5$ok, F5$end + 1L, 1L), ifelse(R5$ok, R5$end + 1L, 1L))
  en <- ifelse(fwd, ifelse(R3$ok, L - w + R3$start - 1L, L), ifelse(F3$ok, L - w + F3$start - 1L, L))
  keep <- keep & !is.na(st) & !is.na(en) & en - st >= 20
  r <- reads[keep]; st <- st[keep]; en <- en[keep]; fw <- fwd[keep]
  r <- ShortRead::narrow(r, start = st, end = en)
  if (any(!fw)) {
    sq <- ShortRead::sread(r)
    qq <- Biostrings::quality(Biostrings::quality(r))           # BStringSet
    sq[!fw] <- Biostrings::reverseComplement(sq[!fw])
    qq[!fw] <- Biostrings::reverse(qq[!fw])
    r <- ShortRead::ShortReadQ(sq, ShortRead::FastqQuality(qq), ShortRead::id(r))
  }
  list(reads = r, strand = ifelse(fw, "+", "-"), n_in = length(reads), n_out = length(r))
}

#' Filtre + orientation pour un échantillon ; écrit le FASTQ obtenu
#' @param p list(min_len, max_len, min_q, primersF, primersR, primer_max_err,
#'               window, require_both, max_reads)
nano_prepare_sample <- function(sample, files, outdir, p) {
  r <- .read_fq(.files_of(files))
  n_in <- length(r)
  L <- Biostrings::width(r); q <- .mean_q(r)
  ok <- L >= p$min_len & L <= p$max_len & q >= p$min_q
  r <- r[ok]
  n_filt <- length(r)
  strand <- rep("+", length(r))
  pF <- parse_primers(p$primersF); pR <- parse_primers(p$primersR)
  if (length(pF) || length(pR)) {
    if (!length(pF) || !length(pR)) stop("Indique les amorces forward ET reverse.")
    o <- nano_orient_trim(r, pF, pR, p$primer_max_err %||% 0.2, p$window %||% 150, isTRUE(p$require_both))
    r <- o$reads; strand <- o$strand
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(outdir, paste0(sample, "_prep.fastq.gz"))
  unlink(out)
  ids <- paste0(sample, "_", seq_along(r), ";strand=", strand)
  r <- ShortRead::ShortReadQ(ShortRead::sread(r), Biostrings::quality(r), Biostrings::BStringSet(ids))
  ShortRead::writeFastq(r, out, compress = TRUE)
  message(sprintf("%s : %d reads -> %d après filtre longueur/qualité -> %d orientés et découpés",
                  sample, n_in, n_filt, length(r)))
  data.frame(sample = sample, input = n_in, filtered = n_filt, oriented = length(r),
             prep_file = out, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# 3. Consensus « poli » à partir d'un alignement multiple
# -----------------------------------------------------------------------------
#' @param msa DNAStringSet aligné (avec "-") ; @param quals liste de vecteurs de
#'   qualité (sans gaps), un par read ; @param strand "+"/"-" par read
polish_msa <- function(msa, quals, strand, hp_correct = TRUE, min_strand_reads = 3) {
  m <- as.matrix(msa)
  n <- nrow(m); nc <- ncol(m)
  # qualités alignées sur les colonnes ; extrémités non couvertes = NA
  W <- matrix(NA_real_, n, nc)
  for (i in seq_len(n)) {
    pos <- which(m[i, ] != "-")
    if (!length(pos)) next
    W[i, pos] <- pmin(quals[[i]][seq_along(pos)], 40) / 40
    inner <- seq(min(pos), max(pos))
    gaps <- setdiff(inner, pos)
    W[i, gaps] <- mean(W[i, pos], na.rm = TRUE)     # poids d'un gap interne
  }
  bases <- c("A", "C", "G", "T", "-")
  strands <- unique(strand)
  usable <- vapply(strands, function(s) sum(strand == s) >= min_strand_reads, TRUE)
  if (!any(usable)) usable[] <- TRUE
  strands <- strands[usable]
  # votes par colonne, normalisés dans chaque brin puis moyennés entre brins
  vote <- matrix(0, length(bases), nc, dimnames = list(bases, NULL))
  for (s in strands) {
    rows <- which(strand == s)
    vs <- sapply(bases, function(b) colSums((m[rows, , drop = FALSE] == b) * W[rows, , drop = FALSE], na.rm = TRUE))
    vs <- t(vs)
    tot <- colSums(vs); tot[tot == 0] <- 1
    vote <- vote + sweep(vs, 2, tot, "/")
  }
  win <- bases[apply(vote, 2, which.max)]
  win[colSums(vote) == 0] <- "-"
  # extrémités : une position couverte par moins de la moitié des reads (bases
  # résiduelles d'amorce ou d'adaptateur sur quelques reads) n'entre pas dans le consensus
  win[colMeans(!is.na(W)) < 0.5] <- "-"
  cols <- which(win != "-")
  cons <- win[cols]
  if (!hp_correct || length(cons) < 3) return(paste(cons, collapse = ""))

  # correction des homopolymères : longueur la plus probable d'après les reads
  r <- rle(cons)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  new_len <- r$lengths
  for (k in which(r$lengths >= 3)) {
    b <- r$values[k]
    c1 <- cols[starts[k]]; c2 <- cols[ends[k]]
    lo <- if (starts[k] > 1) cols[starts[k] - 1] + 1 else 1
    hi <- if (ends[k] < length(cols)) cols[ends[k] + 1] - 1 else nc
    region <- lo:hi
    covered <- rowSums(!is.na(W[, c(c1, c2), drop = FALSE])) == 2
    if (sum(covered) < 3) next
    cnt <- rowSums(m[covered, region, drop = FALSE] == b)
    st_c <- strand[covered]
    dist <- 0
    lens <- sort(unique(cnt))
    for (s in unique(st_c)) {
      tab <- table(factor(cnt[st_c == s], levels = lens))
      dist <- dist + as.numeric(tab) / sum(tab)
    }
    new_len[k] <- lens[which.max(dist)]
  }
  paste(rep(r$values, new_len), collapse = "")
}

# Positions informatives : deux allèles fréquents qui ne sont ni des gaps ni
# dans un homopolymère (où les erreurs nanopore se concentrent)
.informative_cols <- function(m, min_freq) {
  nc <- ncol(m)
  cons <- apply(m, 2, function(x) { t <- table(x[x != "-"]); if (length(t)) names(t)[which.max(t)] else "-" })
  hp <- rep(FALSE, nc)
  run <- rle(cons[cons != "-"])
  hp_nogap <- rep(run$lengths >= 3, run$lengths)
  hp[which(cons != "-")] <- hp_nogap
  which(vapply(seq_len(nc), function(j) {
    if (hp[j]) return(FALSE)
    x <- m[, j]
    if (mean(x == "-") > 0.2) return(FALSE)
    t <- sort(table(x[x != "-"]), decreasing = TRUE)
    length(t) >= 2 && t[2] / sum(t) >= min_freq
  }, TRUE))
}

#' Sépare un groupe de reads alignés en sous-groupes si des variantes co-varient
.split_cluster <- function(m, min_freq, min_sites, min_reads, depth = 0, max_depth = 3) {
  idx <- seq_len(nrow(m))
  if (depth >= max_depth || nrow(m) < 2 * min_reads) return(list(idx))
  inf <- .informative_cols(m, min_freq)
  if (length(inf) < min_sites) return(list(idx))
  sub <- m[, inf, drop = FALSE]
  d <- stats::as.dist(outer(idx, idx, Vectorize(function(a, b) mean(sub[a, ] != sub[b, ]))))
  g <- stats::cutree(stats::hclust(d, "average"), 2)
  if (min(table(g)) < min_reads) return(list(idx))
  # les deux sous-groupes doivent différer sur au moins min_sites positions
  maj <- function(rows) apply(sub[rows, , drop = FALSE], 2, function(x) names(which.max(table(x))))
  if (sum(maj(which(g == 1)) != maj(which(g == 2))) < min_sites) return(list(idx))
  out <- list()
  for (k in 1:2) {
    rows <- which(g == k)
    for (s in .split_cluster(m[rows, , drop = FALSE], min_freq, min_sites, min_reads, depth + 1, max_depth))
      out[[length(out) + 1]] <- rows[s]
  }
  out
}

# -----------------------------------------------------------------------------
# 4. Dorado polish (optionnel, nécessite dorado et Rsamtools)
# -----------------------------------------------------------------------------
dorado_polish <- function(reads, draft, dorado = "dorado", model = "", threads = 1, workdir = tempfile("dorado_")) {
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  fq <- file.path(workdir, "reads.fastq"); fa <- file.path(workdir, "draft.fasta")
  ShortRead::writeFastq(reads, fq, compress = FALSE)
  writeLines(c(">draft", draft), fa)
  bam <- file.path(workdir, "aln.bam")
  st <- system2(dorado, c("aligner", "-t", threads, shQuote(fa), shQuote(fq)), stdout = bam, stderr = FALSE)
  if (!identical(as.integer(st), 0L)) stop("dorado aligner a échoué.")
  sorted <- Rsamtools::sortBam(bam, file.path(workdir, "aln.sorted"))
  Rsamtools::indexBam(sorted)
  out <- file.path(workdir, "polished.fasta")
  args <- c("polish", shQuote(sorted), shQuote(fa), "--device", "cpu", "--threads", threads,
            if (nzchar(model)) c("--model-override", shQuote(model)))
  st <- system2(dorado, args, stdout = out, stderr = FALSE)
  if (!identical(as.integer(st), 0L) || !file.exists(out)) stop("dorado polish a échoué.")
  s <- Biostrings::readDNAStringSet(out)
  if (!length(s)) stop("dorado polish n'a produit aucune séquence.")
  as.character(s[[1]])
}

# -----------------------------------------------------------------------------
# 5. Consensus par échantillon
# -----------------------------------------------------------------------------
#' @param p list(mode = "metabarcoding"|"barcoding", cluster_cutoff, min_reads,
#'   min_frac, max_reads, max_msa, split, split_freq, split_min_sites,
#'   hp_correct, polisher = "R"|"dorado", dorado, dorado_model, threads, seed)
nano_consensus_sample <- function(sample, prep_file, p) {
  if (!requireNamespace("DECIPHER", quietly = TRUE)) stop("Le package DECIPHER est nécessaire.")
  r <- ShortRead::readFastq(prep_file)
  n_all <- length(r)
  if (n_all < (p$min_reads %||% 5)) {
    message(sample, " : trop peu de reads (", n_all, ")")
    return(list(consensus = NULL, clusters = NULL))
  }
  set.seed(p$seed %||% 1)
  if (n_all > p$max_reads) r <- r[sort(sample.int(n_all, p$max_reads))]
  scale <- n_all / length(r)
  sq <- ShortRead::sread(r)
  strand <- sub(".*strand=", "", as.character(ShortRead::id(r)))
  quals <- .qint(r)

  cl <- suppressWarnings(DECIPHER::Clusterize(sq, cutoff = p$cluster_cutoff, processors = p$threads %||% 1, verbose = FALSE))
  grp <- cl$cluster
  sizes <- sort(table(grp), decreasing = TRUE)
  min_n <- if (identical(p$mode, "barcoding")) max(p$min_reads, ceiling(p$min_frac * length(r))) else p$min_reads
  keep <- names(sizes)[sizes >= min_n]
  message(sprintf("%s : %d reads, %d groupes, %d groupe(s) d'au moins %d reads", sample, length(r),
                  length(sizes), length(keep), min_n))

  res <- list(); diag <- list()
  for (g in keep) {
    idx <- which(grp == as.integer(g))
    sub_idx <- if (length(idx) > p$max_msa) sort(sample(idx, p$max_msa)) else idx
    msa <- DECIPHER::AlignSeqs(sq[sub_idx], processors = p$threads %||% 1, verbose = FALSE)
    m <- as.matrix(msa)
    parts <- if (isTRUE(p$split)) .split_cluster(m, p$split_freq, p$split_min_sites, max(3, p$min_reads)) else list(seq_len(nrow(m)))
    drafts <- vapply(parts, function(rows)
      polish_msa(msa[rows], quals[sub_idx[rows]], strand[sub_idx[rows]], isTRUE(p$hp_correct)), "")
    # répartit tous les reads du groupe entre les sous-groupes (si séparation)
    assign <- if (length(parts) == 1) rep(1L, length(idx)) else {
      pa <- .pw("pairwiseAlignment")
      sc <- sapply(drafts, function(d) .pw("score")(pa(sq[idx], Biostrings::DNAString(d), type = "global",
                                                         substitutionMatrix = .SM(), gapOpening = 2, gapExtension = 1)))
      max.col(matrix(sc, ncol = length(drafts)), ties.method = "first")
    }
    for (k in seq_along(parts)) {
      members <- idx[assign == k]
      if (length(members) < p$min_reads) next
      cons <- drafts[k]
      if (identical(p$polisher, "dorado")) {
        cons <- tryCatch(dorado_polish(r[members][seq_len(min(length(members), p$max_msa))], cons,
                                       p$dorado %||% "dorado", p$dorado_model %||% "", p$threads %||% 1),
                         error = function(e) { message("Dorado polish : ", conditionMessage(e), " (consensus R gardé)"); cons })
      }
      id <- paste0(sample, "_c", g, if (length(parts) > 1) letters[k])
      res[[id]] <- data.frame(sample = sample, consensus_id = id, reads = round(length(members) * scale),
                              reads_used = length(members), strand_plus = round(mean(strand[members] == "+"), 2),
                              length = nchar(cons), sequence = cons, stringsAsFactors = FALSE)
    }
    if (length(parts) > 1) message(sprintf("  groupe %s séparé en %d variantes", g, length(parts)))
  }
  cons <- if (length(res)) do.call(rbind, res) else NULL
  if (!is.null(cons)) rownames(cons) <- NULL
  list(consensus = cons,
       clusters = data.frame(sample = sample, cluster = names(sizes), size = as.vector(sizes),
                             kept = names(sizes) %in% keep, stringsAsFactors = FALSE),
       n_reads = n_all, scale = scale)
}

# -----------------------------------------------------------------------------
# 6. Table échantillons x séquences
# -----------------------------------------------------------------------------
#' Fusionne les consensus identiques ou quasi identiques entre échantillons
#' @param merge_cutoff distance maximale pour fusionner (0 = identiques seulement)
#' @param chim_method "pooled" (défaut : un parent présent dans un autre échantillon suffit)
#'   ou "consensus" (décision échantillon par échantillon, comme dans DADA2)
#' @param allow_one_off tolère une différence avec le modèle à deux parents : utile car
#'   le consensus d'une chimère rare repose sur peu de reads et garde quelques erreurs
nano_build_table <- function(cons, merge_cutoff = 0.01, remove_chimeras = TRUE, threads = 1,
                             chim_method = "pooled", allow_one_off = TRUE) {
  if (is.null(cons) || !nrow(cons)) stop("Aucun consensus à assembler.")
  u <- unique(cons$sequence)
  rep_of <- setNames(u, u)
  if (merge_cutoff > 0 && length(u) > 1) {
    cl <- suppressWarnings(DECIPHER::Clusterize(Biostrings::DNAStringSet(u), cutoff = merge_cutoff,
                                                processors = threads, verbose = FALSE))$cluster
    tot <- tapply(cons$reads, cons$sequence, sum)[u]
    for (g in unique(cl)) {
      mem <- u[cl == g]
      rep_of[mem] <- mem[which.max(tot[mem])]
    }
  }
  cons$representative <- unname(rep_of[cons$sequence])
  tab <- tapply(cons$reads, list(cons$sample, cons$representative), sum)
  tab[is.na(tab)] <- 0
  tab <- tab[, order(-colSums(tab)), drop = FALSE]
  storage.mode(tab) <- "integer"
  message(sprintf("Table : %d échantillons x %d séquences (%d consensus fusionnés)",
                  nrow(tab), ncol(tab), length(u) - ncol(tab)))
  chim <- character(0)
  if (remove_chimeras && ncol(tab) > 2) {
    nc <- dada2::removeBimeraDenovo(tab, method = chim_method, allowOneOff = allow_one_off, verbose = FALSE)
    chim <- setdiff(colnames(tab), colnames(nc))
    if (length(chim)) message(length(chim), " consensus chimérique(s) retiré(s)")
    tab <- nc
  }
  list(table = tab, consensus = cons, chimeras = chim)
}

# -----------------------------------------------------------------------------
# Traitement de tous les échantillons (appelé par l'appli en arrière-plan)
# -----------------------------------------------------------------------------
nano_prepare_all <- function(samples, outdir, p) {
  d <- file.path(outdir, "prepared")
  do.call(rbind, lapply(seq_len(nrow(samples)), function(i)
    nano_prepare_sample(samples$sample[i], samples$file[i], d, p)))
}

nano_consensus_all <- function(prep, p) {
  res <- lapply(seq_len(nrow(prep)), function(i) {
    if (prep$oriented[i] < 1) return(NULL)
    nano_consensus_sample(prep$sample[i], prep$prep_file[i], p)
  })
  res <- Filter(Negate(is.null), res)
  cons <- do.call(rbind, lapply(res, `[[`, "consensus"))
  clus <- do.call(rbind, lapply(res, `[[`, "clusters"))
  message(sprintf("Consensus : %d séquences sur %d échantillons", if (is.null(cons)) 0 else nrow(cons), length(res)))
  list(consensus = cons, clusters = clus)
}

#' Script R reproductible (source les fonctions du dépôt)
nano_script <- function(st, repo_dir) {
  P <- st$params
  dp <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  L <- c("# Pipeline nanopore — script généré par l'appli Nanopore Shiny",
         paste0("# ", format(Sys.time(), "%Y-%m-%d %H:%M"), " ; ", R.version.string),
         paste0("repo <- ", dp(repo_dir), "  # dossier contenant pipeline/"),
         "for (f in c(\"pipeline.R\", \"curation.R\", \"nanopore.R\")) source(file.path(repo, \"pipeline\", f))",
         paste0("outdir <- ", dp(P$data$outdir)),
         paste0("samples <- data.frame(sample = ", dp(st$samples$sample), ",\n                      file = ", dp(st$samples$file), ")"))
  if (!is.null(P$prep)) L <- c(L, paste0("prep <- nano_prepare_all(samples, outdir, ", dp(P$prep), ")"))
  if (!is.null(P$consensus)) L <- c(L, paste0("cons <- nano_consensus_all(prep, ", dp(P$consensus), ")"))
  if (!is.null(P$table)) L <- c(L, paste0("tb <- nano_build_table(cons$consensus, merge_cutoff = ", P$table$merge_cutoff,
                                          ", remove_chimeras = ", P$table$remove_chimeras, ", chim_method = ", dp(P$table$chim_method),
                                          ", allow_one_off = ", P$table$allow_one_off, ")"),
                                   "tab <- tb$table")
  if (!is.null(P$taxonomy)) L <- c(L, paste0("taxa <- run_taxonomy(colnames(tab), ", dp(P$taxonomy), ")$taxa"))
  if (!is.null(P$taxfilter)) L <- c(L, paste0("tab <- run_taxfilter(tab, taxa, ", dp(P$taxfilter), ")$final_tax"))
  L <- c(L, "ids <- paste0(\"SEQ\", match(colnames(tab), colnames(tb$table)))  # mêmes identifiants que l'appli",
         "writeLines(paste0(\">\", ids, \"\\n\", colnames(tab)), file.path(outdir, \"sequences.fasta\"))",
         "write.csv(data.frame(ID = ids, sequence = colnames(tab), t(tab), check.names = FALSE), file.path(outdir, \"table.csv\"), row.names = FALSE)")
  paste(L, collapse = "\n")
}

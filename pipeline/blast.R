# =============================================================================
#  blast.R — BLAST à distance sur les serveurs du NCBI (« Common URL API »)
#
#  Les séquences sont envoyées en une seule requête (multi-FASTA), puis le
#  statut est interrogé jusqu'à ce que le résultat soit prêt, et le résultat est
#  récupéré au format JSON (JSON2_S) et mis en tableau.
#
#  Règles d'usage du NCBI respectées : au moins 10 s entre deux requêtes, pas
#  plus d'une interrogation par minute pour une même recherche, paramètres
#  TOOL et EMAIL renseignés, nombre de séquences limité par recherche.
#
#  Dépendances : curl, jsonlite
#  Test hors ligne : si la variable d'environnement DADA2SHINY_BLAST_MOCK
#  contient le chemin d'un fichier JSON2_S, aucune requête n'est envoyée et ce
#  fichier sert de réponse (utilisé pour les tests automatiques).
# =============================================================================

NCBI_BLAST_URL <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
BLAST_MAX_SEQS <- 100

.blast_env <- new.env()
.blast_env$last <- 0

ncbi_request <- function(params, post = FALSE) {
  mock <- Sys.getenv("DADA2SHINY_BLAST_MOCK")
  if (nzchar(mock)) {
    return(switch(params$CMD,
      Put = "<!--QBlastInfoBegin\n    RID = MOCKRID123\n    RTOE = 1\nQBlastInfoEnd\n-->",
      Get = if (identical(params$FORMAT_OBJECT, "SearchInfo"))
              "QBlastInfoBegin\n\tStatus=READY\nQBlastInfoEnd\nQBlastInfoBegin\n\tThereAreHits=yes\nQBlastInfoEnd"
            else paste(readLines(mock, warn = FALSE), collapse = "\n")))
  }
  if (!requireNamespace("curl", quietly = TRUE))
    stop("Le package curl est nécessaire pour BLAST (install.packages('curl')).")
  # au moins 10 s entre deux requêtes au NCBI
  wait <- 10 - (as.numeric(Sys.time()) - .blast_env$last)
  if (wait > 0) Sys.sleep(wait)
  .blast_env$last <- as.numeric(Sys.time())

  params <- params[!vapply(params, function(x) is.null(x) || !nzchar(x), TRUE)]
  enc <- paste(names(params), vapply(params, function(x) curl::curl_escape(as.character(x)), ""),
               sep = "=", collapse = "&")
  h <- curl::new_handle(useragent = "dada2shiny (R, curl)")
  url <- NCBI_BLAST_URL
  if (post) {
    curl::handle_setopt(h, copypostfields = enc)
    curl::handle_setheaders(h, "Content-Type" = "application/x-www-form-urlencoded")
  } else url <- paste0(url, "?", enc)
  res <- curl::curl_fetch_memory(url, handle = h)
  if (res$status_code >= 400) stop("Le serveur BLAST du NCBI a répondu avec l'erreur HTTP ", res$status_code, ".")
  rawToChar(res$content)
}

blast_entrez_query <- function(entrez = "", exclude_uncultured = TRUE) {
  eq <- trimws(entrez %||% "")
  if (!isTRUE(exclude_uncultured)) return(eq)
  ex <- 'NOT uncultured[Title] NOT "environmental samples"[Organism]'
  if (nzchar(eq)) paste0("(", eq, ") ", ex) else paste("all[filter]", ex)
}

#' Met en tableau un résultat JSON2_S (un rapport par séquence)
parse_blast_json2 <- function(txt, query_ids = NULL) {
  j <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  reps <- j$BlastOutput2
  if (is.null(reps)) stop("Réponse du NCBI inattendue (pas de BlastOutput2).")
  if (!is.null(reps$report)) reps <- list(reps)           # une seule séquence
  get <- function(x, k, default = NA) if (is.null(x[[k]])) default else x[[k]]
  rows <- list(); queries <- character(0)
  for (i in seq_along(reps)) {
    s <- reps[[i]]$report$results$search
    q <- get(s, "query_title", "")
    if (!nzchar(q) && !is.null(query_ids) && i <= length(query_ids)) q <- query_ids[i]
    q <- sub("\\s.*$", "", q)                               # garde l'identifiant (ASV12)
    queries <- c(queries, q)
    qlen <- get(s, "query_len", NA)
    for (h in s$hits %||% list()) {
      d <- h$description[[1]]; hs <- h$hsps[[1]]
      rows[[length(rows) + 1]] <- data.frame(
        ASV = q, rang = get(h, "num"),
        espece = get(d, "sciname"), accession = get(d, "accession"),
        identite = round(100 * get(hs, "identity") / get(hs, "align_len"), 2),
        couverture = round(100 * (abs(get(hs, "query_to") - get(hs, "query_from")) + 1) / qlen, 1),
        longueur_align = get(hs, "align_len"), gaps = get(hs, "gaps", 0),
        evalue = get(hs, "evalue"), bitscore = get(hs, "bit_score"),
        taxid = get(d, "taxid"), titre = get(d, "title"),
        stringsAsFactors = FALSE)
    }
  }
  hits <- if (length(rows)) do.call(rbind, rows) else
    data.frame(ASV = character(0), rang = integer(0), espece = character(0), accession = character(0),
               identite = numeric(0), couverture = numeric(0), longueur_align = integer(0),
               gaps = integer(0), evalue = numeric(0), bitscore = numeric(0), taxid = integer(0),
               titre = character(0))
  list(hits = hits, queries = unique(queries))
}

#' Synthèse : meilleur hit par ASV et espèces ex aequo (même score)
blast_summary <- function(hits, asv_ids) {
  do.call(rbind, lapply(asv_ids, function(a) {
    h <- hits[hits$ASV == a, , drop = FALSE]
    if (!nrow(h)) return(data.frame(ASV = a, meilleur_hit = "Aucun hit", identite = NA, couverture = NA,
                                    evalue = NA, accession = NA, ex_aequo = "", n_hits = 0))
    best <- h[which.max(h$bitscore), ]
    tie <- h[h$bitscore >= max(h$bitscore) - 1e-6, , drop = FALSE]
    sp <- unique(stats::na.omit(tie$espece))
    data.frame(ASV = a, meilleur_hit = best$espece, identite = best$identite,
               couverture = best$couverture, evalue = best$evalue, accession = best$accession,
               ex_aequo = if (length(sp) > 1) paste0(length(sp), " espèces : ",
                                                     paste(utils::head(sp, 6), collapse = " / "),
                                                     if (length(sp) > 6) " / …") else "",
               n_hits = nrow(h))
  }))
}

#' Lance un BLAST NCBI et attend le résultat
#' @param seqs vecteur de séquences nommé par les identifiants d'ASV
#' @param p list(database, megablast, max_hits, entrez, exclude_uncultured, email, timeout_min)
run_ncbi_blast <- function(seqs, p) {
  if (!length(seqs)) stop("Aucune séquence sélectionnée.")
  if (length(seqs) > BLAST_MAX_SEQS)
    stop("Trop de séquences (", length(seqs), ") : ", BLAST_MAX_SEQS, " au maximum par recherche.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Le package jsonlite est nécessaire.")
  fasta <- paste0(">", names(seqs), "\n", seqs, collapse = "\n")
  eq <- blast_entrez_query(p$entrez, p$exclude_uncultured)

  message(sprintf("BLAST NCBI : envoi de %d séquence(s) — base %s, %s, %d hits max%s",
                  length(seqs), p$database, if (isTRUE(p$megablast)) "megablast" else "blastn",
                  p$max_hits, if (nzchar(eq)) paste0(", filtre « ", eq, " »") else ""))
  put <- ncbi_request(list(
    CMD = "Put", PROGRAM = "blastn", DATABASE = p$database, QUERY = fasta,
    MEGABLAST = if (isTRUE(p$megablast)) "on", HITLIST_SIZE = p$max_hits,
    ENTREZ_QUERY = eq, TOOL = "dada2shiny", EMAIL = p$email), post = TRUE)
  rid <- regmatches(put, regexec("RID = ([A-Z0-9]+)", put))[[1]][2]
  rtoe <- suppressWarnings(as.numeric(regmatches(put, regexec("RTOE = ([0-9]+)", put))[[1]][2]))
  if (is.na(rid)) {
    err <- regmatches(put, regexec('class="error"[^>]*>([^<]+)', put))[[1]][2]
    stop("Le NCBI n'a pas accepté la recherche", if (!is.na(err)) paste0(" : ", trimws(err)), ".")
  }
  message("BLAST NCBI : RID ", rid, " (estimation du NCBI : ", rtoe %||% "?", " s). ",
          "Résultats consultables aussi sur ", NCBI_BLAST_URL, "?CMD=Get&RID=", rid)

  mock <- nzchar(Sys.getenv("DADA2SHINY_BLAST_MOCK"))
  t0 <- Sys.time(); timeout <- 60 * (p$timeout_min %||% 30)
  # première interrogation après l'estimation du NCBI, puis au plus une par minute
  Sys.sleep(if (mock) 0 else max(60, rtoe %||% 60))
  repeat {
    info <- ncbi_request(list(CMD = "Get", FORMAT_OBJECT = "SearchInfo", RID = rid))
    status <- regmatches(info, regexec("Status=([A-Z]+)", info))[[1]][2]
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
    if (identical(status, "READY")) break
    if (status %in% c("FAILED", "UNKNOWN"))
      stop("La recherche BLAST a échoué côté NCBI (statut ", status, ", RID ", rid, ").")
    if (el > timeout) stop("Délai dépassé (", p$timeout_min, " min). La recherche continue peut-être sur le site du NCBI : RID ", rid, ".")
    message(sprintf("BLAST NCBI : en cours (%d min %02d s)…", el %/% 60, el %% 60))
    Sys.sleep(60)
  }
  if (grepl("ThereAreHits=no", info)) {
    message("BLAST NCBI : aucun hit.")
    res <- list(hits = parse_blast_json2('{"BlastOutput2": []}')$hits)
  } else {
    txt <- ncbi_request(list(CMD = "Get", RID = rid, FORMAT_TYPE = "JSON2_S"))
    res <- parse_blast_json2(txt, names(seqs))
  }
  message(sprintf("BLAST NCBI : %d hits pour %d séquence(s).", nrow(res$hits), length(unique(res$hits$ASV))))
  list(hits = res$hits, rid = rid, asv = names(seqs), time = Sys.time())
}

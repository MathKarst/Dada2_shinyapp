# =============================================================================
#  viz.R — abondances relatives et Krona, indépendants de Shiny
#
#  Conventions : `counts` = matrice échantillons x ASV (comptages),
#                `tax`    = matrice ASV x rangs (NA = non assigné) ou NULL,
#                `asv_ids`= identifiants des ASV (même ordre que les colonnes).
# =============================================================================

OTHERS <- "Autres"
UNASSIGNED <- "Non assigné"

# Palette qualitative (Tableau 20 sans les gris, réservés à Autres / Non assigné)
PAL20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D",
           "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A", "#D37295", "#FABFD2",
           "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6", "#2F4B7C", "#A05195")

taxa_palette <- function(levels) {
  core <- setdiff(levels, c(OTHERS, UNASSIGNED))
  n <- length(core)
  cols <- if (n <= length(PAL20)) PAL20[seq_len(n)] else
    c(PAL20, grDevices::hcl.colors(n - length(PAL20), "Dark 3"))
  c(setNames(cols, core), setNames("#D9D9D9", OTHERS), setNames("#7F7F7F", UNASSIGNED))
}

#' Matrice des abondances relatives (%) : taxons (lignes) x échantillons (colonnes)
#' @param level   "ASV" ou un nom de rang (colonne de `tax`)
#' @param label_rank pour level = "ASV" : rang ajouté entre parenthèses, ex. "ASV3 (Navicula)"
#' @param groups  vecteur (un par échantillon) pour moyenner les échantillons par groupe
#' @param top_n   nombre de taxons affichés ; le reste est regroupé dans « Autres »
abund_matrix <- function(counts, tax = NULL, asv_ids = colnames(counts), level = "ASV",
                         label_rank = NULL, groups = NULL, top_n = 15) {
  if (level == "ASV") {
    lab <- asv_ids
    if (!is.null(tax) && !is.null(label_rank) && label_rank %in% colnames(tax)) {
      g <- tax[, label_rank]
      lab <- ifelse(is.na(g), lab, paste0(lab, " (", g, ")"))
    }
  } else {
    if (is.null(tax) || !level %in% colnames(tax)) stop("Rang « ", level, " » indisponible : assigne d'abord la taxonomie.")
    lab <- ifelse(is.na(tax[, level]) | !nzchar(tax[, level]), UNASSIGNED, tax[, level])
  }
  agg <- rowsum(t(counts), lab, reorder = FALSE)                 # taxons x échantillons
  rel <- sweep(agg, 2, pmax(colSums(agg), 1), "/") * 100
  if (!is.null(groups)) {
    g <- as.character(groups); g[is.na(g) | !nzchar(g)] <- "(manquant)"
    s <- rowsum(t(rel), g, reorder = FALSE)                      # groupes x taxons
    rel <- t(s / as.vector(table(g)[rownames(s)]))
  }
  m <- rowMeans(rel)
  cand <- setdiff(names(sort(m, decreasing = TRUE)), UNASSIGNED)
  keep <- utils::head(cand, top_n)
  others <- setdiff(cand, keep)
  out <- rel[keep, , drop = FALSE]
  if (length(others)) {
    out <- rbind(out, colSums(rel[others, , drop = FALSE]))
    rownames(out)[nrow(out)] <- OTHERS
  }
  if (UNASSIGNED %in% rownames(rel)) {
    out <- rbind(out, rel[UNASSIGNED, ])
    rownames(out)[nrow(out)] <- UNASSIGNED
  }
  out
}

# Dissimilarité de Bray-Curtis entre colonnes d'une matrice de proportions
bray_dist <- function(P) {
  P <- sweep(P, 2, pmax(colSums(P), 1e-12), "/")
  n <- ncol(P)
  D <- matrix(0, n, n, dimnames = list(colnames(P), colnames(P)))
  if (n > 1) for (i in seq_len(n - 1)) for (j in (i + 1):n)
    D[i, j] <- D[j, i] <- 1 - sum(pmin(P[, i], P[, j]))
  stats::as.dist(D)
}

#' Ordre des échantillons (colonnes de `full`, matrice taxons x échantillons non tronquée)
sample_order <- function(full, method = c("name", "cluster", "dominant")) {
  method <- match.arg(method)
  s <- colnames(full)
  if (method == "name" || length(s) < 3) return(sort(s))
  if (method == "cluster") {
    hc <- stats::hclust(bray_dist(full), method = "average")
    return(s[hc$order])
  }
  top <- names(sort(rowMeans(full), decreasing = TRUE))[1]
  s[order(-full[top, ])]
}

#' Format long pour ggplot
abund_long <- function(mat, order = colnames(mat)) {
  df <- data.frame(taxon = rep(rownames(mat), times = ncol(mat)),
                   sample = rep(colnames(mat), each = nrow(mat)),
                   pct = as.vector(mat), stringsAsFactors = FALSE)
  df$taxon <- factor(df$taxon, levels = rev(rownames(mat)))
  df$sample <- factor(df$sample, levels = order)
  df
}

#' Graphique d'abondance relative
#' @param facet  data.frame(sample, facet) optionnel
plot_abund <- function(mat, order = colnames(mat), type = c("bar", "heatmap"),
                       facet = NULL, title = NULL, legend_title = "Taxon") {
  type <- match.arg(type)
  df <- abund_long(mat, order)
  if (!is.null(facet)) df$facet <- facet$facet[match(as.character(df$sample), facet$sample)]
  pal <- taxa_palette(rownames(mat))
  if (type == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(sample, pct, fill = taxon)) +
      ggplot2::geom_col(width = 0.9, colour = "white", linewidth = 0.1) +
      ggplot2::scale_fill_manual(values = pal, breaks = rownames(mat), name = legend_title) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
      ggplot2::labs(x = NULL, y = "Abondance relative (%)", title = title)
  } else {
    df <- df[!df$taxon %in% OTHERS, ]
    p <- ggplot2::ggplot(df, ggplot2::aes(sample, taxon, fill = pct)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(name = "%", trans = "sqrt", option = "magma", direction = -1) +
      ggplot2::labs(x = NULL, y = NULL, title = title)
  }
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   panel.grid = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 9))
  if (!is.null(facet)) p <- p + ggplot2::facet_grid(~facet, scales = "free_x", space = "free_x")
  p
}

# -----------------------------------------------------------------------------
# Krona
# -----------------------------------------------------------------------------
KRONA_URL <- "https://marbl.github.io/Krona/"

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

#' Page HTML Krona (format KronaTools 2.x ; moteur JavaScript chargé depuis le site de Krona)
#' @param counts matrice jeux de données x ASV : une ligne par échantillon ou par groupe
#'               (chaque ligne devient un jeu de données sélectionnable dans Krona)
#' @param add_asv ajoute les ASV comme dernier niveau de l'arbre
krona_html <- function(counts, tax = NULL, asv_ids = colnames(counts), add_asv = TRUE,
                       root = "Tous", magnitude_label = "Reads") {
  keep <- colSums(counts) > 0
  counts <- counts[, keep, drop = FALSE]; asv_ids <- asv_ids[keep]
  if (!is.null(tax)) tax <- tax[keep, , drop = FALSE]
  paths <- lapply(seq_len(ncol(counts)), function(j) {
    p <- character(0)
    if (!is.null(tax)) {
      v <- as.character(tax[j, ])
      k <- which(is.na(v) | !nzchar(v))[1]          # on s'arrête au premier rang non assigné
      p <- if (is.na(k)) v else v[seq_len(k - 1)]
    }
    if (add_asv || !length(p)) p <- c(p, asv_ids[j])
    p
  })
  fmt <- function(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 4))
  out <- character(0)
  walk <- function(name, idx, depth, indent) {
    mags <- rowSums(counts[, idx, drop = FALSE])
    out[[length(out) + 1]] <<- paste0(indent, '<node name="', xml_escape(name), '"><magnitude>',
                                      paste0("<val>", fmt(mags), "</val>", collapse = ""), "</magnitude>")
    nxt <- vapply(paths[idx], function(p) if (length(p) > depth) p[depth + 1] else NA_character_, "")
    has <- !is.na(nxt)
    if (any(has)) {
      kids <- unique(nxt[has])
      w <- vapply(kids, function(k) sum(counts[, idx[has][nxt[has] == k]]), 0)
      for (k in kids[order(-w)]) walk(k, idx[has][nxt[has] == k], depth + 1, paste0(indent, " "))
    }
    out[[length(out) + 1]] <<- paste0(indent, "</node>")
  }
  walk(root, seq_len(ncol(counts)), 0, "   ")

  c('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">',
    '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">',
    ' <head>',
    '  <meta charset="utf-8"/>',
    paste0('  <link rel="shortcut icon" href="', KRONA_URL, 'img/favicon.ico"/>'),
    paste0('  <script id="notfound">window.onload=function(){document.body.innerHTML=',
           '"Impossible de charger Krona depuis ', KRONA_URL,
           ' : une connexion internet est nécessaire pour afficher ce graphique."}</script>'),
    paste0('  <script src="', KRONA_URL, 'src/krona-2.0.js" type="text/javascript"></script>'),
    ' </head>',
    ' <body>',
    paste0('  <img id="hiddenImage" src="', KRONA_URL, 'img/hidden.png" style="display:none"/>'),
    '  <noscript>Javascript must be enabled to view this page.</noscript>',
    '  <div style="display:none">',
    '  <krona collapse="false" key="true">',
    '   <attributes magnitude="magnitude">',
    paste0('    <attribute display="', xml_escape(magnitude_label), '">magnitude</attribute>'),
    '   </attributes>',
    '   <datasets>',
    paste0('    <dataset>', xml_escape(rownames(counts)), '</dataset>'),
    '   </datasets>',
    unlist(out),
    '  </krona>',
    '</div></body></html>')
}

#' Prépare les jeux de données Krona
#' @param mode "sample" (un jeu par échantillon), "group" (moyenne des % par groupe), "all" (tout cumulé)
krona_datasets <- function(counts, mode = c("sample", "group", "all"), groups = NULL) {
  mode <- match.arg(mode)
  if (mode == "sample") return(counts)
  if (mode == "all") return(matrix(colSums(counts), nrow = 1, dimnames = list("Tous les échantillons", colnames(counts))))
  if (is.null(groups)) stop("Choisis une variable de regroupement.")
  rel <- sweep(counts, 1, pmax(rowSums(counts), 1), "/") * 100
  g <- as.character(groups); g[is.na(g) | !nzchar(g)] <- "(manquant)"
  s <- rowsum(rel, g, reorder = FALSE)
  s / as.vector(table(g)[rownames(s)])
}

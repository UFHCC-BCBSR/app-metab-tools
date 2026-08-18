# ===========================
# Shared report helpers
# ===========================
# Small HTML/plot helpers used by both the processing report and the
# differential analysis report sections.

# Serialise a plotly figure into a self-contained <div> + Plotly.newPlot call.
# The report <head> loads plotly from the CDN, so no library is inlined here.
embed_plotly <- function(fig, height = "400px") {
  if (is.null(fig)) return("")
  div_id <- paste0("plotly_", paste(sample(letters, 10, replace = TRUE), collapse = ""))
  fig_json <- plotly::plotly_json(fig, jsonedit = FALSE)
  paste0(
    '<div id="', div_id, '" style="width:100%;height:', height, ';"></div>',
    '<script>',
    'Plotly.newPlot("', div_id, '",',
    fig_json,
    ');</script>'
  )
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

# Numbers for display. Data files keep full precision; only reports round.
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(signif(as.numeric(x), digits), format = "g", digits = digits))
}

fmt_pval <- function(x) {
  ifelse(is.na(x), "",
         ifelse(x < 0.001,
                formatC(x, format = "e", digits = 2),
                formatC(x, format = "f", digits = 4)))
}

# Render a data.frame as an HTML table using the report's .file-table styling.
html_table <- function(df, css_class = "file-table") {
  header <- paste0("<tr>", paste0("<th>", html_escape(colnames(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(r) {
    paste0("<tr>", paste0("<td>", html_escape(r), "</td>", collapse = ""), "</tr>")
  })
  paste0('<table class="', css_class, '">', header, paste(rows, collapse = ""), "</table>")
}

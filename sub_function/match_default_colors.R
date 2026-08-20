#==================================================================================
#   Written by Noushin Farnoud, Aug 2026
#----------------------------------------------------------------------------------
#   Matches a default color palette (from heatmap_colors()) to a set of banner/
#   annotation names or data values in a FLEXIBLE way - i.e., ignoring case and
#   common punctuation differences (spaces, dots, dashes, underscores) so that
#   e.g. "Bcr-Abl1", "BCR_ABL1", "bcr.abl1" and "BCR-ABL1" all resolve to the
#   same default color instead of silently falling back to random colors.
#
#   See also heatmap_colors, add_ALL_banners, generate_complex_oncoprint.
#==================================================================================

# Normalizes a string for flexible/fuzzy comparison (case + punctuation-insensitive)
normalize_for_matching <- function(x) {
  toupper(trimws(gsub("[-_.[:space:]]+", "", x)))
}

# palette: a named color vector/list (e.g. list.ht.colors$ALL.SUBTYPE)
# values:  the actual data values (or banner names) you need colors for
# Returns a named color vector, named with the ORIGINAL `values` (not the
# palette names), so it can be plugged directly into HeatmapAnnotation's `col=`.
match_default_colors <- function(palette, values, quiet = FALSE) {

  values <- unique(as.character(values))

  pal.norm <- normalize_for_matching(names(palette))
  val.norm <- normalize_for_matching(values)

  idx <- match(val.norm, pal.norm)

  matched.colors <- setNames(unlist(palette)[idx], values)

  unmatched <- values[is.na(idx)]

  if (!quiet && length(unmatched) > 0) {
    cat(paste0("\n*** No default color match found for: '", paste(unmatched, collapse = "', '"),
               "'. These will receive random colors unless you assign them manually. ***\n"))
  }

  matched.colors[!is.na(matched.colors)]
}

# Merges custom.banner.colors into list.ht.colors, EMBEDDING new labels into an
# existing named-color-vector palette (keeping its other entries, overriding on
# name clashes) instead of wholesale replacing it like modifyList() would. New
# top-level banner names (not already in `base`) are just added as-is.
merge_banner_colors <- function(base, custom) {

  for (nm in names(custom)) {

    can.embed <- nm %in% names(base) &&
      is.atomic(base[[nm]]) && is.atomic(custom[[nm]]) &&
      !is.null(names(base[[nm]])) && !is.null(names(custom[[nm]]))

    if (can.embed) {
      base[[nm]] <- c(base[[nm]][!(names(base[[nm]]) %in% names(custom[[nm]]))], custom[[nm]])
    } else {
      base[[nm]] <- custom[[nm]]
    }
  }

  base
}

#==================================================================================
#   Alias registry: lets ONE palette in heatmap_colors() serve MANY banner names
#   that don't resemble its name (e.g. RNA_SUBTYPE/DNA_SUBTYPE/FINAL_SUBTYPE all
#   want to reuse the `ALL.SUBTYPE` palette). Add an entry here instead of
#   hand-wiring a new lookup for every banner that should share an existing palette.
#==================================================================================
banner_palette_aliases <- list(
  "ALL.SUBTYPE" = c("RNA_SUBTYPE", "DNA_SUBTYPE", "FINAL_SUBTYPE", "SOC_SUBTYPE",
                    "EXT.DX_SUBTYPE", "EXTENDED_DIAGNOSTIC_SUBTYPE", "HIGH_ANY2_TYPE")
)

# Resolves which palette in `list.ht.colors` a given banner name should use:
# 1) direct name match (case/punctuation/COLORS-suffix-insensitive)
# 2) alias lookup (`banner_palette_aliases`) for banners that share another palette's name
# Returns the whole palette (not filtered to data values) or NULL if nothing matches.
resolve_banner_palette <- function(banner_name, list.ht.colors, aliases = banner_palette_aliases) {

  ht.names.norm <- gsub("(COLORS?|COL|PALETTE|PAL)$", "", normalize_for_matching(names(list.ht.colors)))
  banner.norm <- normalize_for_matching(banner_name)

  idx <- which(ht.names.norm == banner.norm)
  if (length(idx) >= 1) return(list.ht.colors[[idx[1]]])

  for (palette_name in names(aliases)) {
    if (banner.norm %in% normalize_for_matching(aliases[[palette_name]])) {
      pal.idx <- which(ht.names.norm == normalize_for_matching(palette_name))
      if (length(pal.idx) >= 1) return(list.ht.colors[[pal.idx[1]]])
    }
  }

  NULL
}

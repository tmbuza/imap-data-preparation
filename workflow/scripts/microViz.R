check_is_phyloseq <- function(x, argName = NULL, allow_psExtra = TRUE) {
  stopif_ps_extra(x, argName = argName, Ncallers = 2)
  isPhyloseq <- is(x, "phyloseq") && (allow_psExtra || !is(x, "psExtra"))

  if (!isPhyloseq) {
    CLASSES <- if (allow_psExtra) '"phyloseq" or "psExtra"' else '"phyloseq"'

    rlang::abort(call = rlang::caller_env(), message = c(
      paste("argument", argName, "must be a", CLASSES, "object"),
      i = paste0("argument is class: ", paste(class(x), collapse = " "))
    ))
  }
}

check_is_psExtra <- function(x, argName = NULL) {
  stopif_ps_extra(x, argName = argName, Ncallers = 2)
  if (!is(x, "psExtra")) {
    rlang::abort(call = rlang::caller_env(), message = c(
      paste("argument", argName, 'must be a "psExtra" object'),
      i = paste0("argument is class: ", paste(class(x), collapse = " "))
    ))
  }
}

#' @export
#' @rdname psExtra-accessors
ps_get <- function(psExtra, ps_extra, counts = FALSE, warn = TRUE) {
  if (!missing(ps_extra)) psExtra <- ps_extra_arg_deprecation_warning(ps_extra)
  check_is_phyloseq(psExtra)
  if (isTRUE(counts)) {
    return(ps_counts(psExtra, warn = warn))
  }
  return(as(psExtra, "phyloseq"))
}
#' @rdname psExtra-accessors
#' @export
dist_get <- function(psExtra, ps_extra) {
  if (!missing(ps_extra)) psExtra <- ps_extra_arg_deprecation_warning(ps_extra)
  check_is_psExtra(psExtra)
  psExtra@dist
}
#' @rdname psExtra-accessors
#' @export
ord_get <- function(psExtra, ps_extra) {
  if (!missing(ps_extra)) psExtra <- ps_extra_arg_deprecation_warning(ps_extra)
  check_is_psExtra(psExtra)
  psExtra@ord
}
#' @rdname psExtra-accessors
#' @export
info_get <- function(psExtra, ps_extra) {
  if (!missing(ps_extra)) psExtra <- ps_extra_arg_deprecation_warning(ps_extra)
  check_is_phyloseq(psExtra)
  if (!methods::is(psExtra, "psExtra")) {
    return(new_psExtraInfo())
  }
  return(psExtra@info)
}
#' @rdname psExtra-accessors
#' @export
perm_get <- function(psExtra, ps_extra) {
  if (!missing(ps_extra)) psExtra <- ps_extra_arg_deprecation_warning(ps_extra)
  check_is_psExtra(psExtra)
  return(psExtra@permanova)
}
#' @rdname psExtra-accessors
#' @export
bdisp_get <- function(psExtra, ps_extra) {
  if (!missing(ps_extra)) psExtra <- ps_extra_arg_deprecation_warning(ps_extra)
  check_is_psExtra(psExtra)
  return(psExtra@bdisp)
}


#' @rdname psExtra-accessors
#' @export
tax_models_get <- function(psExtra) {
  check_is_psExtra(psExtra, argName = "psExtra")
  return(psExtra@tax_models)
}

#' @rdname psExtra-accessors
#' @export
tax_stats_get <- function(psExtra) {
  check_is_psExtra(psExtra, argName = "psExtra")
  return(psExtra@tax_stats)
}

#' @rdname psExtra-accessors
#' @export
taxatree_models_get <- function(psExtra) {
  check_is_psExtra(psExtra, argName = "psExtra")
  return(psExtra@taxatree_models)
}

#' @rdname psExtra-accessors
#' @export
taxatree_stats_get <- function(psExtra) {
  check_is_psExtra(psExtra, argName = "psExtra")
  return(psExtra@taxatree_stats)
}

#' @param data phyloseq or ps_extra
# @return phyloseq otu_table matrix with taxa as columns
#'
#' @param taxa subset of taxa to return, NA for all (default)
#' @param samples subset of samples to return, NA for all (default)
#' @param counts should ps_get or otu_get attempt to return counts? if present in object
#' @param warn
#' if counts = TRUE, should a warning be emitted if counts are not available?
#' set warn = "error" to stop if counts are not available
#'
#' @rdname psExtra-accessors
#' @export
otu_get <- function(data, taxa = NA, samples = NA, counts = FALSE, warn = TRUE) {
  # get otu_table from object
  if (methods::is(data, "otu_table")) {
    if (isTRUE(counts)) warning("data is otu_table: ignoring `counts = TRUE`")
    otu <- data
  } else {
    ps <- if (isTRUE(counts)) ps_counts(data, warn = warn) else ps_get(data)
    otu <- phyloseq::otu_table(ps)
  }
  if (phyloseq::taxa_are_rows(otu)) otu <- phyloseq::t(otu)

  # subset samples and/or taxa if requested, with slightly more helpful errors
  if (!identical(taxa, NA)) {
    stopifnot(is.character(taxa) || is.numeric(taxa) || is.logical(taxa))
    tmp <- try(expr = otu <- otu[, taxa, drop = FALSE], silent = TRUE)
    if (inherits(tmp, "try-error")) {
      if (is.character(taxa)) {
        wrong <- paste(setdiff(taxa, colnames(otu)), collapse = " / ")
        stop("The following taxa were not found in the otu table:\n", wrong)
      } else {
        stop("Invalid taxa selection")
      }
    }
  }
  if (!identical(samples, NA)) {
    stopifnot(is.character(samples) || is.numeric(samples) || is.logical(samples))
    tmp <- try(expr = otu <- otu[samples, , drop = FALSE], silent = TRUE)
    if (inherits(tmp, "try-error")) {
      if (is.character(samples)) {
        wrong <- paste(setdiff(samples, rownames(otu)), collapse = " / ")
        stop("The following samples were not found in the otu table:\n", wrong)
      } else {
        stop("Invalid sample selection")
      }
    }
  }
  return(otu)
}

#' @rdname psExtra-accessors
#' @export
tt_get <- function(data) {
  if (methods::is(data, "taxonomyTable")) {
    return(data)
  }
  tt <- phyloseq::tax_table(ps_get(data))
  return(tt)
}

#' @param data phyloseq or psExtra
# @return phyloseq sample_data as a tibble,
# with sample_names as new first column called .sample_name
#' @param sample_names_col
#' name of column where sample_names are put.
#' if NA, return data.frame with rownames (sample_names)
#' @rdname psExtra-accessors
#' @export
samdat_tbl <- function(data, sample_names_col = ".sample_name") {
  if (is(data, "sample_data") || is(data, "phyloseq")) {
    df <- samdatAsDataframe(data) # also works for psExtra
  } else {
    rlang::abort(message = c(
      "data must be of class 'phyloseq', 'psExtra', or 'sample_data'",
      i = paste("It is class:", paste(class(data), collapse = " "))
    ))
  }
  if (identical(sample_names_col, NA)) {
    return(df)
  } else {
    df <- tibble::rownames_to_column(df, var = sample_names_col)
    return(tibble::as_tibble(df, .name_repair = "check_unique"))
  }
}

# internal helper that get phyloseq sample_data as plain dataframe
# without changing invalid colnames (like microbiome::meta does)
# or losing rownames / sample_names (like data.frame() with defaults does)
samdatAsDataframe <- function(ps) {
  samdat <- phyloseq::sample_data(ps)
  df <- data.frame(samdat, check.names = FALSE, stringsAsFactors = FALSE)
  return(df)
}

# get phyloseq with counts if available
ps_counts <- function(data, warn = TRUE) {
  check_is_phyloseq(data)
  if (!rlang::is_bool(warn) && !rlang::is_string(warn, string = "error")) {
    stop("warn argument must be TRUE, FALSE, or 'error'")
  }
  counts <- NULL # check this later, warn if still NULL

  # always get ps, regardless of psExtra or phyloseq data or counts presence
  ps <- ps_get(data)

  # get counts and use them if they exist,
  # and check regardless if otutab returned will be counts
  if (is(data, "psExtra")) counts <- data@counts

  # maintain existing taxa_are_rows status for consistency
  if (phyloseq::taxa_are_rows(ps) && !is.null(counts)) counts <- phyloseq::t(counts)
  # put non-null counts table in otu table slot
  if (!is.null(counts)) phyloseq::otu_table(ps) <- counts

  if (isFALSE(warn)) {
    return(ps)
  }

  mess <- paste0(
    "otu_table of counts is NOT available!\n",
    "Available otu_table contains non-zero values that are less than 1"
  )

  # lastly check ps otu_table is counts
  test_matrix <- unclass(otu_get(ps))
  if (any(test_matrix < 1 & test_matrix != 0)) {
    if (identical(warn, "error")) stop(mess)
    if (isTRUE(warn)) warning(mess)
  }
  return(ps)
}

#-------------------------------------


tax_fix <- function(ps,
                    min_length = 4,
                    unknowns = NA,
                    suffix_rank = "classified", # or current
                    sep = " ",
                    anon_unique = TRUE,
                    verbose = TRUE) {
  if (is(ps, "psExtra")) stop("ps is a psExtra, run tax_fix BEFORE tax_agg etc")
  if (methods::is(ps, "phyloseq")) {
    tt <- unclass(phyloseq::tax_table(ps))
  } else if (inherits(ps, "taxonomyTable")) {
    tt <- unclass(ps)
  } else {
    stop("ps must be phyloseq or taxonomyTable class, it is class: ", paste(class(ps), collapse = " "))
  }

  if (identical(unknowns, NA)) {
    unknowns <- tax_common_unknowns(min_length = min_length)
  }

  # get rownames to ensure order doesn't change
  original_rownames <- rownames(tt)
  ranknames <- colnames(tt)
  levels <- ranknames # TODO either fix levels functionality or remove completely
  tt[is.na(tt) | nchar(tt) < min_length | tt %in% unknowns] <- ""
  rowLengthOut <- ncol(tt) # save number of cols before adding .rownames.
  tt_extra <- cbind(tt, .rownames. = original_rownames)

  # split the matrix tt (splits as a vector: see https://stackoverflow.com/questions/6819804/how-to-convert-a-matrix-to-a-list-of-column-vectors-in-r/6823557)
  tt_list <- split(tt_extra, row(tt_extra))
  # repair lost names
  names(tt_list) <- original_rownames

  # replace tax_table unknowns
  tt_out <- vapply(
    X = tt_list,
    FUN.VALUE = character(length = rowLengthOut),
    FUN = function(vec) {
      vec <- unlist(vec)
      is_unknown <- vec == "" # all unknowns now replaced first!
      # replace unknowns with nearest known if required and possible
      if (any(is_unknown[1:rowLengthOut])) {
        if (all(is_unknown[1:rowLengthOut])) {
          if (isTRUE(anon_unique)) {
            out <- vec[rowLengthOut + 1]
          } else {
            out <- paste("unclassified", ranknames[[1]])
          }
          if (isTRUE(verbose)) {
            message(
              "Row named: ", vec[[rowLengthOut + 1]],
              "\ncontains no non-unknown values, returning:\n'",
              out, "' for all replaced levels.\n",
              "Consider editing this tax_table entry manually."
            )
          }
          vec <- rep(out, times = rowLengthOut)
        } else {
          # edit each unknown value in this row
          vec <- tax_fix_value(
            vec = vec, is_unknown = is_unknown, ranknames = ranknames,
            levels = levels, sep = sep, suffix_rank = suffix_rank,
            rowLengthOut = rowLengthOut
          )
        }
      }
      return(vec[1:rowLengthOut])
    }
  )

  # make match original tt format (taxa as rows, ranks as cols)
  if (inherits(tt_out, "matrix")) {
    tt_out <- t(tt_out)
  } else {
    # vapply returns vector if tt had only 1 rank
    tt_out <- as.matrix(tt_out) # returns 1-column matrix
  }

  # ensure original row order
  tt_out <- tt_out[original_rownames, , drop = FALSE]
  # repair colnames
  colnames(tt_out) <- ranknames

  # preserve any columns not listed in levels
  preserved <- !ranknames %in% levels
  if (any(preserved)) tt_out[, preserved] <- tt[, preserved]

  # return phyloseq or tax table (same as input)
  if (inherits(ps, "phyloseq")) {
    phyloseq::tax_table(ps) <- tt_out
    return(ps)
  } else {
    return(phyloseq::tax_table(tt_out))
  }
}

#' internal helper for tax_fix
#'
#' @param vec row over which iteration is occuring
#' @param is_unknown logical vector of length length(row)
#' @param ranknames character vector of length length(row)
#' @param rowLengthOut number of ranks, originally
#'
#' @inheritParams tax_fix
#' @return fixed value
#' @noRd
tax_fix_value <- function(vec,
                          is_unknown,
                          ranknames,
                          levels,
                          rowLengthOut,
                          suffix_rank,
                          sep) {
  vec_out <- vapply(
    X = seq_along(vec),
    FUN.VALUE = vec[1],
    FUN = function(i) {
      if (is_unknown[i]) {
        known_above <- !is_unknown[1:(i - 1)]
        if (!any(known_above)) {
          if (ranknames[[i]] %in% levels) {
            stop(
              "Unknown values detected to the left of known values\n",
              "in row named: ", vec[[rowLengthOut + 1]],
              "\nThis should not happen. Check/fix this row:\n",
              paste(vec[1:rowLengthOut], collapse = "; ")
            )
          }
          return(vec[i])
        }
        nearest_known <- max(which(known_above))
        tax <- vec[nearest_known]
        if (identical(suffix_rank, "classified")) {
          level <- ranknames[[nearest_known]]
        } else {
          level <- ranknames[[i]]
        }
        tax <- paste(tax, level, sep = sep)
        return(tax)
      } else {
        return(vec[i])
      }
    }
  )
  return(vec_out)
}

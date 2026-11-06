# Tree handling: per-tree inverseA caching, multiPhylo validation, and the
# deterministic tree-per-imputation selection rule.
#
# Design notes (see .agents/plans/next-features-2026-11.md, feature 1 + 4a):
# - .get_ainv() caches MCMCglmm::inverseA() per tree. The cache is cleared at
#   the entry of each phase (bace_imp, bace_final_imp, phylo_signal_summary) so
#   memory is bounded; mclapply workers fork AFTER entry and inherit the warm
#   cache copy-on-write (under parallel execution the true behaviour is one
#   compute per tree per worker, which is cheap: ~0.05 s at 10,000 tips).
# - Tree selection for the final runs is DETERMINISTIC in the run index and
#   consumes no RNG, so the per-run L'Ecuyer stream contract (seed, run) and
#   published seeds survive unchanged.
# - .validate_trees() runs .get_ainv() on every tree up front: this surfaces
#   every condition inverseA enforces (most importantly ultrametricity under
#   its default scale = TRUE) with the offending tree named, BEFORE hours of
#   fitting, and pre-warms the cache in the parent process.

# Package-local Ainv cache --------------------------------------------------

.ainv_cache <- new.env(parent = emptyenv())

#' Hash an R object, stably across R versions
#'
#' serialize() output embeds the writer's R version (and, for version 3, the
#' native encoding) in its header; hashing raw serialize() output therefore
#' changes across R upgrades for byte-identical objects. Version-2
#' serialization has a fixed 14-byte header ("X\n" + three XDR ints), which is
#' stripped before hashing so the hash depends only on content.
#' @keywords internal
#' @noRd
.hash_obj <- function(x) {
  raw <- serialize(x, connection = NULL, version = 2)
  raw <- raw[-seq_len(14L)]
  tf <- tempfile(fileext = ".bin")
  on.exit(unlink(tf), add = TRUE)
  writeBin(raw, tf)
  unname(as.character(tools::md5sum(tf)))
}

#' Cache key for a phylogeny
#'
#' Covers everything inverseA's result can depend on, including node labels
#' (inverseA names node rows after node labels when present, "Node<k>"
#' otherwise, so two trees differing only in node labels must not share a
#' cache entry) and Nnode.
#' @keywords internal
#' @noRd
.tree_cache_key <- function(tree) {
  .hash_obj(list(tree$edge, tree$edge.length, tree$tip.label,
                 tree$node.label, tree$Nnode))
}

#' @keywords internal
#' @noRd
.compute_ainv <- function(tree) {
  # nodes = "ALL" is more stable and orders of magnitude faster on large trees
  MCMCglmm::inverseA(tree, nodes = "ALL")$Ainv
}

#' Per-tree cached inverseA
#' @keywords internal
#' @noRd
.get_ainv <- function(tree) {
  key <- .tree_cache_key(tree)
  hit <- get0(key, envir = .ainv_cache, inherits = FALSE)
  if (!is.null(hit)) return(hit)
  A <- .compute_ainv(tree)
  assign(key, A, envir = .ainv_cache)
  A
}

#' @keywords internal
#' @noRd
.clear_ainv_cache <- function() {
  rm(list = ls(envir = .ainv_cache, all.names = TRUE), envir = .ainv_cache)
  invisible(NULL)
}

# multiPhylo handling --------------------------------------------------------

#' First tree of a phylo | multiPhylo | list-of-phylo input
#' @keywords internal
#' @noRd
.first_tree <- function(phylo) {
  if (inherits(phylo, "phylo")) phylo else phylo[[1L]]
}

#' Coerce the phylo argument to a plain list of phylo objects
#'
#' multiPhylo is unpacked with `[[` (not unclass) so ape's "compressed"
#' multiPhylo storage (shared tip labels in attr(x, "TipLabel")) round-trips
#' correctly.
#' @keywords internal
#' @noRd
.as_tree_list <- function(phylo, context = "bace") {
  if (inherits(phylo, "phylo")) return(list(phylo))
  if (inherits(phylo, "multiPhylo")) {
    return(lapply(seq_along(phylo), function(i) phylo[[i]]))
  }
  if (is.list(phylo) && length(phylo) > 0L &&
      all(vapply(phylo, inherits, logical(1L), "phylo"))) {
    return(phylo)
  }
  stop(context, ": 'phylo' must be a phylo object, a multiPhylo, or a list ",
       "of phylo objects (ape).", call. = FALSE)
}

#' Validate a tree set and pre-warm the Ainv cache
#'
#' Per tree: class; data species covered by the tip labels; node labels
#' stripped (inverseA hard-stops on duplicated node labels, e.g. support
#' values, and node labels fragment the Ainv cache key and the pooled Sol
#' column names); negative and zero-length branches error (inverseA rejects
#' zero-length edges unconditionally -- the standard fix is replacing 0 with
#' a tiny length such as 1e-8, and the message says so); and
#' .get_ainv() is run so every inverseA rejection (non-ultrametric trees
#' under its default scale = TRUE, etc.) surfaces NOW with the tree index
#' named rather than hours into a fit. Across trees: a root-depth
#' consistency warning (trees must come from one analysis in common units;
#' CV of root depths > 15% suggests mixed scalings).
#'
#' @param phylo phylo | multiPhylo | list of phylo
#' @param species_values the data's cluster column values (species names)
#' @param context character, used in messages ("bace_imp", "bace_final_imp")
#' @param verbose logical
#' @return plain list of validated phylo objects (node labels stripped)
#' @keywords internal
#' @noRd
.validate_trees <- function(phylo, species_values, context = "bace",
                            verbose = TRUE) {
  trees <- .as_tree_list(phylo, context = context)
  n_trees <- length(trees)
  sp <- unique(as.character(species_values))
  n_stripped <- 0L
  depths <- numeric(n_trees)

  for (i in seq_len(n_trees)) {
    tr <- trees[[i]]
    tag <- if (n_trees > 1L) paste0("tree ", i, " of ", n_trees) else "the tree"

    missing_sp <- setdiff(sp, tr$tip.label)
    if (length(missing_sp) > 0L) {
      stop(context, ": ", length(missing_sp), " species in the data are not ",
           "tip labels of ", tag, " (first few: ",
           paste(utils::head(missing_sp, 5L), collapse = ", "), ").",
           call. = FALSE)
    }

    if (!is.null(tr$node.label)) {
      tr$node.label <- NULL
      n_stripped <- n_stripped + 1L
    }

    if (is.null(tr$edge.length)) {
      stop(context, ": ", tag, " has no branch lengths (edge.length is NULL).",
           call. = FALSE)
    }
    if (any(tr$edge.length < 0)) {
      stop(context, ": ", tag, " has negative branch lengths.", call. = FALSE)
    }
    if (any(tr$edge.length == 0)) {
      # MCMCglmm::inverseA() rejects zero-length edges unconditionally
      # ("some phylogeny edge.lengths are zero"); say so here with the tree
      # named and the standard remedy.
      stop(context, ": ", tag, " has ", sum(tr$edge.length == 0),
           " zero-length branch(es), which MCMCglmm::inverseA() rejects. ",
           "Replace them with a tiny length first, e.g. ",
           "tree$edge.length[tree$edge.length == 0] <- 1e-8.", call. = FALSE)
    }

    # Surface every inverseA rejection now, with the tree named; also
    # pre-warms the per-tree Ainv cache before any worker forks.
    tryCatch(.get_ainv(tr), error = function(e) {
      stop(context, ": ", tag, " was rejected by MCMCglmm::inverseA(): ",
           conditionMessage(e),
           " (a common cause is a non-ultrametric tree, including ",
           "ultrametric trees written to file with limited precision -- ",
           "see ape::is.ultrametric and phytools::force.ultrametric).",
           call. = FALSE)
    })

    depths[i] <- max(ape::node.depth.edgelength(tr))
    trees[[i]] <- tr
  }

  if (n_stripped > 0L && verbose) {
    message(context, ": node labels stripped from ", n_stripped, " tree",
            if (n_stripped > 1L) "s" else "",
            " (they fragment internal bookkeeping and MCMCglmm::inverseA ",
            "rejects duplicated ones; tip labels are untouched).")
  }
  if (n_trees > 1L) {
    cv <- stats::sd(depths) / mean(depths)
    if (is.finite(cv) && cv > 0.15) {
      warning(context, ": root depths vary substantially across the supplied ",
              "trees (CV = ", round(cv, 2), "; range ",
              signif(min(depths), 3), "-", signif(max(depths), 3), "). ",
              "Genuine posterior root-height variation is fine, but mixed ",
              "branch-length scalings (e.g. one rescaled tree among raw ",
              "ones) corrupt pooled variance components. Trees must come ",
              "from one analysis in common units.", call. = FALSE)
    }

    tip_sets_differ <- !all(vapply(
      trees[-1L], function(t2) setequal(t2$tip.label, trees[[1L]]$tip.label),
      logical(1L)))
    if (tip_sets_differ) {
      warning(context, ": the supplied trees do not all share the same tip ",
              "set. Every tree covers the data's species (checked), but ",
              "species-level effects for tips absent from some trees cannot ",
              "be pooled across runs.", call. = FALSE)
    }
  }

  trees
}

#' Deterministic tree-per-run selection
#'
#' n_trees >= n_final: an evenly spaced subsequence of the (chain-ordered)
#' posterior sample -- thinning, in effect. n_trees < n_final: balanced
#' recycling. Both consume no RNG: selection depends only on (n_trees,
#' n_final), so the per-run RNG stream contract and published seeds are
#' unaffected.
#'
#' @return integer vector of length n_final: index into the user's tree set
#'   used by each final run.
#' @keywords internal
#' @noRd
.select_tree_indices <- function(n_trees, n_final) {
  n_trees <- as.integer(n_trees)
  n_final <- as.integer(n_final)
  if (n_trees < 1L || n_final < 1L) {
    stop("n_trees and n_final must both be >= 1")
  }
  if (n_trees == 1L) return(rep(1L, n_final))
  if (n_trees >= n_final) {
    as.integer(round(seq(1L, n_trees, length.out = n_final)))
  } else {
    ((seq_len(n_final) - 1L) %% n_trees) + 1L
  }
}

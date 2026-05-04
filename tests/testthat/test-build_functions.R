# Tests for R/build_functions.R: formula construction + symbol-swap
# helpers used by bace_imp's chained-equation generator.

# ---- .build_formula -------------------------------------------------------

test_that(".build_formula converts a string to a formula object", {
  f <- BACE:::.build_formula("y ~ x1 + x2")
  expect_s3_class(f, "formula")
  expect_equal(length(f), 3L)
  expect_equal(as.character(f[[2]]), "y")
})

test_that(".build_formula accepts a formula object passthrough", {
  f0 <- y ~ x
  expect_error(BACE:::.build_formula(f0), NA)
})

test_that(".build_formula rejects non-formula strings", {
  expect_error(BACE:::.build_formula("not a formula"),
               regexp = "specify a formula string")
  expect_error(BACE:::.build_formula("y x"),
               regexp = "specify a formula string")
})

# ---- .build_formula_string -----------------------------------------------

test_that(".build_formula_string returns one formula per variable", {
  out <- BACE:::.build_formula_string("y ~ x1 + x2")
  expect_type(out, "list")
  expect_named(out, c("y", "x1", "x2"))
  for (f in out) expect_s3_class(f, "formula")
})

test_that(".build_formula_string keeps original formula at LHS=y", {
  out <- BACE:::.build_formula_string("y ~ x1 + x2")
  expect_equal(deparse(out$y), deparse(stats::as.formula("y ~ x1 + x2")))
})

test_that(".build_formula_string swaps response with each predictor", {
  out <- BACE:::.build_formula_string("y ~ x1 + x2")
  # When iterating over x1 it becomes the LHS, with y appearing on
  # the RHS in place of x1's original slot.
  expect_equal(as.character(out$x1[[2]]), "x1")
  rhs_vars <- all.vars(out$x1[[3]])
  expect_true("y" %in% rhs_vars)
  expect_true("x2" %in% rhs_vars)
  expect_false("x1" %in% rhs_vars)  # x1 is on LHS now, not RHS
})

test_that(".build_formula_string handles interactions", {
  out <- BACE:::.build_formula_string("y ~ x1 + x2*x3")
  expect_named(out, c("y", "x1", "x2", "x3"))
  for (f in out) expect_s3_class(f, "formula")
})

test_that(".build_formula_string handles colon (:) interactions", {
  out <- BACE:::.build_formula_string("y ~ x1 + x2:x3")
  expect_named(out, c("y", "x1", "x2", "x3"))
})

test_that(".build_formula_string accepts a formula directly", {
  out <- BACE:::.build_formula_string(y ~ x1 + x2)
  expect_named(out, c("y", "x1", "x2"))
})

test_that(".build_formula_string rejects formulas without LHS", {
  expect_error(BACE:::.build_formula_string("~ x1 + x2"),
               regexp = "of the form: y ~")
})

# ---- .swap_symbols / .replace_symbols ------------------------------------

test_that(".replace_symbols replaces variable names in expressions", {
  e <- quote(a + b * c)
  out <- BACE:::.replace_symbols(e, list(a = "z"))
  expect_equal(deparse(out), "z + b * c")
})

test_that(".replace_symbols leaves unmatched names alone", {
  e <- quote(x + y)
  out <- BACE:::.replace_symbols(e, list(z = "w"))
  expect_equal(deparse(out), deparse(e))
})

test_that(".replace_symbols handles numeric atoms (no replacement needed)", {
  out <- BACE:::.replace_symbols(quote(2 + 3), list(a = "b"))
  expect_equal(deparse(out), "2 + 3")
})

test_that(".replace_symbols handles NULL gracefully", {
  expect_null(BACE:::.replace_symbols(NULL, list(a = "b")))
})

test_that(".swap_symbols swaps two names in an expression", {
  out <- BACE:::.swap_symbols(quote(a + b), a = "a", b = "b")
  # Swapping a<->b means a + b becomes b + a (semantically same; deparse matters)
  expect_equal(deparse(out), "b + a")
})

test_that(".swap_symbols is a no-op when a == b", {
  out <- BACE:::.swap_symbols(quote(a + b), a = "a", b = "a")
  expect_equal(deparse(out), "a + b")
})

# ---- .build_formula_string_random ----------------------------------------

test_that(".build_formula_string_random parses '~ 1 | Species' (single formula)", {
  # species = FALSE returns a single formula, not a list.
  out <- BACE:::.build_formula_string_random("~ 1 | Species", species = FALSE)
  expect_s3_class(out, "formula")
  expect_true("Species" %in% all.vars(out))
})

test_that(".build_formula_string_random with species=TRUE returns dual formulas", {
  out <- BACE:::.build_formula_string_random("~ 1 | Species", species = TRUE)
  expect_type(out, "list")
  expect_named(out, c("phylo", "species"))
  expect_s3_class(out$phylo, "formula")
  expect_s3_class(out$species, "formula")
})

test_that(".build_formula_string_random handles us() random-slope structure", {
  # us(1 + temperature):Species should produce a random-slope phylo formula
  out <- BACE:::.build_formula_string_random(
    "~ us(1 + temperature) | Species", species = TRUE)
  expect_named(out, c("phylo", "species"))
  expect_s3_class(out$phylo, "formula")
})

test_that(".build_formula_string_random rejects formulas without |", {
  expect_error(BACE:::.build_formula_string_random("~ 1", species = FALSE),
               regexp = "No random-effects structure")
})

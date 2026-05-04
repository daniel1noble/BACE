# Tests for R/bace_option_setter_getter.r
# Covers bace_options(), bace_option_defaults(), and the internal
# get/set/validate path through every documented option.

# Always reset BACE.* options after each test so we don't leak state
# into other test files. Each test that mutates options uses
# on.exit(reset_bace_options(), add = TRUE) -- avoids the withr
# dependency.
reset_bace_options <- function() {
  bace_options(.reset = TRUE)
}

test_that("bace_option_defaults returns expected keys", {
  d <- bace_option_defaults()
  expect_named(d, c("verbose", "digits", "gelman"))
  expect_type(d$verbose, "logical")
  expect_type(d$digits, "integer")
  expect_type(d$gelman, "double")
})

test_that("bace_options() with no args returns merged defaults", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  eff <- bace_options()
  expect_named(eff, c("verbose", "digits", "gelman"))
})

test_that("bace_options sets and gets verbose", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  bace_options(verbose = FALSE)
  expect_false(bace_options()$verbose)
  bace_options(verbose = TRUE)
  expect_true(bace_options()$verbose)
})

test_that("bace_options sets digits and coerces to integer", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  bace_options(digits = 5)
  expect_identical(bace_options()$digits, 5L)
})

test_that("bace_options sets gelman to allowed values", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  for (g in c(0, 1, 2)) {
    bace_options(gelman = g)
    expect_equal(bace_options()$gelman, g)
  }
})

test_that("bace_options rejects unknown option names", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  expect_error(bace_options(nope = 1L),
               regexp = "Unknown BACE option")
})

test_that("bace_options requires named args", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  expect_error(bace_options(TRUE),
               regexp = "must be named")
})

test_that("bace_options validates verbose as length-1 logical", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  expect_error(bace_options(verbose = c(TRUE, FALSE)),
               regexp = "verbose must be a single")
  expect_error(bace_options(verbose = "yes"),
               regexp = "verbose must be a single")
})

test_that("bace_options validates digits as length-1 numeric", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  expect_error(bace_options(digits = c(1, 2)),
               regexp = "digits must be length-1")
  expect_error(bace_options(digits = "five"),
               regexp = "digits must be length-1")
})

test_that("bace_options validates gelman against allowed set", {
  reset_bace_options()
  on.exit(reset_bace_options(), add = TRUE)
  expect_error(bace_options(gelman = 3),
               regexp = "gelman must be one of")
  expect_error(bace_options(gelman = "a"),
               regexp = "gelman must be a length-1 number")
})

test_that("bace_options(.reset = TRUE) clears user overrides", {
  bace_options(verbose = FALSE, digits = 7)
  bace_options(.reset = TRUE)
  expect_null(getOption("BACE.verbose"))
  expect_null(getOption("BACE.digits"))
  # But effective values revert to defaults via merge
  eff <- bace_options()
  expect_identical(eff$verbose, bace_option_defaults()$verbose)
  expect_identical(eff$digits, bace_option_defaults()$digits)
})

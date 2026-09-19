# Aliases for internal (unexported) helpers used directly by the test suite.
#
# These functions were exported during development but are internal as of the
# CRAN release (v0.1.0): they are implementation details, not API. The tests
# exercise them deliberately (unit tests of the machinery), so we alias them
# here once rather than sprinkling `BACE:::` through every call site.
# testthat sources helper-*.R files before running tests, and `:::` access to
# our own package's internals from its own test suite is standard practice.

.build_formula               <- BACE:::.build_formula
.build_formula_string        <- BACE:::.build_formula_string
.build_formula_string_random <- BACE:::.build_formula_string_random
.check_mcmc_diagnostics      <- BACE:::.check_mcmc_diagnostics
.count_categorical_fixef     <- BACE:::.count_categorical_fixef
.data_prep                   <- BACE:::.data_prep
.extract_gaussian_attrs      <- BACE:::.extract_gaussian_attrs
.get_type                    <- BACE:::.get_type
.get_variables               <- BACE:::.get_variables
.make_prior                  <- BACE:::.make_prior
.model_fit                   <- BACE:::.model_fit
.pred_cat                    <- BACE:::.pred_cat
.pred_cat_forward            <- BACE:::.pred_cat_forward
.pred_cont                   <- BACE:::.pred_cont
.pred_count                  <- BACE:::.pred_count
.pred_threshold              <- BACE:::.pred_threshold
.pred_threshold_forward      <- BACE:::.pred_threshold_forward
.predict_bace                <- BACE:::.predict_bace
.summarise_var_types         <- BACE:::.summarise_var_types

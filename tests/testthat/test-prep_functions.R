# Tests for prep_functions.R

# Test .get_variables() function ----
test_that(".get_variables() extracts fixed effects correctly", {
  # Simple formula
  expect_equal(.get_variables("y ~ x1 + x2", fix = TRUE)$fix, c("y", "x1", "x2"))

  expect_equal(.get_variables("y ~ x1.1 + x2.2", fix = TRUE)$fix, c("y", "x1.1", "x2.2"))

  expect_equal(.get_variables("y ~ x1_1 + x2_2", fix = TRUE)$fix, c("y", "x1_1", "x2_2"))
  
  # Formula with interaction
  expect_equal(.get_variables("y ~ x1.1 * x2.1", fix = TRUE)$fix, c("y", "x1.1", "x2.1"))
  expect_equal(.get_variables("y ~ x1_1 * x2_2", fix = TRUE)$fix, c("y", "x1_1", "x2_2"))
  
  # Formula with explicit interaction
  expect_equal(.get_variables("y ~ x1 + x2 + x1:x2", fix = TRUE)$fix, c("y", "x1", "x2"))
  expect_equal(.get_variables("y ~ x1.1 + x2.2 + x1.1:x2.2", fix = TRUE)$fix, c("y", "x1.1", "x2.2"))
  expect_equal(.get_variables("y ~ x1_1 + x2_2 + x1_1:x2_2", fix = TRUE)$fix, c("y", "x1_1", "x2_2"))
})

test_that(".get_variables() extracts random effects correctly", {
  # Simple random effect
  result <- .get_variables("~ 1 | Species", fix = FALSE)
  expect_equal(result$ran, c("1"))
  expect_equal(result$cluster, "Species")
  
  # Random effect with slope
  result <- .get_variables("~ 1 + x1 | Species", fix = FALSE)
  expect_equal(result$ran, c("1", "x1"))
  expect_equal(result$cluster, "Species")
  
})

test_that(".get_variables() handles edge cases", {
  # Formula with spaces
  result <- .get_variables("~ 1 | Species", fix = FALSE)
  expect_equal(result$cluster, "Species")
  
  # Formula without spaces
  result <- .get_variables("~1|Species", fix = FALSE)
  expect_equal(result$cluster, "Species")
})

test_that(".get_variables() handles variable names with dots and underscores", {
  # Fixed effects with dots and underscores
  result <- .get_variables("response.var ~ predictor_1 + predictor.2", fix = TRUE)
  expect_equal(result$fix, c("response.var", "predictor_1", "predictor.2"))
  
  # More complex variable names
  result <- .get_variables("my.response_var ~ pred.one_two + pred_three.four", fix = TRUE)
  expect_equal(result$fix, c("my.response_var", "pred.one_two", "pred_three.four"))
  
  # Random effects with dots and underscores
  result <- .get_variables("~ 1 + var.one | cluster_name.group", fix = FALSE)
  expect_equal(result$ran, c("1", "var.one"))
  expect_equal(result$cluster, "cluster_name.group")
  
  # Formula object with special characters
  formula_obj <- as.formula("resp_var.value ~ pred.one + pred_two")
  result <- .get_variables(formula_obj, fix = TRUE)
  expect_equal(result$fix, c("resp_var.value", "pred.one", "pred_two"))
})



# Test .get_type() function ----
test_that(".get_type() correctly identifies binary variables", {
  data <- data.frame(
    bin_unordered = factor(c("A", "B", "A", "B", "A")),
    bin_ordered = ordered(c("Low", "High", "Low", "High", "Low"))
  )
  
  summary_df <- .summarise_var_types(data)
  
  type_unordered <- .get_type(summary_df, "bin_unordered", data)
  type_ordered <- .get_type(summary_df, "bin_ordered", data)
  
  expect_equal(type_unordered, "threshold")
  expect_equal(type_ordered, "threshold")
})

test_that(".get_type() correctly identifies continuous variables", {
  data <- data.frame(
    continuous = c(1.5, 2.7, 3.2, 4.8, 5.1)
  )
  
  summary_df <- .summarise_var_types(data)
  type <- .get_type(summary_df, "continuous", data)
  
  expect_equal(type, "gaussian")
})

test_that(".get_type() correctly identifies count variables", {
  # Create a clearly Poisson-distributed count variable
  set.seed(123)
  data <- data.frame(
    count = rpois(100, lambda = 5)
  )
  
  summary_df <- .summarise_var_types(data)
  type <- .get_type(summary_df, "count", data)
  
  # Should be either poisson or gaussian depending on AIC comparison
  expect_true(type %in% c("poisson", "gaussian"))
})

test_that(".get_type() correctly identifies categorical variables", {
  data <- data.frame(
    cat_unordered = factor(c("A", "B", "C", "A", "B")),
    cat_ordered = ordered(c("Low", "Med", "High", "Low", "Med"))
  )
  
  summary_df <- .summarise_var_types(data)
  
  type_unordered <- .get_type(summary_df, "cat_unordered", data)
  type_ordered <- .get_type(summary_df, "cat_ordered", data)
  
  expect_equal(type_unordered, "categorical")
  expect_equal(type_ordered, "threshold")
})


# Test .data_prep() function ----
test_that(".data_prep() handles missing data in predictors", {
  data <- data.frame(
    y = c(1, 2, 3, 4, 5),
    x1 = c(10, 20, NA, 40, 50),
    x2 = factor(c("A", "B", NA, "A", "B")),
    cluster = c(1, 1, 2, 2, 2)
  )
  
  formula <- as.formula("y ~ x1 + x2")
  types <- list(y = "gaussian", x1 = "gaussian", x2 = "categorical", cluster = "categorical")
  
  result <- .data_prep(formula, data, types, "cluster")
  
  # Check that missing values are imputed
  expect_false(any(is.na(result[[1]]$x1)))
  expect_false(any(is.na(result[[1]]$x2)))
})

test_that(".data_prep() standardizes gaussian variables", {
  data <- data.frame(
    y = c(100, 200, 300, 400, 500),
    x1 = c(10, 20, 30, 40, 50),
    cluster = c(1, 1, 2, 2, 2)
  )
  
  formula <- as.formula("y ~ x1")
  types <- list(y = "gaussian", x1 = "gaussian", cluster = "categorical")
  
  result <- .data_prep(formula, data, types, "cluster")
  
  # Check that gaussian variables are standardized
  expect_equal(mean(result[[1]]$y, na.rm = TRUE), 0, tolerance = 1e-10)
  expect_equal(sd(result[[1]]$y, na.rm = TRUE), 1, tolerance = 1e-10)
  expect_equal(mean(result[[1]]$x1, na.rm = TRUE), 0, tolerance = 1e-10)
  expect_equal(sd(result[[1]]$x1, na.rm = TRUE), 1, tolerance = 1e-10)
})

test_that(".data_prep() returns attributes for gaussian variables", {
  data <- data.frame(
    y = c(100, 200, 300, 400, 500),
    x1 = c(10, 20, 30, 40, 50),
    cluster = c(1, 1, 2, 2, 2)
  )
  
  formula <- as.formula("y ~ x1")
  types <- list(y = "gaussian", x1 = "gaussian", cluster = "categorical")
  
  result <- .data_prep(formula, data, types, "cluster")
  
  # Check that attributes are stored
  expect_true(!is.null(result[[2]]))
  expect_true(!is.null(result[[2]]$y))
  expect_true(!is.null(result[[2]]$x1))
})

test_that(".data_prep() handles count variables with missing data", {
  data <- data.frame(
    y = c(1, 2, 3, 4, 5),
    x1 = c(10, 20, NA, 40, 50),
    cluster = c(1, 1, 2, 2, 2)
  )
  
  formula <- as.formula("y ~ x1")
  types <- list(y = "poisson", x1 = "poisson", cluster = "categorical")
  
  result <- .data_prep(formula, data, types, "cluster")
  
  # Check that missing values are imputed with rounded mean
  expect_false(any(is.na(result[[1]]$x1)))
  expect_true(all(result[[1]]$x1 == round(result[[1]]$x1)))
})


# Test .extract_gaussian_attrs() function ----
test_that(".extract_gaussian_attrs() extracts mean and sd correctly", {
  data <- data.frame(
    y = c(100, 200, 300, 400, 500),
    x1 = c(10, 20, 30, 40, 50),
    x2 = factor(c("A", "B", "C", "A", "B"))
  )
  
  types <- list(y = "gaussian", x1 = "gaussian", x2 = "categorical")
  
  attrs <- .extract_gaussian_attrs(data, types)
  
  expect_equal(attrs$y$mean, 300)
  expect_equal(attrs$y$sd, sd(c(100, 200, 300, 400, 500)))
  expect_equal(attrs$x1$mean, 30)
  expect_equal(attrs$x1$sd, sd(c(10, 20, 30, 40, 50)))
  expect_null(attrs$x2)
})

test_that(".extract_gaussian_attrs() handles missing values", {
  data <- data.frame(
    y = c(100, 200, NA, 400, 500),
    x1 = c(10, 20, 30, 40, 50)
  )
  
  types <- list(y = "gaussian", x1 = "gaussian")
  
  attrs <- .extract_gaussian_attrs(data, types)
  
  # Should compute mean and sd excluding NAs
  expect_equal(attrs$y$mean, mean(c(100, 200, 400, 500)))
  expect_equal(attrs$y$sd, sd(c(100, 200, 400, 500)))
})

test_that(".extract_gaussian_attrs() returns NULL for non-gaussian variables", {
  data <- data.frame(
    y = c(1, 2, 3, 4, 5),
    x1 = factor(c("A", "B", "C", "A", "B"))
  )
  
  types <- list(y = "poisson", x1 = "categorical")
  
  attrs <- .extract_gaussian_attrs(data, types)
  
  expect_null(attrs$y)
  expect_null(attrs$x1)
})


# Test .summarise_var_types() function ----
test_that(".summarise_var_types() correctly identifies variable characteristics", {
  data <- data.frame(
    numeric_var = c(1.5, 2.7, 3.2),
    integer_var = c(1L, 2L, 3L),
    factor_var = factor(c("A", "B", "C")),
    ordered_var = ordered(c("Low", "Med", "High")),
    character_var = c("a", "b", "c"),
    logical_var = c(TRUE, FALSE, TRUE)
  )
  
  summary <- .summarise_var_types(data)
  
  # Check that all variables are present
  expect_equal(nrow(summary), 6)
  expect_equal(summary$variable, names(data))
  
  # Check numeric variable
  expect_true(summary[summary$variable == "numeric_var", "is_numeric"])
  expect_false(summary[summary$variable == "numeric_var", "is_integer"])
  
  # Check integer variable
  expect_true(summary[summary$variable == "integer_var", "is_numeric"])
  expect_true(summary[summary$variable == "integer_var", "is_integer"])
  
  # Check factor variable
  expect_true(summary[summary$variable == "factor_var", "is_factor"])
  expect_false(summary[summary$variable == "factor_var", "is_ordered"])
  expect_equal(summary[summary$variable == "factor_var", "n_levels"], 3)
  
  # Check ordered variable
  expect_true(summary[summary$variable == "ordered_var", "is_factor"])
  expect_true(summary[summary$variable == "ordered_var", "is_ordered"])
  
  # Check character variable
  expect_true(summary[summary$variable == "character_var", "is_character"])
  
  # Check logical variable
  expect_true(summary[summary$variable == "logical_var", "is_logical"])
})

test_that(".summarise_var_types() detects missing values", {
  data <- data.frame(
    with_na = c(1, 2, NA, 4),
    without_na = c(1, 2, 3, 4)
  )
  
  summary <- .summarise_var_types(data)
  
  expect_true(summary[summary$variable == "with_na", "has_na"])
  expect_false(summary[summary$variable == "without_na", "has_na"])
})

test_that(".summarise_var_types() counts unique values correctly", {
  data <- data.frame(
    var1 = c(1, 2, 2, 3, 3, 3),
    var2 = c("A", "A", "A", "A", "A", "A")
  )
  
  summary <- .summarise_var_types(data)
  
  expect_equal(summary[summary$variable == "var1", "n_unique"], 3)
  expect_equal(summary[summary$variable == "var2", "n_unique"], 1)
})

test_that(".summarise_var_types() stores levels for factors", {
  data <- data.frame(
    factor_var = factor(c("A", "B", "C", "A"))
  )
  
  summary <- .summarise_var_types(data, store_levels = TRUE)
  
  expect_equal(summary$levels[[1]], c("A", "B", "C"))
})

test_that(".summarise_var_types() respects max_levels_store", {
  # Create a factor with many levels
  many_levels <- factor(paste0("Level", 1:300))
  data <- data.frame(many_var = many_levels)
  
  summary <- .summarise_var_types(data, store_levels = TRUE, max_levels_store = 200)
  
  # Levels should not be stored because it exceeds max_levels_store
  expect_null(summary$levels[[1]])
})


# Test .get_levels() function ----
test_that(".get_levels() extracts factor levels correctly", {
  x <- factor(c("A", "B", "C", "A"))
  levels_result <- .get_levels(x, store_levels = TRUE, max_levels_store = 200)
  
  expect_equal(levels_result, c("A", "B", "C"))
})

test_that(".get_levels() extracts character unique values correctly", {
  x <- c("apple", "banana", "cherry", "apple")
  levels_result <- .get_levels(x, store_levels = TRUE, max_levels_store = 200)
  
  expect_equal(levels_result, c("apple", "banana", "cherry"))
})

test_that(".get_levels() returns NULL when store_levels is FALSE", {
  x <- factor(c("A", "B", "C"))
  levels_result <- .get_levels(x, store_levels = FALSE, max_levels_store = 200)
  
  expect_null(levels_result)
})

test_that(".get_levels() returns NULL for too many levels", {
  x <- factor(paste0("Level", 1:300))
  levels_result <- .get_levels(x, store_levels = TRUE, max_levels_store = 200)
  
  expect_null(levels_result)
})

test_that(".get_levels() returns NULL for Date and POSIXct", {
  x <- as.Date(c("2024-01-01", "2024-01-02"))
  levels_result <- .get_levels(x, store_levels = TRUE, max_levels_store = 200)
  
  expect_null(levels_result)
  
  x <- as.POSIXct(c("2024-01-01 12:00:00", "2024-01-02 12:00:00"))
  levels_result <- .get_levels(x, store_levels = TRUE, max_levels_store = 200)
  
  expect_null(levels_result)
})

test_that(".get_levels() handles NAs correctly", {
  x <- c("A", "B", NA, "C", NA)
  levels_result <- .get_levels(x, store_levels = TRUE, max_levels_store = 200)
  
  # Should not include NA in unique values
  expect_equal(levels_result, c("A", "B", "C"))
})

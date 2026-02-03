#!/usr/bin/env Rscript

# Test the updated .build_formula_string_random function
source('R/build_functions.R')

cat("=== Testing .build_formula_string_random() ===\n\n")

# Test 1: Backward compatibility
cat("1. Backward compatible - lme4 style:\n")
f1 <- .build_formula_string_random("~ 1 | Species")
cat("   Result:", class(f1), "\n")
cat("   Formula:", deparse(f1), "\n")
cat("   all.vars():", all.vars(f1), "\n\n")

# Test 2: MCMCglmm single
cat("2. MCMCglmm single random effect:\n")
f2 <- .build_formula_string_random("~ us(1):Species")
cat("   Formula:", deparse(f2), "\n")
cat("   all.vars():", all.vars(f2), "\n\n")

# Test 3: MCMCglmm multiple
cat("3. MCMCglmm multiple random effects:\n")
f3 <- .build_formula_string_random("~ us(1):Species + us(1):Genus")
cat("   Formula:", deparse(f3), "\n")
cat("   all.vars():", paste(all.vars(f3), collapse=", "), "\n\n")

# Test 4: Detailed output
cat("4. Detailed structure (return_list=TRUE):\n")
f4 <- .build_formula_string_random("~ us(1):Species + us(1):Genus", return_list = TRUE)
cat("   Clusters:", paste(f4$clusters, collapse=", "), "\n")
cat("   Number of terms:", length(f4$terms), "\n")
cat("   Term 1 - Structure:", f4$terms[[1]]$structure, "| Cluster:", f4$terms[[1]]$cluster, "\n")
cat("   Term 2 - Structure:", f4$terms[[2]]$structure, "| Cluster:", f4$terms[[2]]$cluster, "\n\n")

# Test 5: Different structures
cat("5. Different variance structures:\n")
f5 <- .build_formula_string_random("~ idh(trait):Species", return_list = TRUE)
cat("   Structure:", f5$terms[[1]]$structure, "\n")
cat("   Variables:", paste(f5$terms[[1]]$variables, collapse=", "), "\n\n")

cat("All tests passed!\n")

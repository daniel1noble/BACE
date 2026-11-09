# =================================================================================================
# rubin_missing.R -- which datasets of the campaign have not been produced yet?
#
# WHY THIS EXISTS: a task can die without writing anything (a segfault, a timeout, a node going
# away). The campaign is resumable, so the cure is simply to run the missing ones again -- but
# first you have to know which they are.
#
# HOW IT WORKS, and why this way. The old version of this script kept its OWN copy of the design
# (n = 100/300/1000, lambda = 0.3/0.7/1, rho = 0/0.5). When the design changed, this copy did not,
# and it would have reported the new cells as missing forever while quietly ignoring the old ones.
#
# So it no longer keeps a copy. It reads the task list the launcher wrote -- tasks.txt, one
# rubin_cell.R command per dataset -- works out the filename each of those commands would produce,
# and checks whether that file is there. One definition of the design, in one place, and this
# script cannot drift from it.
#
# USAGE
#   Rscript code/rubin_missing.R <results_dir>
#
# OUTPUT
#   a short summary on screen, and <results_dir>/missing_tasks.txt -- the subset of tasks.txt that
#   still needs running. Feed it straight back to the runner:
#       xargs -P 7 -L 1 Rscript code/rubin_cell.R < <results_dir>/missing_tasks.txt
# =================================================================================================

a <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(a)) a[1] else stop("usage: Rscript code/rubin_missing.R <results_dir>")

## STEP 1. Find the task list the launcher wrote. Without it we cannot know what was meant to run.
tasks_file <- file.path(results_dir, "tasks.txt")
if (!file.exists(tasks_file)) {
  stop("no tasks.txt in ", results_dir, " -- run rubin_campaign_local.sh (DRY_RUN=1 is enough) first")
}
tasks <- readLines(tasks_file)
tasks <- tasks[nzchar(trimws(tasks))]

## STEP 2. Work out the filename a given task line would write. The name itself is built by
##         cell_tag() in rubin_lib.R -- the SAME function rubin_cell.R uses to decide what to
##         write -- so this script cannot disagree with the runner about what a dataset is called.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "rubin_lib.R"))

flag <- function(line, name, default) {
  m <- regmatches(line, regexpr(paste0("--", name, "\\s+[^ ]+"), line))
  if (!length(m)) return(default)
  sub(paste0("--", name, "\\s+"), "", m)
}
expected_file <- function(line) {
  n      <- flag(line, "n", "100")
  seed   <- flag(line, "seed", "1")
  lambda <- flag(line, "lambda", "0.7")
  rho    <- flag(line, "rho", "0")
  miss   <- flag(line, "miss", "mcar")
  frac   <- flag(line, "frac", "0.3")
  M      <- flag(line, "M", "20")
  n_cov  <- flag(line, "n_cov", "0")
  rho_cov<- flag(line, "rho_cov", "0.8")
  out    <- flag(line, "out", results_dir)
  file.path(out, paste0(cell_tag(lambda, rho, miss, frac, n, M, seed,
                                 n_cov = n_cov, rho_cov = rho_cov), ".rds"))
}

## STEP 3. Check each expected file. This is the whole job.
want <- vapply(tasks, expected_file, character(1), USE.NAMES = FALSE)
have <- file.exists(want)

## STEP 4. Report. Broken down by tree size, because that is how the work is split up and how you
##         would decide what to relaunch.
cat(sprintf("campaign : %s\n", results_dir))
cat(sprintf("expected : %d datasets\n", length(want)))
cat(sprintf("present  : %d\n", sum(have)))
cat(sprintf("missing  : %d  (%.1f%%)\n\n", sum(!have), 100 * mean(!have)))

n_of <- vapply(tasks, function(l) flag(l, "n", "?"), character(1), USE.NAMES = FALSE)
for (nn in unique(n_of)) {
  i <- n_of == nn
  cat(sprintf("  n = %-5s expected %6d  present %6d  missing %6d\n", nn, sum(i), sum(have[i]), sum(!have[i])))
}

## STEP 5. Write the missing ones back out in the SAME format as tasks.txt, so relaunching them is
##         one command and needs no editing.
out_file <- file.path(results_dir, "missing_tasks.txt")
if (any(!have)) {
  writeLines(tasks[!have], out_file)
  cat(sprintf("\nwrote %s\n", out_file))
  cat(sprintf("rerun with:  xargs -P 7 -L 1 Rscript code/rubin_cell.R < %s\n", out_file))
} else {
  if (file.exists(out_file)) unlink(out_file)
  cat("\nnothing missing: the campaign is complete.\n")
}

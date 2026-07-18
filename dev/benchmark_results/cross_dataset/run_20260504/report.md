# BACE cross-dataset benchmark — 20260504

- BACE commit: `730eefb2`
- BACE version: 0.0.0.9000
- Datasets: globtherm, pantheria

> **Warning:** Missing output for: avonet, amphibio, bien, leptraits. Check the workflow run logs.

## Per-dataset summary (BACE)

|dataset   | n_traits| mean_NRMSE| mean_cor| mean_cov95| mean_accuracy| mean_bal_acc| mean_brier|
|:---------|--------:|----------:|--------:|----------:|-------------:|------------:|----------:|
|globtherm |       10|      0.924|    0.490|      0.243|            NA|           NA|         NA|
|pantheria |        8|      2.713|    0.755|      0.057|         0.786|        0.613|      0.316|

## BACE vs mean-baseline head-to-head

Pigauto-style reference. `lift_pct` = improvement vs the column-mean / modal-class baseline. Positive = BACE better. For RMSE / NRMSE: `100 * (baseline - bace) / baseline`. For accuracy: `100 * (bace - baseline)`.

### Continuous / count traits

|dataset   |trait         |type       | bace_rmse| base_rmse| lift_pct| bace_cor|
|:---------|:-------------|:----------|---------:|---------:|--------:|--------:|
|globtherm |elevation_max |continuous |        NA|   953.989|       NA|    0.406|
|globtherm |elevation_max |continuous |   827.066|   953.989|     13.3|    0.430|
|globtherm |lat_max       |continuous |        NA|    28.376|       NA|    0.225|
|globtherm |lat_max       |continuous |    27.803|    28.376|      2.0|    0.426|
|globtherm |long_max      |continuous |        NA|    93.604|       NA|    0.622|
|globtherm |long_max      |continuous |   100.231|    93.604|     -7.1|    0.292|
|globtherm |Tmax          |continuous |        NA|     7.757|       NA|    0.728|
|globtherm |Tmax          |continuous |     7.112|     7.757|      8.3|    0.387|
|globtherm |Tmin          |continuous |        NA|    12.852|       NA|    0.589|
|globtherm |Tmin          |continuous |     8.326|    12.852|     35.2|    0.795|

## BACE continuous + count traits (per-trait detail)

|dataset   |trait               |type       |scale | n_hidden|    rmse| nrmse| mae_fit|   mae_raw| correlation| coverage95|
|:---------|:-------------------|:----------|:-----|--------:|-------:|-----:|-------:|---------:|-----------:|----------:|
|globtherm |Tmax                |continuous |raw   |       11|      NA| 0.625|   3.689|     3.318|       0.728|      0.273|
|globtherm |Tmin                |continuous |raw   |       10|      NA| 0.902|   8.748|     8.664|       0.589|      0.100|
|globtherm |lat_max             |continuous |raw   |       11|      NA| 1.191|  34.521|    29.998|       0.225|      0.182|
|globtherm |long_max            |continuous |raw   |       11|      NA| 0.987|  77.602|    72.855|       0.622|      0.182|
|globtherm |elevation_max       |continuous |raw   |        8|      NA| 0.754| 507.666|   457.542|       0.406|      0.375|
|globtherm |Tmax                |continuous |raw   |       33|   7.112| 1.101|   5.637|     4.816|       0.387|      0.333|
|globtherm |Tmin                |continuous |raw   |       31|   8.326| 0.534|   7.326|     7.146|       0.795|      0.161|
|globtherm |lat_max             |continuous |raw   |       33|  27.803| 0.868|  25.849|    24.283|       0.426|      0.212|
|globtherm |long_max            |continuous |raw   |       33| 100.231| 1.043|  88.753|    84.027|       0.292|      0.152|
|globtherm |elevation_max       |continuous |raw   |       26| 827.066| 1.235| 611.637|   504.195|       0.430|      0.462|
|pantheria |body_mass_g         |continuous |log   |       20|      NA| 1.983|   4.606| 24059.268|       0.906|      0.000|
|pantheria |head_body_length_mm |continuous |log   |       11|      NA| 3.547|   3.621|   309.031|       0.857|      0.000|
|pantheria |gestation_d         |continuous |log   |        7|      NA| 2.951|   2.789|   100.492|       0.875|      0.000|
|pantheria |max_longevity_m     |continuous |log   |        5|      NA| 3.739|   3.523|   235.992|       0.570|      0.000|
|pantheria |litter_size         |count      |raw   |       14|      NA| 1.345|   1.750|     1.536|       0.569|      0.286|

## BACE categorical / binary / ordinal traits (per-trait detail)

|dataset   |trait           |type    | n_hidden| accuracy| balanced_accuracy| brier| mae_level|
|:---------|:---------------|:-------|--------:|--------:|-----------------:|-----:|---------:|
|pantheria |diet_breadth    |ordinal |       11|    0.545|             0.375| 0.636|     0.636|
|pantheria |habitat_breadth |ordinal |       16|    0.812|             0.464| 0.312|     0.188|
|pantheria |terrestriality  |binary  |       15|    1.000|             1.000| 0.000|        NA|

## Pre-imputation phylogenetic signal

Reference values per trait, computed before BACE runs. Pagel λ and Blomberg K for continuous; Fritz-Purvis D (OVR mean) for categorical.

|dataset   |trait               |type       | lambda|     K|D  |
|:---------|:-------------------|:----------|------:|-----:|:--|
|globtherm |Tmax                |continuous |  0.874| 0.036|NA |
|globtherm |Tmin                |continuous |  0.914| 0.127|NA |
|globtherm |lat_max             |continuous |  0.274| 0.013|NA |
|globtherm |long_max            |continuous |  0.619| 0.046|NA |
|globtherm |elevation_max       |continuous |  0.197| 0.018|NA |
|globtherm |Tmax                |continuous |  0.716| 0.057|NA |
|globtherm |Tmin                |continuous |  0.882| 0.151|NA |
|globtherm |lat_max             |continuous |  0.517| 0.030|NA |
|globtherm |long_max            |continuous |  0.284| 0.014|NA |
|globtherm |elevation_max       |continuous |  0.355| 0.038|NA |
|pantheria |body_mass_g         |continuous |  0.782| 0.118|NA |
|pantheria |head_body_length_mm |continuous |  0.879| 0.276|NA |
|pantheria |gestation_d         |continuous |  0.944| 0.513|NA |
|pantheria |max_longevity_m     |continuous |  0.897| 0.327|NA |
|pantheria |litter_size         |count      |  0.671| 0.042|NA |
|pantheria |diet_breadth        |ordinal    |     NA|    NA|NA |
|pantheria |habitat_breadth     |ordinal    |     NA|    NA|NA |
|pantheria |terrestriality      |binary     |     NA|    NA|NA |

## Runtime + MCMC config

|dataset   | n_species| n_traits| runtime_min|converged | n_attempts| nitt| thin| burnin| runs| n_final|
|:---------|---------:|--------:|-----------:|:---------|----------:|----:|----:|------:|----:|-------:|
|globtherm |       150|        5|         0.1|FALSE     |          1| 1500|    5|    500|    2|       2|
|globtherm |       150|        5|         0.1|FALSE     |          1| 1500|    5|    500|    2|       2|
|pantheria |       300|        8|         0.2|FALSE     |          1| 1500|    5|    500|    2|       2|


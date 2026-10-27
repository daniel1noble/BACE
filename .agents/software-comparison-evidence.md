# Software-comparison evidence file (ms Table: imputation tools in R)

**Purpose.** Every cell of the manuscript's software-comparison table
(`ms/ms.qmd`, `#tbl-software`) is backed here by page-cited evidence from the
package's OWN reference documentation. Nothing in the table is asserted from
memory or reputation. This file is the audit trail for co-authors and for
referee responses; it can become a supplementary file at submission.

**Method (2026-10-27).** One reader agent per package downloaded the CRAN
reference-manual PDF (or, for GitHub-only packages, README/DESCRIPTION/man
pages) and extracted, with page-cited quotes: approach; proper multiple
imputation (multiple imputed datasets or an integrated Bayesian posterior,
with between-imputation uncertainty usable downstream); multilevel support;
phylogeny support; and per-type imputation support (continuous, count/Poisson,
binary, ordinal-as-ordered, nominal). Anything the documentation does not
state is recorded as "not stated" rather than inferred. The decisive cells
(mice count; miceadds nominal/ordinal; JointAI model types and random-effect
covariance structures; smcfcs method list; Rphylopars prediction output; hmi
CRAN archival) were independently re-read from the PDFs by the lead session.
Manuals are re-downloadable from CRAN (`.../web/packages/<pkg>/<pkg>.pdf`).

**Conventions.** "ordinal" = ordered categorical imputed by an order-aware
model (cumulative/proportional-odds/ordinal-probit), not merely a factor
treated as unordered. "count" = an explicit count/Poisson-type imputation
model, not counts routed through PMM. "proper MI" per Rubin (1987): multiple
completed datasets (or full-Bayes integration) so that downstream analyses
reflect between-imputation uncertainty.

---

## BACE (this package)

- approach: chained equations (FCS); every conditional model is a Bayesian
  phylogenetic mixed model (MCMCglmm)
- proper_mi: yes - n_final independent posterior-predictive imputed datasets
  (`bace_final_imp`, default 50); pooling via stacked posteriors
  (`pool_posteriors`) or Rubin's rules (`with_imputations` + `pool_mi`)
- multilevel: partial - phylogenetic random effect plus an optional
  non-phylogenetic species effect for replicated taxa (`species = TRUE`);
  arbitrary additional grouping factors are NOT supported (honest limit)
- phylogeny: yes - the tree enters every conditional model as a covariance
  structure (`MCMCglmm` ginverse via `ape`/`inverseA`)
- continuous: yes (gaussian); count: yes (poisson, log link); binary: yes
  (probit threshold); ordinal: yes (ordered threshold/probit with cutpoints);
  nominal: yes (multinomial probit; one-vs-rest optional)
- evidence: this repository - `R/prep_functions.R` (`.detect_type`),
  `R/model_functions.R` (per-family prediction), `tests/` (1418 tests),
  `dev/simulation_results/SIMULATION_REPORT.qmd` (all five types validated
  under MCAR/MAR)


## mice (v3.19.0, 2025-12-09, CRAN)

- approach: FCS/chained equations (MICE algorithm): univariate conditional models per variable (pmm, GLM-type, trees)
- proper_mi: yes — returns mids with m imputed datasets; Rubin's-rules pool()
- multilevel: yes — two-level normal and logistic methods (2l.norm/2l.pan/2l.lmer/2l.bin)
- phylogeny: no — never mentioned; verified absence, no tree-based covariance method
- continuous: yes — pmm (default for numeric) and Bayesian linear regression (norm)
- count: no — no dedicated Poisson/count method documented anywhere in the manual; defaultMethod covers only numeric/binary/unordered/ordered types ("poisson": 0 full-text hits)
- binary: yes — Bayesian logistic regression (logreg), default for 2-level factors
- ordinal: yes — proportional odds / cumulative link model (polr), default for ordered factors
- nominal: yes — Bayesian polytomous (multinomial) regression via nnet::multinom (polyreg), default for unordered factors
- reference (as stated in the docs): van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. doi:10.18637/jss.v045.i03 (cited on manual p.1 Description and throughout References sections). CRAN: https://CRAN.R-project.org/package=mice; GitHub: https://github.com/amices/mice

**Evidence (page-cited):**

- approach -> DESCRIPTION, PDF page 1: 'Multiple imputation using Fully Conditional Specification (FCS) implemented by the MICE algorithm as described in Van Buuren and Groothuis-Oudshoorn (2011)... Each variable has its own imputation model.'
- proper_mi -> pool, PDF page 160: 'The pool() function combines the estimates from m repeated complete data analyses... 1. Impute the missing data by the mice() function, resulting in a multiple imputed data set (class mids)'; PDF page 161: 'computes the total variance over the repeated analyses by Rubin's rules (Rubin, 1987, p. 76)'
- proper_mi -> mids class, PDF pages 138-139: 'The mids object is an S3 class that represents a multiply imputed data set'; argument m: 'Number of multiple imputations. The default is m=5.'
- multilevel -> DESCRIPTION, PDF page 1: 'MICE can also impute continuous two-level data (normal model, pan, second-level variables)'
- multilevel -> mice.impute.2l.*, PDF pages 80-85: 2l.bin 'two-level logistic model using lme4::glmer()' (p.80); 2l.lmer 'two-level normal model using lme4::lmer()' (p.81); 2l.norm 'Implements the Gibbs sampler for the linear multilevel model with heterogeneous with-class variance (Kasim and Raudenbush, 1998)' (p.83); 2l.pan 'two-level normal model with homogeneous within group variances' via pan (p.84)
- phylogeny -> whole manual: pdftotext full-text search, 'phylo' occurs 0 times; complete mice.impute.* function index (PDF pages 198-201) contains no phylogenetic/tree-covariance method
- continuous -> mice.impute.norm, PDF page 115: 'Calculates imputations for univariate missing data by Bayesian linear regression, also known as the normal model' (Rubin 1987, p.167 procedure); mice.impute.pmm, PDF page 123: 'Imputation by predictive mean matching'
- continuous (default) -> make.method defaultMethod, PDF page 69: 'By default, the method uses pmm, predictive mean matching (numeric data)...'
- count -> make.method defaultMethod, PDF page 69: defaults cover only '1) numeric data, 2) factor data with 2 levels, 3) factor data with > 2 unordered levels, and 4) factor data with > 2 ordered levels' — no count type; no Poisson/count method in the complete function index (PDF pages 198-201); 'poisson' = 0 hits in full manual text; the pmm entry (pp.123-125) never mentions count data
- binary -> mice.impute.logreg, PDF page 104: 'Imputes univariate missing data using logistic regression'; Details: 'Imputation for binary response variables by the Bayesian logistic regression model (Rubin 1987, p. 169-170)'
- ordinal -> mice.impute.polr, PDF pages 126-127: title 'Imputation of ordered data by polytomous regression'; Details: 'imputes for ordered categorical response variables by the proportional odds logistic regression (polr) model... also known as the cumulative link model'; 'By default, ordered factors with more than two levels are imputed by mice.impute.polr' (uses MASS::polr, multinom fallback if polr fails)
- nominal -> mice.impute.polyreg, PDF pages 128-129: title 'Imputation of unordered data by polytomous regression'; Details: 'imputes categorical response variables by the Bayesian polytomous regression model'; 'By default, unordered factors with more than two levels are imputed by mice.impute.polyreg()'; 'uses the function multinom() from the nnet package'
- reference -> DESCRIPTION, PDF page 1: doi:10.18637/jss.v045.i03; References (e.g. PDF pages 83, 105, 125): 'Van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67.'
- package/version -> PDF pages 1-2: Package 'mice', Version 3.19.0, Date 2025-12-09, URL https://github.com/amices/mice, Repository CRAN

**Not stated / unclear in the docs:**

- count: the 'no' is verified absence, not an explicit disclaimer — the manual never discusses count/Poisson data at all (0 full-text hits for 'poisson'; 'count' appears only incidentally, e.g. a dataset's 'Trauma count (1-5)' variable); whether counts should be routed through pmm is not stated
- phylogeny: verified absence — 'phylo' never appears anywhere in the manual; the docs contain no explicit statement for or against phylogenetic data
- the main mice() function help entry is absent from this PDF build (not in the TOC pages 2-5 or as a bold Index entry pages 198-201), so per-type default methods were cited from make.method's defaultMethod argument (PDF page 69) and MI structure from the mids class (pp.138-139) and pool (pp.160-161) pages
- ordinal caveat stated by the docs themselves (mice.impute.polr Note, PDF page 127): in mice 3.0.0-3.6.6 polr 'could always fail silently', in which case requests for polr were in fact handled by multinom (nominal treatment); fixed per https://github.com/amices/mice/issues/206

---

## miceadds (v3.20-10, 2026-05-28; Robitzsch, Grund & Henke)

- approach: FCS add-on methods for mice; lme4/blme mixed-model imputation (ml.lmer, 2l.*)
- proper_mi: yes — methods run inside mice::mice (m datasets); pool_mi pools
- multilevel: yes — hierarchical/cross-classified, arbitrary levels, random slopes (lme4)
- phylogeny: no — never mentioned; clustering only via categorical cluster/level identifiers
- continuous: yes — ml.lmer model="continuous" (normally distributed); 2l.continuous linear mixed model (lme4::lmer/blme::blmer)
- count: no — no Poisson/count model documented anywhere in the manual; ml.lmer model arg has no count option and smcfcs dep_type offers only norm/lognorm/yj/bc/logistic
- binary: yes — ml.lmer model="binary" = logistic mixed effects model; 2l.binary two-level logistic (lme4::glmer/blme::bglmer)
- ordinal: partial — ml.lmer accepts "ordinal data (model=\"pmm\")", i.e. handled by predictive mean matching, not a cumulative/proportional-odds/ordinal-probit imputation model; ml_mcmc fits multilevel ordinal probit but is documented only as model estimation, not an imputation method
- nominal: partial — multilevel ml.lmer explicitly cannot handle nominal ("Nominal variables with missing values cannot (yet) be handled"); only single-level mice.impute.catpmm imputes a categorical variable via dummy-coded multivariate PMM (L1-distance matching, not a multinomial model; no cluster argument). So no method handles unordered multinomial at level 1 of a multilevel model
- reference (as stated in the docs): Robitzsch, A., & Grund, S. (2026). miceadds: Some Additional Multiple Imputation Functions, Especially for 'mice'. R package version 3.20-10. CRAN: https://CRAN.R-project.org/package=miceadds; GitHub: https://github.com/alexanderrobitzsch/miceadds. Methodological reference stated in docs: Grund, Luedtke & Robitzsch (2018), Organizational Research Methods, 21(1), 111-149

**Evidence (page-cited):**

- proper_mi -> package DESCRIPTION/title, PDF page 1: "Some Additional Multiple Imputation Functions, Especially for 'mice'... Contains functions for multiple imputation which complements existing functionality in R. In particular, several imputation methods for the mice package... are implemented"; mice.1chain, PDF page 63: "modifies the mice::mice function to multiply impute a dataset using a long chain" (nimp = number of imputations); pool_mi, PDF page 143: "Statistical inference for multiply imputed datasets"; examples call mice::mice(dat1, method=method, m=5) (PDF page 78)
- multilevel -> miceadds-package Details, PDF page 5: "mice.impute.ml.lmer is a general function for imputing multilevel data with hierarchical or cross-classified structures for variables at an arbitrary level. This imputation method uses the lme4::lmer function"; mice.impute.ml.lmer, PDF page 82: imputation model "y = X beta + sum_v Z_v u_v for V multivariate random effects... random slopes can be selected with random_slopes"
- phylogeny -> entire manual, PDF pages 1-192: zero occurrences of 'phylog' (verified by full-text search of the decompressed PDF); structure is specified only via cluster identifiers, e.g. levels_id "Specification of the level identifiers" (PDF page 83) and cluster_var "Cluster identifier for Level 2 units" (PDF page 73) — no tree-based covariance anywhere
- continuous -> mice.impute.ml.lmer model argument, PDF page 83: "Can be 'continuous' for normally distributed data"; mice_imputation_2l_lmer, PDF page 110: "mice.impute.2l.continuous imputes values of continuous variables with a linear mixed effects model using lme4::lmer or blme::blmer"
- count -> mice.impute.ml.lmer model argument, PDF page 83 lists only "continuous", "binary", "pmm" (no count/Poisson option); mice.impute.smcfcs dep_type, PDF page 101: "Possible choices are 'norm'..., 'lognorm'..., 'yj'..., 'bc'..., 'logistic'"; full-text search of the whole manual: 0 occurrences of 'poisson' or 'count data'
- binary -> mice.impute.ml.lmer model argument, PDF page 83: "'binary' for dichotomous data specifying a logistic mixed effects model"; mice_imputation_2l_lmer, PDF page 110: "Binary variables can be imputed from a two-level logistic regression model fitted with the lme4::glmer or blme::bglmer function"
- ordinal -> mice.impute.ml.lmer Description, PDF page 82: "The function now only allows continuous data (model='continuous'), ordinal data (model='pmm') or binary data (model='pmm' or model='binary')" — ordinal only via predictive mean matching; ml_mcmc, PDF page 123: "Fits a mixed effects model via MCMC. The outcome can be normally distributed or ordinal" (outcome 'normal' or 'probit') but documented as estimation only — full-text search shows no mice.impute.* function references ml_mcmc
- nominal -> mice.impute.ml.lmer Description, PDF page 82: "Nominal variables with missing values cannot (yet) be handled."; mice.impute.catpmm, PDF page 77: "Imputes a categorical variable using multivariate predictive mean matching... The categorical outcome variable is recoded as a vector of dummy variables. A multivariate linear regression is specified for computing predicted values. The L1 distance... is utilized for predictive mean matching" — usage (y, ry, x, donors, ridge) has no cluster argument, i.e. single-level; the only 'polyreg'/'polr' mention is mice.1chain's pass-through defaultMethod=c('pmm','logreg','polyreg','polr') from mice itself (PDF page 63)

**Not stated / unclear in the docs:**

- count imputation: not stated — no Poisson/count-type method appears anywhere in the manual (absence verified by full-text search, so reported as 'no explicit count method' rather than inferred capability)
- phylogeny: not stated — the manual never mentions phylogenetic trees or tree-based covariance; 'no' reflects absence from the docs, not an explicit negative statement
- whether ml_mcmc (multilevel ordinal probit MCMC, p.123) can serve as an imputation engine: not stated — documented only as 'MCMC Estimation for Mixed Effects Model' with coef/vcov/summary methods; no imputation method in the manual references it
- whether mice.impute.catpmm respects multilevel/cluster structure: not stated — its arguments (y, ry, x, donors, ridge) include no cluster identifier
- random-slopes support in cross-classified ml.lmer models beyond the stated iter_re note (p.83: 'Using iter_re>0 is necessary for cross-classified models with not fully balanced designs', stated for random intercept models with continuous outcomes): not stated in further detail

---

## micemd 1.11.0 (CRAN, 2026-09-28; Audigier & Resche-Rigon)

- approach: mice FCS add-on: two-level univariate methods (two-stage estimators, Bayesian GLMMs, FCS-jomo)
- proper_mi: yes — m imputed datasets (default m=5), pooled downstream
- multilevel: yes — two-level random-effect models via cluster indicator
- phylogeny: no — phylogeny/tree-based covariance never mentioned in manual
- continuous: yes — 2l.2stage.norm, 2l.glm.norm, 2l.jomo, 2l.2stage.pmm (heteroscedastic/homoscedastic normal mixed models)
- count: yes — explicit Poisson: 2l.2stage.pois (two-stage Poisson) and 2l.glm.pois (Bayesian Poisson mixed model)
- binary: yes — 2l.2stage.bin (two-stage logistic), 2l.glm.bin (Bayesian logistic mixed), 2l.jomo (probit/latent variable); 2l.2stage.heckman (binary or continuous, MNAR)
- ordinal: partial — no multilevel 2l.* ordinal method; only mice's single-level polr (proportional odds) as a documented mice.par default
- nominal: partial — no multilevel 2l.* nominal method; only mice's single-level polyreg (polytomous regression) as a documented mice.par default
- reference (as stated in the docs): Audigier V., White I., Jolani S., Debray T., Quartagno M., Carpenter J., van Buuren S., Resche-Rigon M. (2018). Multiple imputation for multilevel data with continuous and binary variables. Statistical Science, doi:10.1214/18-STS646 (cited in DESCRIPTION and package help); for the 2stage methods also Resche-Rigon M. & White I.R. (2016). Statistical Methods in Medical Research 27(6):1634-1649, doi:10.1177/0962280216666564. CRAN: https://cran.r-project.org/package=micemd

**Evidence (page-cited):**

- proper_mi -> mice.par argument m, PDF page 23: "Number of multiple imputations. The default is m=5."; package title p.1: "Multiple Imputation by Chained Equations with Multilevel Data"; package example p.3 pools with pool(ana)
- multilevel -> Description, PDF page 1: "Addons for the 'mice' package to perform multiple imputation using chained equations with two-level data"; 2l.2stage.norm p.13: cluster-wise fits "combined according to a linear random effect model"
- phylogeny -> whole manual (pdftotext grep), 0 hits for "phylo"; only "tree" is "decision tree" (find.defaultMethod, PDF page 6); clustering is a cluster-indicator variable: "The class variable (only one is allowed) is coded as '-2'" (p.9)
- continuous -> mice.impute.2l.2stage.norm, PDF page 13: "Imputes univariate two-level continuous variable from a heteroscedastic normal model"; 2l.glm.norm p.19 (Bayesian linear mixed model); 2l.jomo p.21: "continuous or binary incomplete variable"
- count -> mice.impute.2l.2stage.pois, PDF page 16: "Imputes univariate two-level count variable from a Poisson model"; mice.impute.2l.glm.pois p.20: "Bayesian mixed model (Poisson regression)... dedicated to a count outcome"; find.defaultMethod p.6 labels both "(integer)"
- binary -> mice.impute.2l.2stage.bin, PDF page 9: "Imputes univariate two-level binary variable from a logistic model"; 2l.glm.bin p.17: "Bayesian logistic mixed model... dedicated to a binary outcome"; 2l.jomo p.21: "For a binary variable, a probit link and a latent variables framework are used"
- ordinal -> mice.par defaultMethod, PDF page 24: "polr, proportional odds model for (ordered, >= 2 levels)" (mice's single-level default passed through the wrapper); no 2l.* ordinal method in TOC p.2 or find.defaultMethod p.6 (rows: continuous, binary, integer only)
- nominal -> mice.par defaultMethod, PDF page 24: "polyreg, polytomous regression imputation for unordered categorical data (factor >= 2 levels)" (mice's single-level default); no 2l.* nominal method anywhere (TOC p.2; Description p.1: "Imputation of continuous, binary or count variables are available")
- package scope -> Description, PDF page 1: "Imputation of continuous, binary or count variables are available"; find.defaultMethod p.6: "Returned methods can be: 2l.stage.bin (binary), 2l.stage.norm (continuous), 2l.stage.pois (integer), 2l.glm.bin (binary), 2l.glm.norm (continuous), 2l.glm.pois (integer), 2l.jomo (continuous or binary)"

**Not stated / unclear in the docs:**

- Multilevel ordinal imputation: not stated — no 2l.* ordinal method documented; polr mention is only mice's single-level default inside mice.par's defaultMethod argument
- Multilevel nominal imputation: not stated — no 2l.* nominal/multinomial method documented; polyreg mention is only mice's single-level default inside mice.par
- Phylogenetic structure: not stated — no mention of phylogeny or tree-based covariance anywhere in the 32-page manual
- Overdispersion handling for the Poisson count methods: not stated
- Doc quirk: p.9 2l.2stage.bin Details repeats the 2stage.norm 'continuous... heteroscedastic normal model' text (apparent copy-paste error); its title/Description clearly state binary logistic

---

## pan

- approach: Gibbs sampler, joint multivariate linear mixed model (two-level panel/clustered MI)
- proper_mi: yes: multiple posterior-predictive draws via repeated chain restarts (one draw per call)
- multilevel: yes: two-level subject/cluster random effects
- phylogeny: no: unstructured psi only; phylogeny never mentioned
- continuous: yes: responses modeled multivariate normal (linear mixed model)
- count: no: no count/Poisson imputation method documented anywhere in manual
- binary: no: no binary/threshold method documented; all responses normal
- ordinal: no: no ordinal/cumulative/probit method documented anywhere in manual
- nominal: no: categorical variables never mentioned anywhere in manual
- reference (as stated in the docs): Schafer JL, Yucel RM (2002). Computational strategies for multivariate linear mixed-effects models with missing values. Journal of Computational and Graphical Statistics 11:437-457 (listed in pan() References, p.9, alongside Schafer JL 1997, Tech report 97-04, Penn State, and Schafer JL 2001, Multiple imputation with PAN, APA). CRAN: https://cran.r-project.org/package=pan

**Evidence (page-cited):**

- proper_mi -> pan() Description, PDF page 7: 'This function will typically be used to produce multiple imputations of missing data values in multivariate panel data or clustered data'; Value component y, PDF page 9: 'If "iter" has been set large enough ... this is a proper draw from the posterior predictive distribution of the complete data' - one completed dataset per call; the worked example (PDF page 11) imputes 'm=10 times taking 100 steps between imputations' by restarting the chain (start=result$last) and combines the 10 datasets with mitools::MIcombine
- multilevel -> pan() Description, PDF page 7: model 'yi = Xi%*%beta + Zi%*%bi + ei' with 'bi = (q x r) matrix of coefficients specific to subject or cluster i (random effects)'; subj argument, PDF page 8: 'vector of length nrow(y) giving the subject (or cluster) indicators i for the rows of y' - exactly one clustering level (two-level data)
- phylogeny -> pan() Description, PDF page 7: random effects 'normally distributed with mean zero and unstructured covariance matrix psi'; no phylogenetic tree or tree-based covariance is mentioned anywhere in the 14-page manual (pan.bd, PDF page 12, only restricts psi to block-diagonal)
- continuous -> pan() Description, PDF page 7: 'Gibbs sampler for the multivariate linear mixed model with incomplete data described by Schafer (1997)'; 'the rows of ei are assumed to be independently normal with mean zero and unstructured covariance matrix sigma. Missing values may appear in yi in any pattern'
- count -> whole manual (Contents, PDF page 1: only bitest, ecme, marijuana, pan, pan.bd, YWC.data): no count/Poisson-type imputation model is documented anywhere; the only imputation model stated is the normal multivariate linear mixed model (pan(), PDF page 7)
- binary -> whole manual: no binary/threshold/logistic imputation method documented anywhere; pan() (PDF page 7) models all response columns as jointly normal; y is 'matrix of responses' with no factor/binary handling stated (PDF page 8)
- ordinal -> whole manual: no ordinal/cumulative-link/ordinal-probit method documented on any page; ordered categorical data never mentioned
- nominal -> whole manual: no multinomial/unordered-categorical method documented on any page; the words categorical/factor/multinomial do not appear in any function entry (PDF pages 1-13 read in full)
- approach -> package DESCRIPTION, PDF page 1: Title 'Multiple Imputation for Multivariate Panel or Clustered Data'; Description: 'Gibbs sampler for multivariate linear mixed models with incomplete data, as described in Schafer JL (1997) Imputation of missing covariates under a multivariate linear mixed model'
- reference -> pan() References, PDF page 9: 'Schafer JL, Yucel RM (2002). Computational strategies for multivariate linear mixed-effects models with missing values. Journal of Computational and Graphical Statistics. 11:437-457'

**Not stated / unclear in the docs:**

- Explicit 'continuous only' statement: not stated - the manual never explicitly says categorical/count variables are unsupported; the continuous-only restriction is carried entirely by the fully specified normal linear mixed model (pan(), PDF page 7), so the categorical cells are 'no method documented' rather than an explicit prohibition
- Clustering beyond two levels (e.g. three-level nesting): not stated - only a single subject/cluster indicator (subj) is documented
- Whether non-numeric (factor) input to y is rejected or coerced: not stated
- Default/recommended number of imputations m: not stated - the m=10 in the example is illustrative ('With a sample this small we can afford to be cautious', PDF page 11), not a package default

---

## jomo (v2.7-6, 2023-04-13; Matteo Quartagno, James Carpenter)

- approach: Joint-modelling MI: multivariate latent-normal joint model, Gibbs sampler with Metropolis-Hastings covariance update (Carpenter & Kenward 2013 framework)
- proper_mi: yes — nimp imputations (default 5) stacked long; Rubin-rules pooling shown
- multilevel: yes — clustered (jomo1ran*) and 2-level (jomo2*) imputation, random effects
- phylogeny: no — phylogeny never mentioned; clusters are exchangeable indicator-defined groups
- continuous: yes — jomo1con and relatives impute continuous (gaussian joint-normal) outcomes
- count: not-stated — no Poisson/count method; family limited to gaussian/binomial
- binary: yes — binary handled via latent normal variables (jomo1cat; binomial probit in SMC functions)
- ordinal: partial — ordered treatment only as SMC substantive-model outcome (jomo.polr proportional-odds ordinal probit; jomo.clmm cumulative link); general joint model never distinguishes ordered from unordered ("ordered" appears nowhere; all categorical = n-1 latent normals)
- nominal: yes — n-category variables imputed via n-1 latent normals (jomo1cat/jomo1mix/jomo1rancat)
- reference (as stated in the docs): As stated in the package docs: Author Matteo Quartagno, James Carpenter (maintainer m.quartagno@ucl.ac.uk), CRAN manual p.1; methodological reference cited throughout the manual: Carpenter J.R., Kenward M.G. (2013), Multiple Imputation and its Application, Wiley, ISBN 978-0-470-74052-1 (DESCRIPTION cites doi:10.1002/9781119942283); CRAN: https://cran.r-project.org/package=jomo. NOTE: the Quartagno & Carpenter JSS paper is NOT cited anywhere in the CRAN reference manual itself.

**Evidence (page-cited):**

- proper_mi -> jomo(), PDF page 6: "nimp — Number of Imputations. Default is 5." and Value: "The only argument returned is the imputed dataset in long format. Column 'Imputation' indexes the imputations. Imputation number 0 are the original data."; page 7 example: "# Finally we show how to fit the model and combine estimate with Rubin's rules # Here we use mitml, other options are available in mice, mitools, etc"
- multilevel -> DESCRIPTION, PDF page 1: Title "Multilevel Joint Modelling Multiple Imputation"; jomo1ran, page 58: "JM Imputation of clustered data — a wrapper function linking the six 2-level JM Imputation functions" (args clus = cluster indicator, Z = covariates associated to random effects, u.start = random effects, l2cov = level-2 covariance matrix); jomo2, page 91: "JM Imputation of 2-level data" with Y (level-1 outcomes) and Y2 ("the partially observed level 2 variables", page 92)
- multilevel (cluster-specific covariances) -> jomo(), page 6, meth argument: common / "fixed" / "random" study-specific level-1 covariance matrices (jomo1ranconhr, jomo1rancathr, jomo1ranmixhr); jomo1ranconhr Details, page 78, cites Yucel R.M. (2011), Statistical Modelling 11(4), 351-370
- phylogeny -> whole-manual text search of the decompressed PDF: 0 occurrences of "phylo"; no function in the full contents list (pages 1-2) relates to phylogeny; clusters in examples are cities/schools (cldata page 3: "city — The cluster indicator vector")
- continuous -> jomo1con, PDF page 51: "JM Imputation of single level data with continuous variables only — Impute a single level dataset with continuous outcomes only. A joint multivariate model for partially observed data is assumed and imputations are generated through the use of a Gibbs sampler."
- count -> jomo.glm, PDF page 16, family argument: "One of either 'gaussian' or 'binomial'. For binomial family, a probit link is assumed." (same on page 20 for jomo.glmer and page 41 for jomo.smc: model "can currently be one of lm, glm (binomial), polr, coxph, lmer, clmm or glmer (binomial)"); whole-manual search: 0 occurrences of "poisson" or "count" — no count method stated
- binary -> DESCRIPTION, PDF page 1: "Novel aspects of 'jomo' are the possibility of handling binary and categorical data through latent normal variables"; jomo1cat, page 48: Y.cat = "categorical (or binary) responses", Y.numcat = "number of categories in each categorical (or binary) variable"
- ordinal -> jomo.polr, PDF page 36: title "Joint Modelling Imputation Compatible with Proportional Odds Ordinal Probit Regression"; Description: "substantive model compatible JM imputation, when the substantive model of interest is a simple ordinal regression model"; jomo.clmm, page 7: "Joint Modelling Imputation Compatible with Cumulative Link Mixed Model"; by contrast "ordered" occurs nowhere in the manual and jomo1cat (page 48) describes categorical imputation only via n-1 latent normals with no ordering
- nominal -> jomo1cat, PDF page 48: "Impute a single level dataset with categorical variables as outcomes... imputations are generated through the use of a Gibbs sampler where the covariance matrix is updated with a Metropolis-Hastings step"; beta.start: "For each n-category variable we have a fixed effect parameter for each of the n-1 latent normals"; Details page 49: algorithm "described in detail in Chapter 5 of Carpenter and Kenward (2013)"; mixed types: jomo1mix page 55 (Y.con + Y.cat)
- reference -> PDF page 1: "Author Matteo Quartagno, James Carpenter"; References sections throughout (e.g. pages 6, 49, 52, 56, 78): "Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Wiley, ISBN: 978-0-470-74052-1" (chapter-specific: Ch.3 for jomo1con p.52, Ch.5 for jomo1cat p.49, Ch.9 for jomo1ranconhr p.78); whole-manual search: "Quartagno" appears only in the Author/Maintainer field, no JSS citation present

**Not stated / unclear in the docs:**

- Quartagno & Carpenter JSS-paper citation: not stated in the CRAN reference manual (only the author field and the Carpenter & Kenward 2013 book are given); cite the JSS paper only from a source actually read
- count/Poisson support: not stated (zero mentions of poisson/count; families documented as gaussian/binomial only)
- phylogeny: not stated (zero mentions)
- whether an incomplete ORDERED factor in the general joint imputation model (jomo/jomo1cat/jomo1ran) receives ordered treatment: not stated — manual describes all categorical variables as factors imputed via n-1 latent normals, with ordinal-specific handling documented only for the substantive-model outcome in jomo.polr/jomo.clmm
- the label 'multinomial probit' for the latent-normal categorical model: not used by the manual (it says only 'latent normal variables' / 'n-1 latent normals'; Chapter 5 of Carpenter & Kenward 2013 is cited for the algorithm)
- link function for jomo.polr: title says 'Proportional Odds Ordinal Probit Regression' (probit, not logit) — manual gives no further distributional detail

---

## mitml 0.4-5 (manual dated 2023-03-08)

- approach: front end to pan/jomo: joint multivariate multilevel model MI via MCMC (not FCS)
- proper_mi: yes — m imputed datasets (default 10); Rubin's-rules pooling (testEstimates)
- multilevel: yes — two-level models, cluster indicator plus random effects
- phylogeny: no — never mentioned; only cluster-based random effects described
- continuous: yes — panImpute: continuous multilevel data, missing at level 1; jomoImpute: (mixed) continuous data at levels 1 and 2
- count: not-stated — no count/Poisson imputation method appears anywhere in the manual; stated type support is continuous (pan) and continuous+categorical (jomo) only
- binary: yes — only as 'categorical' data formatted as factor variables (jomoImpute); binary is never named separately
- ordinal: not-stated — manual says only 'categorical' (factor); never distinguishes ordered factors or names an ordinal/cumulative/probit treatment
- nominal: yes — jomoImpute imputes 'categorical' variables formatted as factor variables; the multinomial/latent mechanism itself is not described (deferred to jomo)
- reference (as stated in the docs): Grund, S. [aut,cre], Robitzsch, A. [aut] & Luedtke, O. [aut] (2023). mitml: Tools for Multiple Imputation in Multilevel Modeling. R package version 0.4-5. CRAN: https://CRAN.R-project.org/package=mitml; GitHub: https://github.com/simongrund1/mitml (BugReports, manual p.1). Engine references as cited in the manual: Schafer & Yucel (2002) J. Computational and Graphical Statistics 11:437-457 (pan, p.28); Goldstein, Carpenter, Kenward & Levin (2009) Statistical Modelling 9:173-197 and Carpenter & Kenward (2013) Multiple Imputation and its Application, Wiley (jomo, p.16)

**Evidence (page-cited):**

- proper_mi -> panImpute Description, PDF page 25: 'Performs multiple imputation of multilevel data using the pan package (Schafer & Yucel, 2002)'; argument m, page 13 (jomoImpute) / page 25 (panImpute): 'The number of imputed data sets to generate. Default is 10.'
- proper_mi (downstream) -> mitmlComplete, PDF pages 22-23: 'This function extracts imputed data sets from mitml class objects as produced by panImpute and jomoImpute' ... Value: 'Usually a list of imputed data with class mitml.list'
- proper_mi (pooling) -> testEstimates, PDF pages 40-41: 'Computes final parameter estimates and inferences from multiply imputed data sets' ... 'calculates pooled parameter estimates and inferences as suggested by Rubin (1987, "Rubin's rules") for each parameter of the fitted model'
- multilevel -> package title, PDF page 1: 'Tools for Multiple Imputation in Multilevel Modeling'; jomoImpute Description, page 12: 'Performs single- and multilevel imputation... Supports imputation of missing data at level 1 and 2'; formula operator '|' 'denotes cluster-specific random effects and specifies the cluster indicator (e.g., 1|ID)', pages 14 and 26
- phylogeny -> whole manual: no phylogeny/tree term appears anywhere (function index pages 1-2 and all imputation-function pages checked); the only covariance structures described are cluster random effects and random level-1 residual covariance matrices (jomoImpute random.L1, pages 13-14)
- continuous -> panImpute Description, PDF page 25: 'Supports imputation of continuous multilevel data with missing values at level 1'; Details, page 26: 'This function serves as an interface to the pan package and supports imputation of continuous multilevel data at level 1 (Schafer & Yucel, 2002)'
- count -> panImpute (page 25) and jomoImpute (page 12) Descriptions enumerate supported types as 'continuous' and '(mixed) continuous and categorical' respectively; no count/Poisson model is mentioned anywhere in the manual
- binary -> jomoImpute Details, PDF page 13: 'supports imputation of single-level and multilevel continuous and categorical data at both level 1 and 2... In order for categorical variables to be detected correctly, these must be formatted as a factor variables' (binary never named as its own type)
- ordinal -> jomoImpute, PDF pages 12-14: only the word 'categorical' is used; no mention of ordered factors, proportional-odds/cumulative models, or ordinal probit anywhere in the manual
- nominal -> jomoImpute Description, PDF page 12: 'Performs single- and multilevel imputation for (mixed) continuous and categorical data using the jomo package'; page 13: categorical variables 'must be formatted as a factor variables'; example page 16: 'all variables (mixed continuous and categorical data with missing values at level 1 and 2)'
- approach/delegation -> mitml-package, PDF page 2: 'includes a user-friendly interface to the algorithms implemented in the R packages pan and jomo'; jomoImpute Details, page 13: 'This function serves as an interface to the jomo package'; page 14: 'Missing data are imputed jointly at both levels'; MCMC: keep.chains 'which chains of the MCMC algorithm to store' (pages 13, 26) with n.burn burn-in and n.iter between-imputation iterations (pages 13, 25)

**Not stated / unclear in the docs:**

- count: not stated — no count/Poisson-type imputation model is mentioned anywhere in the manual (type lists name only continuous and categorical)
- ordinal: not stated — the manual never distinguishes ordered from unordered categorical and names no ordinal model
- binary: not named separately — supported only insofar as it is a 'categorical' variable formatted as a factor (jomoImpute)
- phylogeny: not stated — the word/concept never appears in the manual
- latent model used for categorical imputation (e.g., latent-normal/probit) is not described in the mitml manual; it defers to jomo via the cited references (Goldstein et al. 2009; Carpenter & Kenward 2013)
- whether jomoImpute's single-level mode changes type support is not elaborated beyond 'single-level imputation will be performed' when the cluster indicator is omitted (page 14)

---

## hmi (v1.0.0, 2020-10-02 — final version; removed from CRAN 2022-05-10)

- approach: mice-like FCS; per-variable single-level (mice) or multilevel (MCMCglmm) models
- proper_mi: yes — m=5 imputations default; returns mice mids object
- multilevel: yes — random effects via lmer/glmer formula, cluster ID
- phylogeny: no — never mentioned; random effects via cluster ID only
- continuous: yes — dedicated continuous routines (mice "norm" / multilevel)
- count: yes — dedicated count-data routines; distribution not named (no "Poisson")
- binary: yes — dedicated binary routines, single and hierarchical
- ordinal: yes — proportional odds (polr); ordered type detected; multilevel variant
- nominal: yes — multilevel multinomial model; single-level mice regression trees
- reference (as stated in the docs): Matthias Speidel, Joerg Drechsler and Shahab Jolani (2020): "The R Package hmi: A Convenient Tool for Hierarchical Multiple Imputation and Beyond", Journal of Statistical Software, Vol. 95, No. 9, p. 1-48, doi:10.18637/jss.v095.i09 (as stated on PDF p.18). GitHub: https://github.com/matthiasspeidel/hmi (PDF p.1). CRAN: https://CRAN.R-project.org/package=hmi — archived 2022-05-10 ("requires archived package 'linLIR'"); manual PDF recovered from Wayback Machine snapshot (2022-01-20) of cran.r-project.org/web/packages/hmi/hmi.pdf

**Evidence (page-cited):**

- proper_mi -> hmi(), PDF p.16-18: Usage 'hmi(data, model_formula, ..., m = 5, ...)'; 'm An integer defining the number of imputations that should be made' (p.17); 'Value: The function returns a mids object. See mice for further information.' (p.18)
- proper_mi (pooling caveat) -> hmi_pool, PDF p.18-19: 'results from every dataset are pooled by simply averaging... Currently variance estimates for the averaged results are not implemented'; example (p.19): 'if you are interested in fixed effects only, consider pool from mice: mice::pool(with(data = dat_imputed, ...))'
- multilevel -> hmi, PDF p.16-17: title 'hmi: Hierarchical Multilevel Imputation'; model_formula 'designed to handle formula used in lm, glm, lmer and glmer'; package Description (p.1): 'Runs single level and multilevel imputation models'; example formula 'written ~ 1 + gender + coursework + (1 + gender|school)' (p.18)
- multilevel (random effects) -> imp_binary_multi, PDF p.22: 'Z_imp A data.frame with the random effects variables. clID A vector with the cluster ID. nitt An integer defining number of MCMC iterations (see MCMCglmm).'
- phylogeny -> whole manual: no occurrence of 'phylogen*', 'ginverse', 'pedigree' or 'animal model' anywhere (full extracted text searched); random-effect structure enterable only as Z_imp + clID cluster ID (e.g. imp_binary_multi p.22, imp_cont_multi p.25, imp_count_multi p.27)
- continuous -> get_type, PDF p.15: 'continuous (numeric values, or integers with more than 20 different values)'; imp_cont_single p.26: 'It uses mice with the method "norm"'; imp_cont_multi p.25: 'The function for hierarchical imputation of continuous variables'
- count -> get_type, PDF p.15: 'count data (integers)'; imp_count_multi p.27: 'The function for hierarchical imputation of variables with count data' (MCMCglmm nitt argument); imp_count_single p.28; the word 'Poisson' appears nowhere in the manual
- binary -> get_type, PDF p.15: 'binary (two different values - like 0s and 1s or "m" and "f")'; imp_binary_multi p.22: 'The function for hierarchical imputation of binary variables'; imp_binary_single p.23
- ordinal -> imp_orderedcat_single, PDF p.30: 'The function uses the proportional odds logistic regression (polr) approach, implemented in mice'; imp_orderedcat_multi p.29-30 (Z_imp, clID, nitt 'see MCMCglmm'); get_type p.15: 'ordered categorical (the categorical variable is ordered)'
- nominal -> imp_cat_multi, PDF p.23: 'relies on MCMCglmm... here we run a multilevel multinomial model'; imp_cat_single p.24: 'The function uses regression trees for imputation implemented in mice' (donor from same leaf); get_type p.15: 'categorical (the variable is a factor or has more than 3 different values)'
- approach -> package Description, PDF p.1: 'translates this analysis model into commands to impute the data according to it with functions from mice, MCMCglmm or routines build for this package'; imputationcycle p.20: 'cycles through every variable sequentially imputing the values, that are NA... determines the type of the variable and calls the suitable imputation function'
- other supported types (for the table footnote) -> get_type, PDF p.15: semicontinuous, rounded continuous, intercept; imp_interval p.29: interval data ('3000;4000', censored '-Inf;x')

**Not stated / unclear in the docs:**

- Count model distribution: not stated — manual says only 'count data (integers)' (p.15) and 'hierarchical imputation of variables with count data' (p.27); the word 'Poisson' never appears in the manual
- imp_count_single's doc title (p.28) mis-states 'The function for imputation of binary variables' — apparent copy-paste error in the package docs; cited the correctly-titled imp_count_multi instead
- Engine per multilevel routine: 'relies on MCMCglmm' is stated explicitly only for imp_cat_multi (p.23); the other *_multi routines reference MCMCglmm only via their nitt argument ('number of MCMC iterations (see MCMCglmm)')
- Exact statistical model for binary/count/ordinal multilevel routines (link functions, priors): not stated in the manual; the JSS paper is cited for method details (p.18)

---

## mdmb (Model Based Treatment of Missing Data), v1.9-22, Robitzsch & Luedtke

- approach: Factored (sequential) regression: ML/EM (frm_em) or fully Bayesian MCMC (frm_fb); substantive-model-compatible
- proper_mi: yes — frm_fb fully Bayesian; Nimp imputed datasets via frm2datlist, pooled downstream
- multilevel: yes — Bayesian frm_fb only (model="mlreg"); normal and ordinal-probit outcomes
- phylogeny: no — phylogeny/tree never mentioned anywhere in manual
- continuous: yes — linreg (normal); yjtreg (Yeo-Johnson scaled-t, incl. bounded (0,1)); bctreg (Box-Cox, positive)
- count: no — no Poisson/count model in exhaustive model-options list; "poisson"/"count" absent from entire manual
- binary: yes — "logistic" model: logistic regression for dichotomous variables/data
- ordinal: yes — "oprobit": ordinal probit regression for ordinal polytomous data (ordered, threshold/probit)
- nominal: no — manual states nominal-covariate handling "will be included in future mdmb package versions"
- reference (as stated in the docs): Package docs (DESCRIPTION, p.1) cite: Ibrahim, Chen, Lipsitz & Herring (2005) JASA 100:332-346; Luedtke, Robitzsch & West (2020a) Multivariate Behavioral Research 55(3):361-381; Luedtke, Robitzsch & West (2020b) Psychological Methods 25(2):157-181. Authors: Alexander Robitzsch [aut, cre], Oliver Luedtke [aut]. CRAN: https://cran.r-project.org/package=mdmb; GitHub: https://github.com/alexanderrobitzsch/mdmb

**Evidence (page-cited):**

- proper_mi -> frm (frm_fb), PDF page 8: "The function frm_fb employs a fully Bayesian approach with noninformative prior distribution. This function imputes missing values in the models from the posterior distributions. Imputed datasets can be extracted by the function frm2datlist."; page 9 usage: frm_fb(..., Nimp=10, ...) and "frm2datlist(object, as_mids=FALSE) # create list of imputed datasets"; page 10: Nimp = "Number of imputed datasets", as_mids = "whether multiply imputed datasets should be converted into objects of class mids"
- proper_mi (pooling workflow) -> frm Examples, PDF page 14: "create list of multiply imputed datasets ... datlist <- mdmb::frm2datlist(mod2) ... imp2 <- miceadds::datlist2mids(datlist) ... summary( mice::pool(mod2c) )"; also mdmb-package, page 2: "The function frm_fb allows substantive model compatible multiple imputation."
- multilevel -> frm Description, PDF page 8: "Multilevel models (using model=\"mlreg\") for normally distributed (outcome=\"normal\") and ordinal variables (outcome=\"probit\") as well as variables at higher levels (using argument variable_level) are accommodated. Note that multilevel models can only be fit with frm_fb and not with frm_em."; page 9: "mlreg" is a "wrapper to miceadds::ml_mcmc"; frm Example ~page 19-20 shows random-slope formula ( 1 + cwc(x, idcluster) | idcluster )
- phylogeny -> whole manual: full-text search (pdftotext) for "phylo" and "tree" returns zero hits across all 40 pages; no function in the Contents (page 2) relates to trees or phylogenetic covariance
- continuous -> frm dep$model options, PDF page 9: "Options are \"linreg\" (wrapper to stats::lm ...), ... \"yjtreg\" for continuous variables and bounded variables on (0,1) ..., \"bctreg\" for positive continuous variables"; mdmb_regression, page 24: "yjt_regression performs a linear regression in which the response variable is transformed according to the Yeo-Johnson transformation ... residuals are distributed following the scaled t distribution ... bct_regression has same functionality ... but employs a Box-Cox transformation"
- count -> frm dep$model options, PDF page 9: the exhaustive options list ("linreg", "logistic", "oprobit", "yjtreg", "bctreg", "mlreg") contains no count/Poisson model; full-text search for "poisson" and "count" returns zero hits in the entire manual
- binary -> mdmb_regression, PDF page 24: "The function logistic_regression performs logistic regression for dichotomous data."; frm, page 9: "\"logistic\" for dichotomous variables (wrapper to logistic_regression)"; page 8: supported covariates include "dichotomous covariates"
- ordinal -> mdmb_regression, PDF page 24: "The function oprobit_regression performs ordinal probit regression for ordinal polytomous data."; frm, page 9: "\"oprobit\" for ordinal variables (wrapper to oprobit_regression)"; page 8: "ordinal covariates" supported, and multilevel ordinal via outcome=\"probit\"
- nominal -> frm Description, PDF page 8: "The handling of nominal covariates will be included in future mdmb package versions." (only occurrence of "nominal" in the manual, confirmed by full-text search)
- scope of covariate types -> frm Description, PDF page 8: "The current functionality only support missing values on continuous covariates (accommodating skewness and only positive values), dichotomous covariates and ordinal covariates."

**Not stated / unclear in the docs:**

- phylogeny: strictly not stated — the manual never mentions phylogeny or trees (verified by full-text search); reported as 'no' because the capability is absent from the complete function reference, not because the docs explicitly deny it
- ML pathway (frm_em) produces model estimates, not imputed datasets — the manual ties imputed-dataset extraction (frm2datlist) only to frm_fb; whether frm_em output can yield imputations is not stated
- default Nimp=10 shown in usage (page 9); the manual does not state a recommended number of imputations

---

## JointAI (CRAN, v1.1.0)

- approach: Bayesian joint analysis+imputation via JAGS MCMC (automatic covariate sub-models)
- proper_mi: yes — integrated full-Bayes posterior; get_MIdat() also exports m imputed datasets
- multilevel: yes — nlme/lme4-style random effects, multiple grouping levels
- phylogeny: no — not mentioned; rd_vcov only full/blockdiag/indep
- continuous: yes — gaussian (lm default), lognorm, Gamma, beta models
- count: yes — poisson family; glm_poisson_log/identity and glmm_poisson_* model types
- binary: yes — glm_logit default for binary; binomial with logit/probit/log/cloglog
- ordinal: yes — clm/clmm cumulative logit, proportional odds by default (nonprop optional)
- nominal: yes — mlogit/mlogitmm multinomial logit, default for unordered factors
- reference (as stated in the docs): Erler NS, Rizopoulos D, Lesaffre EMEH (2021). "JointAI: Joint Analysis and Imputation of Incomplete Data in R." Journal of Statistical Software, 100(20), 1-56, doi:10.18637/jss.v100.i20 (cited on PDF p.16). CRAN: https://cran.r-project.org/package=JointAI; GitHub: https://github.com/nerler/JointAI

**Evidence (page-cited):**

- proper_mi (integrated Bayes) -> JointAI package overview, PDF page 14: "performs simultaneous imputation and inference for incomplete or complete data under the Bayesian framework"; covariate models are "specified automatically and modelled jointly with the analysis model. MCMC sampling is performed in 'JAGS' via the R package rjags"
- proper_mi (multiple imputed datasets) -> get_MIdat, PDF page 10: "Extract multiple imputed datasets from an object of class JointAI ... returns a dataset containing multiple imputed datasets stacked onto each other (i.e., long format)"; argument m = "number of imputed datasets" (default 10); minspace prevents "strong correlation between imputed datasets"
- proper_mi (caveat) -> get_MIdat Note, PDF page 11: "the imputed values must have been monitored, i.e., imps = TRUE had to be specified in the argument monitor_params"
- multilevel -> JointAI overview, PDF page 14: "lme_imp or lmer_imp for linear mixed models", "glme_imp or glmer_imp for generalized linear mixed models", "clmm_imp for (ordinal) cumulative logit mixed models"
- multilevel -> model_imp / Model formulas, PDF page 28: "It is possible to specify multi-level models as it is done in the package nlme ... or as it is done in the package lme4"; "Multiple levels of grouping: For multiple levels of grouping the specification using formula should be used" (nested and crossed random effects, pp.28-29)
- phylogeny -> model_imp / Random effects variance-covariance structure, PDF page 29: "Implemented structures are full, blockdiag and indep (all off-diagonal elements are zero)" — these govern correlation among random effects of different outcomes; no phylogenetic or user-supplied level covariance matrix is mentioned anywhere in the manual
- continuous -> model_imp / Distribution families, PDF page 30: "gaussian with links: identity, log; Gamma with links: inverse, identity, log"; model types table pages 30-31: lm "linear (normal) model with identity link ...; default for con[tinuous]", lognorm "log-normal model for skewed continuous data", beta "beta model (with logit link) for skewed continuous data in (0, 1)"
- count -> model_imp / Distribution families, PDF page 30: "poisson with links: log, identity"; model types, PDF page 31: glm_poisson_log "Poisson model with log link for count data", glm_poisson_identity "Poisson model with identity link for count data", glmm_poisson_log "Poisson mixed model with log link for count data"
- binary -> model_imp model types, PDF page 31: glm_logit "logistic model for binary data (alternatively: glm_binomial_logit); default for binary variables", glm_probit "probit model for binary data"; families page 30: "binomial with links: logit, probit, log, cloglog"
- ordinal -> JointAI overview, PDF page 14: "clm_imp for (ordinal) cumulative logit models", "clmm_imp for (ordinal) cumulative logit mixed models"; model types, PDF page 31: clm "cumulative logit model for ordered categorical variables; default for ordered factors"; PDF page 33: "By default, proportional odds are assumed and only the intercepts differ per category of the ordinal response" (nonprop argument allows non-proportional effects)
- nominal -> JointAI overview, PDF page 14: "mlogit_imp for multinomial logit models"; model types, PDF page 31: mlogit "multinomial logit model for unordered categorical variables; default for unordered factors with >2 [categories]", mlogitmm "multinomial logit mixed model for unordered categorical variables"; mlogitmm_imp usage shown on PDF page 25
- survival (context) -> JointAI overview, PDF page 14: "survreg_imp for parametric (Weibull) survival models", "coxph_imp for (Cox) proportional hazard models", "JM_imp for joint models of longitudinal and survival data"

**Not stated / unclear in the docs:**

- phylogeny: not stated — the manual never mentions phylogeny or supplying a known covariance matrix among grouping-factor levels; whether the experimental 'custom' JAGS sub-model argument (pages 27 and 33, "feature is experimental and has not been fully tested yet") could encode one is not stated
- downstream pooling of get_MIdat datasets (e.g. Rubin's rules) is not described in the manual; the package's own inference is via the joint MCMC posterior (summary/coef/confint, page 14)
- negative-binomial/overdispersed counts: not listed — the stated families are only gaussian, binomial, Gamma, poisson (page 30)

---

## smcfcs (R, CRAN v2.0.2, 2026-05-01; Bartlett, Keogh, Bonneville, van der Burg, Ekstrom)

- approach: Substantive-model-compatible FCS (Bartlett 2015): chained equations with rejection sampling, each covariate imputation model made compatible with a user-specified substantive model (smtype: lm, logistic, brlogistic, poisson, weibull, coxph, compet; plus design variants casecohort, dtsam discrete-time, finegray, flexsurv Royston-Parmar, nestedcc). Imputes covariates; outcomes auto-imputed only for lm/logistic/poisson substantive models.
- proper_mi: yes — returns m imputed datasets (default m=5); Rubin's-rules pooling via mitools
- multilevel: no — no random-effects/clustered option in documented smtype/method lists
- phylogeny: no — no phylogenetic or tree-based covariance mentioned anywhere
- continuous: yes — method "norm" (normal linear regression); also "latnorm" latent normal for error-prone variables
- count: yes — method "poisson" (Poisson regression), an explicit count imputation model
- binary: yes — methods "logreg" (logistic regression) and "brlogreg" (bias-reduced logistic regression)
- ordinal: yes — method "podds", proportional odds regression for ordered categorical variables
- nominal: yes — method "mlogit", multinomial logistic regression for unordered categorical variables
- reference (as stated in the docs): Bartlett JW, Seaman SR, White IR, Carpenter JR. Multiple imputation of covariates by fully conditional specification: accommodating the substantive model. Statistical Methods in Medical Research 2015; 24(4): 462-487. doi:10.1177/0962280214521348 (cited in manual, smcfcs References). CRAN: https://cran.r-project.org/package=smcfcs; GitHub: https://github.com/jwb133/smcfcs

**Evidence (page-cited):**

- proper_mi -> smcfcs Description + m argument, PDF pages 9-10: 'Multiply imputes missing covariate values using substantive model compatible fully conditional specification'; m = 'The number of imputed datasets to generate. The default is 5.'
- proper_mi (downstream pooling) -> smcfcs Details, PDF page 11: 'impDataset ... is a list of the imputed datasets. Models (e.g. the substantive model) can be fitted to each and results combined using Rubin's rules using the mitools package, as illustrated in the examples'
- multilevel -> smcfcs smtype and method arguments, PDF page 10: the exhaustive 'Possible values are' lists (smtype: lm, logistic, brlogistic, poisson, weibull, coxph, compet; method: norm, logreg, brlogreg, poisson, podds, mlogit, latnorm, passive expressions) contain no random-effects or clustered-data option; Contents (PDF page 2) lists no multilevel function
- phylogeny -> full function documentation, Contents PDF page 2 and function pages 9-24: no phylogenetic, tree, or comparative-data term appears anywhere in the manual
- continuous -> smcfcs method argument, PDF page 10: '"norm" (normal linear regression)'; also '"latnorm" indicates the variable is a latent normal variable which is measured with error'
- count -> smcfcs method argument, PDF page 10: '"poisson" (Poisson regression)'
- binary -> smcfcs method argument, PDF page 10: '"logreg" (logistic regression), "brlogreg" (bias reduced logistic regression)'
- ordinal -> smcfcs method argument, PDF page 10: '"podds" (proportional odds regression for ordered categorical variables)'
- nominal -> smcfcs method argument, PDF page 10: '"mlogit" (multinomial logistic regression for unordered categorical variables)'
- approach/engine -> package DESCRIPTION, PDF page 1: 'a modification of the popular FCS/chained equations multiple imputation approach, and allows imputation of missing covariate values from models which are compatible with the user specified substantive model'; rejection sampling per rjlimit argument, PDF page 10: 'attempts which should be made when using rejection sampling to draw from imputation models'
- substantive models -> smcfcs Details, PDF page 11: 'Imputation is supported for linear regression ("lm"), logistic regression ("logistic"), bias reduced logistic regression ("brlogistic"), Poisson regression ("poisson"), Weibull ("weibull") and Cox regression for time to event data ("coxph"), and Cox models for competing risks data ("compet")'; design variants: smcfcs.casecohort p14, smcfcs.dtsam (discrete-time survival) p15, smcfcs.finegray (Fine-Gray subdistribution hazard) p17, smcfcs.flexsurv (Royston-Parmar flexible parametric) p19, smcfcs.nestedcc (nested case-control) p21, smcfcs.parallel (parallel wrapper) p23
- scope of imputation -> smcfcs Details, PDF page 11: 'If the substantive model is linear, logistic or Poisson regression, smcfcs will automatically impute missing outcomes, if present, using the specified substantive model' (otherwise covariates only)
- reference -> smcfcs References, PDF page 12: 'Bartlett JW, Seaman SR, White IR, Carpenter JR. Multiple imputation of covariates by fully conditional specification: accommodating the substantive model. Statistical Methods in Medical Research 2015; 24(4): 462-487. doi:10.1177/0962280214521348'

**Not stated / unclear in the docs:**

- multilevel: the manual never explicitly discusses clustered/multilevel data; the 'no' rests on the documented exhaustive option lists (p10) containing no such option — no 'single-level only' sentence exists to quote
- phylogeny: never mentioned anywhere; 'no' is by omission from the exhaustive option lists and the full function index (p2)
- overdispersed counts (negative binomial or similar): not stated — only Poisson is documented
- whether podds/mlogit imputation assumes anything beyond standard proportional-odds/multinomial logistic (e.g. probit variants): not stated
- pages 4-8 were not read in full; the Contents (p2) shows they document only example datasets (ex_coxquad ... ex_poisson) and the start of plot.smcfcs, so all capability claims come from the function pages 9-24, page 1, and page 2

---

## Rphylopars (v0.3.10, CRAN manual dated 2024-01-17)

- approach: REML/ML joint multivariate Gaussian trait-evolution model (BM/mvOU/OU/EB/lambda/kappa/delta/star) with within-species variance; point prediction of missing species means + prediction variance (not FCS, not Bayesian)
- proper_mi: no — single predicted values with prediction variances; multiple imputed datasets never mentioned
- multilevel: partial — within-species replicates (intraspecific variation/measurement error) only, no general random effects
- phylogeny: yes — imputation performed on fitted phylogenetic evolutionary model
- continuous: yes — Gaussian trait-evolution models (BM/OU/EB/lambda/kappa/delta), variance-covariance framework
- count: no — not stated; no count/Poisson method anywhere in manual
- binary: no — not stated; binary/threshold support never mentioned in manual
- ordinal: no — not stated; ordinal/ordered-categorical never mentioned in manual
- nominal: no — not stated; nominal/multinomial/factor traits never mentioned in manual
- reference (as stated in the docs): As stated in manual (Rphylopars-package References, PDF p. 2): Goolsby EW, Ane C, Bruggeman J. 2017. "Rphylopars: Fast Multivariate Phylogenetic Comparative Methods for Missing Data and Within-Species Variation." Methods in Ecology & Evolution 8:22-27; original method: Bruggeman J, Heringa J, Brandt BW (2009) "PhyloPars: estimation of missing parameter values using phylogeny", Nucleic Acids Research 37:W179-W184. CRAN: https://cran.r-project.org/package=Rphylopars; GitHub: https://github.com/ericgoolsby/Rphylopars/wiki (URL field, PDF p. 1)

**Evidence (page-cited):**

- proper_mi -> phylopars Description, PDF page 7: 'Reconstructed ancestral states and predicted species means (i.e., for missing data), along with prediction variances, are also provided' — one completed prediction per missing value plus variance; no multiple-dataset or posterior-draw functionality documented anywhere in the manual (pp. 1-18)
- proper_mi (supporting) -> phylopars Value, PDF page 9: 'anc_recon: Reconstructed ancestral states and species means'; 'anc_var: Variance of reconstructed ancestral estimates and imputed species means'; also Examples, PDF page 3: 'PPE$anc_recon # Ancestral state reconstruction and species mean prediction; PPE$anc_var # Prediction variance'
- multilevel -> Rphylopars-package Description, PDF page 2: 'datasets with multiple observations per species (intraspecific variation or measurement error)'; phylopars trait_data arg, page 7: 'multiple observations for species are allowed'; pheno_error arg, page 7: 'parameters are estimated assuming within-species variation'; phenocov_list, page 9: 'species-specific within-species covariance matrices to be held fixed, as in Ives et al (2007)'
- phylogeny -> Rphylopars-package Description, PDF page 2: 'Performs ancestral state reconstruction and missing data imputation on the estimated evolutionary model, which can be specified as Brownian Motion, Ornstein-Uhlenbeck, Early-Burst, Pagel's lambda, kappa, or delta, or a star phylogeny'; Index keyword 'Phylogenetic imputation', page 18; tree arg 'An object of class phylo', page 7
- continuous -> phylopars, PDF page 7: 'estimates parameters for the phylogenetic and phenotypic variance-covariance matrices'; model arg: 'Default is BM ... "mvOU" (for the multivariate Ornstein-Uhlenbeck), or univariate tree transformations: "OU" "lambda", "kappa", "delta", "EB", "star"' — all Gaussian trait-evolution models; numeric example data x <- setNames(rnorm(1e4),...) (anc.recon Examples, page 4)
- count -> whole manual read (pp. 1-18): no count/Poisson-type method, argument, or function documented; model list (phylopars, page 7) contains only continuous trait-evolution models
- binary -> whole manual read (pp. 1-18): binary/threshold/probit support never mentioned; Index (page 18) lists no categorical functionality
- ordinal -> whole manual read (pp. 1-18): ordinal/ordered-categorical handling never mentioned in any function page
- nominal -> whole manual read (pp. 1-18): nominal/multinomial/factor traits never mentioned; trait_data described only as 'one column per trait ... Missing data should be represented with NA' (phylopars, page 7)
- reference -> Rphylopars-package References, PDF page 2: 'Goolsby EW, Ane C, Bruggeman J. 2017. Rphylopars: Fast Multivariate Phylogenetic Comparative Methods for Missing Data and Within-Species Variation. Methods in Ecology & Evolution. 2017. 8:22-27.' and 'Bruggeman J, Heringa J and Brandt BW. (2009) PhyloPars: estimation of missing parameter values using phylogeny. Nucleic Acids Research 37: W179-W184.'
- approach/engine -> phylopars REML arg, PDF page 7: 'If TRUE (default), the algorithm will return REML estimates. If FALSE, maximum likelihood estimates will be returned'; skip_EM/skip_optim args, page 8: EM then BFGS optimization

**Not stated / unclear in the docs:**

- Trait-type restriction is never stated explicitly: the word 'continuous' does not appear; continuous-only support is conveyed by the Gaussian trait-evolution model list (BM/OU/EB/lambda/kappa/delta, p. 7) and numeric examples, not by an explicit statement
- count/binary/ordinal/nominal: manual contains no statement either way (no support documented, no explicit exclusion) — 'no' rests on absence from the complete 18-page manual, i.e. not stated as supported
- Whether users could manually resample from anc_var to construct multiple imputations is not addressed in the docs
- Whether prediction variances account for parameter-estimation uncertainty (vs conditional-on-parameters variance only) is not stated in the manual

---

## missForest (CRAN v1.6.1, manual dated 2026-05-08)

- approach: ranger/randomForest engine; iterative per-variable (chained) random-forest single imputation
- proper_mi: no - single imputed dataset; OOB error is accuracy, not uncertainty
- multilevel: no - random effects/clustering never mentioned in 12-page manual
- phylogeny: no - phylogeny never mentioned anywhere in manual
- continuous: yes - numeric variables via regression forests; NRMSE error reported
- count: no - no count/Poisson type; only numeric vs factor documented
- binary: yes - only as generic factor/classification; 'binary' never named
- ordinal: no - ordered factors never mentioned; no as-ordered handling stated
- nominal: yes - factors via classification forests, PFC error; 3-class iris example
- reference (as stated in the docs): Stekhoven DJ, Buehlmann P (2012). "MissForest — nonparametric missing value imputation for mixed-type data." Bioinformatics, 28(1), 112-118. doi:10.1093/bioinformatics/btr597 (cited on PDF pp.2, 5, 7). CRAN: https://cran.r-project.org/package=missForest ; GitHub: https://github.com/stekhoven/missForest (URL field, PDF p.1)

**Evidence (page-cited):**

- approach -> missForest() Details/Algorithm, PDF p.4: 'The method iteratively imputes each variable with missing values by fitting a random forest on the observed part of that variable and the current imputations of all other variables'; p.3: 'offers two backends: ranger (default) and randomForest (legacy/compatibility)'
- proper_mi -> missForest() Value, PDF p.5: 'ximp — Imputed data matrix (same classes as xmis)' — one completed dataset; no multiple-imputation, pooling, or uncertainty-propagation facility appears on any of the 12 pages (Index p.12 has no such entries)
- proper_mi (OOB is not MI uncertainty) -> DESCRIPTION, PDF p.1: 'It yields an out-of-bag (OOB) imputation error estimate without the need of a test set or elaborate cross-validation'; Value p.5: 'OOBerror — Estimated OOB imputation error. For numeric variables, the normalized root mean squared error (NRMSE); for factors, the proportion falsely classified (PFC)' — an error estimate of the imputation, not between-imputation variance
- multilevel -> entire manual pp.1-12: random effects/clustering/hierarchical data never mentioned; closest argument is strata, PDF p.4: 'List of (factor) variables used for stratified sampling (legacy randomForest semantics). Ignored by ranger' — tree-level sampling, not a multilevel model
- phylogeny -> entire manual pp.1-12: no mention of phylogeny, trees-as-covariance, or comparative data in any function entry or the Index (p.12)
- continuous -> DESCRIPTION, PDF p.1: 'It can be used to impute continuous and/or categorical data'; missForest() Value, PDF p.5: 'For numeric variables, the normalized root mean squared error (NRMSE)'
- count -> varClass() Value, PDF p.10: 'Entries are "numeric" for continuous variables and "factor" for categorical variables' — the documented type system has no count/Poisson type; counts/integers are never mentioned in the manual
- binary -> missForest() Description, PDF p.3: 'imputes missing values for mixed-type data (numeric and categorical)'; classwt arg p.3: 'List of class priors for the categorical variables' — binary is covered only as a generic factor; no binary-specific method is stated
- ordinal -> xmis argument, PDF p.3: 'All columns must be numeric or factor (character columns should be converted to factors beforehand)'; varClass p.10 lists only numeric/factor — ordered factors are never mentioned anywhere, so no proportional-odds/cumulative/ordered treatment is documented
- nominal -> missForest-package, PDF p.2: 'nonparametric missing-value imputation for mixed-type data (continuous and categorical)... models each variable with missingness using random forests'; Value p.5: 'for factors, the proportion falsely classified (PFC)'; cutoff arg p.4: 'List of per-class cutoff vectors for each categorical variable'; examples impute iris incl. 3-level Species (pp.5, 7; p.11: 'four continuous and one categorical variable')

**Not stated / unclear in the docs:**

- Behavior when an ordered factor is supplied: not stated (ordered factors never mentioned; type system documented as numeric vs factor only)
- How integer/count columns are classed or modelled: not stated (no count type documented; manual never discusses counts)
- The word 'binary' never appears; 2-level support is implied only by the generic factor/categorical handling, not stated separately
- Missingness-mechanism assumptions (MAR/MCAR) of the imputation method itself: not stated (MCAR appears only in prodNA, the benchmarking utility that generates missingness)

---

## phylomice 0.0.0.9000 (GitHub-only: pdrhlik/phylomice; not on CRAN, last commit 2016-07-09)

- approach: mice FCS add-on: univariate phylogenetic-GLS imputation (phnorm) using tree-derived covariance matrix psi
- proper_mi: yes — runs inside mice(); docs state "multiple imputation"
- multilevel: no — no random-effects/clustered method among documented methods
- phylogeny: yes — tree-derived covariance matrix psi in GLS imputation
- continuous: yes — phnorm: univariate continuous via phylogenetic GLS
- count: no — no count/Poisson method documented
- binary: no — no binary/logistic/threshold method documented
- ordinal: no — no ordinal/cumulative method documented
- nominal: no — no multinomial/polytomous method documented
- reference (as stated in the docs): No package publication stated in docs. Authors: Patrik Drhlik (aut/cre) and Simon P. Blomberg ("Theoretical background"), 2016. Method docs cite Garland & Ives (2000) The American Naturalist 155(3):346-364 (GLS basis) and van Buuren & Groothuis-Oudshoorn (2011) Journal of Statistical Software 45(3) (mice). GitHub: https://github.com/pdrhlik/phylomice (not on CRAN; install via devtools::install_github)

**Evidence (page-cited):**

- proper_mi -> DESCRIPTION file (no PDF manual exists; GitHub text docs), Title/Description: "Multiple imputation with phylogenetic information." / "methods that allow doing multiple imputation with phylogenetic information"; README 'phnorm' section shows methods run inside mice(): "imp <- mice(data, method = 'phnorm', psi = psi, psiinv = psiinv)"
- multilevel -> README 'Available methods' ("This is a list of methods that are available in the package": phnorm, phpmm) + phylomice.Rd \details ("Currently contained methods: phnorm; phpmm - under development"): no method mentions random effects or clustered data
- phylogeny -> README 'Available methods': "psi is a covariance matrix created from the phylogenetic tree"; mice.impute.phnorm.Rd \note: "based on the mice.impute.norm method from mice package and the GLM approach by Garland and Ives (2000)"; precomputePsi.Rd: "Precomputes psi and psiinv given a tree"
- continuous -> mice.impute.phnorm.Rd \description: "Imputes univariate continuous missing data using the generalized least square approach." (mice.impute.phpp.Rd carries the identical description)
- count -> documented absence: full man/ listing (6 Rd files: phnorm, phpmm, phpp, phppm, phylomice, precomputePsi) and README method list contain no count/Poisson method
- binary -> documented absence: same complete method list contains no binary/logistic/threshold method
- ordinal -> documented absence: same complete method list contains no ordinal/cumulative/proportional-odds method
- nominal -> documented absence: same complete method list contains no multinomial/unordered-categorical method
- maintained/CRAN -> DESCRIPTION: Version 0.0.0.9000 (development-only version number); README installation only via devtools::install_github; GitHub commits API: latest commit 2016-07-09 ("finish phpp method and add phppm method"); CRAN package page returns HTTP 404 and crandb reports not_found

**Not stated / unclear in the docs:**

- phpmm: algorithm and variable types not stated (usage-only stub Rd; README and phylomice.Rd say 'under development')
- phppm: purpose and variable type not stated (stub Rd gives only the usage line and 'nn Number of nearest neighbors')
- phpp: how it differs from phnorm not stated (its Rd title/description duplicate phnorm's; phpp absent from README and package overview Rd)
- Package-level citation/publication: not stated in any doc
- phnorm draw mechanics (whether parameters are drawn from a posterior per imputation): not stated beyond 'based on mice.impute.norm'
- Number of imputations/pooling: delegated to mice, not described in phylomice docs

---

## BHPMF 1.0 (CRAN; archived 2017-11-14, unmaintained)

- approach: Bayesian hierarchical probabilistic matrix factorization (Gibbs sampling) with taxonomic side information
- proper_mi: no — single mean-imputed matrix plus per-cell SD; no multiple datasets or posterior draws exposed
- multilevel: yes — nested taxonomic hierarchy levels are the model core
- phylogeny: partial — taxonomic hierarchy only (species/genus/family/"phylogenetic groups"); no tree or branch lengths
- continuous: yes — numeric data matrix; documented example traits log-transformed and z-scored
- count: not-stated — no count/Poisson-type method mentioned anywhere in the docs
- binary: not-stated — binary variables never mentioned in the docs
- ordinal: not-stated — ordinal/ordered-categorical handling never mentioned in the docs
- nominal: not-stated — nominal imputation never mentioned; taxonomy enters only as side information, not as an imputed variable
- reference (as stated in the docs): Schrodt F., Kattge J., Shan H., Fazayeli F., Karpatne A., Banerjee A., Reichstein M., Boenisch M., Diaz S., Dickie J., Gillison A., Kumar V., Lavorel S., Leadley P.W., Wirth C., Wright I., Wright S.J., Reich P.B. (2015) "BHPMF - a hierarchical Bayesian approach to gap-filling and trait prediction for macroecology and functional biogeography", Global Ecology and Biogeography 24: 1510-1521 (cited in every man page); methods paper: Fazayeli F., Banerjee A., Kattge J., Schrodt F., Reich P.B. (2014) "Uncertainty quantified matrix completion using Bayesian Hierarchical Matrix factorization", ICMLA. CRAN archive: https://cran.r-project.org/src/contrib/Archive/BHPMF/

**Evidence (page-cited):**

- archived status -> CRAN index page (cran_index.html, no PDF manual exists): 'Package BHPMF was removed from the CRAN repository. ... Archived on 2017-11-14 as check problems were not corrected despite reminders.'
- approach -> BHPMF-package.Rd \description (no PDF; Rd file p.1 equivalent): 'filling gaps in a data matrix with hierarchical information (e.g. taxonomy) using Gibbs sampling (BHPMF)'; DESCRIPTION Title: 'Uncertainty Quantified Matrix Completion using Bayesian Hierarchical Matrix Factorization'
- proper_mi -> GapFilling, man/gap_filling.Rd \description: 'Will save the mean and standard deviations (uncertainties) of the predictions for both missing and observed values into the given file paths'; BHPMF-package.Rd \details: 'Output is a gap filled matrix (mean of the predicted distribution) ... along with another matrix containing the standard deviation ... The standard deviation will be used for uncertainity quantification [sic]' — one completed dataset + SD matrix; multiple imputations or draw retrieval appear nowhere
- multilevel -> GapFilling, man/gap_filling.Rd \arguments: 'used.num.hierarchy.levels: Number of hierarchy levels used for gap filling'; DESCRIPTION: 'Fills the gaps of a matrix incorporating a hierarchical side information while providing uncertainty quantification'
- phylogeny -> man/hierarchy.info.Rd \description: 'individual plants are categorized to species, species to genera, genera to families, and families to phylogenetic groups ... the hierarchy as a taxa table'; gap_filling.Rd \arguments: hierarchy.info is 'A matrix containing the hierarchical information' — a taxa table; no phylogenetic tree, phylo object, or branch lengths mentioned in any doc file
- continuous -> man/trait.info.Rd \description: 'Plant traits have log-normal distribution, so we first took the logarithm for each entries in the plant x trait matrix, and then calculate the z-score for each trait'; gap_filling.Rd \arguments X: 'The data matrix (including missing values) ... Missing values must be filled by NA'
- count -> absence across the complete docs (DESCRIPTION + all 6 man/*.Rd files read): no occurrence of count/Poisson or any non-Gaussian data-type method
- binary -> absence across the complete docs: binary/dichotomous variables never mentioned in any man page or DESCRIPTION
- ordinal -> absence across the complete docs: ordinal/ordered-factor handling never mentioned in any man page or DESCRIPTION
- nominal -> absence across the complete docs: no multinomial/categorical imputation mentioned; the only categorical content is the taxonomy predictor table (man/hierarchy.info.Rd), which is side information, not an imputed variable

**Not stated / unclear in the docs:**

- count imputation: not stated (no mention anywhere in DESCRIPTION or the 6 man pages)
- binary imputation: not stated
- ordinal imputation: not stated
- nominal imputation: not stated (taxonomy appears only as hierarchical side information)
- whether input X is restricted to continuous/Gaussian data: not stated explicitly — docs describe only a numeric 'data matrix' and a log-transformed, z-scored continuous trait example
- likelihood/noise distribution of the factorization model: not stated in package docs (only 'mean of the predicted distribution' and 'standard deviation')
- access to full posterior draws (beyond the mean and SD output files): not stated; GapFilling's documented return value covers only tuning outputs
- MAR/MNAR missingness assumptions: not stated anywhere in the docs

---
## TDIP 1.0 (GitHub-only: Matgend/TDIP; not on CRAN; DESCRIPTION dated 2022-10-25, last commits 2024-10-11)

- approach: wrapper/benchmarking toolkit: phylogenetic imputation (Rphylopars
  for continuous; corHMM Markov model for discrete), machine learning
  (missForest, kNN, mice-PMM), deep learning (GAIN via Python/tensorflow);
  each ML/DL method under 3 strategies (no phylogeny / phylogenetic
  eigenvectors / "2-step" with a phylogenetically imputed "hint"); plus a
  hard_voting() majority-vote ensemble
- proper_mi: no - every imputation function's documented return is a single
  completed data frame; mice_phylo takes nbrMI but still documents one data
  frame out; no pooling/Rubin/between-imputation variance anywhere in docs;
  hard_voting() collapses several imputed frames into one consensus frame;
  evaluation is point accuracy (RMSE/PFC)
- multilevel: not stated - no random effects / within-species replicates /
  clustering anywhere in the docs
- phylogeny: yes - three routes: direct (Rphylopars continuous; corHMM
  discrete), tree eigenvectors as predictors (get_eigenvec,
  variance_fraction), and 2-step "hint" concatenation
- continuous: yes (pi_continuous_traits "applying the phylopars approach";
  mice-PMM; missForest)
- count: not stated - no count/Poisson model anywhere in the docs
- binary: yes (mice_phylo: "if binary the algorithm uses a logistic
  regression")
- ordinal: not stated as ordered - "ordinal" appears only as a SIMULATED
  trait class (data_simulator); mice_phylo documents "polytomous regression
  for unordered categorical traits"; corHMM imputation described generically
  for "discrete traits"
- nominal: yes (polytomous regression >2 states; corHMM for factors)
- reference (as stated in the docs): none - no paper cited anywhere in the
  package docs; only a Zenodo DOI badge (10.5281/zenodo.10890000). The ms
  cites @gendre2024 for TDIP on our own authority (the package's authors'
  benchmark paper), not on the package docs'.

**Evidence (page-cited):** README.md "Missing value imputation" section
(method list and the 3 strategies); man/mice_phylo.Rd (nbrMI argument;
single-data-frame Value; PMM continuous / logistic binary / polytomous
unordered); man/pi_continuous_traits.Rd, man/pi_categorical_traits.Rd,
man/imputeOneDiscreteTrait.Rd (corHMM); man/get_eigenvec.Rd and
man/gain_phylo.Rd (eigenvectors / variance_fraction); man/hard_voting.Rd
("selecting the most recurring value for each row"); man/missForest_phylo.Rd
(hint = "dataframe already imputed by a comparative methods"); README
"Error calculation" (RMSE/PFC); DESCRIPTION (Version 1.0, Date 2022-10-25);
GitHub commits API (last commits 2024-10-11); CRAN package page HTTP 404;
README installation via devtools::install_github + NADIA archived-from-CRAN
note.

**Not stated / unclear in the docs:** paper reference; whether mice_phylo can
return multiple datasets internally; ordered handling of ordinal factors at
imputation time; count variables; multilevel; which backend kNN uses;
README/vignette vs man/ function-name inconsistencies
(pi_continuous_traits_phylo vs pi_continuous_traits).

---

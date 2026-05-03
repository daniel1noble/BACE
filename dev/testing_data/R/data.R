#' phyloTraitData: Tidied phylogenetic trait datasets for imputation
#'   benchmarking
#'
#' A bundle of real-world species-by-trait datasets paired with
#' phylogenies, ready for direct use in phylogenetic imputation
#' benchmarks for the sister packages \pkg{pigauto} and \pkg{BACE}.
#'
#' Each newly tidied dataset exposes two objects: \code{<name>_traits}
#' (a \code{data.frame} keyed by \code{rownames} matching the tree
#' tips) and \code{<name>_tree} (an \pkg{ape} \code{phylo}).  Where
#' the source ships higher taxonomy (Order / Family / Genus), those
#' columns are appended at the end of the trait \code{data.frame}
#' under \code{order_name} / \code{family_name} / \code{genus_name}
#' so they can be used as taxonomic random effects in BACE workflows.
#'
#' Five additional bundled datasets shipped from \pkg{pigauto}'s own
#' \code{data/} directory are passed through verbatim
#' (\code{delhey5809}, \code{tree_delhey}, \code{ctmax_sim},
#' \code{tree300}, \code{trees300}) so existing code referring to
#' them keeps working after \code{library(phyloTraitData)}.
#'
#' AVONET is also bundled here as \code{avonet_traits} +
#' \code{avonet_tree}, built from the original AVONET 1.0 release
#' (Tobias et al. 2022) intersected with the Hackett et al. 2008
#' bird tree.  See \code{data-raw/make_avonet.R} for the construction
#' script.
#'
#' @name phyloTraitData-package
#' @aliases phyloTraitData
#' @importFrom ape Ntip
"_PACKAGE"


# ---------------------------------------------------------------------------
# 1. PanTHERIA mammals --------------------------------------------------------
# ---------------------------------------------------------------------------

#' PanTHERIA mammal traits, tidied for imputation benchmarks
#'
#' A wide species-by-trait matrix derived from PanTHERIA 1.0 (Jones
#' et al. 2009 \emph{Ecology} 90:2648), pre-aligned to the
#' \code{\link{pantheria_tree}} tip set.  Eight mixed-type trait
#' columns plus three higher-taxonomy columns.  Continuous values are
#' on the original (raw) scale; users typically log-transform body
#' mass / length / gestation / longevity at use-site.
#'
#' @format A \code{data.frame} with 4,629 rows (one row per species)
#'   and 11 columns.  Row names are species keys ("Genus_species")
#'   matching \code{pantheria_tree$tip.label}.
#' \describe{
#'   \item{body_mass_g}{Adult body mass in g (numeric, raw scale)}
#'   \item{head_body_length_mm}{Head-body length in mm (numeric)}
#'   \item{gestation_d}{Gestation length in days (numeric)}
#'   \item{max_longevity_m}{Max longevity in months (numeric)}
#'   \item{litter_size}{Litter size, integer count}
#'   \item{diet_breadth}{Diet breadth, ordered factor with K=5 levels (1..5)}
#'   \item{habitat_breadth}{Habitat breadth, ordered factor with K=3 levels}
#'   \item{terrestriality}{Terrestrial vs ground-dwelling, factor (1, 2)}
#'   \item{order_name}{MSW93 mammal order, factor}
#'   \item{family_name}{MSW93 mammal family, factor}
#'   \item{genus_name}{MSW93 mammal genus, factor}
#' }
#' @source Jones et al. 2009 \emph{Ecology} 90:2648.
#' @seealso \code{\link{pantheria_tree}}
#' @examples
#' data(pantheria_traits)
#' data(pantheria_tree)
#' stopifnot(setequal(rownames(pantheria_traits), pantheria_tree$tip.label))
"pantheria_traits"

#' Mammal phylogeny aligned to PanTHERIA 1.0
#'
#' A pruned mammal supertree on the PanTHERIA 1.0 species set.  Tip
#' labels are "Genus_species" matching \code{rownames(pantheria_traits)}.
#'
#' @format An \code{ape::phylo} object with 4,629 tips and resolved
#'   branch lengths.
#' @source Bininda-Emonds et al. 2007 \emph{Nature} 446:507 mammal
#'   supertree, intersected with PanTHERIA 1.0 species coverage.
#' @seealso \code{\link{pantheria_traits}}
"pantheria_tree"


# ---------------------------------------------------------------------------
# 2. AmphiBIO amphibians ------------------------------------------------------
# ---------------------------------------------------------------------------

#' AmphiBIO amphibian traits, tidied for imputation benchmarks
#'
#' A wide species-by-trait matrix derived from AmphiBIO v1 (Oliveira
#' et al. 2017 \emph{Sci. Data} 4:170123), pre-aligned to the
#' \code{\link{amphibio_tree}} tip set.  Six mixed-type trait columns
#' plus three higher-taxonomy columns.
#'
#' Note: \code{diurnal} and \code{nocturnal} are encoded
#' presence-only in AmphiBIO (1 = recorded; NA = no record).  Treat
#' NA as truly missing rather than as 0; otherwise class imbalance
#' produces artefactual results.
#'
#' @format A \code{data.frame} with 6,776 rows and 9 columns.  Row
#'   names are species keys ("Genus_species") matching
#'   \code{amphibio_tree$tip.label}.
#' \describe{
#'   \item{body_size_mm}{SVL or total length in mm (numeric)}
#'   \item{body_mass_g}{Body mass in g (numeric)}
#'   \item{diurnal}{Recorded as diurnal: factor "no" / "yes" (presence-only)}
#'   \item{nocturnal}{Recorded as nocturnal: factor "no" / "yes" (presence-only)}
#'   \item{diet_breadth}{Diet breadth (1..5), ordered factor}
#'   \item{habitat}{Habitat preference Fos/Ter/Aqu/Arb, factor (K=4)}
#'   \item{order_name}{Amphibian order, factor}
#'   \item{family_name}{Amphibian family, factor}
#'   \item{genus_name}{Amphibian genus, factor}
#' }
#' @source Oliveira et al. 2017 \emph{Sci. Data} 4:170123.
#' @seealso \code{\link{amphibio_tree}}
"amphibio_traits"

#' Amphibian phylogeny aligned to AmphiBIO v1
#'
#' A taxonomic tree built from AmphiBIO's Order / Family / Genus /
#' Species columns with Grafen rank-based branch lengths and pruned
#' to the species in \code{amphibio_traits}.  This is a taxonomic
#' tree (not a molecular phylogeny); for analyses that need divergence
#' times, use a published amphibian supertree instead.
#'
#' @format An \code{ape::phylo} object with 6,776 tips.
#' @source Oliveira et al. 2017 (taxonomy); Grafen 1989 \emph{Phil.
#'   Trans. R. Soc. B} 326:119 (branch length method).
#' @seealso \code{\link{amphibio_traits}}
"amphibio_tree"


# ---------------------------------------------------------------------------
# 3. BIEN plants --------------------------------------------------------------
# ---------------------------------------------------------------------------

#' BIEN plant trait means, tidied for imputation benchmarks
#'
#' Species-level mean trait values from the Botanical Information and
#' Ecology Network (BIEN) v4 database, intersected with the
#' \code{V.PhyloMaker2} backbone tree (Jin & Qian 2022).  Five
#' continuous trait columns; values are on the raw scale and
#' typically log-transformed at use-site.
#'
#' @format A \code{data.frame} with 19,109 rows and 5 columns.  Row
#'   names are species keys ("Genus_species") matching
#'   \code{bien_tree$tip.label}.
#' \describe{
#'   \item{height_m}{Maximum whole-plant height (m), continuous}
#'   \item{leaf_area}{Leaf area (mm^2), continuous}
#'   \item{sla}{Specific leaf area, leaf area / leaf dry mass (mm^2/mg)}
#'   \item{seed_mass}{Seed mass (mg), continuous}
#'   \item{wood_density}{Stem wood density (g/cm^3), continuous}
#' }
#' @source BIEN v4 \url{https://bien.nceas.ucsb.edu/bien/}; trait data
#'   curated and exposed via the \pkg{BIEN} R package.
#' @seealso \code{\link{bien_tree}}
"bien_traits"

#' Plant phylogeny aligned to BIEN trait means
#'
#' A megaphylogeny-backed tree of seed plants built with
#' \code{V.PhyloMaker2::phylo.maker(scenario = "S3")} on the
#' Smith & Brown 2018 backbone, pruned to the species in
#' \code{bien_traits}.
#'
#' @format An \code{ape::phylo} object with 19,109 tips.
#' @source Jin & Qian 2022 \emph{J. Plant Ecology}; Smith & Brown
#'   2018 \emph{Am. J. Bot.} 105:302.
#' @seealso \code{\link{bien_traits}}
"bien_tree"


# ---------------------------------------------------------------------------
# 4. GlobTherm ectotherms -----------------------------------------------------
# ---------------------------------------------------------------------------

#' GlobTherm ectotherm thermal limits, tidied for imputation benchmarks
#'
#' A multi-class ectotherm dataset of upper / lower thermal limits
#' (Tmax / Tmin) plus geographic covariates and full higher taxonomy
#' for use as taxonomic random effects.
#'
#' @format A \code{data.frame} with 1,969 rows and 9 columns.  Row
#'   names are species keys ("Genus_species") matching
#'   \code{globtherm_tree$tip.label}.
#' \describe{
#'   \item{Tmax}{Upper thermal limit (degC), continuous}
#'   \item{Tmin}{Lower thermal limit (degC), continuous}
#'   \item{lat_max}{Latitude of Tmax measurement (degrees), numeric}
#'   \item{long_max}{Longitude of Tmax measurement, numeric}
#'   \item{elevation_max}{Elevation of Tmax measurement (m), numeric}
#'   \item{class_name}{Class (Insecta / Reptilia / Actinopteri / ...), factor}
#'   \item{order_name}{Order, factor}
#'   \item{family_name}{Family, factor}
#'   \item{genus_name}{Genus, factor}
#' }
#' @source Bennett et al. 2018 \emph{Sci. Data} 5:180022 (GlobTherm
#'   v1.0).
#' @seealso \code{\link{globtherm_tree}}
"globtherm_traits"

#' Taxonomic phylogeny aligned to GlobTherm
#'
#' A taxonomic tree built from GlobTherm's Class / Order / Family /
#' Genus / Species columns with Grafen rank-based branch lengths and
#' pruned to the species in \code{globtherm_traits}.
#'
#' @format An \code{ape::phylo} object with 1,969 tips.
#' @source Bennett et al. 2018 (taxonomy); Grafen 1989 (branch
#'   lengths).
#' @seealso \code{\link{globtherm_traits}}
"globtherm_tree"


# ---------------------------------------------------------------------------
# 5. LepTraits butterflies ---------------------------------------------------
# ---------------------------------------------------------------------------

#' LepTraits butterfly traits, tidied for imputation benchmarks
#'
#' A wide species-by-trait matrix derived from LepTraits 1.0 (Shirey
#' et al. 2022), with four continuous / count traits, twelve monthly
#' flight-presence indicators (covariates), and Family / Genus
#' taxonomy.
#'
#' @format A \code{data.frame} with 8,219 rows and 18 columns.  Row
#'   names are species keys matching \code{leptraits_tree$tip.label}.
#' \describe{
#'   \item{wingspan_lower}{Wingspan lower bound (mm), numeric}
#'   \item{forewing_length_lower}{Forewing length lower bound (mm)}
#'   \item{flight_duration}{Number of flight months per year (count)}
#'   \item{n_hostplant_families}{Number of host-plant families
#'     (count, log-transform at use-site)}
#'   \item{Jan, Feb, ..., Dec}{Twelve binary indicators: 1 if the
#'     species flies in that calendar month, 0 otherwise.  Useful as
#'     phenology covariates.}
#'   \item{family_name}{Lepidopteran family, factor}
#'   \item{genus_name}{Lepidopteran genus, factor}
#' }
#' @source Shirey et al. 2022 \emph{Sci. Data} 9:382.
#' @seealso \code{\link{leptraits_tree}}
"leptraits_traits"

#' Taxonomic phylogeny aligned to LepTraits
#'
#' A taxonomic tree built from LepTraits' Family / Genus / Species
#' columns with Grafen rank-based branch lengths and pruned to the
#' species in \code{leptraits_traits}.
#'
#' @format An \code{ape::phylo} object with 8,219 tips.
#' @source Shirey et al. 2022 (taxonomy); Grafen 1989 (branch
#'   lengths).
#' @seealso \code{\link{leptraits_traits}}
"leptraits_tree"


# ---------------------------------------------------------------------------
# 6. Bundled-from-pigauto datasets (pass-through) ----------------------------
# ---------------------------------------------------------------------------

#' Delhey 5,809-species bird plumage colour dataset (real)
#'
#' A real comparative dataset of 5,809 bird species with plumage
#' lightness (male / female), six climate covariates, and family
#' taxonomy, originally compiled by Delhey 2019 \emph{Glob. Ecol.
#'   Biogeogr.} 28:860 and shipped with \pkg{pigauto}.  Bundled here
#' so it ships alongside the BACE-side benchmarks too.
#'
#' @format A \code{data.frame} with 5,809 rows and 10 columns:
#' \describe{
#'   \item{Species_Key}{Species key matching \code{tree_delhey$tip.label}}
#'   \item{family}{Bird family, character}
#'   \item{annual_mean_temperature, annual_precipitation,
#'         percent_tree_cover, mean_temperature_of_warmest_quarter,
#'         precipitation_of_warmest_quarter, midLatitude}{Continuous
#'         climate / geography covariates}
#'   \item{lightness_male, lightness_female}{Plumage lightness scores,
#'         continuous}
#' }
#' @source Delhey 2019 \emph{Glob. Ecol. Biogeogr.} 28:860.
#' @seealso \code{\link{tree_delhey}}
"delhey5809"

#' Bird phylogeny aligned to delhey5809
#'
#' A bird supertree pruned to the 5,809 species in
#' \code{\link{delhey5809}}.  Tip labels match
#' \code{delhey5809$Species_Key}.
#'
#' @format An \code{ape::phylo} object with 5,809 tips.
#' @source Jetz et al. 2012 \emph{Nature} 491:444 bird supertree,
#'   subset to the Delhey 2019 species set.
#' @seealso \code{\link{delhey5809}}
"tree_delhey"

#' Simulated multi-observation CTmax dataset
#'
#' A simulated multi-observation dataset for testing pigauto's
#' observation-level covariate refinement (\code{obs_refine} MLP).
#' 1,464 observations across 300 species with two covariates per
#' observation (acclimation temperature and CTmax).
#'
#' @format A \code{data.frame} with 1,464 rows and 3 columns:
#' \describe{
#'   \item{species}{Species key matching \code{tree300$tip.label}.}
#'   \item{acclim_temp}{Acclimation temperature (degC), numeric.}
#'   \item{CTmax}{Critical thermal maximum (degC), numeric (the
#'     simulated trait of interest).}
#' }
#' @source Simulated as part of \pkg{pigauto}.
#' @seealso \code{\link{tree300}}
"ctmax_sim"

#' 300-tip bird tree shipped with pigauto
#'
#' A 300-tip bird tree, shared between AVONET-300 (in pigauto) and
#' \code{\link{ctmax_sim}}.  Bundled here for the multi-obs CTmax
#' benchmarks.
#'
#' @format An \code{ape::phylo} object with 300 tips.
#' @source Jetz et al. 2012 bird supertree, subset to a 300-species
#'   sample.
#' @seealso \code{\link{ctmax_sim}}, \code{\link{trees300}}
"tree300"

#' Posterior bird-tree sample for the 300-tip CTmax / AVONET-300 set
#'
#' A 50-tree posterior sample over the 300 species in
#' \code{\link{tree300}}, for tree-uncertainty workflows
#' (Nakagawa & de Villemereuil 2019 \emph{Syst. Biol.} 68:632).
#'
#' @format A \code{multiPhylo} of length 50, each element an
#'   \code{ape::phylo} with 300 tips matching \code{tree300}.
#' @source Posterior sampled from Jetz et al. 2012 bird-supertree
#'   posterior; subset to the 300-species sample matching
#'   \code{tree300}.
#' @seealso \code{\link{tree300}}, \code{\link{ctmax_sim}}
"trees300"

# ---------------------------------------------------------------------------
# 7. AVONET birds ------------------------------------------------------------
# ---------------------------------------------------------------------------

#' AVONET bird traits, tidied for imputation benchmarks
#'
#' A wide species-by-trait matrix derived from AVONET 1.0 (Tobias et
#' al. 2022 \emph{Ecol. Lett.} 25:581), pre-aligned to the
#' \code{\link{avonet_tree}} tip set.  Eight continuous morphometric /
#' range / centroid traits, two ordinals (\code{habitat_density},
#' \code{migration}), three categoricals (\code{trophic_level},
#' \code{primary_lifestyle}, \code{habitat}), plus three higher-
#' taxonomy columns.  Continuous values are on the original (raw)
#' scale; users typically log-transform mass / wing / beak / tarsus /
#' tail / range size at use-site.
#'
#' @format A \code{data.frame} with 9,993 rows (one row per species)
#'   and 16 columns.  Row names are species keys ("Genus_species")
#'   matching \code{avonet_tree$tip.label}.
#' \describe{
#'   \item{mass_g}{Body mass in g (numeric, raw scale)}
#'   \item{wing_length_mm}{Wing length in mm (numeric)}
#'   \item{beak_length_culmen_mm}{Beak length to culmen base in mm}
#'   \item{tarsus_length_mm}{Tarsus length in mm (numeric)}
#'   \item{tail_length_mm}{Tail length in mm (numeric)}
#'   \item{range_size_km2}{Breeding range size in km^2 (numeric)}
#'   \item{centroid_lat}{Range-centroid latitude in decimal degrees}
#'   \item{centroid_lon}{Range-centroid longitude in decimal degrees}
#'   \item{habitat_density}{Habitat density (1=dense, 2=semi-open,
#'     3=open), ordered factor}
#'   \item{migration}{Migration class (1=sedentary, 2=partial,
#'     3=full migrant), ordered factor}
#'   \item{trophic_level}{Carnivore / Herbivore / Omnivore /
#'     Scavenger, factor (K=4)}
#'   \item{primary_lifestyle}{Aerial / Aquatic / Generalist /
#'     Insessorial / Terrestrial, factor (K=5)}
#'   \item{habitat}{Coastal / Desert / Forest / Grassland / Human
#'     Modified / Marine / Riverine / Rock / Shrubland / Wetland /
#'     Woodland, factor (K=11)}
#'   \item{order_name}{AVONET-listed bird order, factor}
#'   \item{family_name}{AVONET-listed bird family, factor}
#'   \item{genus_name}{Genus parsed from species key, factor}
#' }
#' @source Tobias et al. 2022 \emph{Ecol. Lett.} 25:581.  Constructed
#'   by \code{data-raw/make_avonet.R} from the AVONET 1.0 trait CSV.
#' @seealso \code{\link{avonet_tree}}
"avonet_traits"

#' Hackett bird phylogeny aligned to AVONET 1.0
#'
#' A 9,993-tip bird phylogeny matching the AVONET 1.0 species set.
#' Tip labels use the underscore form ("Genus_species") and align
#' one-to-one with \code{rownames(avonet_traits)}.
#'
#' @format An \code{ape::phylo} object with 9,993 tips and resolved
#'   ultrametric branch lengths in millions of years.
#' @source Hackett et al. 2008 \emph{Science} 320:1763 bird
#'   phylogeny, intersected with AVONET 1.0 species coverage.
#' @seealso \code{\link{avonet_traits}}
"avonet_tree"


#' @title A dataset of reef fish growth trajectories
#'
#' @description This dataset was compiled by Morais and Bellwood (2018), and the original can be accessed in the Supporting Information - DS1 of the article with references, locality names and complete names of the categorical traits levels. This is a streamlined version with variables formatted to be used directly within the function \code{\link{predictKmax}}. The only change is that the variable "Position" is here renamed as "Movimentation", which is obviously a less adequate name.
#' 
#' @format A data frame with 1921 rows and 19 variables:
#' \describe{
#'   \item{Family, Taxonomic family}
#'   \item{Species, Species name}
#'   \item{SpecCode, Species code from FishBase}
#'   \item{MaxSizeTL, Maximum recorded size for the species, referring to Total Length in cm}
#'   \item{Diet, Dietary category, in seven levels (see article for explanation)}
#'   \item{Schooling, Schooling behaviour or gregariousness, in six levels (see article for explanation)}
#'   \item{Movimentation, Position relative to the reef, combining horizontal and vertical components, in six levels (see article for explanation)}
#'   \item{a, Length-weight regression parameter ‘a’, estimated from the Bayesian Hierarchical Model from Froese et al. (2014) and FishBase}
#'   \item{b, Length-weight regression parameter ‘b’, estimated from the Bayesian Hierarchical Model from Froese et al. (2014) and FishBase}
#'   \item{FormFactor, Body shape factor (or Form Factor in Froese (2006), reference in article) measuring the extent to which a fish is elongated or deep-bodied. It can be perceived as the ‘a’ parameter value a fish species should have if its b = 3}
#'   \item{Linf, Population asymptotic length in cm as reported by the growth study}
#'   \item{LinfType, Type of measure used to derive Linf. SL = Standard Length; FL = Fork Length; TL = Total Length; NG = Not Given (conservatively assumed to be TL)}
#'   \item{LinfTL, Population asymptotic length in cm, converted to total length if Linf was reported in a measure other than this (see LinfType)}
#'   \item{K, The Von Bertalanffy Growth coefficient K, as reported by the growth study}
#'   \item{Kmax, The standardized growth coefficient Kmax (see article for explanation)}
#'   \item{O, The standardized growth coefficient Ø, also termed Growth Performance Index (see article for explanation)}
#'   \item{lon, Longitudinal geographic coordinate of the population studied (see article for its derivation)}
#'   \item{lat, Latitudinal geographic coordinate of the population studied (see article for its derivation)}
#'   \item{sstmean, Mean sea surface temperature from the geographic coordinate of the population studied, obtained from Bio-ORACLE (see article for explanation)}
#'   \item{pelnpp, Mean pelagic net primary productivity from the geographic coordinate of the population studied, modelled from chlorophyll concentration and photosynthetic active radiation data}
#'
#' @seealso \code{\link{predictKmax}}
#'
#' @source \url{https://research.jcu.edu.au/researchdata/default/detail/a1f460789431ac7aa431212439a7ccb2/}
#'
#' @references <Morais, R.A., and Bellwood, D.R. (2018). Global drivers of reef fish growth. Fish and Fisheries. 19, 874–889. doi:10.1111/faf.12297>
#' @references <Froese, R. (2006). Cube law, condition factor and weight-length relationships: History, meta-analysis and recommendations. Journal of Applied Ichthyology, 22, 241–253. https://doi.org/10.1111/j.1439-0426.2006.00805.x>
#' @references <Froese, R., Thorson, J. T., & Reyes, R. B. (2014). A Bayesian approach for estimating length-weight relationships in fishes. Journal of Applied Ichthyology, 30, 78–85. https://doi.org/10.1111/jai.12299>

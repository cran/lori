#' Alpine plant communities in Aravo, France: Abundance data and covariates
#'
#' Originally published in Choler, P. 2005.
#' Consistent shifts in Alpine plant traits along a
#' mesotopographical gradient. Arctic, Antarctic, and
#' Alpine Research 37: 444–453.
#'
#' Analysed in Dray, S., Choler, P., Dolédec, S.,
#' Peres-Neto, P.R., Thuiler, W., Pavoine, S. & ter Braak,
#' C.J.F. 2014. Combining the fourth-corner and the RLQ
#' methods for assessing trait responses to environmental
#' variation. Ecology 95: 14-21
#'
#' Description from Dray et al. (2014): Community composition of vascular
#' plants was determined in 75 5 × 5 m plots. Each site was
#' described by six environmental variables: mean snowmelt
#' date over the period 1997–1999, slope inclination, aspect,
#' index of microscale landform, index of physical disturbance
#' due to cryoturbation and solifluction, and an index of
#' zoogenic disturbance due to trampling and burrowing activities
#' of the Alpine marmot. All variables are quantitative except
#' the landform and zoogenic disturbance indices that are
#' categorical variables with five and three categories,
#' respectively. Eight quantitative functional traits (i.e.,
#' vegetative height, lateral spread, leaf elevation angle,
#' leaf area, leaf thickness, specific leaf area, mass-based
#' leaf nitrogen content, and seed mass) were measured on the
#'  82 most abundant plant species (out of a total of 132
#'  recorded species).
#'
#' @docType data
#'
#' @usage data(aravo)
#'
#' @format A list with 4 attributes:
#' \describe{
#'   \item{spe}{abundance table of 82 species in 75 environments}
#'   \item{env}{a matrix of 6 covariates for the 75 environments}
#'   \item{traits}{a matrix of 8 covariates for the 82 species}
#'   \item{spe.names}{a vector of 82 species names}
#' }
#' @source \url{http://pbil.univ-lyon1.fr/ade4/ade4-html/aravo.html}
"aravo"

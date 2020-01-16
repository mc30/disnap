###############################################
# Visualisation
###############################################


#' @title Plot PPP objects
#' @description Plot all PPP objects from the provided list.
#'
#' @param lst List of PPP objects.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
plotPPPlist <- function(lst) {
  old.par <- par(mfrow = c(1, length(lst)))
  for (i in 1:length(lst))
    plot(lst[[i]], main = names(lst)[i])
  par(old.par)
}


#' @title Plot density of PPP objects
#' @description Plots density of all PPP objects from the provided list.
#'
#' @param lst List of PPP objects.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
densityPPPlist <- function(lst) {
  # require(spatstat)
  
  old.par <- par(mfrow = c(2, length(lst)))
  
  for (i in 1:length(lst))
    plot(density(lst[[i]]), main = "All tested herds")
  
  for (i in 1:length(lst)) {
    pitem <- lst[[i]]
    plot(density(pitem[which(pitem$marks > 0)]), main = paste0("Positive herds in: ", names(lst)[i])) # TODO: correct for factors
  }
  
  par(old.par)
}
###############################################
# Helper functions
###############################################


#' @title Make a simple plot with a centered text
#' @description Makes an empty plot with a provided text in the middle.
#'
#' @param text Text for the plot.
#' @param ... Additional graphic parameters.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
textPlot <- function(text = "", ...) {
  old.par <- par(mar = c(0, 0, 0, 0))
  
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  
  text(x = 0.5, y = 0.5, text, ...)
  
  par(old.par)
}


#' @title Plot map
#' @description Plots a map for a specified timepoint.
#'
#' @param coords Coordinates of units.
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param column Column (i.e. timepoints) to analyse. Only one column will be used.
#' @param bLegend Whether to add legend with statuses.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
plotMap <- function(coords, statusMat, column = 1, bLegend = FALSE) {
  plot(coords$x, coords$y, col = statusMat[, column] + 1, pch = 19, asp = 1, 
       main = colnames(statusMat)[column], cex = 0.4, xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  
  # Plot points with NA status
  points(coords$x[which(is.na(statusMat[, column]))], coords$y[which(is.na(statusMat[, column]))], cex = 0.4)
  
  types <- sort(unique(as.vector(statusMat)))
  if (bLegend)
    legend("topright", legend = types, pch = 19, col = types + 1)
}
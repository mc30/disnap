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






#' @title Plot strain distribution
#' @description Plots a barplot with strain distributions for considered timepoints.
#'
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param colors Colors to be used.
#' @param exclude Strains that should be excluded.
#' @param ... Additional graphic parameters.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
plotStrainDistribution <- function(statusMat, colors = 1:ncol(statusMat), exclude = c(), ...) {
  statusMat[statusMat %in% exclude] <- NA # Exclude particular STs
  
  # Merge tables for each timepoint
  for (i in 1:ncol(statusMat)) {
    tab <- table(statusMat[statusMat[, i] != 0, i])
    
    if (i == 1)
      temp <- tab
    else
      temp <- merge(temp, tab, by = 'Var1', all = TRUE)
  }
  
  temp[is.na(temp)] <- 0 # Replace NAs with zeros
  mat <- as.matrix(temp)
  mat <- mat[, -1] # Remove first column (STs)
  mat <- apply(mat, 1, as.numeric) # Make STs numeric
  rownames(mat) <- as.character(colnames(statusMat))
  colnames(mat) <- as.numeric(as.character(temp$Var1)) # STs as rownames
  mat <- mat[, order(as.numeric(colnames(mat)))] # Sort by STs
  
  barplot(mat, beside = TRUE, axisnames = FALSE, col = colors, ...)
  
  axis(1, at = (ncol(statusMat) + 1) * c(1:length(colnames(mat))) - ncol(statusMat) / 2,
       labels = colnames(mat), 
       cex.axis = 0.8,
       col.axis = "black",
       las = 2
  )
  
  legend("topright", border = "true", rownames(mat), cex = 1.0, col = colors, pch = 15, bty = 2)
}
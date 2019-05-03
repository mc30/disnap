###############################################
# Infection status analysis
###############################################


#' @title Print probabilities of changing infection status.
#' @description Probabilities of switching status between two timepoints are calculated and printed.
#'
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param cols Columns (i.e. timepoints) to analyse. Only the first two will be used.
#' 
#' @return Returns the matrix with probabilities.
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
confMat <- function(statusMat, cols = c(1, 2)) {
  cm <- matrix(c(length(which(statusMat[, cols[1]] > 0 & statusMat[, cols[2]] > 0)),
                 length(which(statusMat[, cols[1]] > 0 & statusMat[, cols[2]] == 0)),
                 length(which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] > 0)),
                 length(which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] == 0))), nrow = 2, byrow = T)
  
  # print(cm)
  
  
  if (sum(cm[1, ]) != 0)
    cm[1, ] <- cm[1, ] / sum(cm[1, ])
  
  if (sum(cm[2, ]) != 0)
    cm[2, ] <- cm[2, ] / sum(cm[2, ])
  
  # print(cm)
  
  print(paste(colnames(statusMat)[cols], collapse = " -> "))
  
  print(paste0("Probability to stay infected: ", signif(cm[1, 1], 2)))
  print(paste0("Probability to recover: ", signif(cm[1, 2], 2)))
  print(paste0("Probability to become infected: ", signif(cm[2, 1], 2)))
  print(paste0("Probability to stay negative: ", signif(cm[2, 2], 2)))
  
  return(cm)
}















# Distance distribution for case-pairs

colors <- c("red", "blue", "black", "darkgreen", "pink")
lineTypes <- c(2, 4, 1, 3, 3)


#' @title Plot histograms
#' @description .
#'
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param distMat Distance matrix.
#' @param cols Column (i.e. timepoints) to analyse. Only one column will be used.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
plotHists <- function(statusMat, distMat, cols = c(1, 2)) {
  distr1 <- distMat[which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] > 0), ]
  distr2 <- distMat[which(statusMat[, cols[1]] > 0 & statusMat[, cols[2]] > 0), ]
  distr3 <- distMat
  distr4 <- distMat[which(statusMat[, cols[1]] > 0), ]
  
  old.par <- par(mfrow = c(2, 2))
  
  maxDist <- max(distr1, max(distr2, max(distr3, distr4)))
  
  hist(distr1, breaks = seq(0, maxDist, length.out = 100), xlab = "Distance (km)", main = "NI to all", col = colors[1])
  hist(distr2, breaks = seq(0, maxDist, length.out = 100), xlab = "Distance (km)", main = "PS to all", col = colors[2])
  hist(distr3, breaks = seq(0, maxDist, length.out = 100), xlab = "Distance (km)", main = "All to all", col = colors[3])
  hist(distr4, breaks = seq(0, maxDist, length.out = 100), xlab = "Distance (km)", main = "Spring positives to all", col = colors[4])
  
  par(old.par)
}


#' @title Plot distributions
#' @description .
#'
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param distMat Distance matrix.
#' @param bTiff Whether to generate TIFF file.
#' @param cols Column (i.e. timepoints) to analyse. Only one column will be used.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
plotDistr <- function(statusMat, distMat, bTiff = FALSE, cols = c(1, 2)) {
  distr1 <- distMat[which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] > 0), ]
  distr2 <- distMat[which(statusMat[, cols[1]] > 0 & statusMat[, cols[2]] > 0), ]
  distr3 <- distMat
  distr4 <- distMat[which(statusMat[, cols[1]] > 0), ]
  distr5 <- distMat[which(statusMat[, cols[2]] > 0), ]
  
  
  if (bTiff)
    tiff(paste0("graphs/dist_distr", ".tiff"), width = 1200 * 2, height = 800 * 2, res = 72 * 2 *2, compression = "lzw+p")
  
  lwd <- 2.5
  
  plot(density(distr1), col = colors[1], xlab = "Distance (km)",
       ylim = c(0, max(density(distr1)$y) * 1.5),
       lwd = lwd, lty = lineTypes[1],
       main = "")
  
  if (length(distr2) > 0)
    lines(density(distr2), col = colors[2], lty = lineTypes[2], lwd = lwd)
  if (length(distr3) > 0)
    lines(density(distr3), col = colors[3], lty = lineTypes[3], lwd = lwd)
  if (length(distr4) > 0)
    lines(density(distr4), col = colors[4], lty = lineTypes[4], lwd = lwd)
  if (length(distr5) > 0)
    lines(density(distr5), col = colors[5], lty = lineTypes[5], lwd = lwd)
  
  
  legend("topright", 
         legend = c("All to all", "NP to all", "NN to all", 
                    "P* to all",
                    "*P to all"),
         col = colors[c(3, 1, 2, 4, 5)], lty = lineTypes[c(3, 1, 2, 4, 5)], lwd = lwd)
  
  if (bTiff)
    dev.off()
}
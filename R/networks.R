###############################################
# Network data
###############################################

#' @title Plot network data.
#' @description Makes a plot of a provided network.
#'
#' @param coords Coordinates of units.
#' @param adjMat Adjacency matrix for the network.
#' @param color Color of nodes (can be a vector).
#' @param ... Additional graphic parameters.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#' 
#' @export
plotNetwork <- function(coords, adjMat, color = "black", ...) {
  plot(coords$x, coords$y, 
       col = color,
       pch = 19, asp = 1,
       cex = 0.4, xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  
  links <- which(adjMat > 0, arr.ind = T)
  
  from <- links[, 1]
  to <- links[, 2]
  
  arrows(coords$x[from], coords$y[from], coords$x[to], coords$y[to], lwd = 3 * adjMat[which(adjMat > 0)] / max(adjMat), 
         length = 0.1, angle = 30)
}


#' @title Create a random network
#' @description Creates a random network with the specified number of nodes and links.
#'
#' @param nv Number of nodes (vertices).
#' @param ne Number of edges (links).
#' @param weights Weight of links (vector of length \code{ne}).
#' @param bDirected Whether to make the matrix symmetrical (i.e. biderectional links).
#' 
#' @return Adjacency matrix. Non-zero element \code{adjMat[i, j]} indicates the weight of the link between nodes \code{i} and \code{j}.
# @details .
#'
#' @author Mikhail Churakov
#' 
#' @export
createRandomNetwork <- function(nv, ne, weights = rep(1, ne), bDirected = FALSE) {
  adjMat <- matrix(rep(0, nv ^ 2), nrow = nv)
  
  src <- sample(1:nv, ne, replace = F)
  dest <- sample(1:nv, ne, replace = F)
  
  for (i in 1:ne) {
    adjMat[src[i], dest[i]] <- weights[i]
    
    if (bDirected) 
      adjMat[dest[i], src[i]] <- -weights[i]
    else
      adjMat[dest[i], src[i]] <- weights[i]
  }
  
  return(adjMat)
}


#' @title Permute links.
#' @description Permutes links in the network.
#'
#' @param adjMat Adjacency matrix for the network.
#' 
#' @return Adjacency matrix with permuted links.
#' @details Permutation is equivalent to creating a new adjustment matrix with the same number of nodes and edges.
#'
#' @author Mikhail Churakov
#' 
#' @export
permuteLinks <- function(adjMat) {
  nv <- nrow(adjMat)
  ne <- length(which(adjMat > 0))
  
  if (isSymmetric(adjMat))
    newMat <- createRandomNetwork(nv, ne / 2, bDirected = FALSE) # undirected
  else
    newMat <- createRandomNetwork(nv, ne, bDirected = TRUE) # directed
  return(newMat)
}


#' @title Calculate node exposure.
#' @description Calculates incoming links for a node.
#'
#' @param node Index of the specified node.
#' @param adjMat Adjacency matrix for the network.
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param column Column to analyze.
#' 
#' @return Statuses of incoming links for the specified node.
# @details .
#'
#' @author Mikhail Churakov
#' 
#' @export
calculateNodeExposure <- function(node, adjMat, statusMat, column = 1) {
  links <- which(adjMat > 0, arr.ind = T)
  
  sel <- which(links[, 1] == node)
  
  # print(length(sel)) # Incoming links
  # print(links[sel, 2]) # Sources
  
  return(statusMat[links[sel, 2], column])
}


#' @title Calculate mean exposure.
#' @description Calculates mean exposure from incoming links for a set of node.
#'
#' @param selNodes Indices of the specified nodes.
#' @param adjMat Adjacency matrix for the network.
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param column Column to analyze.
#' @param bLog Whether to print infective links.
#' 
#' @return A vector with two elements: number of incoming links, number of incoming links from infectious units (status > 0).
#' @details Implemented via matrix multiplications.
#'
#' @author Mikhail Churakov
#' 
#' @export
calculateMeanExposure <- function(selNodes, adjMat, statusMat, column = 1, bLog = FALSE) {
  vec <- rep(0, nrow(statusMat))
  vec[selNodes] <- 1
  rownames(statusMat)[which(vec > 0)]
  
  res <- adjMat %*% vec # from nodes
  
  res[which(res > 0)] # Should be 1
  rownames(statusMat)[which(res > 0)]
  
  links <- length(which(res > 0))
  infLinks <- length(which(statusMat[which(res > 0), column] > 0))
  
  sel <- which(res > 0)
  # sel <- which(statusMat[which(res > 0), column] > 0)
  
  if (bLog & length(which(res > 0)) > 0)
    print(paste0("Infective links from: ", paste(rownames(statusMat)[sel], collapse = ", ")))
  
  return(c(infLinks, links))
}





#' @title Calculate mean exposure.
#' @description Calculates mean exposure for NI and PS units.
#'
#' @param adjMat Adjacency matrix for the network.
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param cols Column (i.e. timepoints) to analyse. Only one column will be used.
#' @param nsim Number of simulations.
#' @param bLog Whether to print extra info.
#' 
#' @return List of observed and simulated mean exposures for NI and PS units: obsNI, obsPS, mexpNI, mexpPS.
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
meanExpFun <- function(adjMat, statusMat, cols = 1:2, nsim = 1000, bLog = FALSE) {
  if (bLog)
    print("Observed mean exposure")
  
  selNodes <- which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] > 0)
  if (bLog)
    print(paste0(length(selNodes), " NI herds: ", paste(rownames(statusMat)[selNodes], collapse = ", ")))
  else
    print(paste0(length(selNodes), " NI herds"))
  res <- calculateMeanExposure(selNodes, adjMat, statusMat, column = cols[1], bLog = bLog)
  # if (bLog)
  print(paste0("Mean exposure of NI: ", signif(res[1] / res[2], 3), " (", res[1], " out of ", res[2], ")"))
  obsNI <- res[1] / res[2]
  
  
  selNodes <- which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] == 0)
  if (bLog)
    print(paste0(length(selNodes), " PS herds: ", paste(rownames(statusMat)[selNodes], collapse = ", ")))
  else
    print(paste0(length(selNodes), " PS herds"))
  res <- calculateMeanExposure(selNodes, adjMat, statusMat, column = cols[1], bLog = bLog)
  # if (bLog)
  print(paste0("Mean exposure of PS: ", signif(res[1] / res[2], 3), " (", res[1], " out of ", res[2], ")"))
  obsPS <- res[1] / res[2]
  
  
  mexpNI <- vector(mode = "numeric", length = nsim)
  mexpPS <- vector(mode = "numeric", length = nsim)
  
  pb <- txtProgressBar(style = 3)
  for (i in 1:nsim) {
    # print(paste0(i, " / ", nsim))
    
    num <- length(which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] > 0))
    
    # repeat {
    selNodes <- sample(1:nrow(statusMat), num, replace = F) # Permute nodes
    res <- calculateMeanExposure(selNodes, adjMat, statusMat, column = cols[1], bLog = bLog)
    if (bLog)
      print(paste0("Mean exposure of NI: ", signif(res[1] / res[2], 3), " (", res[1], " out of ", res[2], ")"))
    mexpNI[i] <- res[1] / res[2]
    
    #   if (!is.nan(mexpNI[i]) & !mexpNI[i] == 0) 
    #     break
    # }
    
    num <- length(which(statusMat[, cols[1]] == 0 & statusMat[, cols[2]] == 0))
    
    # repeat {
    selNodes <- sample(1:nrow(statusMat), num, replace = F) # Permute nodes
    res <- calculateMeanExposure(selNodes, adjMat, statusMat, column = cols[1], bLog = bLog)
    if (bLog)
      print(paste0("Mean exposure of PS: ", signif(res[1] / res[2], 3), " (", res[1], " out of ", res[2], ")"))
    mexpPS[i] <- res[1] / res[2]
    
    #   if (!is.nan(mexpPS[i]) & !mexpPS[i] == 0)
    #     break
    # }
    
    setTxtProgressBar(pb, i / nsim)
  }
  
  mexpNI[which(is.nan(mexpNI))] <- 0 # In case of NaN from 0 / 0
  mexpPS[which(is.nan(mexpPS))] <- 0 # In case of NaN from 0 / 0
  
  if (bLog) {
    print(mexpNI)
    print(mexpPS)
  }
  
  return(list(obsNI = obsNI, obsPS = obsPS, mexpNI = mexpNI, mexpPS = mexpPS))
}



#' @title Plot mean exposure.
#' @description Plot mean exposure for NI and PS units.
#'
#' @param obsNI Observed mean exposure for NI units.
#' @param obsPS Observed mean exposure for PS units.
#' @param mexpNI Simulated mean exposures for NI units.
#' @param mexpPS Simulated mean exposures for PS units.
#' @param bEllipse Whether to show ellipse containing a fraction of points.
#' @param elLevel Fraction of points that should be inside the ellipse (if bEllipse is TRUE).
#' @param quantileVals Vector of quantiles used to define significance levels.
#' @param bShowPerc Whether to show percentages of simulated values below the observed exposure.
#' @param ... Additional graphic parameters.
#' 
# @return List of observed and simulated mean exposures for NI and PS units: obsNI, obsPS, mexpNI, mexpPS.
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
plotMeanExp <- function(obsNI, obsPS, mexpNI, mexpPS, 
                        bEllipse = FALSE, elLevel = 0.95, 
                        quantileVals = c(0.05, 0.95), bShowPerc = TRUE, ...) {
  plot(mexpNI, mexpPS, pch = 19, cex = 0.5, col = "darkgray", ...)
  abline(0, 1, lty = 2)
  
  if (bEllipse) {
    require(car)
    el <- dataEllipse(mexpNI, mexpPS, levels = elLevel, draw = F)
    polygon(el, col = adjustcolor("darkgray", alpha.f = 0.5))
  }
  
  
  # abline(v = mean(mexpNI) + c(-2 * sd(mexpNI), 2 * sd(mexpNI)), lty = 2, col = "darkgray")
  # abline(h = mean(mexpPS) + c(-2 * sd(mexpPS), 2 * sd(mexpPS)), lty = 2, col = "darkgray")
  
  abline(v = quantile(mexpNI, quantileVals), lty = 1, col = "darkgray")
  abline(v = quantile(mexpNI, c(0, 1)), lty = 3, col = "lightgray")
  
  abline(h = quantile(mexpPS, quantileVals), lty = 1, col = "darkgray")
  abline(h = quantile(mexpPS, c(0, 1)), lty = 3, col = "lightgray")
  
  points(mean(mexpNI), mean(mexpPS), col = "black", pch = 19, lwd = 2) # Mean simulated
  points(obsNI[1], obsPS[1], col = "red", pch = 19, lwd = 2) # Observed
  
  
  x <- ecdf(res$mexpNI)(res$obsNI)
  y <- ecdf(res$mexpPS)(res$obsPS)
  
  # Show percentages of values below the observed ones
  if (bShowPerc) {
    mtext(paste0("P(NI) = ", signif(x, 3)), side = 3, at = obsNI[1], line = 0.05, col = "red", las = 1, cex = 0.7)
    mtext(paste0("P(PS) = ", signif(y, 3)), side = 4, at = obsPS[1], line = 0.05, col = "red", las = 0, cex = 0.7)
  }
  
  # mtext(signif(obsNI[1], 2), side = 3, at = obsNI[1], line = 0.5, col = "red", las = 2, cex = 0.7)
  # mtext(signif(mean(mexpNI), 2), side = 3, at = mean(mexpNI), line = 0.5, las = 2, cex = 0.7)
  
  # mtext(signif(obsPS[1], 2), side = 4, at = obsPS[1], line = 0.5, col = "red", las = 1, cex = 0.7)
  # mtext(signif(mean(mexpPS), 2), side = 4, at = mean(mexpPS), line = 0.5, las = 1, cex = 0.7)
}
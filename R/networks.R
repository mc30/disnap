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
#' @param nodes Indices of the specified nodes.
#' @param adjMat Adjacency matrix for the network.
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param column Column to analyze.
#' 
#' @return A vector with two elements: number of incoming links, number of incoming links from infectious units (status > 0).
# @details .
#'
#' @author Mikhail Churakov
#' 
#' @export
calculateMeanExposure <- function(nodes, adjMat, statusMat, column = 1) {
  links <- 0
  infLinks <- 0
  for (node in nodes) {
    contacts <- calculateNodeExposure(node, adjMat, statusMat, column)
    print(paste0(node, ": ", paste0(contacts, collapse = ", ")))
    
    links <- links + length(contacts)
    infLinks <- infLinks + length(which(contacts > 0))
  }
  
  return(c(infLinks, links))
}
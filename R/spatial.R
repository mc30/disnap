###############################################
# Spatial analysis
###############################################


#' @title Create a list of PPP objects
#' @description Status matrix and corresponding coordinates are used to create a list of PPP objects for further spatial analysis.
#' 
#' @param statusMat A matrix with statuses (NA, 0, any positive integer) of each unit. Rownames indicate location names, colnames indicate names of timepoints (e.g. years).
#' @param coords Coordinates of units.
#' @param cols Columns (i.e. timepoints) to analyse.
#' @param bPosOnly Whether to keep only locations of positive units. 
#' 
#' @return List of PPP objects.
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
createPPPList <- function(statusMat, coords, cols = c(1, 2), bPosOnly = FALSE) {
  # library(spatstat)
  # library(maptools)
  # library(raster)
  
  # Create a list of ppp
  lst <- vector("list", length(cols))
  for (i in cols) {
    if (bPosOnly)
      sel <- which(statusMat[ , cols[i]] > 0)
    else
      sel <- 1:nrow(statusMat)
    pitem <- ppp(x = coords$x[sel], y = coords$y[sel],
                 window = owin(range(coords$x[sel]), range(coords$y[sel])),
                 # marks = as.factor(statusMat[, cols[1]])
                 marks = statusMat[sel, cols[i]]
    )
    
    lst[[i]] <- pitem
  }
  names(lst) <- colnames(statusMat)[cols]
  
  return(lst)
}


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


#' @title Perform SaTScan test on PPP objects
#' @description Runs SaTScan test (implemented in smacpod package) on each PPP object from the provided list.
#'
#' @param lst List of PPP objects.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
satscanPPPlist <- function(lst) {
  # require(smacpod)
  
  old.par <- par(mfrow = c(1, length(lst)))
  
  for (i in 1:length(lst)) {
    if (length(levels(lst[[i]]$marks)) == 0) { # Empty plot if no cases
      textPlot("Could not run spscan.test")
      next
    }
    
    out <- spscan.test(lst[[i]], case = 2, alpha = 0.1)
    plot(out, chars = c(1, 20), main = paste0("Most likely cluster: ", names(lst)[i]))
  }
  
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


#' @title Plot log-ratio for PPP objects
#' @description Plots log-ratio [...] for each PPP object from the provided list.
#'
#' @param lst List of PPP objects.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
logrrPPPlist <- function(lst) {
  # require(spatstat)
  
  old.par <- par(mfrow = c(1, length(lst)))
  
  for (i in 1:length(lst)) {
    if (length(levels(lst[[i]]$marks)) == 0) { # Empty plot if no cases
      textPlot("Could not run logrr")
      next
    }
    
    r2 <- logrr(lst[[i]], sigma = spatstat::bw.scott)
    plot(r2, main = "")
    title(main = paste0("Log-ratio for: ", names(lst)[i]), cex.main = 0.75)
  }
  
  par(old.par)
}



#' @title Perform K-function test
#' @description Plots results of K-function spatial clustering test.
#'
#' @param lst List of PPP objects.
#' @param nsim Number of simulations. Default is 10.
#' 
# @return .
#' @details K-function is measured for the observed and simulated spatial dictributions.
#'
#' @author Mikhail Churakov (\email{mikhail.churakov@@gmail.com}).
#'
#' @export
spatclustKfun <- function(lst, nsim = 10) {
  old.par <- par(mfrow = c(1, length(lst)))
  
  for (i in 1:length(lst)) {
    if (length(levels(lst[[i]]$marks)) == 0) { # Empty plot if no cases
      textPlot("Could not run K-function test")
      next
    }
    
    # K function
    kd1 = kdest(lst[[i]])
    # plot(kd1, iso ~ r, ylab = "difference", legend = F, main = "")
    kd2 = kdest(lst[[i]], nsim = nsim, level = 0.8)
    plot(kd2, legend = F)
    
    # kdplus.test performs a global test of clustering for comparing cases and controls using the method
    # of Diggle and Chetwynd (1991). It relies on the difference in estimated K functions.
    kdplus.test(kd2)
    
    # plot(r2, main = "")
    title(main = paste0("K function: ", names(lst)[i]), cex.main = 0.75)
  }
  
  par(old.par)
}



#' @title Perform Cuzick-Edwards' kNN test
#' @description Plots results of Cuzick-Edwards' kNN clustering test.
#'
#' @param lst List of PPP objects.
#' @param nn Number of nearest neighbours. Default is 20.
#' @param nsim Number of simulations. Default is 500.
#' @param signifLevel Significance level (p-value cut-off). Default is 0.05.
#' 
# @return .
# @details .
#'
#' @author Mikhail Churakov
#'
#' @export
spatclustkNN <- function(lst, nn = 20, nsim = 500, signifLevel = 0.05) {
  # require(smacpod)
  
  old.par <- par(mfrow = c(1, length(lst)))
  
  for (i in 1:length(lst)) {
    if (length(levels(lst[[i]]$marks)) == 0) { # Empty plot if no cases
      textPlot("Could not run kNN test")
      next
    }
    
    # kNN
    x <- qnn.test(lst[[i]], nsim = nsim, q = 1:min(nn, lst[[i]]$n))
    
    # plot(x$qsum$q, x$qsum$Tq)
    
    plot(x$qsum$q, x$qsum$pvalue, 
         # ylim = c(0, max(x$qsum$pvalue)), 
         ylim = c(0, 1), 
         type = "l", xlab = "Number of neighbours", ylab = "p-value")
    points(x$qsum$q, x$qsum$pvalue, col = ifelse(x$qsum$pvalue < signifLevel, "red", "black"))
    abline(h = signifLevel, col = "red")
    
    title(main = paste0("kNN: ", names(lst)[i]), cex.main = 0.75)
  }
  
  par(old.par)
}
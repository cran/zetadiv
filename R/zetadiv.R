.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is zetadiv 1.2.1")
  packageStartupMessage("Package \"zetadiv\" was built under R.4.1.3")
}


##################
##MAIN FUNCTIONS##
##################


#' Zeta diversity decline using Monte Carlo sampling
#'
#' Computes zeta diversity, the number of species shared by multiple assemblages, for a range of orders (number of assemblages or sites), using combinations of sampled sites, and fits the decline to an exponential and a power law relationship.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param xy Site coordinates. This is only used if \code{NON} = TRUE or \code{DIR} = TRUE.
#' @param orders  Range of number of assemblages or sites for which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta diversity is computed for each number of assemblages or sites.
#' @param sd.correct Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) or not (using the number of site combinations as the denominator).
#' @param sd.correct.adapt Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) if \code{sam} is higher than the number of possible combinations, or not (using the number of site combinations as the denominator) if \code{sam} is lower than the number of possible combinations. If \code{sd.correct.adapt = TRUE}, it takes precedence over \code{sd.correct}.
#' @param confint.level  Percentage for the confidence intervals of the coefficients from the regressions.
#' @param sd.plot  Boolean value (TRUE or FALSE) indicating if the standard deviation of each zeta diversity value must be plotted.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1. Has no effect if \code{normalize} != \code{FALSE}.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param NON Boolean value (TRUE or FALSE) indicating if the number of species in common should only be counted for the nearest neighbours.
#' @param FPO A vector with the coordinates of the fixed point origin from which the zeta diversity will be computed (overrides NON). In that case, \eqn{\zeta_1} is the number of species in the closest site to the FPO, \eqn{\zeta_2} is the number of species shared by the 2 closest sites, etc.
#' @param DIR Boolean value (TRUE or FALSE) indicating if zeta diversity must be computed using a directed nearest neighbour scheme in the direction away from the FPO, starting from any site.
#' @param empty.row Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity, and computation errors for the normalized version of zeta due to divisions by 0. Options are "\code{empty}" to let the data untreated, "\code{remove}" to remove the empty rows, 0 to set the normalized zeta to 0 when zeta is divided by 0 during normalization (sites share no species, so are completely dissimilar), and 1 to set the normalized zeta to 1 when zeta is divided by 0 during normalization (i.e. sites are perfectly similar).
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @param silent  Boolean value (TRUE or FALSE) indicating if messages must be printed.
#' @details If the number of combinations of sites is lower than the value of the parameter \code{sam}, all the combinations are used and an exact solution is computed. In that case, using the number of site combinations as the denominator may be appropriate to compute the standard deviation, if all sites were sampled and the zeta values. This can be adjusted with parameters \code{sd.correct} and \code{sd.correct.adapt}.
#' @details \code{Zeta.decline.mc} is faster than \code{\link{Zeta.decline.ex}} to compute the exact value of zeta diversity when the number of species is higher than \eqn{C^N_{i}}, where \emph{N} is the total number of sites and \emph{i} is the order of zeta.
#' @details The exponential and the power law fit are performed using linear regressions on log-transformed data (only the zeta values are log-transformed for the exponential fit, and both the orders and the zeta values are log-transformed for the power law fit).
#' @details \code{Zeta.decline.mc} enables accomodating richness heterogeneity by setting \code{normalize = "Jaccard"}, \code{normalize = "Sorensen"} or \code{normalize = "Simpson"}. This cannot be performed by \cr \code{\link{Zeta.decline.ex}}.
#' @return \code{Zeta.decline.mc} returns a list containing the following components:
#' @return \item{zeta.order}{The number of assemblages or sites for which the zeta diversity was computed.}
#' @return \item{combinations}{The number of possible combinations of sites for the chosen orders.}
#' @return \item{zeta.val}{The zeta diversity values.}
#' @return \item{zeta.val.sd}{The zeta diversity standard deviation values.}
#' @return \item{zeta.ratio}{The ratio of zeta diversity values by the zeta diversity values at the lower order \eqn{\zeta_i / \zeta_{i-1}}.}
#' @return \item{zeta.exp}{Object of class "\code{lm}", containing the output of the exponential regression.}
#' @return \item{zeta.exp.confint}{The confidence intervals of the coefficients of the exponential regression.}
#' @return \item{zeta.pl}{Object of class "\code{lm}", containing the output of the power law regression.}
#' @return \item{zeta.pl.confint}{The confidence intervals of the coefficients of the power law regression.}
#' @return \item{aic}{AIC values for \code{zeta.exp} and \code{zeta.pl}.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Zeta.order.mc}}, \code{\link{Plot.zeta.decline}}
#' @examples
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[,1:2]
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' dev.new(width = 12, height = 4)
#' zeta.bird <- Zeta.decline.mc(data.spec.bird, xy.bird, orders = 1:5, sam = 100,
#'    NON = TRUE)
#' zeta.bird
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' dev.new(width = 12, height = 4)
#' zeta.marion <- Zeta.decline.mc(data.spec.marion, orders = 1:5, sam = 100,
#'      normalize = "Jaccard")
#' zeta.marion
#'
#' @export
#'
Zeta.decline.mc <- function(data.spec, xy = NULL, orders = 1:10, sam = 1000, sd.correct = TRUE, sd.correct.adapt = FALSE, confint.level = 0.95, sd.plot = TRUE, rescale = FALSE, normalize = FALSE, NON = FALSE, FPO = NULL, DIR = FALSE, empty.row = "empty", plot = TRUE, silent=TRUE){

  if (!inherits(data.spec,"data.frame")){
    stop(paste("Error: ",deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.",  sep = ""))
  }
  if(max(orders)>dim(data.spec)[1]){
    stop("Error: rrong value for \"orders\": the maximum value must be equal or lower than the number of sites.")
  }
  if(NON == TRUE && is.null(xy)){
    stop("Error: if NON = TRUE, xy must be non null.")
  }
  if(NON == TRUE && nrow(data.spec) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }

  if(empty.row == "remove"){
    if(length(which(rowSums(data.spec)))>0){
      data.spec <- data.spec[-which(rowSums(data.spec)==0),]
    }
  }

  x <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()

  if(is.null(FPO)){
    if(NON == FALSE){
      for(j in orders){
        if (j == 1){
          zeta.val[j]<-mean(rowSums(data.spec))
          if(sd.correct == TRUE & sd.correct.adapt == FALSE){
            zeta.val.sd[j] <- stats::sd(rowSums(data.spec))
          }else{
            zeta.val.sd[j] <- stats::sd(rowSums(data.spec))*nrow(data.spec)
          }
          if(rescale == TRUE || normalize != FALSE){
            zeta.val[j] <- 1
            zeta.val.sd[j] <- zeta.val.sd[j]/mean(rowSums(data.spec))
          }
        }else{
          if(choose(x, j)>sam){
            if(silent==FALSE){
              print(paste("Monte Carlo sampling for order",j))
            }
            u <- rep(NA, sam)
            for(z in 1:sam){
              samp <- sample(1:x, j, replace = FALSE)
              u[z] <- sum(apply(data.spec[samp, ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }
            }
          }else{
            if(silent==FALSE){
              print(paste("Exact solution for order",j))
            }
            u <- rep(NA, choose(x, j))
            samp <- utils::combn(1:x, j)
            for(z in 1:dim(samp)[2]){
              u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp[, z], ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[samp[, z], ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[samp[, z], ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }
            }
          }
          if(rescale == TRUE){
            z1 <- mean(rowSums(data.spec))
            u <- u / z1
          }

          zeta.val[j]<-mean(u)

          if(sd.correct.adapt == FALSE){
            if(sd.correct == TRUE){
              zeta.val.sd[j] <- stats::sd(u)
            }else{
              zeta.val.sd[j] <- stats::sd(u)*sqrt((length(u)-1)/length(u))
            }
          }else{
            if(x>sam){
              zeta.val.sd[j] <- stats::sd(u)*sqrt((length(u)-1)/length(u))
            }else{
              zeta.val.sd[j] <- stats::sd(u)
            }
          }
        }
      }
    }else{
      for(j in orders){
        if (j == 1){
          zeta.val[j]<-mean(rowSums(data.spec))
          if(sd.correct == TRUE & sd.correct.adapt == FALSE){
            zeta.val.sd[j] <- stats::sd(rowSums(data.spec))
          }else{
            zeta.val.sd[j] <- stats::sd(rowSums(data.spec))*nrow(data.spec)
          }
          if(rescale == TRUE || normalize != FALSE){
            zeta.val[j] <- 1
            zeta.val.sd[j] <- zeta.val.sd[j]/mean(rowSums(data.spec))
          }
          if(rescale == TRUE || normalize != FALSE){
            zeta.val[j] <- 1
            zeta.val.sd[j] <- zeta.val.sd[j]/mean(rowSums(data.spec))
          }
        }else{
          if(x>sam){
            u <- rep(NA, sam)
            samps <- sample(1:x, sam, replace = FALSE)
            for(z in 1:sam){
              samp <- samps[z]
              xy.dist <- (xy[,1]-xy[samp,1])^2+(xy[,2]-xy[samp,2])^2
              samp <- c(samp,order(xy.dist)[2:j])
              u[z] <- sum(apply(data.spec[samp, ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }
            }
          }else{
            u <- rep(NA, x)
            samps <- 1:x
            for(z in 1:x){
              samp <- samps[z]
              xy.dist <- (xy[,1]-xy[samp,1])^2+(xy[,2]-xy[samp,2])^2
              samp <- c(samp,order(xy.dist)[2:j])
              u[z] <- sum(apply(data.spec[samp, ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }
            }
          }

          if(rescale == TRUE & normalize == FALSE){
            z1 <- mean(rowSums(data.spec))
            u <- u / z1
          }

          zeta.val[j]<-mean(u)

          if(sd.correct.adapt == FALSE){
            if(sd.correct == TRUE){
              zeta.val.sd[j] <- stats::sd(u)
            }else{
              zeta.val.sd[j] <- stats::sd(u)*sqrt((length(u)-1)/length(u))
            }
          }else{
            if(x>sam){
              zeta.val.sd[j] <- stats::sd(u)*sqrt((length(u)-1)/length(u))
            }else{
              zeta.val.sd[j] <- stats::sd(u)
            }
          }

        }
      }
    }
  }else{
    if(DIR == FALSE){
      xy.dist <- (FPO[1]-xy[,1])^2+(FPO[2]-xy[,2])^2
      for(j in orders){
        samp <- order(xy.dist)[1:j]
        zeta.val[j] <- sum(apply(data.spec[samp, ], 2, prod))
        if (normalize == "Jaccard"){
          toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
          if(toto==0){
            if(empty.row == 0){
              zeta.val[j] <- 0
            }else if(empty.row == 1){
              zeta.val[j] <- 1
            }
          }else
            zeta.val[j] <- zeta.val[j] / toto
        }else if (normalize == "Sorensen"){
          toto <- (mean(apply(data.spec[samp, ], 1, sum)))
          if(toto==0){
            if(empty.row == 0){
              zeta.val[j] <- 0
            }else if(empty.row == 1){
              zeta.val[j] <- 1
            }
          }else
            zeta.val[j] <- zeta.val[j] / toto
        }else if (normalize == "Simpson"){
          toto <- (min(apply(data.spec[samp, ], 1, sum)))
          if(toto==0){
            if(empty.row == 0){
              zeta.val[j] <- 0
            }else if(empty.row == 1){
              zeta.val[j] <- 1
            }
          }else
            zeta.val[j] <- zeta.val[j] / toto
        }
        zeta.val.sd[j] <- 0
      }
    }else{
      for(j in orders){
        if(j==1){
          zeta.val[j]<-mean(rowSums(data.spec))
          if(sd.correct == TRUE & sd.correct.adapt == FALSE){
            zeta.val.sd[j] <- stats::sd(rowSums(data.spec))
          }else{
            zeta.val.sd[j] <- stats::sd(rowSums(data.spec))*nrow(data.spec)
          }
          if(rescale[j] == TRUE || normalize != FALSE){
            zeta.val[j] <- 1
            zeta.val.sd[j] <- zeta.val.sd/mean(rowSums(data.spec))
          }
        }else{
          xy.FPO <- as.matrix(xy-FPO)
          if(x>sam){
            u <- rep(NA, sam)
            samps <- sample(1:x, sam, replace = FALSE)
            for(z in 1:sam){
              samp <- samps[z]

              xy0 <- xy.FPO[samp,]
              no <- sqrt(sum(xy0^2))
              R <- matrix(c(xy0[2]/no,xy0[1]/no,-xy0[1]/no,xy0[2]/no),2,2)
              xy.FPO.tr <- apply(xy.FPO,1,function(xy.FPO,R,xy0,no){R %*% matrix(xy.FPO,2,1) - c(0,no)},R,xy0,no)
              xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
              xy.FPO.tr[samp,] <- c(0,0)
              xy.FPO.tr[which(xy.FPO.tr[,2]<=0),] <- NA

              xy.dist <- xy.FPO.tr[,1]^2+xy.FPO.tr[,2]^2

              if(length(which(!is.na(xy.dist)))>=(j-1)){
                samp <- c(samp,order(xy.dist)[1:(j-1)])
                u[z] <- sum(apply(data.spec[samp, ], 2, prod))
                if (normalize == "Jaccard"){
                  toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[z] <- 0
                    }else if(empty.row == 1){
                      u[z] <- 1
                    }
                  }else
                    u[z] <- u[z] / toto
                }else if (normalize == "Sorensen"){
                  toto <- (mean(apply(data.spec[samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[z] <- 0
                    }else if(empty.row == 1){
                      u[z] <- 1
                    }
                  }else
                    u[z] <- u[z] / toto
                }else if (normalize == "Simpson"){
                  toto <- (min(apply(data.spec[samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[z] <- 0
                    }else if(empty.row == 1){
                      u[z] <- 1
                    }
                  }else
                    u[z] <- u[z] / toto
                }
              }else{
                if(silent==FALSE){
                  print(paste("warning: number of sites away from the FPO too low to compute zeta for order",j))
                }
                u[z] <- NA
              }
            }
          }else{
            u <- rep(NA, x)
            samps <- 1:x
            for(z in 1:(x-j+1)){
              samp <- samps[z]

              xy0 <- xy.FPO[samp,]
              no <- sqrt(sum(xy0^2))
              R <- matrix(c(xy0[2]/no,xy0[1]/no,-xy0[1]/no,xy0[2]/no),2,2)
              xy.FPO.tr <- apply(xy.FPO,1,function(xy.FPO,R,xy0,no){R %*% matrix(xy.FPO,2,1) - c(0,no)},R,xy0,no)
              xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
              xy.FPO.tr[samp,] <- c(0,0)
              xy.FPO.tr[which(xy.FPO.tr[,2]<=0),] <- NA

              xy.dist <- xy.FPO.tr[,1]^2+xy.FPO.tr[,2]^2

              if(length(which(!is.na(xy.dist)))>=(j-1)){
                samp <- c(samp,order(xy.dist)[1:(j-1)])
                u[z] <- sum(apply(data.spec[samp, ], 2, prod))
                if (normalize == "Jaccard"){
                  toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[z] <- 0
                    }else if(empty.row == 1){
                      u[z] <- 1
                    }
                  }else
                    u[z] <- u[z] / toto
                }else if (normalize == "Sorensen"){
                  toto <- (mean(apply(data.spec[samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[z] <- 0
                    }else if(empty.row == 1){
                      u[z] <- 1
                    }
                  }else
                    u[z] <- u[z] / toto
                }else if (normalize == "Simpson"){
                  toto <- (min(apply(data.spec[samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[z] <- 0
                    }else if(empty.row == 1){
                      u[z] <- 1
                    }
                  }else
                    u[z] <- u[z] / toto
                }
              }else{
                if(silent==FALSE){
                  print(paste("warning: number of sites away from the FPO too low to compute zeta for order",j))
                }
                u[z] <- NA
              }
            }
          }
          zeta.val[j] <- mean(u,na.rm=TRUE)
          if(sd.correct.adapt == FALSE){
            if(sd.correct == TRUE){
              zeta.val.sd[j] <- stats::sd(u,na.rm=TRUE)
            }else{
              zeta.val.sd[j] <- stats::sd(u,na.rm=TRUE)*sqrt((length(which(!is.na(u)))-1)/length(which(!is.na(u))))
            }
          }else{
            if(x>sam){
              zeta.val.sd[j] <- stats::sd(u,na.rm=TRUE)*sqrt((length(which(!is.na(u)))-1)/length(which(!is.na(u))))
            }else{
              zeta.val.sd[j] <- stats::sd(u,na.rm=TRUE)
            }
          }
        }
      }
    }
  }

  ##create a single list for output
  zeta <- list()
  zeta$zeta.order <- orders
  zeta$combinations <- choose(x, orders)
  zeta$zeta.val <- zeta.val
  zeta$zeta.val.sd <- zeta.val.sd
  zeta$ratio <- zeta.val[2:length(zeta.val)]/zeta.val[1:(length(zeta.val)-1)]

  ##regression - exponential
  zeta.val.log <- log10(zeta.val)
  zeta.val.log[which(is.infinite(zeta.val.log))] <- NA
  zeta.exp <- stats::lm(zeta.val.log ~ c(orders), na.action = stats::na.omit)
  zeta$zeta.exp <- zeta.exp
  zeta$zeta.exp.confint <- suppressMessages(stats::confint(zeta.exp,level=confint.level))


  ##regression - power law
  zeta.pl <- stats::lm(zeta.val.log ~ log10(c(orders)), na.action = stats::na.omit)
  zeta$zeta.pl <- zeta.pl
  zeta$zeta.pl.confint <- suppressMessages(stats::confint(zeta.pl,level=confint.level))

  zeta$aic <- stats::AIC(zeta$zeta.exp, zeta$zeta.pl)


  ##Plot zeta and regressions
  if(plot == TRUE){
    Plot.zeta.decline(zeta, sd.plot=sd.plot)
  }

  return(zeta)

}





#' Zeta diversity for a specific number of assemblages or sites using Monte Carlo sampling
#'
#' Computes zeta diversity, the number of species shared by multiple assemblages, for a specific order (number of assemblages or sites).
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param xy Site coordinates. This is only used if \code{NON} = TRUE or \code{DIR} = TRUE.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta diversity is computed.
#' @param sd.correct Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) or not (using the number of site combinations as the denominator).
#' @param sd.correct.adapt Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) if \code{sam} is higher than the number of possible combinations, or not (using the number of site combinations as the denominator) if \code{sam} is lower than the number of possible combinations. If \code{sd.correct.adapt == TRUE}, it takes precedence over \code{sd.correct}.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param NON Boolean value (TRUE or FALSE) indicating if the number of species in common should only be counted for the nearest neighbours.
#' @param FPO A vector with the coordinates of the fixed point origin from which the zeta diversity will be computed (overrides NON). In that case, \eqn{\zeta_1} is the number of species in the closest site to the FPO, \eqn{\zeta_2} is the number of species shared by the 2 closest sites, etc.
#' @param DIR Boolean value (TRUE or FALSE) indicating if zeta diversity must be computed using a directed nearest neighbour scheme in the direction away from the FPO, starting from any site.
#' @param empty.row Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity, and computation errors for the normalized version of zeta due to divisions by 0. Options are "\code{empty}" to let the data untreated, "\code{remove}" to remove the empty rows, 0 to set the normalized zeta to 0 when zeta is divided by 0 during normalization (sites share no species, so are completely dissimilar), and 1 to set the normalized zeta to 1 when zeta is divided by 0 during normalization (i.e. sites are perfectly similar).
#' @param silent  Boolean value (TRUE or FALSE) indicating if messages must be printed.
#' @details If the number of combinations of sites is lower than the value of the parameter \code{sam}, all the combinations are used and an exact solution is computed. In that case, using the number of site combinations as the denominator may be appropriate to compute the standard deviation, if all sites were sampled and the zeta values. This can be adjusted with parameters \code{sd.correct} and \code{sd.correct.adapt}.
#' @details \code{Zeta.order.mc} is faster than \code{\link{Zeta.order.ex}} to compute the exact value of zeta diversity when the number of species is higher than \eqn{C^N_{i}}, where \emph{N} is the total number of sites and \emph{i} is the order of zeta.
#' @details \code{Zeta.order.mc} enables accomodating richness heterogeneity by setting \code{normalize = "Jaccard"}, \code{normalize = "Sorensen"} or \code{normalize = "Simpson"}. This cannot be performed by \cr \code{\link{Zeta.order.ex}}.
#' @return \code{Zeta.order.mc}  returns a list containing the following components:
#' @return \item{zeta.order}{The number of assemblages or sites for which the zeta diversity was computed.}
#' @return \item{combinations}{The number of possible combinations of sites for the chosen order.}
#' @return \item{zeta.val}{The zeta diversity values.}
#' @return \item{zeta.val.sd}{The standard deviation of zeta diversity.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}
#' @examples
#'
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[,1:2]
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' zeta.bird <- Zeta.order.mc(data.spec.bird, order = 3, sam=100)
#' zeta.bird
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' zeta.marion <- Zeta.order.mc(data.spec.marion, xy.marion, order = 3, sam = 100,
#'    NON = TRUE)
#' zeta.marion
#'
#' @export
Zeta.order.mc <- function(data.spec, xy=NULL, order = 1, sam = 1000, sd.correct = TRUE, sd.correct.adapt = FALSE, rescale = FALSE, normalize = FALSE, NON = FALSE, FPO = NULL, DIR = FALSE, empty.row = "empty", silent=TRUE){

  if (!inherits(data.spec,"data.frame")){
    stop("Error: ",paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  if((NON == TRUE || !is.null(FPO)) & is.null(xy)){
    stop("Error: if NON = TRUE or !is.null(FPO), xy must be non null.")
  }
  if((NON == TRUE || !is.null(FPO)) && nrow(data.spec) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }

  if(empty.row == "remove"){
    if(length(which(rowSums(data.spec)))>0){
      data.spec <- data.spec[-which(rowSums(data.spec)==0),]
    }
  }

  x <- dim(data.spec)[1]

  if (is.null(FPO)){
    if(order==1){
      zeta.val<-mean(rowSums(data.spec))
      if(sd.correct == TRUE & sd.correct.adapt == FALSE){
        zeta.val.sd <- stats::sd(rowSums(data.spec))
      }else{
        zeta.val.sd <- stats::sd(rowSums(data.spec))*nrow(data.spec)
      }
      if(rescale == TRUE || normalize != FALSE){
        zeta.val <- 1
        zeta.val.sd <- zeta.val.sd/mean(rowSums(data.spec))
      }
    }else{
      if(NON == FALSE){
        if(choose(x, order)>sam){
          if(silent==FALSE){
            print(paste("Monte Carlo sampling for order",order))
          }
          u <- rep(NA, sam)
          for(z in 1:sam){
            samp <- sample(1:x, order, replace = FALSE)
            u[z] <- sum(apply(data.spec[samp, ], 2, prod))
            if (normalize == "Jaccard"){
              toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Sorensen"){
              toto <- (mean(apply(data.spec[samp, ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Simpson"){
              toto <- (min(apply(data.spec[samp, ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }
          }
        }else{
          if(silent==FALSE){
            print(paste("Exact solution for order",order))
          }
          u <- rep(NA, choose(x, order))
          samp <- utils::combn(1:x, order)
          for(z in 1:dim(samp)[2]){
            u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
            if (normalize == "Jaccard"){
              toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp[, z], ]), 2, prod)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Sorensen"){
              toto <- (mean(apply(data.spec[samp[, z], ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Simpson"){
              toto <- (min(apply(data.spec[samp[, z], ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }
          }
        }
      }else{
        if(x>sam){
          u <- rep(NA, sam)
          samps <- sample(1:x, sam, replace = FALSE)
          for(z in 1:sam){
            samp <- samps[z]
            xy.dist <- (xy[,1]-xy[samp,1])^2+(xy[,2]-xy[samp,2])^2
            samp <- c(samp,order(xy.dist)[2:order])
            u[z] <- sum(apply(data.spec[samp, ], 2, prod))

            if (normalize == "Jaccard"){
              toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Sorensen"){
              toto <- (mean(apply(data.spec[samp, ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Simpson"){
              toto <- (min(apply(data.spec[samp, ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }
          }
        }else{
          u <- rep(NA, x)
          samps <- 1:x
          for(z in 1:x){
            samp <- samps[z]
            xy.dist <- (xy[,1]-xy[samp,1])^2+(xy[,2]-xy[samp,2])^2
            samp <- c(samp,order(xy.dist)[2:order])
            u[z] <- sum(apply(data.spec[samp, ], 2, prod))
            if (normalize == "Jaccard"){
              toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Sorensen"){
              toto <- (mean(apply(data.spec[samp, ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }else if (normalize == "Simpson"){
              toto <- (min(apply(data.spec[samp, ], 1, sum)))
              if(toto==0){
                if(empty.row == 0){
                  u[z] <- 0
                }else if(empty.row == 1){
                  u[z] <- 1
                }
              }else
                u[z] <- u[z] / toto
            }
          }
        }
      }

      if(rescale == TRUE & normalize == FALSE){
        #z1 <- mean(rowSums(data.spec))
        #u <- u / z1
        u <- u / ncol(data.spec)
      }

      zeta.val <- mean(u)
      if(sd.correct.adapt == FALSE){
        if(sd.correct == TRUE){
          zeta.val.sd <- stats::sd(u)
        }else{
          zeta.val.sd <- stats::sd(u)*sqrt((length(u)-1)/length(u))
        }
      }else{
        if(x>sam){
          zeta.val.sd <- stats::sd(u)*sqrt((length(u)-1)/length(u))
        }else{
          zeta.val.sd <- stats::sd(u)
        }
      }

    }
  }else{
    if(DIR == FALSE){
      xy.dist <- (FPO[1]-xy[,1])^2+(FPO[2]-xy[,2])^2
      samp <- order(xy.dist)[1:order]
      u <- sum(apply(data.spec[samp, ], 2, prod))
      if (normalize == "Jaccard"){
        toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
        if(toto==0){
          if(empty.row == 0){
            u <- 0
          }else if(empty.row == 1){
            u <- 1
          }
        }else
          u <- u / toto
      }else if (normalize == "Sorensen"){
        toto <- (mean(apply(data.spec[samp, ], 1, sum)))
        if(toto==0){
          if(empty.row == 0){
            u <- 0
          }else if(empty.row == 1){
            u <- 1
          }
        }else
          u <- u / toto
      }else if (normalize == "Simpson"){
        toto <- (min(apply(data.spec[samp, ], 1, sum)))
        if(toto==0){
          if(empty.row == 0){
            u <- 0
          }else if(empty.row == 1){
            u <- 1
          }
        }else
          u <- u / toto
      }
      #zeta.val.sd <- 0
    }else{
      if(order==1){
        zeta.val<-mean(rowSums(data.spec))
        if(sd.correct == TRUE & sd.correct.adapt == FALSE){
          zeta.val.sd <- stats::sd(rowSums(data.spec))
        }else{
          zeta.val.sd <- stats::sd(rowSums(data.spec))*nrow(data.spec)
        }
        if(rescale == TRUE || normalize != FALSE){
          zeta.val <- 1
          zeta.val.sd <- zeta.val.sd/mean(rowSums(data.spec))
        }
      }else{
        xy.FPO <- as.matrix(xy-FPO)
        if(x>sam){
          u <- rep(NA, sam)
          samps <- sample(1:x, sam, replace = FALSE)
          for(z in 1:sam){
            samp <- samps[z]

            xy0 <- xy.FPO[samp,]
            no <- sqrt(sum(xy0^2))
            R <- matrix(c(xy0[2]/no,xy0[1]/no,-xy0[1]/no,xy0[2]/no),2,2)
            xy.FPO.tr <- apply(xy.FPO,1,function(xy.FPO,R,xy0,no){R %*% matrix(xy.FPO,2,1) - c(0,no)},R,xy0,no)
            xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
            xy.FPO.tr[samp,] <- c(0,0)
            xy.FPO.tr[which(xy.FPO.tr[,2]<=0),] <- NA

            xy.dist <- xy.FPO.tr[,1]^2+xy.FPO.tr[,2]^2

            if(length(which(!is.na(xy.dist)))>=(order-1)){
              samp <- c(samp,order(xy.dist)[1:(order-1)])
              u[z] <- sum(apply(data.spec[samp, ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }
            }else{
              if(silent==FALSE){
                print(paste("warning: number of sites away from the FPO too low to compute zeta for order",order))
              }
              u[z] <- NA
            }
          }
        }else{
          u <- rep(NA, x)
          samps <- 1:x
          for(z in 1:(x-order+1)){
            samp <- samps[z]

            xy0 <- xy.FPO[samp,]
            no <- sqrt(sum(xy0^2))
            R <- matrix(c(xy0[2]/no,xy0[1]/no,-xy0[1]/no,xy0[2]/no),2,2)
            xy.FPO.tr <- apply(xy.FPO,1,function(xy.FPO,R,xy0,no){R %*% matrix(xy.FPO,2,1) - c(0,no)},R,xy0,no)
            xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
            xy.FPO.tr[samp,] <- c(0,0)
            xy.FPO.tr[which(xy.FPO.tr[,2]<=0),] <- NA

            xy.dist <- xy.FPO.tr[,1]^2+xy.FPO.tr[,2]^2

            if(length(which(!is.na(xy.dist)))>=(order-1)){
              samp <- c(samp,order(xy.dist)[1:(order-1)])
              u[z] <- sum(apply(data.spec[samp, ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[z] <- 0
                  }else if(empty.row == 1){
                    u[z] <- 1
                  }
                }else
                  u[z] <- u[z] / toto
              }
            }else{
              if(silent==FALSE){
                print(paste("warning: number of sites away from the FPO too low to compute zeta for order",order))
              }
              u[z] <- NA
            }
          }
        }
      }
    }

    zeta.val <- mean(u,na.rm=TRUE)
    if(sd.correct.adapt == FALSE){
      if(sd.correct == TRUE){
        zeta.val.sd <- stats::sd(u,na.rm=TRUE)
      }else{
        zeta.val.sd <- stats::sd(u,na.rm=TRUE)*sqrt((length(which(!is.na(u)))-1)/length(which(!is.na(u))))
      }
    }else{
      if(x>sam){
        zeta.val.sd <- stats::sd(u,na.rm=TRUE)*sqrt((length(which(!is.na(u)))-1)/length(which(!is.na(u))))
      }else{
        zeta.val.sd <- stats::sd(u,na.rm=TRUE)
      }
    }
  }

  zeta.order <- list()
  zeta.order$zeta.order <- order
  zeta.order$combinations <- choose(x, order)
  zeta.order$zeta.val <- zeta.val
  zeta.order$zeta.val.sd <- zeta.val.sd

  return(zeta.order)

}




#' Number of species in common between a specific number of assemblages or sites using Monte Carlo sampling, for multiple combinations and several groups of taxa
#'
#' Computes the number of species shared by multiple assemblages, for a specific order (number of assemblages or sites), for multiple combinations and several groups of taxa.
#' @param data.spec  A list of site-by-species presence-absence data frames, with sites as rows and species as columns.
#' @param xy Site coordinates. This is only used if \code{NON} = TRUE or \code{DIR} = TRUE.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta diversity is computed.
#' @param sd.correct Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) or not (using the number of site combinations as the denominator).
#' @param sd.correct.adapt Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) if \code{sam} is higher than the number of possible combinations, or not (using the number of site combinations as the denominator) if \code{sam} is lower than the number of possible combinations. If \code{sd.correct.adapt == TRUE}, it takes precedence over \code{sd.correct}.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param NON Boolean value (TRUE or FALSE) indicating if the number of species in common should only be counted for the nearest neighbours.
#' @param FPO A vector with the coordinates of the fixed point origin from which the zeta diversity will be computed (overrides NON). In that case, \eqn{\zeta_1} is the number of species in the closest site to the FPO, \eqn{\zeta_2} is the number of species shared by the 2 closest sites, etc.
#' @param DIR Boolean value (TRUE or FALSE) indicating if zeta diversity must be computed using a directed nearest neighbour scheme in the direction away from the FPO, starting from any site.
#' @param empty.row Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity, and computation errors for the normalized version of zeta due to divisions by 0. Options are "\code{empty}" to let the data untreated, "\code{remove}" to remove the empty rows, 0 to set the normalized zeta to 0 when zeta is divided by 0 during normalization (sites share no species, so are completely dissimilar), and 1 to set the normalized zeta to 1 when zeta is divided by 0 during normalization (i.e. sites are perfectly similar).
#' @param silent  Boolean value (TRUE or FALSE) indicating if messages must be printed.
#' @details Contrary to \code{Zeta.order.mc}, the number of species shared by the different combinations of assemblages are not averaged, but returned as is. This is useful to then compare local zeta diversity for different groups of taxa.
#' @details As for \code{Zeta.order.mc}, if the number of combinations of sites is lower than the value of the parameter \code{sam}, all the combinations are used and an exact solution is computed. In that case, using the number of site combinations as the denominator may be appropriate to compute the standard deviation, if all sites were sampled and the zeta values. This can be adjusted with parameters \code{sd.correct} and \code{sd.correct.adapt}.
#' @return \code{Zeta.order.mc.mult}  returns a list containing the following components:
#' @return \item{zeta.order}{The number of assemblages or sites for which the zeta diversity was computed.}
#' #' @return \item{sites}{A matrix in which each row contains the indices of a given combination, i.e. of the specific  \code{sam} assemblages.}
#' @return \item{zeta.val}{A data frame in which each column is the number of species shared by the assemblages.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}
#' @examples
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[1:2]
#' 
#' data.spec.marion <- Marion.species[3:33]
#' 
#' ##random other communities
#' data.spec.marion2a <- data.spec.marion
#' data.spec.marion2a[which(data.spec.marion2a==1,arr.ind=TRUE)] <- 0
#' for(i in 1:ncol(data.spec.marion2a))
#'   data.spec.marion2a[sample(nrow(data.spec.marion2a),8),i] <- 1

#' data.spec.marion2b <- data.spec.marion
#' data.spec.marion2b[which(data.spec.marion2b==1,arr.ind=TRUE)] <- 0
#' for(i in 1:ncol(data.spec.marion2b))
#' data.spec.marion2b[sample(nrow(data.spec.marion2b),8),i] <- 1
#' 
#' dat.spec.tot <- list(data.spec.marion,data.spec.marion2a,data.spec.marion2b)
#' zeta.tot <- Zeta.order.mc.mult(data.spec=dat.spec.tot,order=3,sam=200)
#'
#' @export
Zeta.order.mc.mult <- function(data.spec, xy=NULL, order = 1, sam = 1000, sd.correct = TRUE, sd.correct.adapt = FALSE, rescale = FALSE, normalize = FALSE, NON = FALSE, FPO = NULL, DIR = FALSE, empty.row = "empty", silent=TRUE){
  
  if(!inherits(data.spec,"list")){
    stop("Error: ",paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a list of ata frames.", sep = ""))
  }
  
  if((NON == TRUE || !is.null(FPO)) & is.null(xy)){
    stop("Error: if NON = TRUE or !is.null(FPO), xy must be non null.")
  }
  if((NON == TRUE || !is.null(FPO)) && nrow(data.spec[[1]]) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }
  
  x <- dim(data.spec[[1]])[1]
  
  zeta.val <- list()
  u <- list()
  #sites <- matrix(NA,sam,order)
  
  if (is.null(FPO)){
    if(order==1){
      for(i in 1:length(data.spec))
        zeta.val[[i]]<-(rowSums(data.spec[[i]]))
      sites <- matrix(1:length(data.spec),length(data.spec),1)
    }else{
      if(NON == FALSE){
        if(choose(x, order)>sam){
          if(silent==FALSE){
            print(paste("Monte Carlo sampling for order",order))
          }
          for(i in 1:length(data.spec)){
            u[[i]] <- rep(NA, sam)
          }
          sites <- matrix(NA,sam,order)
          for(z in 1:sam){
            samp <- sample(1:x, order, replace = FALSE)
            sites[z,] <- samp
            for(i in 1:length(data.spec)){
              u[[i]][z] <- sum(apply(data.spec[[i]][samp, ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec[[i]])-sum(apply((1-data.spec[[i]][samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[[i]][samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[[i]][samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }
            }
          }
        }else{
          if(silent==FALSE){
            print(paste("Exact solution for order",order))
          }
          u <- list()
          for(i in 1:length(data.spec)){
            u[[i]] <- rep(NA, choose(x, order))
          }
          samp <- utils::combn(1:x, order)
          sites <- t(samp)
          for(i in 1:length(data.spec)){
            for(z in 1:dim(samp)[2]){
              u[[i]][z] <- sum(apply(data.spec[[i]][samp[, z], ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec[[i]])-sum(apply((1-data.spec[[i]][samp[, z], ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[[i]][samp[, z], ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[[i]][samp[, z], ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }
            }
          }
        }
      }else{
        if(x>sam){
          u <- list()
          for(i in 1:length(data.spec)){
            u[[i]] <- rep(NA, sam)
          }
          samps <- sample(1:x, sam, replace = FALSE)
          sites <- matrix(NA,sam,order)
          for(z in 1:sam){
            samp <- samps[z]
            xy.dist <- (xy[,1]-xy[samp,1])^2+(xy[,2]-xy[samp,2])^2
            samp <- c(samp,order(xy.dist)[2:order])
            sites[z,] <- samp
            for(i in 1:length(data.spec)){
              u[[i]][z] <- sum(apply(data.spec[[i]][samp, ], 2, prod))
              
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec[[i]])-sum(apply((1-data.spec[[i]][samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[[i]][samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[[i]][samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }
            }
          }
        }else{
          u <- list()
          for(i in 1:length(data.spec)){
            u[[i]] <- rep(NA, x)
          }
          samps <- 1:x
          sites <- matrix(NA,x,order)
          for(z in 1:x){
            samp <- samps[z]
            xy.dist <- (xy[,1]-xy[samp,1])^2+(xy[,2]-xy[samp,2])^2
            samp <- c(samp,order(xy.dist)[2:order])
            sites[z,] <- samp
            for(i in 1:length(data.spec)){
              u[[i]][z] <- sum(apply(data.spec[[i]][samp, ], 2, prod))
              if (normalize == "Jaccard"){
                toto <- (ncol(data.spec[[i]])-sum(apply((1-data.spec[[i]][samp, ]), 2, prod)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Sorensen"){
                toto <- (mean(apply(data.spec[[i]][samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }else if (normalize == "Simpson"){
                toto <- (min(apply(data.spec[[i]][samp, ], 1, sum)))
                if(toto==0){
                  if(empty.row == 0){
                    u[[i]][z] <- 0
                  }else if(empty.row == 1){
                    u[[i]][z] <- 1
                  }
                }else
                  u[[i]][z] <- u[[i]][z] / toto
              }
            }
          }
        }
      }
      
      if(rescale == TRUE & normalize == FALSE){
        for(i in 1:length(data.spec)){
          u[[i]] <- u[[i]] / ncol(data.spec[[i]][[1]])
        }
      }
      
      zeta.val <- u
      
    }
  }else{
    if(DIR == FALSE){
      xy.dist <- (FPO[1]-xy[,1])^2+(FPO[2]-xy[,2])^2
      samp <- order(xy.dist)[1:order]
      sites=samp
      for(i in 1:length(data.spec)){
        zeta.val[[i]] <- sum(apply(data.spec[[i]][samp, ], 2, prod))
        if (normalize == "Jaccard"){
          toto <- (ncol(data.spec[[i]])-sum(apply((1-data.spec[[i]][samp, ]), 2, prod)))
          if(toto==0){
            if(empty.row == 0){
              zeta.val[[i]] <- 0
            }else if(empty.row == 1){
              zeta.val[[i]] <- 1
            }
          }else
            zeta.val[[i]] <- zeta.val[[i]] / toto
        }else if (normalize == "Sorensen"){
          toto <- (mean(apply(data.spec[[i]][samp, ], 1, sum)))
          if(toto==0){
            if(empty.row == 0){
              zeta.val[[i]] <- 0
            }else if(empty.row == 1){
              zeta.val[[i]] <- 1
            }
          }else
            zeta.val[[i]] <- zeta.val[[i]] / toto
        }else if (normalize == "Simpson"){
          toto <- (min(apply(data.spec[[i]][samp, ], 1, sum)))
          if(toto==0){
            if(empty.row == 0){
              zeta.val[[i]] <- 0
            }else if(empty.row == 1){
              zeta.val[[i]] <- 1
            }
          }else
            zeta.val[[i]] <- zeta.val[[i]] / toto
        }
      }
    }else{
      if(order==1){
        for(i in 1:length(data.spec)){
          zeta.val[[i]]<-rowSums(data.spec[[i]])
        }
      }else{
        xy.FPO <- as.matrix(xy-FPO)
        if(x>sam){
          for(i in 1:length(data.spec)){
            u[[i]] <- rep(NA, sam)
          }
          samps <- sample(1:x, sam, replace = FALSE)
          sites <- matrix(NA,sam,order)
          for(z in 1:sam){
            samp <- samps[z]
            xy0 <- xy.FPO[samp,]
            no <- sqrt(sum(xy0^2))
            R <- matrix(c(xy0[2]/no,xy0[1]/no,-xy0[1]/no,xy0[2]/no),2,2)
            xy.FPO.tr <- apply(xy.FPO,1,function(xy.FPO,R,xy0,no){R %*% matrix(xy.FPO,2,1) - c(0,no)},R,xy0,no)
            xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
            xy.FPO.tr[samp,] <- c(0,0)
            xy.FPO.tr[which(xy.FPO.tr[,2]<=0),] <- NA
            
            xy.dist <- xy.FPO.tr[,1]^2+xy.FPO.tr[,2]^2
            
            if(length(which(!is.na(xy.dist)))>=(order-1)){
              samp <- c(samp,order(xy.dist)[1:(order-1)])
              sites[z,] <- samp
              for(i in 1:length(data.spec)){
                u[[i]][z] <- sum(apply(data.spec[[i]][samp, ], 2, prod))
                if (normalize == "Jaccard"){
                  toto <- (ncol(data.spec[[i]])-sum(apply((1-data.spec[[i]][samp, ]), 2, prod)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[[i]][z] <- 0
                    }else if(empty.row == 1){
                      u[[i]][z] <- 1
                    }
                  }else
                    u[[i]][z] <- u[[i]][z] / toto
                }else if (normalize == "Sorensen"){
                  toto <- (mean(apply(data.spec[[i]][samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[[i]][z] <- 0
                    }else if(empty.row == 1){
                      u[[i]][z] <- 1
                    }
                  }else
                    u[[i]][z] <- u[[i]][z] / toto
                }else if (normalize == "Simpson"){
                  toto <- (min(apply(data.spec[[i]][samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[[i]][z] <- 0
                    }else if(empty.row == 1){
                      u[[i]][z] <- 1
                    }
                  }else
                    u[[i]][z] <- u[[i]][z] / toto
                }
              }
            }else{
              if(silent==FALSE){
                print(paste("warning: number of sites away from the FPO too low to compute zeta for order",order))
              }
              for(i in 1:length(data.spec)){
                u[[i]][z] <- NA
              }
            }
          }
        }else{
          for(i in 1:length(data.spec)){
            u[[i]] <- rep(NA, x)
          }
          samps <- 1:x
          sites <- matrix(NA,x,order)
          for(z in 1:(x-order+1)){
            samp <- samps[z]
            
            xy0 <- xy.FPO[samp,]
            no <- sqrt(sum(xy0^2))
            R <- matrix(c(xy0[2]/no,xy0[1]/no,-xy0[1]/no,xy0[2]/no),2,2)
            xy.FPO.tr <- apply(xy.FPO,1,function(xy.FPO,R,xy0,no){R %*% matrix(xy.FPO,2,1) - c(0,no)},R,xy0,no)
            xy.FPO.tr <- as.data.frame(t(xy.FPO.tr))
            xy.FPO.tr[samp,] <- c(0,0)
            xy.FPO.tr[which(xy.FPO.tr[,2]<=0),] <- NA
            
            xy.dist <- xy.FPO.tr[,1]^2+xy.FPO.tr[,2]^2
            
            if(length(which(!is.na(xy.dist)))>=(order-1)){
              samp <- c(samp,order(xy.dist)[1:(order-1)])
              sites[z,] <- samp
              for(i in 1:length(data.spec)){
                u[[i]][z] <- sum(apply(data.spec[[i]][samp, ], 2, prod))
                if (normalize == "Jaccard"){
                  toto <- (ncol(data.spec[[i]])-sum(apply((1-data.spec[[i]][samp, ]), 2, prod)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[[i]][z] <- 0
                    }else if(empty.row == 1){
                      u[[i]][z] <- 1
                    }
                  }else
                    u[[i]][z] <- u[[i]][z] / toto
                }else if (normalize == "Sorensen"){
                  toto <- (mean(apply(data.spec[[i]][samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[[i]][z] <- 0
                    }else if(empty.row == 1){
                      u[[i]][z] <- 1
                    }
                  }else
                    u[[i]][z] <- u[[i]][z] / toto
                }else if (normalize == "Simpson"){
                  toto <- (min(apply(data.spec[[i]][samp, ], 1, sum)))
                  if(toto==0){
                    if(empty.row == 0){
                      u[[i]][z] <- 0
                    }else if(empty.row == 1){
                      u[[i]][z] <- 1
                    }
                  }else
                    u[[i]][z] <- u[[i]][z] / toto
                }
              }
            }else{
              if(silent==FALSE){
                print(paste("warning: number of sites away from the FPO too low to compute zeta for order",order))
              }
              for(i in 1:length(data.spec)){
                u[[i]][z] <- NA
              }
            }
          }
        }
      }
      zeta.val <- u
    }
    
    
    
  }
  
  zeta.order.mult <- list()
  zeta.order.mult$zeta.order <- order
  zeta.order.mult$sites <- sites
  zeta.order.mult$zeta.val <- data.frame(t(matrix(unlist(zeta.val), nrow=length(zeta.val), byrow=T)))
  
  return(zeta.order.mult)
}





#' Expectation of zeta diversity decline
#'
#' Computes the expectation of zeta diversity, the number of species shared by multiple assemblages for a range of orders (number of assemblages or sites), using a formula based on the occupancy of the species, and fits the decline to an exponential and a power law relationship.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param orders  Range of number of assemblages or sites for which zeta diversity is computed.
#' @param sd.correct Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) or not (using the number of site combinations as the denominator).
#' @param confint.level  Percentage for the confidence intervals of the coefficients from the regressions.
#' @param sd.plot  Boolean value (TRUE or FALSE) indicating if the standard deviation of each zeta diversity value must be plotted.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1.
#' @param empty.row Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity. Options are  "\code{empty}" to let the data untreated or "\code{remove}" to remove the empty rows.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @details \code{Zeta.decline.ex} is much faster than \code{\link{Zeta.decline.mc}} to compute the exact value of zeta diversity when the number of species is lower than \eqn{C^N_{i}}, where \emph{N} is the total number of sites and \emph{i} is the order of zeta.
#' @details \code{sd.correct} should be set to \code{TRUE} if the assemblages represent a subsample of the whole system. It can be set to \code{FALSE} if the sampling is exhaustive, for example in case of a continuous regular grid covering the whole study area.
#' @details The exponential and the power law fit are performed using linear regressions on log-transformed data (only the zeta values are log-transformed for the exponential fit, and both the orders and the zeta values are log-transformed for the power law fit).
#' @return \code{Zeta.decline.ex} returns a list containing the following components:
#' @return \item{zeta.order}{The number of assemblages or sites for which the zeta diversity was computed.}
#' @return \item{combinations}{The number of possible combinations of sites for the chosen orders.}
#' @return \item{zeta.val}{The zeta diversity values.}
#' @return \item{zeta.val.sd}{The zeta diversity standard deviation values.}
#' @return \item{zeta.ratio}{The ratio of zeta diversity values by the zeta diversity values at the lower order \eqn{\zeta_i / \zeta_{i-1}}.}
#' @return \item{zeta.exp}{Object of class "\code{lm}", containing the output of the exponential regression.}
#' @return \item{zeta.exp.confint}{The confidence intervals of the coefficients of the exponential regression.}
#' @return \item{zeta.pl}{Object of class "\code{lm}", containing the output of the power law regression.}
#' @return \item{zeta.pl.confint}{The confidence intervals of the coefficients of the power law regression.}
#' @return \item{aic}{AIC values for \code{zeta.exp} and \code{zeta.pl}.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references McGeoch M. A., Latombe G., Andrew N. R., Nakagawa S., Nipperess D. A., Roige M., Marzinelli E. M., Campbell A. H., Verges A., Thomas T., Steinberg P. D., Selwood K. E., Henriksen M. V. & Hui C. (2019). Measuring continuous compositional change using decline and decay in zeta diversity. \emph{Ecology}, 100(11), e02832.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.order.ex}}, \code{\link{Plot.zeta.decline}}
#' @examples
#' utils::data(bird.spec.coarse)
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' dev.new(width = 12, height = 4)
#' zeta.bird <- Zeta.decline.ex(data.spec.bird, orders = 1:5)
#' zeta.bird
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' dev.new(width = 12, height = 4)
#' zeta.marion <- Zeta.decline.ex(data.spec.marion, orders = 1:5)
#' zeta.marion
#'
#' @export
#'
Zeta.decline.ex <- function(data.spec, orders = 1:10, sd.correct = TRUE, confint.level = 0.95, sd.plot = TRUE, rescale = FALSE, empty.row = "empty", plot = TRUE){

  if(max(orders)>dim(data.spec)[1]){
    stop("Error: wrong value for \"orders\": the maximum value must be equal or lower than the number of sites.")
  }

  if(empty.row == "remove"){
    data.spec <- data.spec[-which(rowSums(data.spec)==0),]
  }

  data.spec <- as.matrix(data.spec)
  intercept_mat <- t(data.spec) %*% data.spec
  occupancy <- colSums(data.spec)

  zeta.val <- numeric(length=length(orders))
  zeta.val.sd <- numeric(length=length(orders))

  for(i in orders) {
    p <- exp(lchoose(occupancy,i)-lchoose(nrow(data.spec),i))
    zeta.val[i] <- sum(p)
    varmat <- exp(lchoose(intercept_mat,i)-lchoose(nrow(data.spec),i))
    for(j in 1:length(occupancy)) {
      for(k in 1:length(occupancy)) {
        varmat[j,k] <- varmat[j,k] - p[j]*p[k]
      }
    }
    if(sd.correct == TRUE){
      zeta.val.sd[i] <- sqrt(sum(varmat)*choose(nrow(data.spec),i)/(choose(nrow(data.spec),i)-1))
    }else{
      zeta.val.sd[i] <- sqrt(sum(varmat))
    }
  }

  if(rescale == TRUE){
    z1 <- mean(rowSums(data.spec))
    zeta.val <- zeta.val / z1
    zeta.val.sd <- zeta.val.sd / z1
  }


  ##create a single list for output
  zeta <- list()
  zeta$zeta.order <- orders
  zeta$combinations <- choose(x <- dim(data.spec)[1], orders)
  zeta$zeta.val <- zeta.val
  zeta$zeta.val.sd <- zeta.val.sd
  zeta$ratio <- zeta.val[2:length(zeta.val)]/zeta.val[1:(length(zeta.val)-1)]

  ##regression - exponential
  zeta.val.log <- log10(zeta.val)
  zeta.val.log[which(is.infinite(zeta.val.log))] <- NA
  zeta.exp <- stats::lm(zeta.val.log ~ c(orders), na.action = stats::na.omit)
  zeta$zeta.exp <- zeta.exp
  zeta$zeta.exp.confint <- suppressMessages(stats::confint(zeta.exp,level=confint.level))


  ##regression - power law
  zeta.pl <- stats::lm(zeta.val.log ~ log10(c(orders)), na.action = stats::na.omit)
  zeta$zeta.pl <- zeta.pl
  zeta$zeta.pl.confint <- suppressMessages(stats::confint(zeta.pl,level=confint.level))

  zeta$aic <- stats::AIC(zeta$zeta.exp, zeta$zeta.pl)

  ##Plot zeta and regressions
  if(plot == TRUE){
    Plot.zeta.decline(zeta, sd.plot=sd.plot)
  }

  return(zeta)

}





#' Expectation of zeta diversity for a specific number of assemblages or sites
#'
#' Computes the expectation of zeta diversity, the number of species shared by multiple assemblages, for a specific order (number of assemblages or sites) using a formula based on the occupancy of the species.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sd.correct Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) or not (using the number of site combinations as the denominator).
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1.
#' @param empty.row Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity. Options are  "empty" to let the data untreated or "remove" to remove the empty rows.
#' @details \code{Zeta.order.ex} is much faster than \code{\link{Zeta.order.mc}} to compute the exact value of zeta diversity when the number of species is lower than \eqn{C^N_{i}}, where \emph{N} is the total number of sites and \emph{i} is the order of zeta.
#' @details \code{sd.correct} should be set to \code{TRUE} if the assemblages represent a subsample of the whole system. It can be set to \code{FALSE} if the sampling is exhaustive, for example in case of a continuous regular grid covering the whole study area.
#' @return \code{zeta.order.ex}  returns a list containing the following components:
#' @return \item{zeta.order}{The number of assemblages or sites for which the zeta diversity was computed.}
#' @return \item{combinations}{The number of possible combinations of sites for the chosen order.}
#' @return \item{zeta.val}{The zeta diversity values.}
#' @return \item{zeta.val.sd}{The standard deviation of zeta diversity.}
#' @references Hui C. & McGeoch M.A. (2014). zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references McGeoch M. A., Latombe G., Andrew N. R., Nakagawa S., Nipperess D. A., Roige M., Marzinelli E. M., Campbell A. H., Verges A., Thomas T., Steinberg P. D., Selwood K. E., Henriksen M. V. & Hui C. (2019). Measuring continuous compositional change using decline and decay in zeta diversity. \emph{Ecology}, 100(11), e02832.
#' @seealso \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.decline.mc}}
#' @examples
#'
#' utils::data(bird.spec.coarse)
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' zeta.bird <- Zeta.order.ex(data.spec.bird, order = 3)
#' zeta.bird
#'
#' ##########
#'s
#' utils::data(Marion.species)
#' data.spec.marion <- Marion.species[,3:33]
#'
#' zeta.marion <- Zeta.order.ex(data.spec.marion, order = 3)
#' zeta.marion
#'
#' @export
Zeta.order.ex <- function(data.spec, order = 1, sd.correct = TRUE, rescale = FALSE, empty.row="empty"){

  if(empty.row == "remove"){
    data.spec <- data.spec[-which(rowSums(data.spec)==0),]
  }

  if(order>dim(data.spec)[1]){
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }

  data.spec <- as.matrix(data.spec)
  intercept_mat <- t(data.spec) %*% data.spec
  occupancy <- colSums(data.spec)

  p <- exp(lchoose(occupancy,order)-lchoose(nrow(data.spec),order))
  zeta.val <- sum(p)
  varmat <- exp(lchoose(intercept_mat,order)-lchoose(nrow(data.spec),order))
  for(j in 1:length(occupancy)) {
    for(k in 1:length(occupancy)) {
      varmat[j,k] <- varmat[j,k] - p[j]*p[k]
    }
  }
  if(sd.correct == TRUE){
    zeta.val.sd <- sqrt(sum(varmat)*choose(nrow(data.spec),order)/(choose(nrow(data.spec),order)-1))
  }else{
    zeta.val.sd <- sqrt(sum(varmat))
  }

  zeta.order <- list()
  zeta.order$zeta.order <- order
  zeta.order$combinations <- choose(x <- dim(data.spec)[1], order)
  zeta.order$zeta.val <- zeta.val
  zeta.order$zeta.val.sd <- zeta.val.sd

  return(zeta.order)

}




#' Zeta diversity decline plotting
#'
#' Plots the output of the functions \code{Zeta.decline.mc} and \code{Zeta.decline.ex}.
#' @param zeta  A list produced by the function \code{Zeta.decline.mc} or \code{Zeta.decline.ex}.
#' @param sd.plot Boolean value (TRUE or FALSE) indicating if the standard deviation of each zeta diversity value must be plotted.
#' @param arrange.plots Boolean value (TRUE or FALSE) indicating if the graphics device must be divided into 4 subplots.
#' @return A plot of the zeta decline with 4 subplots displaying (i) the raw decline, (ii) the ratios of the zeta values (computed as \eqn{\zeta_i / \zeta_{i-1}}), (iii) the fit in a log plot and (iv) the fit in a log-log plot.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}
#' @examples
#'
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[1:2]
#' data.spec.bird <- bird.spec.coarse[3:193]
#'
#' dev.new(width = 12, height = 4)
#' zeta.bird <- Zeta.decline.mc(data.spec.bird, orders = 1:5, sam=100, plot = FALSE)
#' Plot.zeta.decline(zeta.bird)
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[1:2]
#' data.spec.marion <- Marion.species[3:33]
#'
#' dev.new(width = 12, height = 4)
#' zeta.marion <- Zeta.decline.ex(data.spec.marion, orders = 1:5, plot = FALSE)
#' Plot.zeta.decline(zeta.marion)
#'
#' @export
#'
Plot.zeta.decline <- function(zeta, sd.plot = TRUE, arrange.plots = TRUE){

  if(arrange.plots == TRUE){
    graphics::par(mfrow = c(1, 4))
  }
  if (sd.plot == TRUE){
    graphics::plot(zeta$zeta.order, zeta$zeta.val, xlab = "Zeta order", ylab = "Zeta diversity", pch = 20, ylim = c(0, zeta$zeta.val[1] + zeta$zeta.val.sd[1]), main = "Zeta diversity decline")
    graphics::lines(zeta$zeta.order, zeta$zeta.val)
    ##sd of zeta as error bars
    graphics::lines(zeta$zeta.order,zeta$zeta.val + zeta$zeta.val.sd,lty=2)
    graphics::lines(zeta$zeta.order,zeta$zeta.val - zeta$zeta.val.sd,lty=2)
  }else{
    graphics::plot(zeta$zeta.order, zeta$zeta.val, xlab = "Zeta order", ylab = "Zeta diversity", pch = 20, ylim = c(0, zeta$zeta.val[1]), main = "Zeta diversity decline")
    graphics::lines(zeta$zeta.order, zeta$zeta.val)
  }
  graphics::plot(zeta$zeta.order[1:(length(zeta$zeta.order)-1)],zeta$ratio,pch=20,xlab = "Zeta order", ylab = "Zeta ratio", main = "Ratio of zeta diversity decline")
  graphics::lines(zeta$zeta.order[1:(length(zeta$zeta.order)-1)],zeta$ratio)
  graphics::plot(zeta$zeta.order, zeta$zeta.val, log = "y", pch = 20, xlab = "Zeta order", ylab = "Zeta diversity", main = "Exponential regression")
  graphics::lines(zeta$zeta.order, 10^stats::predict.lm(zeta$zeta.exp, data.frame(zeta$zeta.order)))
  graphics::plot(zeta$zeta.order, zeta$zeta.val, log = "xy", pch = 20, xlab = "Zeta order", ylab = "Zeta diversity", main = "Power law regression")
  graphics::lines(zeta$zeta.order, 10^stats::predict.lm(zeta$zeta.pl, data.frame(zeta$zeta.order)))
}




#' Sensitivity analysis for the sample size of zeta
#'
#' Computes zeta diversity for a given order (number of assemblages or sites) for a range of sample sizes, to assess the sensitivity to this parameter.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param xy Site coordinates. This is only used if \code{NON} = TRUE or \code{DIR} = TRUE.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam.seq Sequence of samples for which the zeta diversity is computed.
#' @param reps  Number of replicates of zeta diversity computations for each sample size.
#' @param sd.correct Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) or not (using the number of site combinations as the denominator).
#' @param sd.correct.adapt Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) if \code{sam} is higher than the number of possible combinations, or not (using the number of site combinations as the denominator) if \code{sam} is lower than the number of possible combinations. If \code{sd.correct.adapt == TRUE}, it takes precedence over \code{sd.correct}.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1.
#' @param NON Boolean value (TRUE or FALSE) indicating if the number of species in common should only be counted for the nearest neighbours.
#' @param FPO A vector with the coordinates of the fixed point origin from which the zeta diversity will be computed (overrides NON). In that case, \eqn{\zeta_1} is the number of species in the closest site to the FPO, \eqn{\zeta_2} is the number of species shared by the 2 closest sites, etc.
#' @param DIR Boolean value (TRUE or FALSE) indicating if zeta diversity must be computed using a directed nearest neighbour scheme in the direction away from the FPO, starting from any site.
#' @param display  Boolean value (TRUE or FALSE) indicating if the current value of the sample size must be displayed. Acts as a counter.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted as a boxplot of the zeta diversity distributions for each sample size
#' @param notch  Boolean value (TRUE or FALSE) indicating if the notches must be plotted in the boxplot.
#' @return \code{Zeta.sam.sensitivity} returns a matrix with \code{(sam.max-sam.min)/sam.incr} columns and \code{reps} rows.
#' @details Note that the execution of \code{Zeta.sam.sensitivity} can be quite lengthy, because of the number of replicates needed.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}
#' @examples
#' \donttest{
#' #Note that the sensitivity analyses in the following two examples are quite long to run,
#' #typically around 10 minutes for the first example and 1-2 minutes for the second.
#'
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[1:2]
#' data.spec.bird <- bird.spec.coarse[3:193]
#'
#' dev.new()
#' zeta.sens.bird <- Zeta.sam.sensitivity(data.spec.bird, order = 3,
#'     sam.seq = seq(250,1000,250), reps = 20, display = TRUE, plot = TRUE, notch = TRUE)
#' zeta.sens.bird
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[1:2]
#' data.spec.marion <- Marion.species[3:33]
#'
#' dev.new()
#' zeta.sens.marion <- Zeta.sam.sensitivity(data.spec.marion, order = 3,
#'     sam.seq = seq(50,250,50), reps = 20, plot = TRUE, notch = TRUE)
#' zeta.sens.marion
#' }
#'
#' @export
Zeta.sam.sensitivity <- function(data.spec, xy = NULL, order = 1, sam.seq, reps = 20, sd.correct = TRUE, sd.correct.adapt = FALSE, rescale = FALSE, normalize = FALSE, NON = FALSE, FPO = NULL, DIR = FALSE, display = TRUE, plot = TRUE, notch = TRUE){

  if (!inherits(data.spec,"data.frame")){
    stop("Error: ",paste(deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }

  x <- dim(data.spec)[1]
  zeta.val <- matrix(NA, reps, length(sam.seq))
  zeta.val.sd <- matrix(NA, reps, length(sam.seq))

  i.sam <- 0
  for(sam in sam.seq){
    if(display == TRUE){print(sam)}
    i.sam <- i.sam + 1
    for(i in 1:reps){
      u <- rep(NA, sam)
      for(j in 1:order){
        zeta <- Zeta.order.mc(data.spec = data.spec, xy = xy, order = order, sam = sam, sd.correct = sd.correct, sd.correct.adapt = sd.correct.adapt, rescale = rescale, normalize = normalize, NON = NON, FPO = FPO, DIR = DIR, silent=TRUE)
        zeta.val[i, i.sam] <- zeta$zeta.val
        zeta.val.sd[i, i.sam] <- zeta$zeta.val.sd
      }
    }
  }

  if (plot == TRUE){
    graphics::boxplot(zeta.val, notch = notch, names = sam.seq, xlab = "number of samples", ylab = paste("zeta ", order, sep = ""), main = "Distributions of zeta diversities for different number of samples")
  }

  zeta.sens <- zeta.val

  return(zeta.sens)

}







#' Multi-site generalised dissimilarity modelling for a set of environmental variables and distances
#'
#' Computes a regression model of zeta diversity for a given order (number of assemblages or sites) against a set of environmental variables and distances between sites. The different regression models available are generalised linear models, generalised linear models with negative constraints, generalised additive models, shape constrained additive models, and I-splines.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param data.env  Site-by-variable data frame, with sites as rows and environmental variables as columns.
#' @param xy Site coordinates, to account for distances between sites.
#' @param data.spec.pred Site-by-species presence-absence data frame or list of data frames, with sites as rows and species as columns, for which zeta diversity will be computed and used as a predictor of the zeta diversity of \code{data.spec}.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta diversity is computed.
#' @param reg.type Type of regression used in the multi-site generalised dissimilarity modelling. Options are "\code{glm}" for generalised linear models, "\code{ngls}" for negative linear models, "\code{gam}" for generalised additive models, "\code{scam}" for shape constrained additive models (with monotonic decreasing by default), and "\code{ispline}" for I-spline models (forcing monotonic decline), as recommended in generalised dissimilarity modelling by Ferrier \emph{et al}. (2007).
#' @param family A description of the error distribution and link function to be used in the \code{glm}, \code{gam} and \code{scam} models (see \code{\link[stats]{family}} for details of family functions).
#' @param method.glm Method used in fitting the generalised linear model. The default method \cr "glm.fit.cons" is an adaptation of method \code{glm.fit2} from package \code{glm2} using a constrained least squares regression (default is negative coefficients) in the reweighted least squares. Another option is "glm.fit2", which calls \code{glm.fit2}; see help documentation for glm.fit2 in package \code{glm2}.
#' @param cons type of constraint in the glm if \code{method.glm = "glm.fit.cons"}. Default is -1 for negative coefficients on the predictors. The other option is 1 for positive coefficients on the predictors.
#' @param cons.inter type of constraint for the intercept. Default is 1 for positive intercept, suitable for Gaussian family. The other option is -1 for negative intercept, suitable for binomial family.
#' @param confint.level  Percentage for the confidence intervals of the coefficients from the generalised linear models.
#' @param bs A two-letter character string indicating the (penalized) smoothing basis to use in the scam model. Default is "\code{mpd}" for monotonic decreasing splines. see \code{\link[mgcv]{smooth.terms}} for an overview of what is available.
#' @param kn Number of knots in the GAM and SCAM. Default is -1 for determining kn automatically using Generalized Cross-validation.
#' @param order.ispline Order of the I-spline.
#' @param kn.ispline Number of knots in the I-spline.
#' @param distance.type Method to compute distance. Default is "\code{Euclidean}", for Euclidean distance. The other options are (i) "\code{ortho}" for orthodromic distance, if xy correspond to longitudes and latitudes (orthodromic distance is computed using the \code{geodist} function from package \code{geodist}); and (ii) "\code{custom}", in which case the user must provide a distance matrix for \code{dist.custom}.
#' @param dist.custom Distance matrix provided by the user when \code{distance.type} = \code{"custom"}.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by the total number of species in the dataset, to get a range of values between 0 and 1. Has no effect if \code{normalize} != \code{FALSE}.
#' @param rescale.pred  Boolean value (TRUE or FALSE) indicating if the spatial distances and differences in environmental variables should be rescaled between 0 and 1.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param method  Name of a function (as a string) indicating how to combine the pairwise differences and distances for more than 3 sites. It can be a basic R-function such as "\code{mean}" or "\code{max}", but also a custom function.
#' @param silent  Boolean value (TRUE or FALSE) indicating if warnings must be printed.
#' @param empty.row Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity, and computation errors for the normalized version of zeta due to divisions by 0. Options are "\code{empty}" to let the data untreated, "\code{remove}" to remove the empty rows, 0 to set the normalized zeta to 0 when zeta is divided by 0 during normalization (sites share no species, so are completely dissimilar), and 1 to set the normalized zeta to 1 when zeta is divided by 0 during normalization (i.e. sites are perfectly similar).
#' @param control As for \code{\link{glm}}.
#' @param glm.init Boolean value, indicating if the initial parameters for fitting the glm with constraint on the coefficients signs for \code{reg.type == "ispline"} should be initialised based on the correlation coefficients betwen the zeta values and the environmental difference or distance. \code{glm.init = TRUE} helps preventing the error message: \code{error: cannot find valid starting values:} \cr \code{please specify some}.
#' @return \code{Zeta.msgdm} returns a list whose component vary depending on the regression technique. The list can contain the following components:
#' @return \item{val}{Vector of zeta values used in the MS-GDM.}
#' @return \item{predictors}{Data frame of the predictors used in the MS-GDM.}
#' @return \item{range.min}{Vector containing the minimum values of the numeric variables, used for rescaling the variables between 0 and 1 for I-splines (see Details).}
#' @return \item{range.max}{Vector containing the maximum values of the numeric variables, used for rescaling the variables between 0 and 1 for I-splines (see Details).}
#' @return \item{rescale.factor}{Factor by which the predictors were divided if \code{rescale.pred = TRUE} and \code{order>1}.}
#' @return \item{order.ispline}{The value of the original parameter, to be used in \code{Plot.ispline}.}
#' @return \item{kn.ispline}{The value of the original parameter, to be used in \code{Plot.ispline}.}
#' @return \item{model}{An object whose class depends on the type of regression (\code{glm}, \code{nnnpls}, \code{gam} or \code{scam}; I-splines return and object of class \code{glm}), corresponding to the regression over distance for the number of assemblages or sites specified in \code{order}.}
#' @return \item{confint}{The confidence intervals for the coefficients from generalised linear models with no constraint. \code{confint} is not generated for the other types of regression.}
#' @return \item{vif}{The variance inflation factors for all the variables for the generalised linear regression. \code{vif} is not generated for the other types of regression.}
#' @details The environmental variables can be numeric or factorial.
#' @details If \code{order = 1}, the variables are used as such in the regression, and factorial variables must be dummy for the output of the regression to be interpretable.
#' @details For numeric variables, if \code{order>1} the pairwise difference between sites is computed and combined according to \code{method}. For factorial variables, the distance corresponds to the number of unique values over the number of assemblages of sites specified by \code{order}.
#' @details If \code{xy = NULL}, \code{Zeta.msgdm} only uses environmental variables in the regression. Otherwise, it also computes and uses euclidian distance (average or maximum distance between multiple sites, depending on the parameters \code{method}) as an explanatory variable.
#' @details If \code{rescale.pred = TRUE}, zeta is regressed against the differences of values of the environmental variables divided by the maximum difference for each variable, to be rescaled between 0 and 1. If \code{!is.null(xy)}, distances between sites are also divided by the maximum distance. If \code{order = 1}, the variables are transformed by first subtracting their minimum value, and dividing by the difference of their maximum and minimum values.
#' @details If \code{reg.type = "ispline"}, the variables are rescaled between 0 and 1 prior to computing the I-splines by subtracting their minimum value, and dividing by the difference of their maximum and minimum values.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Ferrier, S., Manion, G., Elith, J., & Richardson, K. (2007). Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity and Distributions}, 13(3), 252-264.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Predict.msgdm}},
#' @seealso \code{\link{Ispline}}
#' @import scam
#' @examples
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[1:2]
#' data.spec.bird <- bird.spec.coarse[3:193]
#' utils::data(bird.env.coarse)
#' data.env.bird <- bird.env.coarse[,3:9]
#'
#' zeta.glm <- Zeta.msgdm(data.spec.bird, data.env.bird, sam = 100, order = 3)
#' zeta.glm
#' dev.new()
#' graphics::plot(zeta.glm$model)
#'
#' zeta.ngls <- Zeta.msgdm(data.spec.bird, data.env.bird, xy.bird, sam = 100, order = 3,
#'     reg.type = "ngls", rescale = TRUE)
#' zeta.ngls
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[1:2]
#' data.spec.marion <- Marion.species[3:33]
#' utils::data(Marion.env)
#' data.env.marion <- Marion.env[3]
#'
#' zeta.gam <- Zeta.msgdm(data.spec.marion, data.env.marion, sam = 100, order = 3,
#'     reg.type = "gam")
#' zeta.gam
#' dev.new()
#' graphics::plot(zeta.gam$model)
#'
#' zeta.ispline <- Zeta.msgdm(data.spec.marion, data.env.marion, xy.marion, sam = 100,
#'     order = 3, normalize = "Jaccard", reg.type = "ispline")
#' zeta.ispline
#' 
#' zeta.ispline.r <- Return.ispline(zeta.ispline, data.env.marion, distance = TRUE)
#' zeta.ispline.r
#' 
#' dev.new()
#' Plot.ispline(isplines = zeta.ispline.r, distance = TRUE)
#' 
#' dev.new()
#' Plot.ispline(msgdm = zeta.ispline, data.env = data.env.marion, distance = TRUE)
#'
#' @export
Zeta.msgdm <- function (data.spec, data.env, xy = NULL, data.spec.pred = NULL, order = 1, sam = 1000, reg.type = "glm", family = stats::gaussian(), method.glm = "glm.fit.cons", cons = -1, cons.inter = 1, confint.level = 0.95, bs = "mpd", kn = -1, order.ispline = 2, kn.ispline = 1, distance.type = "Euclidean", dist.custom = NULL, rescale = FALSE, rescale.pred = TRUE, method = "mean", normalize = FALSE, silent = FALSE, empty.row = 0, control = list(), glm.init = FALSE) {
  
  
  if (nrow(data.spec) != nrow(data.env)) {
    stop("Error: data.spec and data.env must have the same number of rows.")
  }
  if (!is.null(xy)) {
    if (nrow(data.spec) != nrow(xy) || nrow(data.env) != 
        nrow(xy)) {
      stop("Error: data.spec, data.env and xy must have the same number of rows.")
    }
  }
  if (!inherits(data.spec,"data.frame")) {
    stop(paste("Error: ", deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if (!inherits(data.env,"data.frame")) {
    stop(paste("Error: ", deparse(substitute(data.env)), " is a ", class(data.env), ". It must be a data frame.", sep = ""))
  }
  if (order > dim(data.spec)[1]) {
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  if (length(which(sapply(data.env, inherits, c("factor", "numeric"))==0)) > 0) {
    stop("Error: variables must be numeric or factor")
  }
  if (order == 1 & (!is.null(xy) | distance.type == "custom")) {
    stop("Error: cannot include distance for order = 1")
  }
  if (silent == FALSE & order == 1 & sum(sapply(data.env, inherits, "factor")) > 0) {
    warning("factor variables should be dummy for order = 1")
  }
  if (silent == FALSE & !is.null(dist.custom)) {
    if (!isSymmetric(dist.custom)) {
      warning("Distance matrix is not symmetrical")
    }
  }
  if (empty.row == "remove") {
    if (length(which(rowSums(data.spec) == 0)) > 0) {
      data.env <- data.env[-which(rowSums(data.spec) == 0), ]
      if (!is.null(xy)) 
        xy <- xy[-which(rowSums(data.spec) == 0), ]
      if (!is.null(data.spec.pred)) 
        if (inherits(data.spec.pred,"data.frame") )
          data.spec.pred <- data.spec.pred[-which(rowSums(data.spec) == 0), ]
      if (inherits(data.spec.pred,"list")) {
        for (p in 1:length(data.spec.pred)) data.spec.pred[[p]] <- data.spec.pred[[p]][-which(rowSums(data.spec) == 0), ]
      }
      data.spec <- data.spec[-which(rowSums(data.spec) == 0), ]
    }
  }
  if (reg.type == "ispline") {
    num <- which(sapply(data.env, inherits, "numeric"))
    if (length(num) > 1) {
      range.min <- apply(data.env[, num], 2, min)
      range.max <- apply(data.env[, num], 2, max)
    }else {
      range.min <- min(data.env[, num])
      range.max <- max(data.env[, num])
    }
  }
  if (reg.type == "ispline") {
    data.env.num <- as.data.frame(data.env[, which(sapply(data.env, inherits, "numeric"))])
    names(data.env.num) <- names(data.env)[which(sapply(data.env, inherits, "numeric"))]
    ts <- matrix(NA, ncol(data.env.num), 2 * order.ispline + kn.ispline)
    for (i in 1:ncol(data.env.num)) {
      data.env.num[, i] <- (data.env.num[, i] - min(data.env.num[, i]))/(max(data.env.num[, i]) - min(data.env.num[, i]))
      ts[i, ] <- c(rep(0, order.ispline), stats::quantile(data.env.num[, i], probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 1), 1/(kn.ispline + 1))), rep(1, order.ispline))
    }
    IE <- matrix(NA, nrow(data.env.num), (ncol(data.env.num) * (order.ispline + kn.ispline)))
    k = order.ispline
    for (j in 1:ncol(data.env.num)) {
      for (i in 1:(order.ispline + kn.ispline)) {
        xx <- 0
        for (x in data.env.num[, j]) {
          xx <- xx + 1
          if (x == 1) {
            IE[xx, (j - 1) * (order.ispline + kn.ispline) + i] <- 1
          }else {
            IE[xx, (j - 1) * (order.ispline + kn.ispline) + i] <- .Ii(i, k, x, ts[j, ])
          }
        }
      }
    }
    IE <- data.frame(IE)
    for (i in 1:(ncol(IE)/(order.ispline + kn.ispline))) {
      for (j in 1:(order.ispline + kn.ispline)) {
        names(IE)[(i - 1) * (order.ispline + kn.ispline) + j] <- paste(names(data.env.num)[i], j, sep = "")
      }
    }
    Fa <- as.data.frame(data.env[, which(sapply(data.env, inherits, "factor"))])
    names(Fa) <- names(data.env)[which(sapply(data.env, inherits, "factor"))]
    data.env <- cbind(data.frame(IE), Fa)
  }
  x.dim <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
  if (order == 1) {
    if (nrow(data.spec) < sam) {
      zeta.val <- rowSums(data.spec)
      data.var <- data.env
    }
    else {
      samp <- sample(nrow(data.spec), sam)
      zeta.val <- rowSums(data.spec[samp, ])
      data.var <- data.env[samp, ]
    }
    if (rescale == TRUE || normalize != FALSE) {
      zeta.val <- zeta.val/max(zeta.val)
    }
  }else {
    if (choose(x.dim, order) > sam) {
      u <- rep(NA, sam)
      if (!is.null(data.spec.pred)) {
        if (inherits(data.spec.pred,"data.frame") )
          u2 <- rep(NA, sam)
        if (inherits(data.spec.pred,"list")) {
          u2 <- list()
          for (p in 1:length(data.spec.pred)) u2[[p]] <- rep(NA, sam)
        }
      }
      data.var <- as.data.frame(matrix(NA, sam, dim(data.env)[2]))
      distance <- rep(NA, sam)
      for (z in 1:sam) {
        samp <- sample(1:x.dim, order, replace = FALSE)
        u[z] <- sum(apply(data.spec[samp, ], 2, prod))
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred,"data.frame") )
            u2[z] <- sum(apply(data.spec.pred[samp, 
            ], 2, prod))
          if (inherits(data.spec.pred,"list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]][z] <- sum(apply(data.spec.pred[[p]][samp, ], 2, prod))
          }
        }
        if (normalize == "Jaccard") {
          toto <- (ncol(data.spec) - sum(apply((1 - data.spec[samp, ]), 2, prod)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }else if (empty.row == 1) {
              u[z] <- 1
            }
          }else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred,"data.frame") ){
              tata <- (ncol(data.spec.pred) - sum(apply((1 - data.spec.pred[samp, ]), 2, prod)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred,"list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (ncol(data.spec.pred[[p]]) - sum(apply((1 - data.spec.pred[[p]][samp, ]), 2, prod)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }else if (normalize == "Sorensen") {
          toto <- (mean(apply(data.spec[samp, ], 1, 
                              sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }else if (empty.row == 1) {
              u[z] <- 1
            }
          }else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred,"data.frame") ){
              tata <- (mean(apply(data.spec.pred[samp, ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred,"list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (mean(apply(data.spec.pred[[p]][samp, ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }else if (normalize == "Simpson") {
          toto <- (min(apply(data.spec[samp, ], 1, sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }else if (empty.row == 1) {
              u[z] <- 1
            }
          }else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred,"data.frame") ){
              tata <- (min(apply(data.spec.pred[samp, ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred,"list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (min(apply(data.spec.pred[[p]][samp, ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        fac <- which(sapply(data.env, inherits, "factor"))
        num <- which(sapply(data.env, inherits, "numeric"))
        if (order > 2) {
          if (length(num) > 1) {
            toto <- data.env[samp,num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(apply(toto, 2, stats::dist), 2, get(method))
          }else if (length(num) > 0) {
            toto <- data.env[samp,num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(as.matrix(c(stats::dist(toto))), 2, get(method))
          }
        }else {
          if (length(num) > 1) {
            toto <- data.env[samp,num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(toto, 2, stats::dist)
          }else if (length(num) > 0) {
            toto <- data.env[samp,num]
            rownames(toto) <- c()
            data.var[z, num] <- stats::dist(toto)
          }
        }
        if (length(fac) > 1) {
          toto <- data.env[samp,fac]
          rownames(toto) <- c()
          data.var[z, fac] <- apply(toto, 2, function(x) {length(unique(x))}) - 1
        }else if (length(fac) > 0) {
          toto <- data.env[samp,fac]
          rownames(toto) <- c()
          data.var[z, fac] <- length(unique(toto)) - 1
        }
        if (!is.null(xy)) {
          if (distance.type == "Euclidean") {
            distance[z] <- apply(as.matrix(c(stats::dist(xy[samp, ]))), 2, get(method))
          }else if (distance.type == "ortho") {
            distance[z] <- apply(as.matrix(c(.gdist_matrix(xy[samp, ]))), 2, get(method))
          }else {
            stop("Error: invalid distance type")
          }
        }else if (distance.type == "custom") {
          if (is.null(dist.custom)) {
            stop("Error: a distance matrix must be provided if distance.type = 'custom'")
          }
          distance[z] <- apply(as.matrix(dist.custom[t(utils::combn(sort(samp), 2))]), 2, get(method))
        }
      }
      if (rescale == TRUE & normalize == FALSE) {
        u <- u/ncol(data.spec)
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred,"data.frame") )
            u2 <- u2/ncol(data.spec.pred)
          if (inherits(data.spec.pred,"list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]] <- u2[[p]]/ncol(data.spec.pred[[p]])
          }
        }
      }
    }else {
      u <- rep(NA, choose(x.dim, order))
      if (!is.null(data.spec.pred)) {
        if (inherits(data.spec.pred,"data.frame") )
          u2 <- rep(NA, choose(x.dim, order))
        if (inherits(data.spec.pred,"list")) {
          u2 <- list()
          for (p in 1:length(data.spec.pred)) u2[[p]] <- rep(NA, choose(x.dim, order))
        }
      }
      data.var <- as.data.frame(matrix(NA, choose(x.dim, order), dim(data.env)[2]))
      distance <- rep(NA, choose(x.dim, order))
      samp <- utils::combn(1:x.dim, order)
      for (z in 1:dim(samp)[2]) {
        u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred,"data.frame") )
            u2[z] <- sum(apply(data.spec.pred[samp[, z], ], 2, prod))
          if (inherits(data.spec.pred,"list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]][z] <- sum(apply(data.spec.pred[[p]][samp[, z], ], 2, prod))
          }
        }
        if (normalize == "Jaccard") {
          toto <- (ncol(data.spec) - sum(apply((1 - data.spec[samp[, z], ]), 2, prod)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }else if (empty.row == 1) {
              u[z] <- 1
            }
          }else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred,"data.frame") ){
              tata <- (ncol(data.spec.pred) - sum(apply((1 - data.spec.pred[samp[, z], ]), 2, prod)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred,"list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (ncol(data.spec.pred[[p]]) - sum(apply((1 - data.spec.pred[[p]][samp[, z], ]), 2, prod)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }else if (normalize == "Sorensen") {
          toto <- (mean(apply(data.spec[samp[, z], ], 1, sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }else if (empty.row == 1) {
              u[z] <- 1
            }
          }else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred,"data.frame") ){
              tata <- (mean(apply(data.spec.pred[samp[, z], ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred,"list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (mean(apply(data.spec.pred[[p]][samp[, z], ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }else if (normalize == "Simpson") {
          toto <- (min(apply(data.spec[samp[, z], ], 1, sum)))
          if (toto == 0) {
            if (empty.row == 0) {
              u[z] <- 0
            }else if (empty.row == 1) {
              u[z] <- 1
            }
          }else {
            u[z] <- u[z]/toto
          }
          if (!is.null(data.spec.pred)) {
            if (inherits(data.spec.pred,"data.frame") ){
              tata <- (min(apply(data.spec.pred[samp[, z], ], 1, sum)))
              if (tata == 0) {
                if (empty.row == 0) {
                  u2[z] <- 0
                }else if (empty.row == 1) {
                  u2[z] <- 1
                }
              }else {
                u2[z] <- u2[z]/tata
              }
            }
            if (inherits(data.spec.pred,"list")) {
              for (p in 1:length(data.spec.pred)) {
                tata <- (min(apply(data.spec.pred[[p]][samp[, z], ], 1, sum)))
                if (tata == 0) {
                  if (empty.row == 0) {
                    u2[[p]][z] <- 0
                  }else if (empty.row == 1) {
                    u2[[p]][z] <- 1
                  }
                }else {
                  u2[[p]][z] <- u2[[p]][z]/tata
                }
              }
            }
          }
        }
        fac <- which(sapply(data.env, inherits, "factor"))
        num <- which(sapply(data.env, inherits, "numeric"))
        if (order > 2) {
          if (length(num) > 1) {
            toto <- data.env[samp[,z],num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(apply(toto, 2, stats::dist), 2, get(method))
          }else if (length(num) > 0) {
            toto <- data.env[samp[,z],num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(as.matrix(c(stats::dist(toto))), 2, get(method))
          }
        }else {
          if (length(num) > 1) {
            toto <- data.env[samp[,z],num]
            rownames(toto) <- c()
            data.var[z, num] <- apply(toto, 2, stats::dist)
          }else if (length(num) > 0) {
            toto <- data.env[samp[,z],num]
            rownames(toto) <- c()
            data.var[z, num] <- stats::dist(toto)
          }
        }
        if (length(fac) > 1) {
          toto <- data.env[samp[,z],fac]
          rownames(toto) <- c()
          data.var[z, fac] <- apply(toto, 2, function(x) {length(unique(x))}) - 1
        }else if (length(fac) > 0) {
          toto <- data.env[samp[,z],fac]
          rownames(toto) <- c()
          data.var[z, fac] <- length(unique(toto)) - 1
        }
        if (!is.null(xy)) {
          if (distance.type == "Euclidean") {
            distance[z] <- apply(as.matrix(c(stats::dist(xy[samp[, z], ]))), 2, get(method))
          }else if (distance.type == "ortho") {
            distance[z] <- apply(as.matrix(c(.gdist_matrix(xy[samp[, z], ]))), 2, get(method))
          }else {
            stop("Error: invalid distance type")
          }
        }else if (distance.type == "custom") {
          if (is.null(dist.custom)) {
            stop("Error: a distance matrix must be provided if distance.type = 'custom'")
          }
          distance[z] <- apply(as.matrix(dist.custom[t(utils::combn(sort(samp[, z]), 2))]), 2, get(method))
        }
      }
      if (rescale == TRUE & normalize == FALSE) {
        u <- u/ncol(data.spec)
        if (!is.null(data.spec.pred)) {
          if (inherits(data.spec.pred,"data.frame") )
            u2 <- u2/ncol(data.spec.pred)
          if (inherits(data.spec.pred,"list")) {
            for (p in 1:length(data.spec.pred)) u2[[p]] <- u2[[p]]/ncol(data.spec.pred[[p]])
          }
        }
      }
    }
    zeta.val <- u
  }
  if (!is.null(xy) | distance.type == "custom") {
    distance.raw <- distance
    d <- max(distance)
  }
  if (rescale.pred == TRUE) {
    if (order > 1) {
      fac <- apply(data.var, 2, max)
      data.var <- data.var/matrix(rep(apply(data.var, 2, max), min(choose(x.dim, order), sam)), min(choose(x.dim, order), sam), dim(data.env)[2], byrow = T)
    }
    else {
      if (reg.type != "ispline") {
        num <- which(sapply(data.env, inherits, "numeric"))
        if (length(num) > 1) {
          range.min <- apply(data.var[, num], 2, min)
          range.max <- apply(data.var[, num], 2, max)
          data.var[, num] <- (data.var[, num] - matrix(rep(apply(data.var[, num], 2, min), min(choose(x.dim, order), sam)), min(choose(x.dim, order), sam), dim(data.var[, num])[2], , byrow = T))/matrix(rep((apply(data.var[, num], 2, max) - apply(data.var[, num], 2, min)), min(choose(x.dim, order), sam)), min(choose(x.dim, order), sam), dim(data.var[, num])[2], byrow = T)
        }else {
          range.min <- min(data.var[, num])
          range.max <- max(data.var[, num])
          data.var[, num] <- (data.var[, num] - min(data.var[, num]))/(max(data.var[, num]) - min(data.var[, num]))
        }
      }
    }
    if (!is.null(xy) | distance.type == "custom") {
      distance <- distance/max(distance)
    }
  }
  names(data.var) <- names(data.env)
  if ((!is.null(xy) | distance.type == "custom") & reg.type == 
      "ispline") {
    dist2 <- distance/max(distance)
    ts <- c(rep(0, order.ispline), stats::quantile(dist2, probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 1), 1/(kn.ispline + 1))), rep(1, order.ispline))
    IE <- matrix(NA, length(distance), (order.ispline + kn.ispline))
    k = order.ispline
    for (i in 1:(order.ispline + kn.ispline)) {
      xx <- 0
      for (x in dist2) {
        xx <- xx + 1
        if (x == 1) {
          IE[xx, i] <- 1
        }else {
          IE[xx, i] <- .Ii(i, k, x, ts)
        }
      }
    }
    distance <- data.frame(IE)
    for (i in 1:(order.ispline + kn.ispline)) {
      names(distance)[i] <- paste("distance", i, sep = "")
    }
  }
  if (reg.type == "ispline") {
    if (!is.null(data.spec.pred)) {
      if (inherits(data.spec.pred,"data.frame") ){
        sp <- 1 - u2
        spp <- matrix(sp, length(sp), 1)
        ts <- c(rep(0, order.ispline), stats::quantile(sp, probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 1), 1/(kn.ispline + 1))), rep(1, order.ispline))
        IE <- matrix(NA, length(u2), (order.ispline + kn.ispline))
        k = order.ispline
        for (i in 1:(order.ispline + kn.ispline)) {
          xx <- 0
          for (x in sp) {
            xx <- xx + 1
            if (x == 1) {
              IE[xx, i] <- 1
            }else {
              IE[xx, i] <- .Ii(i, k, x, ts)
            }
          }
        }
        sp.prev <- data.frame(IE)
        for (i in 1:(order.ispline + kn.ispline)) {
          names(sp.prev)[i] <- paste("Biotic", i, sep = "")
        }
      }
      if (inherits(data.spec.pred,"list")) {
        for (p in 1:length(data.spec.pred)) {
          sp <- 1 - u2[[p]]
          if (p == 1) {
            spp <- matrix(sp, length(sp), 1)
          }else {
            spp <- cbind(spp, sp)
          }
          ts <- c(rep(0, order.ispline), stats::quantile(sp, probs = seq(1/(kn.ispline + 1), 1 - 1/(kn.ispline + 1), 1/(kn.ispline + 1))), rep(1, order.ispline))
          IE <- matrix(NA, length(u2[[p]]), (order.ispline + kn.ispline))
          k = order.ispline
          for (i in 1:(order.ispline + kn.ispline)) {
            xx <- 0
            for (x in sp) {
              xx <- xx + 1
              if (x == 1) {
                IE[xx, i] <- 1
              }else {
                IE[xx, i] <- .Ii(i, k, x, ts)
              }
            }
          }
          if (p == 1) {
            sp.prev <- data.frame(IE)
          }else {
            sp.prev <- cbind(sp.prev, data.frame(IE))
          }
        }
        pp <- 0
        for (p in 1:length(data.spec.pred)) {
          for (i in 1:(order.ispline + kn.ispline)) {
            pp <- pp + 1
            names(sp.prev)[pp] <- paste("Biotic_", p, "_", i, sep = "")
          }
        }
      }
    }
  }else {
    if (!is.null(data.spec.pred)) {
      if (inherits(data.spec.pred,"data.frame") ){
        sp.prev <- data.frame(1 - u2)
        names(sp.prev) <- c("Biotic")
      }
      if (inherits(data.spec.pred,"list")) {
        for (p in 1:length(data.spec.pred)) {
          if (p == 1) {
            sp.prev <- data.frame(1 - u2[[p]])
          }else {
            sp.prev <- cbind(sp.prev, data.frame(1 - u2[[p]]))
          }
        }
        for (p in 1:length(data.spec.pred)) {
          names(sp.prev)[p] <- paste("Biotic", p, sep = "")
        }
      }
    }
  }
  if (is.null(xy) & distance.type != "custom" & is.null(data.spec.pred)) {
    data.tot <- data.var
  }else if (is.null(xy) & distance.type != "custom" & !is.null(data.spec.pred)) {
    data.tot <- cbind(data.var, sp.prev)
  }else if (!is.null(xy) & distance.type != "custom" & is.null(data.spec.pred)) {
    data.tot <- cbind(data.var, distance)
  }else {
    data.tot <- cbind(data.var, sp.prev, distance)
  }
  zeta.msgdm <- list()
  zeta.msgdm$val <- zeta.val
  zeta.msgdm$predictors <- data.tot
  if (reg.type == "ispline") {
    zeta.msgdm$range.min <- range.min
    zeta.msgdm$range.max <- range.max
    if (!is.null(data.spec.pred)) 
      zeta.msgdm$biotic <- spp
    if (!is.null(xy) | distance.type == "custom") 
      zeta.msgdm$distance <- distance.raw
  }
  if (rescale.pred == TRUE) {
    if (order > 1) {
      if (!is.null(xy) | distance.type == "custom") {
        zeta.msgdm$rescale.factor <- c(fac, d)
      }else {
        zeta.msgdm$rescale.factor <- fac
      }
    }else {
      if (reg.type != "ispline") {
        num <- which(sapply(data.env, inherits, "numeric"))
        if (length(num) > 1) {
          zeta.msgdm$range.min <- range.min
          zeta.msgdm$range.max <- range.max
        }
      }
    }
  }
  zeta.msgdm$my.order <- order
  zeta.msgdm$order.ispline <- order.ispline
  zeta.msgdm$kn.ispline <- kn.ispline
  if (reg.type == "glm") {
    if (method.glm == "glm.fit.cons") {
      zeta.msgdm.model <- glm.cons(zeta.val ~ ., data = data.tot, family = family, method = method.glm, cons = cons, cons.inter = cons.inter, control = control)
    }else {
      zeta.msgdm.model <- glm2::glm2(zeta.val ~ ., data = data.tot, family = family, method = method.glm, control = control)
      zeta.msgdm.confint <- suppressMessages(stats::confint(zeta.msgdm.model, level = confint.level))
    }
    if (dim(data.env)[2] > 1) {
      zeta.msgdm.vif <- car::vif(zeta.msgdm.model)
    }else {
      zeta.msgdm.vif <- NA
    }
    zeta.msgdm$model <- zeta.msgdm.model
    if (method.glm == "glm.fit2") 
      zeta.msgdm$confint <- zeta.msgdm.confint
    zeta.msgdm$vif <- zeta.msgdm.vif
  }else if (reg.type == "ngls") {
    data.tot2 <- cbind(rep(1, nrow(data.tot)), data.tot)
    start <- c(1, rep(-1, ncol(data.tot)))
    zeta.msgdm$model <- nnls::nnnpls(as.matrix(data.tot2), zeta.val, con = start)
  }else if (reg.type == "gam") {
    xnam <- names(data.tot)
    fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ", kn, ") + s(", sep = ""), sep = ""), ", k = ", kn, ")", sep = ""))
    zeta.msgdm$model <- mgcv::gam(fm, data = data.tot, family = family)
  }else if (reg.type == "scam") {
    xnam <- names(data.tot)
    fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ", kn, ", bs = '", bs, "') + s(", sep = ""), sep = ""), ", k = ", kn, ",bs='", bs, "')", sep = ""))
    zeta.msgdm$model <- scam::scam(fm, data = data.tot, family = family)
  }else if (reg.type == "ispline") {
    if (glm.init == TRUE) {
      tutu <- stats::cor(data.tot, zeta.val)
      tutu[which(tutu > 0)] <- 0
      zeta.msgdm$model <- glm.cons(zeta.val ~ ., data = data.tot, family = family, method = "glm.fit.cons", cons = cons, cons.inter = cons.inter, control = control, start = c(-1, tutu))
    }else {
      zeta.msgdm$model <- glm.cons(zeta.val ~ ., data = data.tot, family = family, method = "glm.fit.cons", cons = cons, cons.inter = cons.inter, control = control)
    }
  }else {
    stop("Error: unknown regression type.")
  }
  return(zeta.msgdm)
}





#' Fitting Generalized Linear Models with constraint on the coefficients signs
#'
#' \code{glm.cons} is an adaptation of function \code{glm2} from package \{glm2\} in which the least squares estimation is replaced by a regression with signs constraint on the coefficients using function \code{nnnpls} from package \{nnls\}.
#' @param formula as for \code{\link{glm}}
#' @param family as for \code{\link{glm}}
#' @param data as for \code{\link{glm}}
#' @param weights as for \code{\link{glm}}
#' @param subset as for \code{\link{glm}}
#' @param na.action as for \code{\link{glm}}
#' @param start as for \code{\link{glm}}
#' @param etastart as for \code{\link{glm}}
#' @param mustart as for \code{\link{glm}}
#' @param offset as for \code{\link{glm}}
#' @param control as for \code{\link{glm}}
#' @param model as for \code{\link{glm}}
#' @param method the method used in fitting the model. The default method "\code{glm.fit.cons}" uses function {nnnpls} from package nnls instead of \code{lm.fit} to impose the sign of the coefficients. As in \code{glm}, the alternative method "\code{model.frame}" returns the model frame and does no fitting.
#' @param cons type of constraint. Default is -1 for negative coefficients on the predictors. The other option is 1 for positive coefficients on the predictors.
#' @param cons.inter type of constraint for the intercept. Default is 1 for positive intercept, suitable for Gaussian family. The other option is -1 for negative intercept, suitable for binomial family.
#' @param x as for \code{\link{glm}}
#' @param y as for \code{\link{glm}}
#' @param contrasts as for \code{\link{glm}}
#' @param ... as for \code{\link{glm}}
#' @return The value returned by \code{glm.cons} has exactly the same structure as the value returned by \code{glm} and \code{glm.2}.
#' @references Marschner, I.C. (2011) glm2: Fitting generalized linear models with convergence problems. \emph{The R Journal}, 3(2), 12-15.
#' @seealso \code{\link{glm}}, \code{\link{glm2}}
#' @examples
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts))
#' glm.D93 <- glm.cons(counts ~ outcome + treatment, family = poisson())
#' glm.D93.ngl <- glm.cons(counts ~ outcome + treatment, family = poisson(),
#'    method="glm.fit.cons")
#' summary(glm.D93)
#' summary(glm.D93.ngl)
#' @export
glm.cons <- function (formula, family = stats::gaussian(), data, weights, subset,
                      na.action, start = NULL, etastart, mustart, offset, control = list(...),
                      model = TRUE, method = "glm.fit.cons", cons = -1, cons.inter = 1, x = FALSE, y = TRUE, contrasts = NULL,
                      ...)
{
  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame"))
    return(mf)
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  if (identical(method, "glm.fit.cons"))
    control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- stats::model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!stats::is.empty.model(mt))
    stats::model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(stats::model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  mustart <- stats::model.extract(mf, "mustart")
  etastart <- stats::model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, cons = cons, cons.inter = cons.inter, start = start, etastart = etastart,
                   mustart = mustart, offset = offset, family = family,
                   control = control, intercept = attr(mt, "intercept") >
                     0L))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method,
                      x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights, , cons = cons, cons.inter = cons.inter,
                      offset = offset, family = family, control = control,
                      intercept = TRUE))
    if (!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase maxit?")
    fit$null.deviance <- fit2$deviance
  }
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula, terms = mt,
                     data = data, offset = offset, control = control, method = method, cons = cons, cons.inter = cons.inter,
                     contrasts = attr(X, "contrasts"), xlevels = stats::.getXlevels(mt,
                                                                                    mf)))
  class(fit) <- c(fit$class, c("glm", "lm"))
  fit
}




#' Generalized Linear Models fitting method with negative coefficients constraint
#'
#' \code{glm.fit.cons} is an adaptation of function \code{glm.fit2} from package \{glm2\} in which the least squares estimation is replaced by a non-positive regression using function \code{nnnpls} from package \{nnls\}.
#' @param x as for \code{\link{glm.fit}}
#' @param y as for \code{\link{glm.fit}}
#' @param weights as for \code{\link{glm.fit}}
#' @param cons type of constraint. Default is -1 for negative coefficients on the predictors. The other option is 1 for positive coefficients on the predictors.
#' @param cons.inter type of constraint for the intercept. Default is 1 for positive intercept, suitable for Gaussian family. The other option is -1 for negative intercept, suitable for binomial family.
#' @param start as for \code{\link{glm.fit}}
#' @param etastart as for \code{\link{glm.fit}}
#' @param mustart as for \code{\link{glm.fit}}
#' @param offset as for \code{\link{glm.fit}}
#' @param family as for \code{\link{glm.fit}}
#' @param control as for \code{\link{glm.fit}}
#' @param intercept as for \code{\link{glm.fit}}
#' @return The value returned by \code{glm.fit.cons} has exactly the same structure as the value returned by \code{glm.fit} and \code{glm.fit2}.
#' @references Marschner, I.C. (2011) glm2: Fitting generalized linear models with convergence problems. \emph{The R Journal}, 3(2), 12-15.
#' @seealso \code{\link{glm.fit}}, \code{\link{glm.fit2}}
#' @examples
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts))
#' glm.D93 <- glm.cons(counts ~ outcome + treatment, family = poisson())
#' glm.D93.ngl <- glm.cons(counts ~ outcome + treatment, family = poisson(),
#'    method="glm.fit.cons")
#' summary(glm.D93)
#' summary(glm.D93.ngl)
#' @export
glm.fit.cons <- function (x, y, weights = rep(1, nobs), cons = -1, cons.inter = 1, start = NULL, etastart = NULL,
                          mustart = NULL, offset = rep(0, nobs), family = stats::gaussian(),
                          control = list(), intercept = TRUE)
{
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) rownames(y) else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv))
    stop("'family' argument seems not to be a valid family object",
         call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if (is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model",
           call. = FALSE)
    mu <- linkinv(eta)
    if (!validmu(mu))
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  }
  else {
    coefold <- NULL
    eta <- if (!is.null(etastart))
      etastart
    else if (!is.null(start))
      if (length(start) != nvars)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                      nvars, paste(deparse(xnames), collapse = ", ")),
             domain = NA)
    else {
      coefold <- start
      offset + as.vector(if (NCOL(x) == 1L) x*start else x%*%start)
    }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
      stop("cannot find valid starting values: please specify some",
           call. = FALSE)
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE
    for (iter in 1L:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (any(is.na(varmu)))
        stop("NAs in V(mu)")
      if (any(varmu == 0))
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good])))
        stop("NAs in d(mu)/d(eta)")
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning("no observations informative at iteration ",
                iter)
        break
      }
      z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
      w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
      ngoodobs <- as.integer(nobs - sum(!good))
      fit.lm <- stats::lm.fit(x=x[good, , drop = FALSE]*w, y=z*w, singular.ok=FALSE, tol=min(1e-07, control$epsilon/1000))
      A <- x[good, , drop = FALSE]*w
      b <- z*w
      fit <- nnls::nnnpls(A = A, b = b,con=c(cons.inter,rep(cons,ncol(A)-1)))
      if (any(!is.finite(fit$x))) {
        conv <- FALSE
        warning(gettextf("non-finite coefficients at iteration %d",
                         iter), domain = NA)
        break
      }
      if (nobs < fit.lm$rank)
        stop(gettextf("X matrix has rank %d, but only %d observations",
                      fit.lm$rank, nobs), domain = NA)
      start[fit.lm$qr$pivot] <- fit$x
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace)
        cat("Deviance =", dev, "Iterations -", iter,
            "\n")
      boundary <- FALSE
      if (!is.finite(dev)) {
        if (is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values",
               call. = FALSE)
        warning("step size truncated due to divergence",
                call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit)
            stop("inner loop 1; cannot correct step size",
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace)
          cat("Step halved: new deviance =", dev, "\n")
      }
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values",
               call. = FALSE)
        warning("step size truncated: out of bounds",
                call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit)
            stop("inner loop 2; cannot correct step size",
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace)
          cat("Step halved: new deviance =", dev, "\n")
      }
      if (((dev - devold)/(0.1 + abs(dev)) >= control$epsilon)&(iter>1)) {
        if (is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values",
               call. = FALSE)
        warning("step size truncated due to increasing deviance", call. = FALSE)
        ii <- 1
        while ((dev - devold)/(0.1 + abs(dev)) > -control$epsilon) {
          if (ii > control$maxit) break
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        if (ii > control$maxit) warning("inner loop 3; cannot correct step size")
        else if (control$trace) cat("Step halved: new deviance =", dev, "\n")
      }
      if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      }
      else {
        devold <- dev
        coef <- coefold <- start
      }
    }
    if (!conv)
      warning("glm.fit.cons: algorithm did not converge. Try increasing the maximum iterations", call. = FALSE)
    if (boundary)
      warning("glm.fit.cons: algorithm stopped at boundary value",
              call. = FALSE)
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps))
        warning("glm.fit.cons: fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE)
    }
    if (family$family == "poisson") {
      if (any(mu < eps))
        warning("glm.fit.cons: fitted rates numerically 0 occurred",
                call. = FALSE)
    }
    if (fit.lm$rank < nvars)
      coef[fit.lm$qr$pivot][seq.int(fit.lm$rank + 1, nvars)] <- NA
    xxnames <- xnames[fit.lm$qr$pivot]
    residuals <- (y - mu)/mu.eta(eta)
    fit.lm$qr$qr <- as.matrix(fit.lm$qr$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit.lm$qr$qr[1L:nr, 1L:nvars]
    }
    else Rmat <- fit.lm$qr$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit.lm$qr$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY) names(fit.lm$effects) <- c(xxnames[seq_len(fit.lm$rank)], rep.int("", sum(good) - fit.lm$rank))
  wtdmu <- if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) 0 else fit.lm$rank
  resdf <- n.ok - rank
  aic.model <- aic(y, n.ok, mu, weights, dev) + 2 * rank
  list(coefficients = coef, residuals = residuals, fitted.values = mu, effects = if (!EMPTY) fit.lm$effects, R = if (!EMPTY) Rmat, rank = rank, qr = if (!EMPTY) structure(fit.lm$qr[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"), family = family, linear.predictors = eta, deviance = dev, aic = aic.model, null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, df.residual = resdf, df.null = nulldf, y = y, converged = conv, boundary = boundary)
}




#' Transform data using I-splines
#'
#' Evaluates the I-splines for all variables of a data frame, as performed in \code{Zeta.msgdm}.
#' @param dat A data frame whose columns are variables to be transformed using I-splines.
#' @param order.ispline Order of the I-spline.
#' @param kn.ispline Number of knots in the I-spline.
#' @param rescale Indicates how to rescale the values between 0 and 1. Default is 0, which divides the data by the maximum value. Any other value corresponds to setting the minimum value to 0.
#' @return \code{Ispline} returns a data frame with the same number of rows as dat and
#' @return \code{ncol(dat)} * \code{(order.ispline} + \code{kn.ispline)} columns.
#' @references Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical Science}, 425-441.
#' @references Ferrier, S., Manion, G., Elith, J., & Richardson, K. (2007). Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity and Distributions}, 13(3), 252-264.
#' @seealso \code{\link{Zeta.msgdm}}
#' @examples
#' utils::data(bird.env.coarse)
#' data.env <- bird.env.coarse[,3:9]
#' data.env.splines <- Ispline(data.env)
#' @export
Ispline <- function(dat, order.ispline = 2, kn.ispline = 1, rescale = 0){
  ts <- matrix(NA,ncol(dat),2*order.ispline+ kn.ispline)
  for(i in 1:ncol(dat)){
    if(rescale == 0){
      dat[,i] <- dat[,i]/max(dat[,i])
    }else{
      dat[,i] <- (dat[,i]-min(dat[,i]))/(max(dat[,i])-min(dat[,i]))
    }
    ts[i,] <- c(rep(0,order.ispline),stats::quantile(dat[,i],probs=seq(1/(kn.ispline+1),1-1/(kn.ispline+1),1/(kn.ispline+1))),rep(1,order.ispline))
  }

  IE <- matrix(NA,nrow(dat),(ncol(dat)*(order.ispline+kn.ispline)))
  k=order.ispline
  for(j in 1:ncol(dat)){
    for(i in 1:(order.ispline+kn.ispline)){
      xx <- 0
      for(x in dat[,j]){
        xx <- xx+1
        if(x==1){
          IE[xx,(j-1)*(order.ispline+kn.ispline)+i] <- 1
        }else{
          IE[xx,(j-1)*(order.ispline+kn.ispline)+i] <- .Ii(i,k,x,ts[j,])
        }
      }
    }
  }
  IE <- data.frame(IE)
  for(i in 1:(ncol(IE)/(order.ispline+kn.ispline))){
    for(j in 1:(order.ispline+kn.ispline)){
      names(IE)[(i-1)*(order.ispline+kn.ispline)+j] <- paste(names(dat)[i],j,sep="")
    }
  }
  
  out <- list()
  out$splines <- IE
  out$order.ispline = order.ispline
  out$kn.ispline = kn.ispline

  return(out)
}



#' Perform an I-spline regression
#'
#' Evaluates the I-splines for all variables of a data frame of predictor variables, and perform a generalised linear regression with constraint on the parameters.
#' @param response A vector of numeric values representing the response variable.
#' @param predictor A data frame of numeric variables representing the predictors.
#' @param order.ispline Order of the I-spline.
#' @param kn.ispline Number of knots in the I-spline.
#' @param family A description of the error distribution and link function to be used in the \code{glm} model (see \code{\link[stats]{family}} for details of family functions).
#' @param method.glm Method used in fitting the generalised linear model. The default method \cr "glm.fit.cons" is an adaptation of method \code{glm.fit2} from package \code{glm2} using a constrained least squares regression in the reweighted least squares. Another option is "glm.fit2", which calls \code{glm.fit2}; see help documentation for glm.fit2 in package \code{glm2}.
#' @param cons type of constraint in the glm if \code{method.glm = "glm.fit.cons"}. Default is 1 for positive coefficients on the predictors. The other option is -1 for negative coefficients on the predictors.
#' @param cons.inter type of constraint for the intercept. Default is 1 for positive intercept, suitable for Gaussian family. The other option is -1 for negative intercept, suitable for binomial family.
#' @param Plot Boolean value indicating if the I-splines must be plotted.
#' @param lty Line types to be used in the plotting. If nothing is provided, \code{lty} is a sequence of integers from 1 to the number of variables used for the computation of \code{msgdm}.
#' @param lwd Line width.
#' @param control As for \code{\link{glm}}.
#' @details \code{Reg.ispline} performs a non-linear regression using a combination of GLM and I-splines. It can, for example, be used to compare regression outputs when using MS-GDM with I-splines on environmental variables and biotic variables as in \code{Zetya.msgdm} to the same regression approach without environmental variables.
#' @return \code{Reg.ispline} returns a list of the following elements:
#' @return \item{splines}{A data frame in which each columns contains the value resulting from the transformation of the predictors into individual I-splines. The number of columns of \code{splines} is the number of predictors times the number of splines (determined as the sum of \code{order.ispline} and \code{kn.ispline}).}
#' @return \item{spline}{A data frame in which each columns contains the value resulting from the combinations of the individual I-splines. This combination is obtained by multiplying the coefficients of \code{model} and the values of the individual I-splines \code{splines}}.
#' @return \item{model}{A \code{\link{glm}} model using \code{response} as the response variable, and \code{splines} as the predictors}.
#' @references Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical Science}, 425-441.
#' @seealso \code{\link{Zeta.msgdm}},\code{\link{Ispline}}
#' @examples
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[1:2]
#' 
#' data.spec.marion <- Marion.species[3:33]
#' 
#' ##random other communities
#' data.spec.marion2a <- data.spec.marion
#' data.spec.marion2a[which(data.spec.marion2a==1,arr.ind=TRUE)] <- 0
#' for(i in 1:ncol(data.spec.marion2a))
#'   data.spec.marion2a[sample(nrow(data.spec.marion2a),8),i] <- 1

#' data.spec.marion2b <- data.spec.marion
#' data.spec.marion2b[which(data.spec.marion2b==1,arr.ind=TRUE)] <- 0
#' for(i in 1:ncol(data.spec.marion2b))
#' data.spec.marion2b[sample(nrow(data.spec.marion2b),8),i] <- 1
#' 
#' dat.spec.tot <- list(data.spec.marion,data.spec.marion2a,data.spec.marion2b)
#' zeta.tot <- Zeta.order.mc.mult(data.spec=dat.spec.tot,order=3,sam=200)
#' zeta.splines <- Ispline(zeta.tot$zeta.val[,2:3])
#' data.tot <- data.frame(zeta.val=zeta.tot$zeta.val[,1],zeta.splines$splines)
#' 
#' dev.new()
#' Reg.ispline(response = zeta.tot$zeta.val[,1], predictor = zeta.tot$zeta.val[,2:3], lwd=2, cons=1)
#' @export
Reg.ispline <- function(response, predictor, order.ispline = 2, kn.ispline = 1, family = stats::gaussian(), method.glm = "glm.fit.cons", cons = 1, cons.inter = 1, control = list(), Plot = TRUE, lty = NULL, lwd = 1){
  
  if(is.null(lty)){
    lty=1:ncol(predictor)
  }
  
  n.splines <- order.ispline+kn.ispline
  
  predictor.spline <- Ispline(predictor,order.ispline = order.ispline, kn.ispline = kn.ispline)
  
  data.tot <- data.frame(response=response,predictor=predictor.spline)
  
  reg.model <- glm.cons(response ~ ., data = predictor.spline$splines, family = family, method = method.glm, cons = cons, cons.inter = cons.inter, control = control)
  
  splines <- matrix(NA,nrow(predictor),ncol(predictor))
  for(i in 1:ncol(predictor)){
    splines[,i] <- rowSums(matrix(reg.model$coefficients[(2+(i-1)*n.splines):(i*n.splines+1)],nrow(predictor.spline$splines),3,byrow = TRUE)*predictor.spline$splines[,(1+(i-1)*n.splines):(i*n.splines)])
  }
  
  if(Plot==TRUE){
    graphics::plot(sort(predictor[,1]),splines[order(predictor[,1]),1],type="l",xlim=range(predictor),ylim=range(splines),lwd=lwd,lty=lty[1])
    if(ncol(predictor)>1){
      for(i in 2:ncol(predictor)){
        graphics::lines(sort(predictor[,i]),splines[order(predictor[,i]),2],lwd=lwd,lty=lty[i])
      }
    }
    graphics::legend(x="topleft",legend = names(predictor),lwd=lwd,lty=lty,bty="n")
  }
  
  out.list <- list()
  out.list$splines <- predictor.spline
  out.list$spline.tot <- splines
  out.list$model <- reg.model
  
  return(out.list)
}




#' Predict zeta values for new environmental and distance data
#'
#' Predict the zeta values for new environmental and distance data from the models returned by \code{Zeta.msgdm}.
#' @param model.msgdm A model returned by \code{Zeta.msgdm}. The class of the model depends on the type of regression used in \code{Zeta.msgdm}.
#' @param reg.type Type of regression used in \code{Zeta.msgdm}. Options are "\code{glm}" for generalised linear models, "\code{ngls}" for negative linear models, "\code{gam}" for generalised additive models, "\code{scam}" for shape constrained additive models (with monotonic decreasing by default), and "\code{ispline}" for I-spline models (forcing monotonic decreasing), as recommended in generalised dissimilarity modelling by Ferrier \emph{et al}. (2007).
#' @param newdata A data frame with the new environmental and distance data. The names of the columns must be the same as the names used in the data frame used in \code{Zeta.msgdm}. For I-splines, the data frame must be generated beforehand from the original data by the function \code{\link{Ispline}}.
#' @param type The type of prediction required, as for \code{predict.glm}. The default is on the scale of the response variable; the alternative "link" is on the scale of the linear predictors.
#' @return \code{Predict.msgdm} returns a vector of predicted zeta values.
#' @references Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical Science}, 425-441.
#' @references Ferrier, S., Manion, G., Elith, J., & Richardson, K. (2007). Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity and Distributions}, 13(3), 252-264.
#' @seealso \code{\link{Zeta.msgdm}}
#' @import scam
#' @examples
#' utils::data(bird.spec.fine)
#' xy.bird <- bird.spec.fine[1:500,1:2]
#' data.spec.bird <- bird.spec.fine[1:500,3:192]
#' utils::data(bird.env.fine)
#' data.env.bird <- bird.env.fine[1:500,3:9]
#'
#' zeta.glm <- Zeta.msgdm(data.spec.bird, data.env.bird, sam = 100, order = 3)
#' newdata <- data.frame(matrix(NA,100,ncol(data.env.bird)))
#' names(newdata) <- names(data.env.bird)
#' for(z in 1:100){
#'   samp <- sample(1:104, 3, replace = FALSE)
#'   newdata[z,] <- apply(apply(bird.env.fine[501:604,3:9][samp,], 2,
#'      stats::dist), 2, mean)
#' }
#' ##rescale the data like during MS-GDM
#' newdata <- newdata/matrix(rep(zeta.glm$rescale.factor,100),
#'    100,length(zeta.glm$rescale.factor),byrow=TRUE)
#' new.zeta.glm <- Predict.msgdm(model.msgdm = zeta.glm$model, reg.type = "glm",
#'    newdata = newdata)
#'
#'
#'
#' zeta.ngls <- Zeta.msgdm(data.spec.bird, data.env.bird, sam = 100, order = 3,
#'    reg.type = "ngls", normalize = "Jaccard")
#' newdata <- data.frame(matrix(NA,100,ncol(data.env.bird)))
#' names(newdata) <- names(data.env.bird)
#' for(z in 1:100){
#'   samp <- sample(1:104, 3, replace = FALSE)
#'   newdata[z,] <- apply(apply(bird.env.fine[501:604,3:9][samp,], 2, stats::dist),
#'      2, mean)
#' }
#' ##rescale the data like during MS-GDM
#' newdata <- newdata/matrix(rep(zeta.ngls$rescale.factor,100),
#'    100,length(zeta.ngls$rescale.factor),byrow=TRUE)
#' new.zeta.ngls <- Predict.msgdm(model.msgdm = zeta.ngls$model, reg.type = "ngls",
#'    newdata = newdata)
#' @export
Predict.msgdm <- function(model.msgdm, reg.type, newdata, type = "response"){
  if(reg.type=="glm" || reg.type=="ispline"){
    new.zeta <- stats::predict.glm(object = model.msgdm, newdata = newdata, type = type)
  }else if(reg.type=="ngls"){
    coefs <- stats::coef(model.msgdm)
    newdata <- cbind(rep(1,nrow(newdata)), newdata)
    new.zeta <- c(coefs%*%t(as.matrix(newdata)))
  }else if(reg.type=="gam"){
    new.zeta <- mgcv::predict.gam(object = model.msgdm, newdata = newdata, type = type)
  }else if(reg.type=="scam"){
    new.zeta <- scam::predict.scam(object = model.msgdm, newdata = newdata, type = type)
  }else{
    stop("Error: unknown regression type.")
  }
  return(new.zeta)
}





#' Computing splines coordinates from I-spline-based multi-site generalised dissimilarity modelling
#'
#' Stores the coordinates of the I-splines resulting from \code{Zeta.msgdm} for plotting.
#' @param msgdm  Output of function \code{Zeta.msgdm} computed with \code{reg.type = ispline}.
#' @param data.env  Site-by-variable data frame used for the computation of \code{msgdm}, with sites as rows and environmental variables as columns.
#' @param distance Boolean, indicates is distance was used in the computation of \code{msgdm}.
#' @param biotic Integer, indicates the number of other groups of taxa for which zeta diversity was computed and used in the computation of \code{msgdm}.
#' @details \code{Return.ispline} allows to store the same number of coordinates for all I-splines, to average replicates and obtain confidence intervals.
#' @return \code{Return.ispline} returns a list containing the following components used to plot the I-splines:
#' @return \item{env}{A data frame containing the rescaled environmental (numeric and factor), distance and biotic x-values.}
#' @return \item{Ispline}{A data frame containing the I-spline values corresponding to the rescaled environmental (numeric and factor), distance and biotic x-values.}
#' @seealso \code{\link{Zeta.msgdm}}, \code{\link{Ispline}}
#' @examples
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[1:2]
#' data.spec.marion <- Marion.species[3:33]
#'
#' utils::data(Marion.env)
#' data.env.marion <- Marion.env[3]
#'
#' zeta.ispline <- Zeta.msgdm(data.spec.marion, data.env.marion, xy.marion, sam = 100,
#'     order = 3, normalize = "Jaccard", reg.type = "ispline")
#' zeta.ispline
#' zeta.ispline.r <- Return.ispline(zeta.ispline, data.env.marion, distance = TRUE)
#' zeta.ispline.r
#' 
#' dev.new()
#' Plot.ispline(isplines = zeta.ispline.r, distance = TRUE)
#' 
#' dev.new()
#' Plot.ispline(msgdm = zeta.ispline, data.env = data.env.marion, distance = TRUE)
#'
#' @export
#' 
Return.ispline <- function (msgdm, data.env, distance = FALSE, biotic = 0){
  
  
  ##for retro-compatibility
  if(biotic==FALSE)
    biotic <- 0
  if(biotic==TRUE)
    biotic <- 1
  ##
  
  my.order <- msgdm$my.order
  order.ispline = msgdm$order.ispline
  kn.ispline = msgdm$kn.ispline
  
  num.splines <- order.ispline+kn.ispline
  
  data.env.num <- as.data.frame(data.env[,which(sapply(data.env, inherits, "numeric"))])
  names(data.env.num) <- names(data.env)[which(sapply(data.env, inherits, "numeric"))]
  
  Fa <- as.data.frame(data.env[,which(sapply(data.env, inherits, "factor"))])
  names(Fa) <- names(data.env)[which(sapply(data.env, inherits, "factor"))]
  
  if(distance == FALSE & biotic == 0){
    XX <- as.data.frame(matrix(rep(seq(0,1,0.01),ncol(data.env)),101,ncol(data.env)))
    names(XX) <- names(data.env)
  }else if(distance == TRUE & biotic == 0){
    XX <- as.data.frame(matrix(rep(seq(0,1,0.01),ncol(data.env)+1),101,ncol(data.env)+1))
    names(XX) <- c(names(data.env),"Distance")
  }else if(distance == FALSE & biotic > 0){
    XX <- as.data.frame(matrix(rep(seq(0,1,0.01),ncol(data.env)+biotic),101,ncol(data.env)+biotic))
    names(XX) <- c(names(data.env),paste("Biotic",1:biotic,sep=""))
  }else{
    XX <- as.data.frame(matrix(rep(seq(0,1,0.01),ncol(data.env)+biotic+1),101,ncol(data.env)+biotic+1))
    names(XX) <- c(names(data.env),paste("Biotic",1:biotic,sep=""),"Distance")
  }
  
  
  env.ispline <- Ispline(data.env.num,rescale=1, order.ispline = order.ispline, kn.ispline = kn.ispline)$splines
  subsamp <- 0
  if(length(msgdm$val) < nrow(data.env)){
    ind.sel <- sample(nrow(data.env),length(msgdm$val))
    env.ispline <- env.ispline[ind.sel,]
    data.env <- data.env[ind.sel,]
    subsamp <- 1
  }
  
  
  if(distance == TRUE){
    d.ind <- c(sample(1:length(msgdm$val),min(nrow(data.env)-2,length(msgdm$val)-2)),which.max(msgdm$distance),which.min(msgdm$distance))
    d <- msgdm$distance[d.ind]
    d.spline <- msgdm$predictors[d.ind,(ncol(msgdm$predictors)-(num.splines-1)):ncol(msgdm$predictors)]
    
    if(biotic >0){
      for(b in 1:biotic){
        bio.ind <- c(sample(1:length(msgdm$val),min(nrow(data.env)-2,length(msgdm$val)-2)),which.max(msgdm$biotic[,b]),which.min(msgdm$biotic[,b]))
        if(b==1){
          bio <- msgdm$biotic[bio.ind,b]
          bio <- matrix(bio,length(bio),1)
          bio.spline <- msgdm$predictors[bio.ind,(ncol(msgdm$predictors)-((biotic-(b-2))*3-1)):(ncol(msgdm$predictors)-((biotic-(b-1))*3))]
        }else{
          bio <- cbind(bio,msgdm$biotic[bio.ind,b])
          bio.spline <- cbind(bio.spline,msgdm$predictors[bio.ind,(ncol(msgdm$predictors)-((biotic-(b-2))*3-1)):(ncol(msgdm$predictors)-((biotic-(b-1))*3))])
        }
      }
      
      if(length(which(sapply(data.env, inherits, "factor")))==0){
        X.ispline <- cbind(env.ispline,bio.spline,d.spline)
        Isplines.pred <- matrix(NA,ncol(data.env)+biotic+1,nrow(data.env))
        for(i in 1:(ncol(data.env)+biotic+1)){
          Isplines.pred[i,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(i-1)*num.splines):(i*num.splines+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(i-1)*num.splines):(i*num.splines)])
        }
      }else{
        X.ispline <- cbind(env.ispline,matrix(rep(1,ncol(Fa)*nrow(env.ispline)),nrow(env.ispline),ncol(Fa)),bio.spline,d.spline)
        Isplines.pred <- matrix(NA,ncol(data.env)+biotic+1,nrow(data.env))
        for(i in 1:ncol(data.env.num)){
          Isplines.pred[i,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(i-1)*num.splines):(i*num.splines+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(i-1)*num.splines):(i*num.splines)])
        }
        ii <- 0
        for(i in (ncol(data.env.num)+1):(ncol(data.env.num)+ncol(Fa))){
          ii <- ii+1
          Isplines.pred[i,] <- -stats::coef(msgdm$model)[(1+ncol(data.env.num)*num.splines)+ii]* X.ispline[,ncol(data.env.num)*num.splines+ii]
        }
        for(b in 1:biotic){
          Isplines.pred[ncol(data.env)+b,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(ncol(data.env.num)+b-1)*num.splines+ncol(Fa)):((ncol(data.env.num)+b)*num.splines+ncol(Fa)+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(ncol(data.env.num)+b-1)*num.splines+ncol(Fa)):((ncol(data.env.num)+b)*num.splines+ncol(Fa))])
        }
        Isplines.pred[ncol(data.env)+biotic+1,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(ncol(data.env.num)+biotic-1)*num.splines+ncol(Fa)+num.splines):((ncol(data.env.num)+biotic)*num.splines+ncol(Fa)+1+num.splines)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(ncol(data.env.num)+biotic-1)*num.splines+ncol(Fa)+num.splines):((ncol(data.env.num)+biotic)*num.splines+ncol(Fa)+num.splines)])
      }
    }else{
      if(length(which(sapply(data.env, inherits, "factor")))==0){
        X.ispline <- cbind(env.ispline,d.spline)
        Isplines.pred <- matrix(NA,ncol(data.env)+1,min(nrow(data.env),length(msgdm$val)))
        for(i in 1:(ncol(data.env.num)+1)){
          Isplines.pred[i,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(i-1)*num.splines):(i*num.splines+1)],min(nrow(data.env),length(msgdm$val)),num.splines,byrow=TRUE)* X.ispline[,(1+(i-1)*num.splines):(i*num.splines)])
        }
      }else{
        X.ispline <- cbind(env.ispline,matrix(rep(1,ncol(Fa)*nrow(env.ispline)),nrow(env.ispline),ncol(Fa)),d.spline)
        Isplines.pred <- matrix(NA,ncol(data.env)+1,nrow(data.env))
        for(i in 1:(ncol(data.env.num)+1)){
          Isplines.pred[i,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(i-1)*num.splines):(i*num.splines+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(i-1)*num.splines):(i*num.splines)])
        }
        ii <- 0
        for(i in (ncol(data.env.num)+1):(ncol(data.env.num)+ncol(Fa))){
          ii <- ii+1
          Isplines.pred[i,] <- -stats::coef(msgdm$model)[(1+ncol(data.env.num)*num.splines)+ii]* X.ispline[,ncol(data.env.num)*num.splines+ii]
        }
        Isplines.pred[ncol(data.env)+1,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(ncol(data.env.num))*num.splines+ncol(Fa)):((ncol(data.env.num)+1)*num.splines+ncol(Fa)+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+ncol(data.env.num)*num.splines+ncol(Fa)):((ncol(data.env.num)+1)*num.splines+ncol(Fa))])
      }
    }
  }else{
    if(biotic > 0){
      for(b in 1:biotic){
        bio.ind <- c(sample(1:length(msgdm$val),min(nrow(data.env)-2,length(msgdm$val)-2)),which.max(msgdm$biotic[,b]),which.min(msgdm$biotic[,b]))
        if(b==1){
          bio <- msgdm$biotic[bio.ind,b]
          bio <- matrix(bio,length(bio),1)
          bio.spline <- msgdm$predictors[bio.ind,(ncol(msgdm$predictors)-((biotic-(b-1))*3-1)):(ncol(msgdm$predictors)-((biotic-(b))*3))]
        }else{
          bio <- cbind(bio,msgdm$biotic[bio.ind,b])
          bio.spline <- cbind(bio.spline,msgdm$predictors[bio.ind,(ncol(msgdm$predictors)-((biotic-(b-1))*3-1)):(ncol(msgdm$predictors)-((biotic-(b))*3))])
        }
      }
      
      if(length(which(sapply(data.env, inherits, "factor")))==0){
        X.ispline <- cbind(env.ispline,bio.spline)
        Isplines.pred <- matrix(NA,ncol(data.env)+biotic,nrow(data.env))
        for(i in 1:(ncol(data.env)+biotic)){
          Isplines.pred[i,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(i-1)*num.splines):(i*num.splines+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(i-1)*num.splines):(i*num.splines)])
        }
      }else{
        X.ispline <- cbind(env.ispline,matrix(rep(1,ncol(Fa)*nrow(env.ispline)),nrow(env.ispline),ncol(Fa)),bio.spline)
        Isplines.pred <- matrix(NA,ncol(data.env)+biotic,nrow(data.env))
        for(i in 1:ncol(data.env.num)){
          Isplines.pred[i,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(i-1)*num.splines):(i*num.splines+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(i-1)*num.splines):(i*num.splines)])
        }
        ii <- 0
        for(i in (ncol(data.env.num)+1):(ncol(data.env.num)+ncol(Fa))){
          ii <- ii+1
          Isplines.pred[i,] <- -stats::coef(msgdm$model)[(1+ncol(data.env.num)*num.splines)+ii]* X.ispline[,ncol(data.env.num)*num.splines+ii]
        }
        for(b in 1:biotic){
          Isplines.pred[ncol(data.env)+b,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(ncol(data.env.num)+b-1)*num.splines+ncol(Fa)):((ncol(data.env.num)+b)*num.splines+ncol(Fa)+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(ncol(data.env.num)+b-1)*num.splines+ncol(Fa)):((ncol(data.env.num)+b)*num.splines+ncol(Fa))])
        }
      }
    }else{
      if(length(which(sapply(data.env, inherits, "factor")))==0){
        X.ispline <- env.ispline
      }else{
        X.ispline <- cbind(env.ispline,matrix(rep(1,ncol(Fa)*nrow(env.ispline)),nrow(env.ispline),ncol(Fa)))
      }
      Isplines.pred <- matrix(NA,ncol(data.env),nrow(data.env))
      for(i in 1:ncol(data.env.num)){
        Isplines.pred[i,] <- rowSums(matrix(-stats::coef(msgdm$model)[(2+(i-1)*num.splines):(i*num.splines+1)],nrow(data.env),num.splines,byrow=TRUE)* X.ispline[,(1+(i-1)*num.splines):(i*num.splines)])
      }
      if(length(which(sapply(data.env, inherits, "factor")))>0){
        ii <- 0
        for(i in (ncol(data.env.num)+1):(ncol(data.env.num)+ncol(Fa))){
          ii <- ii+1
          Isplines.pred[i,] <- -stats::coef(msgdm$model)[(1+ncol(data.env.num)*num.splines)+ii]* X.ispline[,ncol(data.env.num)*num.splines+ii]
        }
      }
    }
  }
  
  
  if(subsamp == 1){
    env.resc <- data.env.num[ind.sel,]
  }else{
    env.resc <- data.env.num
  }
  for(i in 1:ncol(env.resc)){
    env.resc[,i] <- (env.resc[i]-min(env.resc[i]))/(max(env.resc[i])-min(env.resc[i]))
  }
  
  ##create output list
  splines.out <- list()
  env.resc.out <- env.resc
  Isplines.pred.out <- Isplines.pred
  for(i in 1:ncol(data.env.num)){
    env.resc.out[,i] <- sort(env.resc[,i])
    Isplines.pred.out[i,] <- Isplines.pred[i,order(env.resc[,i])]
  }
  splines.out$env <- data.frame(env.resc.out)
  
  if(ncol(Fa)>0){
    Fa.out <- matrix(0,ncol(Fa),ncol(Isplines.pred))
    ii <- 0
    for(i in (ncol(data.env.num)+1):(ncol(data.env.num)+ncol(Fa))){
      ii <- ii+1
      Fa.out[ii,] <- seq(0,1,1/(ncol(Isplines.pred)-1))
      Isplines.pred.out[i,] <- seq(0,max(Isplines.pred[i,]),max(Isplines.pred[i,])/(ncol(Isplines.pred)-1))
    }
    splines.out$env <- cbind(splines.out$env,data.frame(t(Fa.out)))
    for(i in (ncol(data.env.num)+1):(ncol(data.env.num)+ncol(Fa))){
      names(splines.out$env)[i] <- paste("Fa",i-ncol(data.env.num),sep="")
    }
  }
  
  if(biotic > 0){
    for(b in 1:biotic){
      i <- i+1
      bio.out <- sort(bio[,b])
      Isplines.pred.out[i,] <- Isplines.pred[i,order(bio[,b])]
      splines.out$env <- cbind(splines.out$env,bio.out)
      names(splines.out$env)[i] <- paste("biotic",b,sep="")
    }
  }
  
  if(distance == TRUE){
    i <- i+1
    d.out <- sort(d/max(d))
    Isplines.pred.out[i,] <- Isplines.pred[i,order(d)]
    splines.out$env <- cbind(splines.out$env,d.out)
    names(splines.out$env)[i] <- "distance"
  }
  
  splines.out$Ispline <- data.frame(t(Isplines.pred.out))
  names(splines.out$Ispline) <- names(splines.out$env)
  
  splines.out$env.num <- data.env.num
  splines.out$Fa <- Fa
  splines.out$distance <- distance
  splines.out$biotic <- biotic
  splines.out$my.order <- my.order
  # if(distance == TRUE)
  
  return(splines.out)
}




#' Plots I-splines for Multi-Site Generalised Dissimilarity Modelling
#'
#'Plots I-splines computed by \code{Return.ispline}, or calls \code{Return.ispline} if the outputs from \code{Zeta.msgdm} are provided before plotting.
#' @param isplines Output of function \code{Return.ispline}.
#' @param msgdm  Output of function \code{Zeta.msgdm} computed with \code{reg.type = ispline}.
#' @param data.env  Site-by-variable data frame used for the computation of \code{msgdm}, with sites as rows and environmental variables as columns.
#' @param distance Boolean, indicates is distance was used in the computation of \code{msgdm}.
#' @param biotic Boolean, indicates is zeta diversity from another community was used in the computation of \code{msgdm}.
#' @param pch Shapes of the points to be used in the plotting. If nothing is provided, \code{pch} is a sequence of integers from 1 to the number of variables used for the computation of \code{msgdm}.
#' @param lty Line types to be used in the plotting. If nothing is provided, \code{lty} is a sequence of integers from 1 to the number of variables used for the computation of \code{msgdm}.
#' @param legend  Boolean, indicates if the legend must be drawn.
#' @param lwd Line width.
#' @param cex Point size.
#' @param num.quantiles Number of points to plot on the I-splines. Default is 11 to plot a point every 10 percents of the range of values.
#' @return \code{Plot.ispline} returns a data frame with the same number of rows as dat and \code{ncol(dat)} * \code{(order.ispline} + \code{kn.ispline)} columns.
#' @references Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical Science}, 425-441.
#' @references Ferrier, S., Manion, G., Elith, J., & Richardson, K. (2007). Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity and Distributions}, 13(3), 252-264.
#' @seealso \code{\link{Zeta.msgdm}}
#' @examples
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[1:2]
#' data.spec.marion <- Marion.species[3:33]
#'
#' utils::data(Marion.env)
#' data.env.marion <- Marion.env[3]
#'
#' zeta.ispline <- Zeta.msgdm(data.spec.marion, data.env.marion, xy.marion, sam = 100,
#'     order = 3, normalize = "Jaccard", reg.type = "ispline")
#' zeta.ispline
#' zeta.ispline.r <- Return.ispline(zeta.ispline, data.env.marion, distance = TRUE)
#' zeta.ispline.r
#' 
#' dev.new()
#' Plot.ispline(isplines = zeta.ispline.r, distance = TRUE)
#' 
#' dev.new()
#' Plot.ispline(msgdm = zeta.ispline, data.env = data.env.marion, distance = TRUE)
#' 
#' 
#' @export
#' 
Plot.ispline <- function (isplines = NULL,msgdm, data.env, distance = FALSE, biotic = 0, pch = NULL, lty = NULL, legend = TRUE, lwd = 1, cex = 1, num.quantiles = 11){
  
  if(is.null(isplines)){
    isplines <- Return.ispline(msgdm=msgdm, data.env=data.env, distance = distance, biotic = biotic)
  }
  
  env.resc <- isplines$env
  Isplines.pred <- isplines$Ispline
  data.env.num <- isplines$env.num
  Fa <- isplines$Fa
  distance <- isplines$distance
  biotic <- isplines$biotic
  my.order <- isplines$my.order
  
  if(is.null(pch)){
    pch <- 1:ncol(env.resc)
  }
  if(is.null(lty)){
    lty <- 1:ncol(env.resc)
  }
  
  graphics::plot(sort(env.resc[,1]),Isplines.pred[order(env.resc[,1]),1], type="l",ylim=c(0,max(Isplines.pred)),main="",xlab="Rescaled range",ylab="I-splines",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,cex=cex,lwd=lwd)
  ind.points <- numeric()
  for(i in 1:num.quantiles){
    ind.points[i] <- which.min(abs(stats::quantile(env.resc[,1],seq(0,1,1/(num.quantiles-1)))[i]-sort(env.resc[,1])))
  }
  graphics::points(sort(env.resc[,1])[ind.points],Isplines.pred[order(env.resc[,1]),1][ind.points],pch=pch[1],cex=cex,lwd=lwd)
  if(ncol(data.env.num) > 1){
    for(i in 2:ncol(data.env.num)){
      graphics::lines(sort(env.resc[,i]),Isplines.pred[order(env.resc[,i]),i],lty=lty[i],lwd=lwd)
      ind.points <- numeric()
      for(j in 1:num.quantiles){
        ind.points[j] <- which.min(abs(stats::quantile(env.resc[,i],seq(0,1,1/(num.quantiles-1)))[j]-sort(env.resc[,i])))
      }
      graphics::points(sort(env.resc[,i])[ind.points],Isplines.pred[order(env.resc[,i]),i][ind.points],pch=pch[i],cex=cex,lwd=lwd)
    }
  }
  
  i <- ncol(data.env.num)
  
  if(!is.null(Fa)){
    if(ncol(Fa)>0){
      for(i in (ncol(data.env.num)+1):(ncol(data.env.num)+ncol(Fa))){
        if(max(Isplines.pred[,i])>0){
          graphics::lines(seq(0,1,1/nrow(Isplines.pred)),seq(0,max(Isplines.pred[,i]),max(Isplines.pred[,i])/nrow(Isplines.pred)),lty=lty[i],col="red",lwd=lwd)
          graphics::points (seq(0,1,1/(my.order-1)),seq(0,max(Isplines.pred[,i]),max(Isplines.pred[,i])/(my.order-1)),pch=pch[i],col="red",lwd=lwd,cex=cex)
        }else{
          graphics::lines(c(0,1),c(0,0),lty=lty[i],col="red",lwd=lwd)
          graphics::points (seq(0,1,1/(my.order-1)),rep(0,my.order),pch=pch[i],col="red",lwd=lwd,cex=cex)
        }
      }
    }
  }
  
  if(biotic >0){
    for(b in 1:biotic){
      i <- i+1
      graphics::lines(env.resc[,i],Isplines.pred[,i],lty=lty[i],lwd=lwd,col="green")
      ind.points <- numeric()
      for(j in 1:num.quantiles){
        ind.points[j] <- which.min(abs(stats::quantile(env.resc[,i],seq(0,1,1/(num.quantiles-1)))[j]-env.resc[,i]))
      }
      graphics::points(env.resc[ind.points,i],Isplines.pred[ind.points,i],pch=pch[i],cex=cex,lwd=lwd,col="green")
    }
  }
  
  if(distance == TRUE){
    i <- i+1
    graphics::lines(env.resc$distance,Isplines.pred$distance,lty=lty[i],lwd=lwd,col="blue")
    ind.points <- numeric()
    for(j in 1:num.quantiles){
      ind.points[j] <- which.min(abs(stats::quantile(env.resc$distance,seq(0,1,1/(num.quantiles-1)))[j]-env.resc$distance))
    }
    graphics::points(env.resc$distance[ind.points],Isplines.pred$distance[ind.points],pch=pch[i],cex=cex,lwd=lwd,col="blue")
  }
  
  if(is.null(Fa) | length(Fa)==0) {
    if(distance == FALSE & biotic == 0){
      legend.text <- c(names(data.env.num))
      col="black"
    }else if(distance == TRUE & biotic == 0){
      legend.text <- c(names(data.env.num),"Distance")
      col=c(rep("black",ncol(data.env.num)),"blue")
    }else if(distance == FALSE & biotic > 0){
      legend.text <- c(names(data.env.num),paste("Biotic",1:biotic))
      col=c(rep("black",ncol(data.env.num)),rep("green",biotic))
    }else{
      legend.text <- c(names(data.env.num),paste("Biotic",1:biotic),"Distance")
      col=c(rep("black",ncol(data.env.num)),rep("green",biotic),"blue")
    }
  }else{
    if(distance == FALSE & biotic == 0){
      legend.text <- c(names(data.env.num),names(Fa))
      col=c(rep("black",ncol(data.env.num)),rep("red",ncol(Fa)))
    }else if(distance == TRUE & biotic == 0){
      legend.text <- c(names(data.env.num),names(Fa),"Distance")
      col=c(rep("black",ncol(data.env.num)),rep("red",ncol(Fa)),"blue")
    }else if(distance == FALSE & biotic > 0){
      legend.text <- c(names(data.env.num),names(Fa),paste("Biotic",1:biotic))
      col=c(rep("black",ncol(data.env.num)),rep("red",ncol(Fa)),rep("green",biotic))
    }else{
      legend.text <- c(names(data.env.num),names(Fa),paste("Biotic",1:biotic),"Distance")
      col=c(rep("black",ncol(data.env.num)),rep("red",ncol(Fa)),rep("green",biotic),"blue")
    }
  }
  if(legend == TRUE){
    #legend("topleft",lty=lty,pch=pch,names(XX),lwd=lwd,cex=cex,bty="n")
    legend("topleft",lty=lty,pch=pch,legend=legend.text,lwd=lwd,cex=cex,bty="n",col=col)
  }
  
  ##create output list
  invisible(isplines)
  
}





#' Zeta distance decay for a specific number of assemblages or sites
#'
#' Computes the distance decay of zeta diversity for a specific order (number of assemblages or sites), using either a generalised linear model with possible constraint on the coefficients, a generalised additive model, or a shape constrained additive model.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta diversity is computed.
#' @param distance.type Method to compute distance. Default is "\code{Euclidean}", for Euclidean distance. The other options are (i) "\code{ortho}" for orthodromic distance, if xy correspond to longitudes and latitudes (orthodromic distance is computed using the \code{geodist} function from package \code{geodist}); and (ii) "\code{custom}", in which case the user must provide a distance matrix for \code{dist.custom}.
#' @param dist.custom Distance matrix provided by the user when \code{distance.type} = \code{"custom"}.
#' @param method Name of a function (as a string) indicating how to combine the pairwise differences and distances for more than 3 sites. It can be a basic R-function such as "\code{mean}" or "\code{max}", but also a custom function.
#' @param method.glm Method used in fitting the generalised linear model. The default method \cr "glm.fit.cons" is an adaptation of method \code{glm.fit2} from package \code{glm2} using a negative least squares regression in the reweighted least squares. Another option is "glm.fit2", which calls \code{glm.fit2}.; see help documentation for \code{glm.fit2} in package \code{glm}.
#' @param cons type of constraint in the glm if \code{method.glm = "glm.fit.cons"}. Default is -1 for negative coefficients on the predictors. The other option is 1 for positive coefficients on the predictors.
#' @param cons.inter type of constraint for the intercept. Default is 1 for positive intercept, suitable for Gaussian family. The other option is -1 for negative intercept, suitable for binomial family.
#' @param reg.type Type of regression. Options are "\code{glm}" for generalised linear models "\code{gam}" for generalised additive models and "\code{scam}" for shape constrained additive models (with monotonic decreasing by default).
#' @param family A description of the error distribution and link function to be used in the \code{glm}, \code{gam} and \code{scam} models (see \code{\link[stats]{family}} for details of family functions).
#' @param confint.level  Percentage for the confidence intervals of the coefficients from the generalised linear models.
#' @param kn Number of knots in the GAM and SCAM. Default is -1 for determining kn automatically using Generalized Cross-validation.
#' @param bs A two-letter character string indicating the (penalized) smoothing basis to use in the scam model. Default is "\code{mpd}" for monotonic decreasing splines. see \code{\link[mgcv]{smooth.terms}} for an overview of what is available.
#' @param trsf Name of a function (as a string) indicating how to transform distance.
#' @param cutoff If specified, maximum distance value for which the linear regression must be performed.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1. Has no effect if \code{normalize} != \code{FALSE}.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param empty.row Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity, and computation errors for the normalized version of zeta due to divisions by 0. Options are "\code{empty}" to let the data untreated, "\code{remove}" to remove the empty rows, 0 to set the normalized zeta to 0 when zeta is divided by 0 during normalization (sites share no species, so are completely dissimilar), and 1 to set the normalized zeta to 1 when zeta is divided by 0 during normalization (i.e. sites are perfectly similar).
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @return \code{Zeta.ddecay} returns a list containing the following components:
#' @return \item{order}{The order of zeta for which the distance decay was computed.}
#' @return \item{reg.type}{A character string indicating the type of regression that was performed.}
#' @return \item{reg}{An object whose class depends on the type of regression (\code{glm}, \code{gam} or \code{scam}), corresponding to the regression over distance for the number of assemblages or sites specified in 'order'.}
#' @return \item{confint}{The confidence intervals for the coefficients from the generalised linear model. \code{confint} is not generated for generalised additive models and shape constrained additive models.}
#' @return \item{zeta.val}{The values of zeta for the sampled sites used in the regression.}
#' @return \item{distance}{The distances for the sampled sites used in the regression.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Zeta.ddecays}},
#' @seealso \code{\link{Plot.zeta.ddecay}}
#' @import scam
#' @examples
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[,1:2]
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' dev.new()
#' zeta.ddecay.bird <- Zeta.ddecay(xy.bird, data.spec.bird, sam = 100, order = 3,
#'     method.glm = "glm.fit2", confint.level = 0.95)
#'     
#' dev.new()
#' zeta.ddecay.bird <- Zeta.ddecay(data.spec=data.spec.bird, distance.type = "custom",
#'     dist.custom = as.matrix(dist(xy.bird)), cutoff = 800000, sam = 100, order = 3,
#'     reg.type = "gam", confint.level = 0.95)
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' dev.new()
#' zeta.ddecay.marion <- Zeta.ddecay(xy.marion, data.spec.marion, sam = 100, order = 3,
#'     method.glm = "glm.fit2", confint.level = 0.95, trsf = "log", normalize = "Jaccard")
#'
#' @export
Zeta.ddecay <- function(xy, data.spec, order = 2, sam = 1000, distance.type = "Euclidean", dist.custom = NULL, method = "mean", reg.type = "glm", family = stats::gaussian(), method.glm = "glm.fit.cons", cons = -1, cons.inter = 1, confint.level = 0.95, kn = -1, bs = "mpd", trsf = "NULL", cutoff = NULL, rescale = FALSE, normalize = FALSE, empty.row = "remove", plot = TRUE){
  
  
  if(distance.type != "custom"){
    if(nrow(data.spec) != nrow(xy)){
      stop("Error: data.spec and xy must have the same number of rows.")
    }
  }
  
  if (!inherits(data.spec,"data.frame")){
    stop(paste("Error: ",deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  
  if(empty.row == "remove" & distance.type != "custom"){
    if(length(which(rowSums(data.spec)==0))>0){
      xy <- xy[-which(rowSums(data.spec)==0),]
      data.spec <- data.spec[-which(rowSums(data.spec)==0),]
    }
  }
  
  x <- dim(data.spec)[1]
  zeta.val <- numeric()
  zeta.val.sd <- numeric()
  
  if(choose(x, order)>sam){
    u <- rep(NA, sam)
    distance <- rep(NA, sam)
    
    for(z in 1:sam){
      samp <- sample(1:x, order, replace = FALSE)
      u[z] <- sum(apply(data.spec[samp, ], 2, prod))
      if (normalize == "Jaccard"){
        toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp, ]), 2, prod)))
        if(toto==0){
          if(empty.row == 0){
            u[z] <- 0
          }else if(empty.row == 1){
            u[z] <- 1
          }
        }else
          u[z] <- u[z] / toto
      }else if (normalize == "Sorensen"){
        toto <- (mean(apply(data.spec[samp, ], 1, sum)))
        if(toto==0){
          if(empty.row == 0){
            u[z] <- 0
          }else if(empty.row == 1){
            u[z] <- 1
          }
        }else
          u[z] <- u[z] / toto
      }else if (normalize == "Simpson"){
        toto <- (min(apply(data.spec[samp, ], 1, sum)))
        if(toto==0){
          if(empty.row == 0){
            u[z] <- 0
          }else if(empty.row == 1){
            u[z] <- 1
          }
        }else
          u[z] <- u[z] / toto
      }
      if(order == 1){
        stop("Error: distance decay cannot be computed for zeta 1")
      }else {
        if(distance.type == "Euclidean"){
          distance[z] <- apply(as.matrix(c(stats::dist(xy[samp, ]))),2,get(method))
        }else if(distance.type == "ortho"){
          distance[z] <- apply(as.matrix(c(.gdist_matrix(xy[samp, ]))),2,get(method))
        }else if(distance.type == "custom"){
          if(is.null(dist.custom))
          {
            stop("Error: a distance matrix must be provided if distance.type = 'custom'")
          }
          distance[z] <- apply(as.matrix(dist.custom[t(utils::combn(sort(samp),2))]),2,get(method))
        }else{
          stop("Error: invalid distance type")
        }
      }
    }
    
  }else{
    u <- rep(NA, choose(x, order))
    distance <- rep(NA, choose(x, order))
    samp <- utils::combn(1:x, order)
    
    for(z in 1:dim(samp)[2]){
      u[z] <- sum(apply(data.spec[samp[, z], ], 2, prod))
      if (normalize == "Jaccard"){
        toto <- (ncol(data.spec)-sum(apply((1-data.spec[samp[, z], ]), 2, prod)))
        if(toto==0){
          if(empty.row == 0){
            u[z] <- 0
          }else if(empty.row == 1){
            u[z] <- 1
          }
        }else
          u[z] <- u[z] / toto
      }else if (normalize == "Sorensen"){
        toto <- (mean(apply(data.spec[samp[, z], ], 1, sum)))
        if(toto==0){
          if(empty.row == 0){
            u[z] <- 0
          }else if(empty.row == 1){
            u[z] <- 1
          }
        }else
          u[z] <- u[z] / toto
      }else if (normalize == "Simpson"){
        toto <- (min(apply(data.spec[samp[, z], ], 1, sum)))
        if(toto==0){
          if(empty.row == 0){
            u[z] <- 0
          }else if(empty.row == 1){
            u[z] <- 1
          }
        }else
          u[z] <- u[z] / toto
      }
      if(order == 1){
        stop("Error: distance decay cannot be computed for zeta 1")
      }else {
        if(distance.type == "Euclidean"){
          distance[z] <- apply(as.matrix(c(stats::dist(xy[samp[, z], ]))),2,get(method))
        }else if(distance.type == "ortho"){
          distance[z] <- apply(as.matrix(c(.gdist_matrix(xy[samp[, z], ]))),2,get(method))
        }else if(distance.type == "custom"){
          if(is.null(dist.custom))
          {
            stop("Error: a distance matrix must be provided if distance.type = 'custom'")
          }
          distance[z] <- apply(as.matrix(dist.custom[t(utils::combn(sort(samp[, z]),2))]),2,get(method))
        }
        else{
          stop("Error: invalid distance type")
        }
      }
    }
  }
  
  if(rescale == TRUE & normalize == FALSE){
    #z1 <- mean(rowSums(data.spec))
    #u <- u / z1
    u <- u / ncol(data.spec)
  }
  zeta.val <- u
  
  distance.reg <- distance
  zeta.val.reg <- zeta.val
  if(!is.null(cutoff)){
    distance.reg <- distance[which(distance <= cutoff)]
    zeta.val.reg <- zeta.val.reg[which(distance <= cutoff)]
  }
  if(trsf != "NULL"){
    distance.reg <- c(apply(as.matrix(distance.reg),2,get(trsf)))
  }
  
  
  zeta.ddecay <- list()
  zeta.ddecay$order <- order
  zeta.ddecay$reg.type <- reg.type
  if(reg.type == "glm"){
    if(method.glm == "glm.fit.cons")
      zeta.ddecay.reg <- glm.cons(zeta.val.reg ~ distance.reg, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
    else
      zeta.ddecay.reg <- glm2::glm2(zeta.val.reg ~ distance.reg, family = family, method = method.glm)
    zeta.ddecay$reg <- zeta.ddecay.reg
    if(method.glm == "glm.fit2")
      zeta.ddecay$confint <- suppressMessages(stats::confint(zeta.ddecay.reg, level = confint.level))
    
    if(plot == TRUE){
      preds <- stats::predict(zeta.ddecay.reg, newdata = data.frame(distance.reg = sort(distance.reg)), type = "link", se.fit = TRUE)
      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit
      fit2 <- zeta.ddecay.reg $family$linkinv(fit)
      upr2 <- zeta.ddecay.reg $family$linkinv(upr)
      lwr2 <- zeta.ddecay.reg $family$linkinv(lwr)
      graphics::plot(distance.reg, zeta.val.reg, xlab = "Distance", ylab = paste("Zeta ", order, sep = ""), pch = 16)
      graphics::lines(sort(distance.reg), fit2, col = "red", lwd = 2)
      graphics::lines(sort(distance.reg), upr2, col = "red", lty = 2, lwd = 2)
      graphics::lines(sort(distance.reg), lwr2, col = "red", lty = 2, lwd = 2)
    }
  }else if(reg.type == "gam"){
    fm <- stats::as.formula(paste("zeta.val.reg ~ s(distance.reg, k = ", kn ,")",sep=""))
    zeta.ddecay.reg <- mgcv::gam(fm, family = family)
    zeta.ddecay$reg <- zeta.ddecay.reg
    if(plot == TRUE){
      preds <- stats::predict(zeta.ddecay.reg, newdata = data.frame(distance.reg = sort(distance.reg)), se.fit = TRUE)
      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit
      graphics::plot(distance.reg, zeta.val.reg, xlab = "Distance", ylab = paste("Zeta ", order, sep = ""), pch = 16)
      graphics::lines(sort(distance.reg), fit, col = "red", lwd = 2)
      graphics::lines(sort(distance.reg), lwr, col = "red", lty = 2, lwd = 2)
      graphics::lines(sort(distance.reg), upr, col = "red", lty = 2, lwd = 2)
    }
  }else if(reg.type == "scam"){
    data.reg <- data.frame(zeta.val.reg,distance.reg)
    fm <- stats::as.formula(paste("zeta.val.reg ~ s(distance.reg, k = ", kn ,", bs = '", bs ,"')",sep=""))
    zeta.ddecay.reg <- scam::scam(fm,data=data.reg, family = family)
    zeta.ddecay$reg <- zeta.ddecay.reg
    if(plot == TRUE){
      preds <- scam::predict.scam(zeta.ddecay.reg, newdata = data.frame(distance.reg = sort(distance.reg)), se.fit = TRUE)
      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit
      graphics::plot(distance.reg, zeta.val.reg, xlab = "Distance", ylab = paste("Zeta ", order, sep = ""), pch = 16)
      graphics::lines(sort(distance.reg), fit, col = "red", lwd = 2)
      graphics::lines(sort(distance.reg), lwr, col = "red", lty = 2, lwd = 2)
      graphics::lines(sort(distance.reg), upr, col = "red", lty = 2, lwd = 2)
    }
  }else{
    stop("Error: unknown regression type.")
  }
  
  
  zeta.ddecay$zeta.val <- zeta.val.reg
  zeta.ddecay$distance <- distance.reg
  
  return(zeta.ddecay)
  
}


#' Zeta distance-decay plotting
#'
#' Plots the output of the function \code{Zeta.ddecay}.
#' @param zeta.ddecay A list produced by the function \code{Zeta.ddecay}.
#' @return A plot of the zeta distance-decay with distance on the x-axis and the value of zeta on the y-axis.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Zeta.ddecay}},
#' @seealso \code{\link{Zeta.ddecays}}
#' @examples
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[,1:2]
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' dev.new()
#' zeta.ddecay.bird <- Zeta.ddecay(xy.bird, data.spec.bird, sam = 100, order = 3,
#'     confint.level = 0.95,plot=FALSE)
#' Plot.zeta.ddecay(zeta.ddecay.bird)
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' zeta.ddecay.marion <- Zeta.ddecay(xy.marion, data.spec.marion, sam = 100, order = 3,
#'     confint.level = 0.95, trsf = "log", normalize = "Jaccard",plot=FALSE)
#' dev.new()
#' Plot.zeta.ddecay(zeta.ddecay.marion)
#'
#' @export
Plot.zeta.ddecay <- function(zeta.ddecay){

  if(zeta.ddecay$reg.type == "glm"){
    preds <- stats::predict(zeta.ddecay$reg, newdata = data.frame(distance.reg = sort(zeta.ddecay$distance)), type = "link", se.fit = TRUE)
    critval <- 1.96 ## approx 95% CI
    upr <- preds$fit + (critval * preds$se.fit)
    lwr <- preds$fit - (critval * preds$se.fit)
    fit <- preds$fit
    fit2 <- zeta.ddecay$reg $family$linkinv(fit)
    upr2 <- zeta.ddecay$reg $family$linkinv(upr)
    lwr2 <- zeta.ddecay$reg $family$linkinv(lwr)
    graphics::plot(zeta.ddecay$distance, zeta.ddecay$zeta.val, xlab = "Distance", ylab = paste("Zeta ", zeta.ddecay$order, sep = ""), pch = 16)
    graphics::lines(sort(zeta.ddecay$distance), fit2, col = "red", lwd = 2)
    graphics::lines(sort(zeta.ddecay$distance), upr2, col = "red", lty = 2, lwd = 2)
    graphics::lines(sort(zeta.ddecay$distance), lwr2, col = "red", lty = 2, lwd = 2)
  }else if(zeta.ddecay$reg.type == "gam"){
    preds <- stats::predict(zeta.ddecay$reg, newdata = data.frame(distance.reg = sort(zeta.ddecay$distance)), se.fit = TRUE)
    critval <- 1.96 ## approx 95% CI
    upr <- preds$fit + (critval * preds$se.fit)
    lwr <- preds$fit - (critval * preds$se.fit)
    fit <- preds$fit
    graphics::plot(zeta.ddecay$distance, zeta.ddecay$zeta.val, xlab = "Distance", ylab = paste("Zeta ", zeta.ddecay$order, sep = ""), pch = 16)
    graphics::lines(sort(zeta.ddecay$distance), fit, col = "red", lwd = 2)
    graphics::lines(sort(zeta.ddecay$distance), lwr, col = "red", lty = 2, lwd = 2)
    graphics::lines(sort(zeta.ddecay$distance), upr, col = "red", lty = 2, lwd = 2)
  }else if(zeta.ddecay$reg.type == "scam"){
    preds <- scam::predict.scam(zeta.ddecay$reg, newdata = data.frame(distance.reg = sort(zeta.ddecay$distance)), se.fit = TRUE)
    critval <- 1.96 ## approx 95% CI
    upr <- preds$fit + (critval * preds$se.fit)
    lwr <- preds$fit - (critval * preds$se.fit)
    fit <- preds$fit
    graphics::plot(zeta.ddecay$distance, zeta.ddecay$zeta.val, xlab = "Distance", ylab = paste("Zeta ", zeta.ddecay$order, sep = ""), pch = 16)
    graphics::lines(sort(zeta.ddecay$distance), fit, col = "red", lwd = 2)
    graphics::lines(sort(zeta.ddecay$distance), lwr, col = "red", lty = 2, lwd = 2)
    graphics::lines(sort(zeta.ddecay$distance), upr, col = "red", lty = 2, lwd = 2)
  }else{
    stop("Error: Unknown regression type.")
  }

}



#' Zeta distance decay for a range of numbers of assemblages or sites
#'
#' Computes the distance decay of zeta diversity for a range of orders (number of assemblages or sites), using generalised linear models.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param orders  Range of number of assemblages or sites at which zeta diversity is computed. All the orders must be striclty greater than 1.
#' @param sam  Number of samples for which the zeta diversity is computed.
#' @param family A description of the error distribution and link function to be used in the generalised linear models (see \code{\link[stats]{family}} for details of family functions).
#' @param confint.level  Percentage for the confidence intervals of the coefficients from the linear regression.
#' @param distance.type Method to compute distance. Default is "\code{Euclidean}", for Euclidean distance. The other options are (i) "\code{ortho}" for orthodromic distance, if xy correspond to longitudes and latitudes (orthodromic distance is computed using the \code{geodist} function from package \code{geodist}); and (ii) "\code{custom}", in which case the user must provide a distance matrix for \code{dist.custom}.
#' @param dist.custom Distance matrix provided by the user when \code{distance.type} = \code{"custom"}.
#' @param method Name of a function (as a string) indicating how to combine the pairwise differences and distances for more than 3 sites. It can be a basic R-function such as "\code{mean}" or "\code{max}", but also a custom function.
#' @param trsf Name of a function (as a string) indicating how to transform distance. Default is "NULL" for the identity transformation.
#' @param cutoff If specified, maximum distance value for which the linear regression must be performed.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1. Has no effect if \code{normalize} != \code{FALSE}.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @return \code{Zeta.ddecays} returns a list containing the following components:
#' @return \item{orders}{Range of number of assemblages or sites at which zeta diversity was computed.}
#' @return \item{coefs}{A vector of the coefficients from the generalised linear models for the numbers of sites specified by \code{orders}.}
#' @return \item{confint}{The confidence intervals for the coefficients from the generalised linear models.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Zeta.ddecay}}
#' @examples
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[,1:2]
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' dev.new()
#' zeta.ddecays.bird <- Zeta.ddecays(xy.bird, data.spec.bird, sam = 100, orders = 2:5,
#'     plot = TRUE, confint.level = 0.95)
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' dev.new()
#' zeta.ddecays.marion <- Zeta.ddecays(xy.marion, data.spec.marion, sam = 100,
#'    orders = 2:5, plot = TRUE, confint.level = 0.95)
#'
#' @export
Zeta.ddecays <- function(xy, data.spec, orders = 2:10, sam = 1000, family = stats::gaussian(), distance.type = "Euclidean", dist.custom = NULL, method = "mean", confint.level = 0.95, trsf = "NULL", cutoff = NULL, rescale = FALSE, normalize = FALSE, plot = TRUE){

  if(nrow(data.spec) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }

  if (!inherits(data.spec,"data.frame")){
    stop(paste("Error: ",deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(max(orders)>dim(data.spec)[1]){
    stop("Error: wrong value for \"orders\": the maximum value must be equal or lower than the number of sites.")
  }

  if(length(which(orders<= 1))>0){stop("Error: orders must be striclty greater than 1")}

  zeta.ddecays.coefs <- rep(NA, length(orders))
  zeta.ddecays.confint <- matrix(NA, length(orders), 2)
  ii <- 0
  for (i in orders){
    print(i)
    ii <- ii + 1
    temp <- Zeta.ddecay(xy, data.spec, sam = sam, order = i, distance.type = distance.type,  dist.custom = dist.custom, method = method, reg.type = "glm", family = family, method.glm = "glm.fit2", confint.level = confint.level, trsf = trsf, cutoff = cutoff, normalize = normalize, plot = FALSE)
    zeta.ddecays.coefs[ii] <- stats::coef(temp$reg)[2]
    zeta.ddecays.confint[ii, ] <- temp$confint[2, ]

  }

  zeta.ddecays <- list()
  zeta.ddecays$orders <- orders
  zeta.ddecays$coefs <- zeta.ddecays.coefs
  zeta.ddecays$confint <- zeta.ddecays.confint

  if (plot == TRUE){
    graphics::plot(orders, zeta.ddecays$coefs, pch = 16, ylim = c(min(0, range(zeta.ddecays$confint)[1]), max(0, range(zeta.ddecays$confint)[2])), xlab = "Order of zeta", ylab = "Slope", main = "Distance decay of zeta diversity")
    graphics::lines(orders, zeta.ddecays$coefs)
    suppressWarnings(graphics::arrows(x0 = c(orders), y0 = zeta.ddecays$coefs, x1 = c(orders), y1 = zeta.ddecays$confint[, 1], angle = 90, length = 0.2))
    suppressWarnings(graphics::arrows(x0 = c(orders), y0 = zeta.ddecays$coefs, x1 = c(orders), y1 = zeta.ddecays$confint[, 2], angle = 90, length = 0.2))
    graphics::lines(0:(max(orders) + 2), rep(0, max(orders) + 3), lty = 2)
  }

  return(zeta.ddecays)

}



#' Zeta distance-decay plotting for multiple orders
#'
#' Plots the output of the function \code{Zeta.ddecays}.
#' @param zeta.ddecays A list produced by the function \code{Zeta.ddecays}.
#' @return A plot of the zeta distance-decay with the orders on the x-axis and the slope of the linear distance-decays on the y-axis.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Zeta.ddecays}},
#' @seealso \code{\link{Zeta.ddecay}}, \code{\link{Plot.zeta.ddecay}}
#' @examples
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[,1:2]
#' data.spec.bird <- bird.spec.coarse[,3:193]
#'
#' dev.new()
#' zeta.ddecays.bird <- Zeta.ddecays(xy.bird, data.spec.bird, sam = 100, orders = 2:5,
#'     plot = FALSE, confint.level = 0.95)
#' Plot.zeta.ddecays(zeta.ddecays.bird)
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' dev.new()
#' zeta.ddecays.marion <- Zeta.ddecays(xy.marion, data.spec.marion, sam = 100,
#'    orders = 2:5, plot = FALSE, confint.level = 0.95)
#' Plot.zeta.ddecays(zeta.ddecays.marion)
#'
#' @export
Plot.zeta.ddecays <- function(zeta.ddecays){

  graphics::plot(zeta.ddecays$orders, zeta.ddecays$coefs, pch = 16, ylim = c(min(0, range(zeta.ddecays$confint)[1]), max(0, range(zeta.ddecays$confint)[2])), xlab = "Order of zeta", ylab = "Slope", main = "Distance decay of zeta diversity")
  graphics::lines(zeta.ddecays$orders, zeta.ddecays$coefs)
  suppressWarnings(graphics::arrows(x0 = c(zeta.ddecays$orders), y0 = zeta.ddecays$coefs, x1 = c(zeta.ddecays$orders), y1 = zeta.ddecays$confint[, 1], angle = 90, length = 0.2))
  suppressWarnings(graphics::arrows(x0 = c(zeta.ddecays$orders), y0 = zeta.ddecays$coefs, x1 = c(zeta.ddecays$orders), y1 = zeta.ddecays$confint[, 2], angle = 90, length = 0.2))
  graphics::lines(0:(max(zeta.ddecays$orders) + 2), rep(0, max(zeta.ddecays$orders) + 3), lty = 2)

}



#' Rescaling of data following a hierarchical increase in grain size
#'
#' Increases grain by hierarchically nesting regularly spaced sites.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param data.env  Site-by-variable data frame, with sites as rows and environmental variables as columns.
#' @param method Name of a function (as a string) indicating how to combine the coordinates and the environmental variables. It can be a basic R-function such as "\code{mean}" or "\code{max}", but also a custom function.
#' @param n  Mapping grain (the number of sites combined to generate data at a coarser grain). Regularly spaced sites are grouped as \code{n} x \code{n} sites.
#' @details The sites (plots or quadrates) are aggregated as nearest neighbouring groups of \code{n} x \code{n} sites, using a nested approach, starting from the lowest x and y, to increase the grain. The sites can be spatially contiguous or discontiguous, as long as they are regularly spaced. This function is not suitable for irregularly spaced sites. If the total number of sites is not a multiple of \code{n} x \code{n}, the extra sites are discarded.
#' @return \code{rescale.regular} returns a data frame with the rescaled data.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Zeta.scale.regular}}, \code{\link{Zeta.scale.min.dist}}, \code{\link{rescale.min.dist}}
#' @examples
#'
#' utils::data(bird.spec.fine)
#' xy.bird <- bird.spec.fine[1:2]
#' data.spec.bird <- bird.spec.fine[3:192]
#'
#' data.rescale <- rescale.regular(xy.bird, data.spec.bird, n = 4)
#' @export
rescale.regular <- function(xy, data.spec, data.env=NULL, method = "mean", n){

  if(nrow(data.spec) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }

  if (!inherits(data.spec,"data.frame")){
    stop(paste("Error: ",deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }

  ##sort data according to the plots coordinates
  xy <- xy[order(xy[,1], xy[,2]), ]
  data.spec <- data.spec[order(xy[,1], xy[,2]), ]

  ##compute the scale dependence for the specified grains
  Ux <- sort(unique(xy[,1]))
  Uy <- sort(unique(xy[,2]))

  max.x <- Ux[length(Ux) - length(Ux)%%(n^2)]
  max.y <- Uy[length(Uy) - length(Uy)%%(n^2)]


  if(length(which(xy[,1]>max.x | xy[,2]>max.y))>0){
    xy2 <- xy[-which(xy[,1]>max.x | xy[,2]>max.y), ]
    data.spec2 <- data.spec[-which(xy[,1]>max.x | xy[,2]>max.y), ]
    if(!is.null(data.env)){
      data.env2 <- data.env[-which(xy[,1]>max.x | xy[,2]>max.y), ]
    }
  }else{
    xy2 <- xy
    data.spec2 <- data.spec
    if(!is.null(data.env)){
      data.env2 <- data.env
    }
  }

  names <- c(names(xy), names(data.spec))
  data2 <- data.frame(stats::setNames(replicate(length(names), numeric(0), simplify = F), names))

  for (i in seq(1, length(unique(xy2$x)), n)){
    for (j in seq(1, length(unique(xy2$y)), n)){
      if(length(which(xy2$x %in% Ux[i:(i + n - 1)] & xy2$y %in% Uy[j:(j + n - 1)]))>0){
        temp.xy <- xy2[which(xy2$x %in% Ux[i:(i + n - 1)] & xy2$y %in% Uy[j:(j + n - 1)]), ]
        temp.xy <- apply(temp.xy, 2, get(method))

        temp.spec <- data.spec2[which(xy2$x %in% Ux[i:(i + n - 1)] & xy2$y %in% Uy[j:(j + n - 1)]), ]
        temp.spec <- (apply(temp.spec, 2, sum)>0) * 1

        if(!is.null(data.env)){
          temp.env <- data.env2[which(xy2$x %in% Ux[i:(i + n - 1)] & xy2$y %in% Uy[j:(j + n - 1)]), ]
          temp.env <- apply(temp.env, 2, get(method))
        }

        ## add the vector to the new coarser dataset
        if(is.null(data.env)){
          data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec)))
        }else{
          data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec, temp.env)))
        }
      }
    }
  }

  if(is.null(data.env)){
    names(data2) <- c(names(xy),names(data.spec))
  }else{
    names(data2) <- c(names(xy),names(data.spec),names(data.env))
  }

  return(data2)

}



#' Rescaling of data based on the minimum distance between sites
#'
#' Combines sites based on the minimum distance between them.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param data.env  Site-by-variable data frame, with sites as rows and environmental variables as columns.
#' @param m  Mapping grain (the number of sites combined to generate data at a coarser grain). The \code{m} closest sites are grouped together.
#' @param distance.type Method to compute distance. Default is "\code{Euclidean}", for Euclidean distance. The other options are (i) "\code{ortho}" for orthodromic distance, if xy correspond to longitudes and latitudes (orthodromic distance is computed using the \code{geodist} function from package \code{geodist}); and (ii) "\code{custom}", in which case the user must provide a distance matrix for \code{dist.custom}.
#' @param dist.custom Distance matrix provided by the user when \code{distance.type} = \code{"custom"}.
#' @param method  Name of a function (as a string) indicating how to combine the coordinates and the environmental variables. It can be a basic R-function such as "\code{mean}" or "\code{max}", but also a custom function.
#' @param shuffle Boolean value (TRUE or FALSE) indicating if the order of the sites must be randomised, which can have an impact on the outputs if some distances are equal.
#' @details The nearest neighbouring sites (plots, quadrates, or areas of varying shapes) are grouped as spatial clusters of 2, 3, 4, etc. sites, based on the minimum distance between them. Since the procedure is based on the relative distance between sites, the site order can have an impact on the output. This function is suitable for both regularly and irregularly spaced sites, contiguous or non contiguous. For regularly spaced sites, the use of \code{\link{rescale.regular}} is recommended.
#' @return \code{rescale.min.dist} returns a data frame with the rescaled data.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}},
#' @seealso \code{\link{Zeta.scale.min.dist}}, \code{\link{Zeta.scale.regular}}, \code{\link{rescale.regular}}
#' @examples
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' data.rescale <- rescale.min.dist(xy.marion, data.spec.marion, m=2)
#'
#' @export
rescale.min.dist <- function(xy, data.spec, data.env = NULL, m,distance.type = "Euclidean", dist.custom = NULL, method = "mean", shuffle =FALSE){

  if(nrow(data.spec) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }

  if (!inherits(data.spec,"data.frame")){
    stop(paste("Error: ",deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(!is.null(dist.custom)){
    if(!isSymmetric(dist.custom)){
      stop("Error: distance matrix is not symmetrical")
    }
  }

  ##randomize the plot orders
  if (shuffle == TRUE){
    xy <- xy[sample(nrow(xy)), ]
    data.spec <- data.spec[sample(nrow(xy)), ]
  }

  #pairwise distance
  if(distance.type == "Euclidean"){
    D <- as.matrix(stats::dist(xy))
  }else if(distance.type == "ortho"){
    D <- as.matrix(.gdist_matrix(stats::dist(xy)))
  }else if(distance.type == "custom"){
    if(is.null(dist.custom))
    {
      stop("Error: a distance matrix must be provided if distance.type = 'custom'")
    }
    D <- dist.custom
  }else{
    stop("Error: invalid distance type")
  }

  D[which(D == 0)] <- NA

  names <- c(names(xy), names(data.spec))
  data2 <- data.frame(stats::setNames(replicate(length(names), numeric(0), simplify = F), names))

  i.plot <- 1
  for(i in 1:floor(dim(xy)[1] / m)){

    while(length(which(!is.na(D[, i.plot]))) == 0){
      i.plot <- i.plot + 1
    }

    ##select the plots in the group
    if(distance.type != "custom"){
      temp.xy <- xy[c(i.plot, order(D[, i.plot])[1:(m - 1)]), ]
    }
    temp.spec <- data.spec[c(i.plot, order(D[, i.plot])[1:(m - 1)]), ]

    if(!is.null(data.env)){
      temp.env <- data.env[c(i.plot, order(D[, i.plot])[1:(m - 1)]), ]
    }

    ##compute the mean coordinates, the unions of species, and the mean of environmental variables
    temp.xy <- apply(temp.xy, 2, get(method))
    temp.spec <- (apply(temp.spec, 2, sum)>0) * 1
    if(!is.null(data.env)){
      temp.env <- apply(temp.env, 2, get(method))
    }

    ## add the vector to the new coarser dataset
    if(is.null(data.env)){
      data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec)))
    }else{
      data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec, temp.env)))
    }

    ##set the column corresponding to the plots to NA in the distance matrix to avoid using them in the following
    DD <- D
    D[, order(DD[, i.plot])[1:(m - 1)]] <- NA
    D[order(DD[, i.plot])[1:(m - 1)], ] <- NA
    D[, i.plot] <- NA
    D[i.plot, ] <- NA
  }

  if(is.null(data.env)){
    names(data2) <- c(names(xy),names(data.spec))
  }else{
    names(data2) <- c(names(xy),names(data.spec),names(data.env))
  }

  return(data2)

}




#' Zeta diversity scaling with sample grain using hierarchical increases in grain size
#'
#' Computes zeta diversity scaling with sample grain for a specific order (number of assemblages or sites), increasing grain by hierarchically nesting of regularly spaced sites.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param n  Vector of mapping grains: regularly spaced sites are grouped as \code{n[i]} x \code{n[i]} sites to generate data at a coarser grain.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param sam  Number of samples for which the zeta diversity is computed.
#' @param method Name of a function (as a string) indicating how to combine the coordinates. It can be a basic R-function such as "\code{mean}" or "\code{max}", but also a custom function.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1. Has no effect if \code{normalize} != \code{FALSE}.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @param zeta.type The function that must be used for the computation of zeta diversity. Default is "\code{exact}" for calling \code{Zeta.order.ex}. Use "\code{monte carlo}" for calling \code{Zeta.order.mc}.
#' @details The sites (plots or quadrates) are incrementally aggregated as nearest neighbouring groups of 4, 9, etc. sites, using a nested approach, starting from the lowest x and y, to increase the grain. The sites can be spatially contiguous or discontiguous, as long as they are regularly spaced (see Scheiner et al., 2011). If the total number of sites is not a multiple of \code{n[i]} x \code{n[i]}, the extra sites are discarded.
#' @return \code{Zeta.scale.regular} returns a list containing the following components:
#' @return \item{order}{The order of zeta.}
#' @return \item{n}{The vector of mapping grains: regularly spaced sites are grouped as \code{n[i]} x \code{n[i]} sites to generate data at a coarser grain.}
#' @return \item{values}{The zeta diversity values for each grain.}
#' @return \item{sd}{The standard deviation of zeta diversity for each grain.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}
#' @seealso \code{\link{Zeta.scale.min.dist}}, \code{\link{rescale.regular}}, \code{\link{rescale.min.dist}}
#' @examples
#' utils::data(bird.spec.fine)
#' xy.bird <- bird.spec.fine[1:400,1:2]
#' data.spec.bird <- bird.spec.fine[1:400,3:192]
#'
#' dev.new()
#' ##sam = 25 is used here for fast execution, but a higher value is advised
#' zeta.scale.reg <- Zeta.scale.regular(xy.bird, data.spec.bird, n = 1:3, order = 3,
#'     sam = 25, normalize = "Jaccard", zeta.type="monte carlo")
#' @export
Zeta.scale.regular <- function(xy, data.spec, n, order = 1, sam = 1000, method = "mean", rescale = FALSE, normalize = FALSE, plot = TRUE, zeta.type="exact"){

  if(nrow(data.spec) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }

  if (!inherits(data.spec,"data.frame")){
    stop(paste("Error: ",deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }

  n <- sort(n)
  n2 <- n

  zeta.val <- rep(NA, length(n))
  zeta.val.sd <- rep(NA, length(n))

  names <- c(names(xy), names(data.spec))

  if(n[1] == 1){
    if(zeta.type=="monte carlo"){
      zeta <- Zeta.order.mc(data.spec, order = order, sam = sam, rescale = rescale, normalize = normalize)
    }else if (zeta.type=="exact"){
      zeta <- Zeta.order.ex(data.spec, order = order, rescale = rescale)
    }else{
      stop("Error: unknown method for the computation of zeta")
    }
    zeta.val[1] <- zeta$zeta.val
    zeta.val.sd[1] <- zeta$zeta.val.sd
    n2 <- n[2:length(n)]
  }



  ##sort data according to the plots coordinates
  xy <- xy[order(xy[,1], xy[,2]), ]
  data.spec <- data.spec[order(xy[,1], xy[,2]), ]

  ##compute the scale dependence for the specified grains
  for (nn in 1:length(n2)){

    Ux <- sort(unique(xy[,1]))
    Uy <- sort(unique(xy[,2]))

    max.x <- Ux[length(Ux) - length(Ux)%%n2[nn]]
    max.y <- Uy[length(Uy) - length(Uy)%%n2[nn]]


    if(length(which(xy[,1]>max.x | xy[,2]>max.y))>0){
      xy2 <- xy[-which(xy[,1]>max.x | xy[,2]>max.y), ]
      data.spec2 <- data.spec[-which(xy[,1]>max.x | xy[,2]>max.y), ]
    }else{
      xy2 <- xy
      data.spec2 <- data.spec
    }

    data2 <- data.frame(stats::setNames(replicate(length(names), numeric(0), simplify = F), names))

    for (i in seq(1, length(unique(xy2$x)), n2[nn])){
      for (j in seq(1, length(unique(xy2$y)), n2[nn])){
        if(length(which(xy2$x %in% Ux[i:(i + n2[nn] - 1)] & xy2$y %in% Uy[j:(j + n2[nn] - 1)]))>0){
          temp.xy <- xy2[which(xy2$x %in% Ux[i:(i + n2[nn] - 1)] & xy2$y %in% Uy[j:(j + n2[nn] - 1)]), ]
          temp.xy <- apply(temp.xy, 2, get(method))

          temp.spec <- data.spec2[which(xy2$x %in% Ux[i:(i + n2[nn] - 1)] & xy2$y %in% Uy[j:(j + n2[nn] - 1)]), ]
          temp.spec <- (apply(temp.spec, 2, sum)>0) * 1

          ## add the vector to the new coarser dataset
          data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec)))
        }
      }
    }

    names(data2) <- names

    ##compute the zeta diversity of the new grain
    if(n[1] == 1){
      if(zeta.type=="monte carlo"){
        zeta <- Zeta.order.mc(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam, rescale = rescale, normalize = normalize)
      }else if (zeta.type=="exact"){
        zeta <- Zeta.order.ex(data2[3:(2 + dim(data.spec)[2])], order = order, rescale = rescale)
      }else{
        stop("Error: unknown method for the computation of zeta")
      }
      zeta.val[nn + 1] <- zeta$zeta.val
      zeta.val.sd[nn + 1] <- zeta$zeta.val.sd
    }else{
      if(zeta.type=="monte carlo"){
        zeta <- Zeta.order.mc(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam, rescale = rescale, normalize = normalize)
      }else if (zeta.type=="exact"){
        zeta <- Zeta.order.ex(data2[3:(2 + dim(data.spec)[2])], order = order, rescale = rescale)
      }else{
        stop("Error: unknown method for the computation of zeta")
      }
      zeta.val[nn] <- zeta$zeta.val
      zeta.val.sd[nn] <- zeta$zeta.val.sd
    }

  }

  if(plot == TRUE){
    graphics::plot(n^2, zeta.val, xlab = "Grain", ylab = paste("Zeta ", order, sep = ""), main = "Zeta-Scale Relationship", pch = 16)
    graphics::lines(n^2, zeta.val)
  }

  zeta.scale.reg <- list()
  zeta.scale.reg$order <- order
  zeta.scale.reg$n <- n
  zeta.scale.reg$values <- zeta.val
  zeta.scale.reg$sd <- zeta.val.sd

  return(zeta.scale.reg)

}



#' Zeta diversity scaling with sample grain dependency based on the minimum distance between sites
#'
#' Computes zeta diversity scaling with sample grain for a specific order (number of assemblages or sites), increasing grain by sequentially adding sites based on the minimum distance between them.
#' @param xy  Site-by-coordinate data frame, with sites as rows and coordinates as columns.
#' @param data.spec  Site-by-species presence-absence data frame, with sites as rows and species as columns.
#' @param m  Vector of mapping grains: \code{m[i]} sites are grouped together to generate data at a coarser grain.
#' @param order  Specific number of assemblages or sites at which zeta diversity is computed.
#' @param reorder  Number of times the sites are rearranged and grouped together for the computation of zeta (see Details).
#' @param shuffle Boolean value (TRUE or FALSE) indicating if the order of the sites must be randomised, which can have an impact on the outputs if some distances are equal.
#' @param sam  Number of samples for which the zeta diversity is computed.
#' @param method Name of a function (as a string) indicating how to combine the coordinates. It can be a basic R-function such as "\code{mean}" or "\code{max}", but also a custom function.
#' @param rescale Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1. Has no effect if \code{normalize} != \code{FALSE}.
#' @param normalize Indicates if the zeta values for each sample should be divided by the total number of species for this specific sample (\code{normalize = "Jaccard"}), by the average number of species per site for this specific sample (\code{normalize = "Sorensen"}), or by the minimum number of species in the sites of this specific sample \cr (\code{normalize = "Simpson"}). Default value is \code{FALSE}, indicating that no normalization is performed.
#' @param plot  Boolean value (TRUE or FALSE) indicating if the outputs must be plotted.
#' @param sd  Boolean value (TRUE or FALSE) indicating if the standard deviation must be plotted for each grain.
#' @param distance.type Method to compute distance. Default is "\code{Euclidean}", for Euclidean distance. The other options are (i) "\code{ortho}" for orthodromic distance, if xy correspond to longitudes and latitudes (orthodromic distance is computed using the \code{geodist} function from package \code{geodist}); and (ii) "\code{custom}", in which case the user must provide a distance matrix for \code{dist.custom}.
#' @param dist.custom Distance matrix provided by the user when \code{distance.type} = \code{"custom"}.
#' @param zeta.type The function that must be used for the computation of zeta diversity. Default is "\code{exact}" for calling \code{Zeta.order.ex}. Use "\code{monte carlo}" for calling \code{Zeta.order.mc}.
#' @details The nearest neighbouring sites (plots, quadrates, or areas of varying shapes) are grouped as spatial clusters of 2, 3, 4, etc. sites, based on the minimum distance between them. Since the procedure is based on the relative distance between sites, the site order can have an impact on the output. The procedure is therefore performed 'reorder' times, for which sites are randomly reordered each time, and the mean zeta is computed. This function is suitable for both regularly and irregularly spaced sites, contiguous or non contiguous (\emph{sensu} Scheiner et al., 2011). For regularly spaced sites, the use of \code{\link{Zeta.scale.regular}} is recommended.
#' @return \code{zeta.scale.min.dist} returns a list containing the following components:
#' @return \item{order}{The order of zeta.}
#' @return \item{m}{The vector of mapping grains: m[i] sites are grouped together to generate data at a coarser grain.}
#' @return \item{values}{A matrix containing the zeta diversity values over the '\code{reorder}' computations, for each grain.}
#' @return \item{sd}{A matrix containing the standard deviation of zeta diversity over the '\code{reorder}' computations, for each grain.}
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}},
#' @seealso \code{\link{Zeta.scale.regular}}, \code{\link{rescale.regular}}
#' @examples
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' dev.new()
#' zeta.scale.irreg.species <- Zeta.scale.min.dist(xy.marion, data.spec.marion, m = 1:3,
#'     order = 3, reorder = 3, sam = 50, normalize = "Jaccard")
#'
#' @export
Zeta.scale.min.dist <- function(xy, data.spec, m, order = 1, reorder = 100, shuffle = TRUE, sam = 1000, method = "mean", rescale = FALSE, normalize = FALSE, plot = TRUE, sd = TRUE, distance.type = "Euclidean", dist.custom = NULL, zeta.type="exact"){

  if(nrow(data.spec) != nrow(xy)){
    stop("Error: data.spec and xy must have the same number of rows.")
  }

  if (!inherits(data.spec,"data.frame")){
    stop(paste("Error: ",deparse(substitute(data.spec)), " is a ", class(data.spec), ". It must be a data frame.", sep = ""))
  }
  if(order>dim(data.spec)[1]){
    stop("Error: wrong value for \"order\": it must be equal or lower than the number of sites.")
  }
  if(!is.null(dist.custom)){
    if(!isSymmetric(dist.custom)){
      stop("Error: distance matrix is not symmetrical")
    }
  }

  m <- sort(m)
  m2 <- m

  zeta.val <- matrix(NA, reorder, length(m))
  zeta.val.sd <- matrix(NA, reorder, length(m))

  names <- c(names(xy), names(data.spec))

  if (m[1] == 1){
    if(zeta.type=="monte carlo"){
      zeta <- Zeta.order.mc(data.spec, order = order, sam = sam, rescale = rescale, normalize = normalize)
    }else if (zeta.type=="exact"){
      zeta <- Zeta.order.ex(data.spec, order = order, rescale = rescale)
    }else{
      stop("Error: unknown method for the computation of zeta")
    }
    zeta.val[, 1] <- rep(zeta$zeta.val, reorder)
    zeta.val.sd[, 1] <- rep(zeta$zeta.val.sd, reorder)
    m2 <- m[2:length(m)]
  }

  ##compute the scale dependence for the specified grains
  for (nn in 1:length(m2)){

    zeta.scale.temp <- rep(NA, reorder)

    ##repeat 'reorder times'
    for(reord in 1:reorder){

      ##randomize the plot orders
      if (shuffle == TRUE){
        xy <- xy[sample(nrow(xy)), ]
        data.spec <- data.spec[sample(nrow(xy)), ]
      }

      #pairwise distance
      if(distance.type == "Euclidean"){
        D <- as.matrix(stats::dist(xy))
      }else if(distance.type == "ortho"){
        D <- as.matrix(.gdist_matrix(stats::dist(xy)))
      }else if(distance.type == "custom"){
        if(is.null(dist.custom))
        {
          stop("Error: a distance matrix must be provided if distance.type = 'custom'")
        }
        D <- dist.custom
      }else{
        stop("Error: invalid distance type")
      }

      D[which(D == 0)] <- NA

      data2 <- data.frame(stats::setNames(replicate(length(names), numeric(0), simplify = F), names))

      i.plot <- 1
      for(i in 1:floor(dim(xy)[1] / m2[nn])){

        while(length(which(!is.na(D[, i.plot]))) == 0){
          i.plot <- i.plot + 1
        }

        ##select the plots in the group
        temp.xy <- xy[c(i.plot, order(D[, i.plot])[1:(m2[nn] - 1)]), ]
        temp.spec <- data.spec[c(i.plot, order(D[, i.plot])[1:(m2[nn] - 1)]), ]

        ##compute the mean coordinates, the unions of species, and the mean of environmental variables
        temp.xy <- apply(temp.xy, 2, get(method))
        temp.spec <- (apply(temp.spec, 2, sum)>0) * 1

        ## add the vector to the new coarser dataset
        data2 <- rbind(data2, as.vector(c(temp.xy, temp.spec)))

        ##set the column corresponding to the plots to NA in the distance matrix to avoid using them in the following
        DD <- D
        D[, order(DD[, i.plot])[1:(m2[nn] - 1)]] <- NA
        D[order(DD[, i.plot])[1:(m2[nn] - 1)], ] <- NA
        D[, i.plot] <- NA
        D[i.plot, ] <- NA

      }

      names(data2) <- names

      ##compute the zeta diversity of the new grain
      if (m[1] == 1){
        if(zeta.type=="monte carlo"){
          zeta <- Zeta.order.mc(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam, rescale = rescale, normalize = normalize)
        }else if (zeta.type=="exact"){
          zeta <- Zeta.order.ex(data2[3:(2 + dim(data.spec)[2])], order = order, rescale = rescale)
        }else{
          stop("Error: unknown method for the computation of zeta")
        }
        zeta.val[reord, nn + 1] <- zeta$zeta.val
        zeta.val.sd[reord, nn + 1] <- zeta$zeta.val.sd
      }else{
        if(zeta.type=="monte carlo"){
          zeta <- Zeta.order.mc(data2[3:(2 + dim(data.spec)[2])], order = order, sam = sam, rescale = rescale, normalize = normalize)
        }else if (zeta.type=="exact"){
          zeta <- Zeta.order.ex(data2[3:(2 + dim(data.spec)[2])], order = order, rescale = rescale)
        }else{
          stop("Error: unknown method for the computation of zeta")
        }
        zeta.val[reord, nn] <- zeta$zeta.val
        zeta.val.sd[reord, nn] <- zeta$zeta.val.sd
      }

    }
  }


  if(plot == TRUE){
    graphics::plot(m, apply(zeta.val, 2, mean), ylim = c(apply(zeta.val, 2, mean)[1] - apply(zeta.val, 2, stats::sd)[1], apply(zeta.val, 2, mean)[length(m)] + apply(zeta.val, 2, stats::sd)[length(m)]), xlab = "Grain", ylab = paste("Zeta ", order, sep = ""), main = "Zeta-Scale Relationship", pch = 16)
    graphics::lines(m, apply(zeta.val, 2, mean))
    if(sd == TRUE){
      for(i in m){
        suppressWarnings(graphics::arrows(m[i], apply(zeta.val, 2, mean)[i], m[i], apply(zeta.val, 2, mean)[i] + apply(zeta.val, 2, stats::sd)[i], angle = 90, length = 0.1))
        suppressWarnings(graphics::arrows(m[i], apply(zeta.val, 2, mean)[i], m[i], apply(zeta.val, 2, mean)[i] - apply(zeta.val, 2, stats::sd)[i], angle = 90, length = 0.1))
      }
    }
  }

  zeta.scale.irreg <- list()
  zeta.scale.irreg$order <- order
  zeta.scale.irreg$m <- m
  zeta.scale.irreg$values <- zeta.val
  zeta.scale.irreg$sd <- zeta.val.sd

  return(zeta.scale.irreg)

}





#' Plotting of zeta diversity scaling with sample grain using hierarchical increases in grain size
#'
#' Plots the output of the function \code{Zeta.scale.regular}.
#' @param zeta.scale.reg  A list generated by the function \code{Zeta.scale.regular}.
#' @param size.init initial Size of the plots before aggregation.
#' @param add Boolean value indicating if the graph must be plotted in a new graphics device or added to the active one.
#' @param ylim Numeric vectors of length 2, giving the range of y values.
#' @param col String indicating the color of the graph.
#' @return A plot of the zeta diversity scaling with the mapping grain \code{n} x \code{n} (the number of sites combined to generate data at a coarser grain) on the x-axis and the value of zeta on the y-axis.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}},
#' @seealso \code{\link{Zeta.scale.regular}}, \code{\link{Zeta.scale.min.dist}}, \code{\link{rescale.regular}},
#' @seealso \code{\link{Plot.zeta.scale.min.dist}}
#' @examples
#' utils::data(bird.spec.fine)
#' xy.bird <- bird.spec.fine[1:400,1:2]
#' data.spec.bird <- bird.spec.fine[1:400,3:192]
#'
#' ##sam = 25 is used here for fast execution, but a higher value is advised
#' zeta.scale.reg <- Zeta.scale.regular(xy.bird, data.spec.bird, n = 1:3, order = 3,
#'     sam = 25, normalize = "Jaccard",plot=FALSE)
#' dev.new()
#' Plot.zeta.scale.regular(zeta.scale.reg)
#' @export
Plot.zeta.scale.regular <- function(zeta.scale.reg, size.init = 1, add = FALSE, ylim = NULL, col = "black"){

  if(is.null(ylim)){
    ylim=c(0,max(zeta.scale.reg$values))
  }
  if(add == FALSE){
    graphics::plot((size.init*zeta.scale.reg$n)^2, zeta.scale.reg$values, xlab = "Grain", ylab = paste("Zeta ", zeta.scale.reg$order, sep = ""), main = "Hierarchical scaling of Zeta", pch = 16,ylim=ylim,col = col)
  }else{
    graphics::points((size.init*zeta.scale.reg$n)^2, zeta.scale.reg$values, pch = 16,ylim=ylim,col = col)
  }
  graphics::lines((size.init*zeta.scale.reg$n)^2, zeta.scale.reg$values,col = col)

}





#' Plotting of zeta diversity scaling with sample grain dependency based on the minimum distance between sites
#'
#' Plots the output of the function \code{Zeta.scale.min.dist}.
#' @param zeta.scale.irreg  A list generated by the function \code{Zeta.scale.min.dist}.
#' @param size.init Initial size of the plots before aggregation.
#' @param add Boolean value indicating if the graph must be plotted in a new graphics device or added to the active one.
#' @param ylim Numeric vectors of length 2, giving the range of y values.
#' @param col String indicating the color of the graph.
#' @return A plot of the zeta diversity scaling with the mapping grain m (the number of sites combined to generate data at a coarser grain) on the x-axis and the value of zeta on the y-axis.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}},
#' @seealso \code{\link{Zeta.scale.min.dist}}, \code{\link{rescale.regular}}, \code{\link{Zeta.scale.regular}}, \code{\link{rescale.regular}},
#' @seealso \code{\link{Plot.zeta.scale.regular}}
#' @examples
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#'
#' zeta.scale.irreg.species <- Zeta.scale.min.dist(xy.marion, data.spec.marion, m = 1:3,
#'     order = 3, reorder = 3, sam = 50, normalize = "Jaccard",plot=FALSE)
#' dev.new()
#' Plot.zeta.scale.min.dist(zeta.scale.irreg.species)
#' @export
Plot.zeta.scale.min.dist <- function(zeta.scale.irreg, size.init = 1, add = FALSE, ylim = NULL, col = "black"){

  if(is.null(ylim)){
    ylim = c(apply(zeta.scale.irreg$values, 2, mean)[1] - apply(zeta.scale.irreg$values, 2, stats::sd)[1], apply(zeta.scale.irreg$values, 2, mean)[length(zeta.scale.irreg$m)] + apply(zeta.scale.irreg$values, 2, stats::sd)[length(zeta.scale.irreg$m)])
  }

  if(add == FALSE){
    graphics::plot(size.init*zeta.scale.irreg$m, apply(zeta.scale.irreg$values, 2, mean), ylim = ylim, xlab = "Grain", ylab = paste("Zeta ", zeta.scale.irreg$order, sep = ""), main = "Hierarchical scaling of Zeta", pch = 16, col = col)
  }else{
    graphics::points(size.init*zeta.scale.irreg$m, apply(zeta.scale.irreg$values, 2, mean), pch = 16, col = col)
  }
  graphics::lines(size.init*zeta.scale.irreg$m, apply(zeta.scale.irreg$values, 2, mean), col = col)
  for(i in zeta.scale.irreg$m){
    suppressWarnings(graphics::arrows(size.init*zeta.scale.irreg$m[i], apply(zeta.scale.irreg$values, 2, mean)[i], size.init*zeta.scale.irreg$m[i], apply(zeta.scale.irreg$values, 2, mean)[i] + apply(zeta.scale.irreg$values, 2, stats::sd)[i], angle = 90, length = 0.1))
    suppressWarnings(graphics::arrows(size.init*zeta.scale.irreg$m[i], apply(zeta.scale.irreg$values, 2, mean)[i], size.init*zeta.scale.irreg$m[i], apply(zeta.scale.irreg$values, 2, mean)[i] - apply(zeta.scale.irreg$values, 2, stats::sd)[i], angle = 90, length = 0.1))
  }
}





#' Variation partitioning for zeta diversity
#'
#' Variation partitioning of zeta diversity for a specific order (number of assemblages or sites) over distance and environmental variables.
#' @param msgdm.mod An object return by function \code{\link{Zeta.msgdm}}.
#' @param num.part Number of partitions of zeta diversity. Can be 2 or 3.
#' @param reg.type Type of regression for the multi-site generalised dissimilarity modelling. Options are "glm" for generalised linear models, "ngls" for negative linear models, "gam" for generalised additive models, "scam" for shape constrained additive models, and "ispline" for I-spline models, as recommended in generalised dissimilarity modelling by Ferrier \emph{et al}. (2007).
#' @param family A description of the error distribution and link function to be used in the \code{glm}, \code{gam} and \code{scam} models (see \code{\link[stats]{family}} for details of family functions).
#' @param method.glm Method used in fitting the generalised linear model. The default method \cr "glm.fit.cons" is an adaptation of method \code{glm.fit2} from package \code{glm2} using a negative least squares regression in the reweighted least squares. Another option is "glm.fit2", which corresponds to method \code{glm.fit2}; see help documentation for \code{glm.fit2} in package \code{glm}.
#' @param cons type of constraint in the glm if \code{method.glm = "glm.fit.cons"}. Default is -1 for negative coefficients on the predictors. The other option is 1 for positive coefficients on the predictors.
#' @param cons.inter type of constraint for the intercept. Default is 1 for positive intercept, suitable for Gaussian family. The other option is -1 for negative intercept, suitable for binomial family.
#' @param kn Number of knots in the GAM and SCAM. Default is -1 for determining kn automatically using Generalized Cross-validation.
#' @param bs A two-letter character string indicating the (penalized) smoothing basis to use in the scam model. Default is "\code{mpd}" for monotonic decreasing splines. see \code{\link[mgcv]{smooth.terms}} for an overview of what is available.
#' @return \code{Zeta.varpart} returns a data frame with one column containing the variation explained by each component \code{a} (the variation explained by distance alone), \code{b} (the variation explained by either distance or the environment), \code{c} (the variation explained by the environment alone) and \code{d} (the unexplained variation).
#' @details Note that, for a given regression, the variation explained is computed as 1-(RSS/TSS)*(v-1)/(v-p-1), where RSS is the residual sum of squares and TSS is the total sum of squares, v is the number of variables used in the regression (which is greater than the original number of variables for I-splines) and p is the number of samples. 1-(RSS/TSS) corresponds to the classical R-squared for linear regression only, and results for non-linear regressions should be interpreted with caution.
#' @details The environmental variables can be numeric or factorial, and \code{order} must be greater than 1.
#' @details For numeric variables, the pairwise difference between sites is computed and combined according to \code{method}. For factorial variables, the distance corresponds to the number of unique values over the number of assemblages of sites specified by \code{order}.
#' @details Zeta is regressed against the differences of values of the environmental variables divided by the maximum difference for each variable, to be rescaled between 0 and 1. If \code{!is.null(xy)}, distances between sites are also divided by the maximum distance.
#' @references Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.
#' @references Borcard, D., Legendre, P. & Drapeau, P. (1992). Partialling out the spatial component of ecological variation. \emph{Ecology} 73, 1045-1055.
#' @references Legendre, P. &  Legendre, L.F. (2012). \emph{Numerical ecology}, 3rd English edition. Elsevier Science BV, Amsterdam.
#' @seealso \code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}}, \code{\link{Zeta.msgdm}}, \code{\link{pie.neg}}
#' @import scam
#' @examples
#' utils::data(bird.spec.coarse)
#' xy.bird <- bird.spec.coarse[,1:2]
#' data.spec.bird <- bird.spec.coarse[,3:193]
#' utils::data(bird.env.coarse)
#' data.env.bird <- bird.env.coarse[,3:9]
#'
#' zeta.bird <- Zeta.msgdm(data.spec.bird, data.env.bird, xy.bird, sam = 100, order = 3)
#' zeta.varpart.bird <- Zeta.varpart(zeta.bird, method.glm = "glm.fit2")
#' zeta.varpart.bird
#' dev.new()
#' pie.neg(zeta.varpart.bird[4:7,1], density = c(4, 0, 8, -1),
#'     angle = c(90, 0, 0, 0),
#'     labels = c("distance", "undistinguishable", "environment", "unexplained"),
#'     radius = 0.9)
#'
#' ##########
#'
#' utils::data(Marion.species)
#' xy.marion <- Marion.species[,1:2]
#' data.spec.marion <- Marion.species[,3:33]
#' utils::data(Marion.env)
#' data.env.marion <- Marion.env[3:4]
#'
#' zeta.marion <- Zeta.msgdm(data.spec.marion, data.env.marion, xy.marion, sam = 100,
#'     order = 3, normalize = "Jaccard")
#' zeta.varpart.marion <- Zeta.varpart(zeta.marion, method.glm = "glm.fit2")
#' zeta.varpart.marion
#' dev.new()
#' pie.neg(zeta.varpart.marion[4:7,1], density = c(4, 0, 8, -1),
#'     angle = c(90, 0, 0, 0),
#'     labels = c("distance", "undistinguishable", "environment", "unexplained"),
#'     radius = 0.9)
#'
#' @export
Zeta.varpart <- function(msgdm.mod, num.part = 2, reg.type = "glm", family = stats::gaussian(), method.glm = "glm.fit.cons", cons = -1, cons.inter = 1, kn = -1, bs="mpd"){

  zeta.val <- msgdm.mod$val
  data.tot <- msgdm.mod$predictors

  if(reg.type == "ispline"){
    if(num.part==2){
      data.var <- data.tot[,1:(ncol(data.tot)-3)]
      distance <- data.tot[,(ncol(data.tot)-2):ncol(data.tot)]
    }else{
      data.var <- data.tot[,1:(ncol(data.tot)-6)]
      sp.pred <- data.tot[,(ncol(data.tot)-5):(ncol(data.tot)-3)]
      distance <- data.tot[,(ncol(data.tot)-2):ncol(data.tot)]
      data.var.sp <- data.tot[,1:(ncol(data.tot)-3)]
      data.var.dist <- data.tot[,c(1:(ncol(data.tot)-6),(ncol(data.tot)-2):ncol(data.tot))]
      data.sp.dist <- data.tot[,(ncol(data.tot)-5):ncol(data.tot)]
    }
  }else{
    if(num.part==2){
      data.var <- data.tot[,1:(ncol(data.tot)-1),drop=FALSE]
      distance <- data.tot[,ncol(data.tot)]
    }else{
      data.var <- data.tot[,1:(ncol(data.tot)-2),drop=FALSE]
      sp.pred <- data.tot[,ncol(data.tot)-1,drop=FALSE]
      distance <- data.tot[,ncol(data.tot)]
      data.var.sp <- data.tot[,1:(ncol(data.tot)-1),drop=FALSE]
      data.var.dist <- data.tot[,c(1:(ncol(data.tot)-2),ncol(data.tot)),drop=FALSE]
      data.sp.dist <- data.tot[,(ncol(data.tot)-1):ncol(data.tot),drop=FALSE]
    }
  }

  if(num.part == 2){
    if(reg.type == "glm"){
      if(method.glm == "glm.fit.cons"){
        abc <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ ., data = data.tot, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        ab <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ distance, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        bc <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ ., data = data.var, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
      }else{
        abc <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ ., data = data.tot, family = family, method = method.glm))$r.squared,0)
        bc <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ ., data = data.var, family = family, method = method.glm))$r.squared,0)
        ab <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ distance, family = family, method = method.glm))$r.squared,0)
      }
      b <- max(ab+bc-abc,0)
      a <- max(ab-b,0)
      c <- max(bc-b,0)
    }else if(reg.type == "gam"){
      ##create formula to be used in gam
      xnam <- names(data.tot)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      abc <- max(summary(mgcv::gam(fm, data = data.tot, family = family))$r.sq,0)
      ab <- max(summary(mgcv::gam(zeta.val ~ s(distance, k = kn), family = family,))$r.sq,0)
      xnam <- names(data.var)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      bc <- max(summary(mgcv::gam(fm, data = data.var, family = family))$r.sq,0)
      b <- max(ab+bc-abc,0)
      a <- max(ab-b,0)
      c <- max(bc-b,0)
    }else if(reg.type == "scam"){
      ##create formula to be used in scam
      xnam <- names(data.tot)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      abc <- max(summary(scam::scam(fm, data = data.tot, family = family))$r.sq,0)
      ab <- max(summary(scam::scam(zeta.val ~ s(distance, k=kn, bs=bs), data=as.data.frame(distance), family = family))$r.sq,0)
      xnam <- names(data.var)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      bc <- max(summary(scam::scam(fm, data = data.var, family = family))$r.sq,0)
      b <- max(ab+bc-abc,0)
      a <- max(ab-b,0)
      c <- max(bc-b,0)
    }else if(reg.type == "ngls"){
      TSS <- sum((zeta.val - mean(zeta.val))^2)
      toto <- rep(NA,ncol(data.tot))
      for(i in 1:ncol(data.tot)){
        toto[i] <- (paste("b",i,"*",names(data.tot)[i],sep=""))
      }
      X <- rep(1,nrow(data.tot))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))

      start <- data.frame(matrix(NA,1,ncol(data.tot)+1))
      start.names <- rep(NA,ncol(data.tot))
      for(i in 1:(ncol(data.tot)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1

      m  <- stats::nls(fm, data=data.tot,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(data.tot))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      abc <- 1-(1-R2)*(nrow(data.tot))/(nrow(data.tot)-ncol(data.tot)-1)

      start <- data.frame(matrix(NA,1,1+1))
      start.names <- rep(NA,1)
      for(i in 1:(1+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(zeta.val ~ b0 * X + b1 * distance, start=start,algorithm="port",upper=c(Inf,0))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      ab <- 1-(1-R2)*(nrow(data.tot))/(nrow(data.tot)-1-1)

      start <- data.frame(matrix(NA,1,ncol(data.var)+1))
      start.names <- rep(NA,ncol(data.var))
      for(i in 1:(ncol(data.var)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      toto <- rep(NA,ncol(data.var))
      for(i in 1:ncol(data.var)){
        toto[i] <- (paste("b",i,"*",names(data.var)[i],sep=""))
      }
      X <- rep(1,nrow(data.var))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      m  <- stats::nls(fm, data=data.tot,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(data.var))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      bc <- 1-(1-R2)*(nrow(data.tot))/(nrow(data.tot)-ncol(data.tot)-1)

      abc <- max(abc,0)
      ab <- max(ab,0)
      bc <- max(bc,0)
      b <- max(ab+bc-abc,0)
      a <- max(ab-b,0)
      c <- max(bc-b,0)
    }else if(reg.type == "ispline"){
      toto <- glm.cons(zeta.val ~ ., data = data.tot, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      abc <- 1 - (toto$deviance/toto$null.deviance)*(nrow(data.tot)-1)/(nrow(data.tot)-ncol(data.tot)-1)
      toto <- glm.cons(zeta.val ~ ., data = distance, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      ab <- 1 - (toto$deviance/toto$null.deviance)*(nrow(distance)-1)/(nrow(distance)-ncol(distance)-1)
      toto <- glm.cons(zeta.val ~ ., data = data.var, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      bc <- 1 - (toto$deviance/toto$null.deviance)*(nrow(data.var)-1)/(nrow(data.var)-ncol(data.var)-1)
      b <- max(ab+bc-abc,0)
      a <- max(ab-b,0)
      c <- max(bc-b,0)
    }else{
      stop("Error: Unknwon regression type.")
    }
  }else{
    if(reg.type == "glm"){
      if(method.glm == "glm.fit.cons"){
        ABC <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ ., data = data.tot, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        AB <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ data.var.sp, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        AC <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ data.var.dist, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        BC <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ data.sp.dist, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        A <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ ., data = data.var, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        B <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ ., data = sp.pred, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        C <- max(vegan::RsquareAdj(glm.cons(zeta.val ~ ., data = distance, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
      }else{
        ABC <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ ., data = data.tot, family = family, method = method.glm, cons = cons, cons.inter = cons.inter))$r.squared,0)
        AB <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ data.var.sp, family = family, method = method.glm))$r.squared,0)
        AC <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ data.var.dist, family = family, method = method.glm))$r.squared,0)
        BC <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ data.sp.dist, family = family, method = method.glm))$r.squared,0)
        A <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ ., data = data.var, family = family, method = method.glm))$r.squared,0)
        B <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ ., data = sp.pred, family = family, method = method.glm))$r.squared,0)
        C <- max(vegan::RsquareAdj(glm2::glm2(zeta.val ~ ., data = distance, family = family, method = method.glm))$r.squared,0)
      }
      g <- max(A+B+C-AB-AC-BC+ABC,0)
      d <- max(A+B-AB-g,0)
      e <- max(B+C-BC,0)
      f <- max(A+C-AC,0)
      a <- max(A-d-f-g,0)
      b <- max(B-d-e-g,0)
      c <- max(C-e-f-g,0)
    }else if(reg.type == "gam"){
      ##create formula to be used in gam
      xnam <- names(data.tot)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      ABC <- max(summary(mgcv::gam(fm, data = data.tot, family = family))$r.sq,0)
      xnam <- names(data.var.sp)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      AB <- max(summary(mgcv::gam(fm, data = data.var.sp, family = family))$r.sq,0)
      xnam <- names(data.var.dist)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      AC <- max(summary(mgcv::gam(fm, data = data.var.dist, family = family))$r.sq,0)
      xnam <- names(data.sp.dist)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      BC <- max(summary(mgcv::gam(fm, data = data.sp.dist, family = family))$r.sq,0)
      xnam <- names(data.var)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      A <- max(summary(mgcv::gam(fm, data = data.var, family = family))$r.sq,0)
      xnam <- names(sp.pred)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      B <- max(summary(mgcv::gam(fm, data = sp.pred, family = family))$r.sq,0)
      xnam <- names(distance)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,") + s(",sep=""),sep=""), ", k = ",kn,")",sep=""))
      C <- max(summary(mgcv::gam(fm, data = distance, family = family))$r.sq,0)
      g <- max(A+B+C-AB-AC-BC+ABC,0)
      d <- max(A+B-AB-g,0)
      e <- max(B+C-BC,0)
      f <- max(A+C-AC,0)
      a <- max(A-d-f-g,0)
      b <- max(B-d-e-g,0)
      c <- max(C-e-f-g,0)
    }else if(reg.type == "scam"){
      ##create formula to be used in scam
      xnam <- names(data.tot)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      ABC <- max(summary(scam::scam(fm, data = data.tot, family = family))$r.sq,0)
      xnam <- names(data.var.sp)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      AB <- max(summary(scam::scam(fm, data = data.var.sp, family = family))$r.sq,0)
      xnam <- names(data.var.dist)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      AC <- max(summary(scam::scam(fm, data = data.var.dist, family = family))$r.sq,0)
      xnam <- names(data.sp.dist)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      BC <- max(summary(scam::scam(fm, data = data.sp.dist, family = family))$r.sq,0)
      xnam <- names(data.var)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      A <- max(summary(scam::scam(fm, data = data.var, family = family))$r.sq,0)
      xnam <- names(sp.pred)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      B <- max(summary(scam::scam(fm, data = sp.pred, family = family))$r.sq,0)
      xnam <- names(distance)
      fm <- stats::as.formula(paste("zeta.val ~ s(", paste(xnam, collapse = paste(", k = ",kn,", bs = '", bs ,"') + s(",sep=""),sep=""),", k = ",kn, ",bs='", bs, "')",sep=""))
      C <- max(summary(scam::scam(fm, data = distance, family = family))$r.sq,0)
      g <- max(A+B+C-AB-AC-BC+ABC,0)
      d <- max(A+B-AB-g,0)
      e <- max(B+C-BC,0)
      f <- max(A+C-AC,0)
      a <- max(A-d-f-g,0)
      b <- max(B-d-e-g,0)
      c <- max(C-e-f-g,0)
    }else if(reg.type == "ngls"){
      TSS <- sum((zeta.val - mean(zeta.val))^2)

      toto <- rep(NA,ncol(data.tot))
      for(i in 1:ncol(data.tot)){
        toto[i] <- (paste("b",i,"*",names(data.tot)[i],sep=""))
      }
      X <- rep(1,nrow(data.tot))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      start <- data.frame(matrix(NA,1,ncol(data.tot)+1))
      start.names <- rep(NA,ncol(data.tot))
      for(i in 1:(ncol(data.tot)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(fm, data=data.tot,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(data.tot))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      ABC <- 1-(1-R2)*(nrow(data.tot))/(nrow(data.tot)-ncol(data.tot)-1)

      toto <- rep(NA,ncol(data.var.sp))
      for(i in 1:ncol(data.var.sp)){
        toto[i] <- (paste("b",i,"*",names(data.var.sp)[i],sep=""))
      }
      X <- rep(1,nrow(data.var.sp))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      start <- data.frame(matrix(NA,1,ncol(data.var.sp)+1))
      start.names <- rep(NA,ncol(data.var.sp))
      for(i in 1:(ncol(data.var.sp)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(fm, data=data.var.sp,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(data.var.sp))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      AB <- 1-(1-R2)*(nrow(data.var.sp))/(nrow(data.var.sp)-ncol(data.var.sp)-1)

      toto <- rep(NA,ncol(data.var.dist))
      for(i in 1:ncol(data.var.dist)){
        toto[i] <- (paste("b",i,"*",names(data.var.dist)[i],sep=""))
      }
      X <- rep(1,nrow(data.var.dist))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      start <- data.frame(matrix(NA,1,ncol(data.var.dist)+1))
      start.names <- rep(NA,ncol(data.var.dist))
      for(i in 1:(ncol(data.var.dist)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(fm, data=data.var.dist,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(data.var.dist))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      AC <- 1-(1-R2)*(nrow(data.var.dist))/(nrow(data.var.dist)-ncol(data.var.dist)-1)

      toto <- rep(NA,ncol(data.sp.dist))
      for(i in 1:ncol(data.sp.dist)){
        toto[i] <- (paste("b",i,"*",names(data.sp.dist)[i],sep=""))
      }
      X <- rep(1,nrow(data.sp.dist))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      start <- data.frame(matrix(NA,1,ncol(data.sp.dist)+1))
      start.names <- rep(NA,ncol(data.sp.dist))
      for(i in 1:(ncol(data.sp.dist)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(fm, data=data.sp.dist,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(data.sp.dist))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      BC <- 1-(1-R2)*(nrow(data.sp.dist))/(nrow(data.sp.dist)-ncol(data.sp.dist)-1)

      toto <- rep(NA,ncol(data.var))
      for(i in 1:ncol(data.var)){
        toto[i] <- (paste("b",i,"*",names(data.var)[i],sep=""))
      }
      X <- rep(1,nrow(data.var))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      start <- data.frame(matrix(NA,1,ncol(data.var)+1))
      start.names <- rep(NA,ncol(data.var))
      for(i in 1:(ncol(data.var)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(fm, data=data.var,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(data.var))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      A <- 1-(1-R2)*(nrow(data.var))/(nrow(data.var)-ncol(data.var)-1)

      toto <- rep(NA,ncol(sp.pred))
      for(i in 1:ncol(sp.pred)){
        toto[i] <- (paste("b",i,"*",names(sp.pred)[i],sep=""))
      }
      X <- rep(1,nrow(sp.pred))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      start <- data.frame(matrix(NA,1,ncol(sp.pred)+1))
      start.names <- rep(NA,ncol(sp.pred))
      for(i in 1:(ncol(sp.pred)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(fm, data=sp.pred,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(sp.pred))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      B <- 1-(1-R2)*(nrow(sp.pred))/(nrow(sp.pred)-ncol(sp.pred)-1)

      toto <- rep(NA,ncol(distance))
      for(i in 1:ncol(distance)){
        toto[i] <- (paste("b",i,"*",names(distance)[i],sep=""))
      }
      X <- rep(1,nrow(distance))
      fm <- stats::as.formula(paste("zeta.val ~ b0*X +",paste(toto,collapse=" + ",sep=""),sep=""))
      start <- data.frame(matrix(NA,1,ncol(distance)+1))
      start.names <- rep(NA,ncol(distance))
      for(i in 1:(ncol(distance)+1)){
        start[1,i] <- -1
        start.names[i] <- (paste("b",i-1,sep=""))
      }
      names(start) <- start.names
      start <- as.list(start)
      start[[1]] <- 1
      m  <- stats::nls(fm, data=distance,start=start,algorithm="port",upper=c(Inf,rep(0,ncol(distance))))
      RSS.p <- sum(stats::residuals(m)^2)
      R2 <- 1 - (RSS.p/TSS)
      C <- 1-(1-R2)*(nrow(distance))/(nrow(distance)-ncol(distance)-1)

      ABC <- max(ABC,0)
      AB <- max(AB,0)
      AC <- max(AC,0)
      BC <- max(BC,0)

      g <- max(A+B+C-AB-AC-BC+ABC,0)
      d <- max(A+B-AB-g,0)
      e <- max(B+C-BC,0)
      f <- max(A+C-AC,0)
      a <- max(A-d-f-g,0)
      b <- max(B-d-e-g,0)
      c <- max(C-e-f-g,0)
    }else if(reg.type == "ispline"){
      toto <- glm.cons(zeta.val ~ ., data = data.tot, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      ABC <- 1 - (toto$deviance/toto$null.deviance)*(nrow(data.tot)-1)/(nrow(data.tot)-ncol(data.tot)-1)
      toto <- glm.cons(zeta.val ~ ., data = data.var.sp, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      AB <- 1 - (toto$deviance/toto$null.deviance)*(nrow(data.var.sp)-1)/(nrow(data.var.sp)-ncol(distance)-1)
      toto <- glm.cons(zeta.val ~ ., data = data.var.dist, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      BC <- 1 - (toto$deviance/toto$null.deviance)*(nrow(data.var.dist)-1)/(nrow(data.var.dist)-ncol(data.var.dist)-1)
      toto <- glm.cons(zeta.val ~ ., data = data.var, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      A <- 1 - (toto$deviance/toto$null.deviance)*(nrow(data.var)-1)/(nrow(data.var)-ncol(data.var)-1)
      toto <- glm.cons(zeta.val ~ ., data = sp.pred, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      B <- 1 - (toto$deviance/toto$null.deviance)*(nrow(sp.pred)-1)/(nrow(sp.pred)-ncol(sp.pred)-1)
      toto <- glm.cons(zeta.val ~ ., data = distance, family = family, method = method.glm, cons = cons, cons.inter = cons.inter)
      C <- 1 - (toto$deviance/toto$null.deviance)*(nrow(distance)-1)/(nrow(distance)-ncol(distance)-1)
      g <- max(A+B+C-AB-AC-BC+ABC,0)
      d <- max(A+B-AB-g,0)
      e <- max(B+C-BC,0)
      f <- max(A+C-AC,0)
      a <- max(A-d-f-g,0)
      b <- max(B-d-e-g,0)
      c <- max(C-e-f-g,0)
    }else{
      stop("Error: Unknwon regression type.")
    }
  }

  if(num.part==2){
    zeta.varpart <- data.frame(c(abc,ab,bc,a,b,c,1-abc))
    row.names(zeta.varpart) <- c("[abc]","[ab]","[bc]","[a]","[b]","[c]","[d]")
    names(zeta.varpart) <- "Adjusted Rsq"
  }else{
    zeta.varpart <- data.frame(c(ABC,AB,AC,BC,A,B,C,a,b,c,d,e,f,g,1-ABC))
    row.names(zeta.varpart) <- c("[abcdefg]","[abdefg]","[acefg]","[bcefg]","[adfg]","[bdeg]","[cefg]","[a]","[b]","[c]","[d]","[e]","[f]","[g]","[h]")
    names(zeta.varpart) <- "Adjusted Rsq"
  }

  return(zeta.varpart)

}




#' Pie Chart, considering negative values as zeros
#'
#' Plots a pie chart, considering negative values as zeros, for the purpose of illustrating variation partitioning.
#' @param x  A vector of non-negative numerical quantities. The values in x are displayed as the areas of pie slices.
#' @param labels One or more expressions or character strings giving names for the slices. Other objects are coerced by \code{\link[grDevices]{as.graphicsAnnot}}. For empty or NA (after coercion to character) labels, no label nor pointing line is drawn.
#' @param edges  The circular outline of the pie is approximated by a polygon with this many edges.
#' @param radius  The pie is drawn centered in a square box whose sides range from -1 to 1. If the character strings labeling the slices are long it may be necessary to use a smaller radius.
#' @param clockwise  Logical indicating if slices are drawn clockwise or counter clockwise (i.e., mathematically positive direction, used by default).
#' @param init.angle  number specifying the starting angle (in degrees) for the slices. Defaults to 0 (i.e., '3 o'clock') unless clockwise is true where init.angle defaults to 90 (degrees), (i.e., '12 o'clock').
#' @param density  The density of shading lines, in lines per inch. The default value of NULL means that no shading lines are drawn. Non-positive values of density also inhibit the drawing of shading lines.
#' @param angle  The slope of shading lines, given as an angle in degrees (counter-clockwise).
#' @param col  A vector of colors to be used in filling or shading the slices. If missing a set of 6 pastel colours is used, unless density is specified when par("fg") is used.
#' @param border,lty  (possibly vectors) arguments passed to polygon which draws each slice.
#' @param main  An overall title for the plot.
#' @param warning Boolean value. Set to FALSE to avoid displaying a warning if some values are negative and set to 0.
#' @param ...  Graphical parameters can be given as arguments to pie. They will affect the main title and labels only.
#' @details This function is identical to the function \code{\link[graphics]{pie}} in \{graphics\}, except that it considers all negative values as zeros, to allow for plotting variation partitioning outputs. The original \code{\link[graphics]{pie}} function returns an error when negative values are present. However, variation partitioning can return negative values, which can then be treated as zeros (Legendre & Legendre, 2008). This function allows direct use of the results from \code{\link{Zeta.varpart}} without editing the data.
#' @seealso \code{\link[graphics]{pie}}, \code{\link{Zeta.varpart}}
#' @references  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988). \emph{The new S language}. Wadsworth & Brooks/Cole.
#' @references  Cleveland, W. S. (1985). \emph{The elements of graphing data}. Wadsworth: Monterey, CA, USA.
#' @references Legendre, P. &  Legendre, L.F. (2012). \emph{Numerical ecology}, 3rd English edition. Elsevier Science BV, Amsterdam.
#' @examples
#' pie.neg(rep(1, 24), col = rainbow(24), radius = 0.9)
#' @export
pie.neg <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE,
                     init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45,
                     col = NULL, border = NULL, lty = NULL, main = NULL, warning = TRUE, ...)
{
  if (!is.numeric(x) || any(is.na(x)))
    stop("Error: 'x' values must be numeric.")
  if (sum(x < 0)>0){
    x[which(x<0)] <- 0
    warning("Negative values set to 0.")
  }

  if (is.null(labels))
    labels <- as.character(seq_along(x))
  else labels <- grDevices::as.graphicsAnnot(labels)
  x <- c(0, cumsum(x) / sum(x))
  dx <- diff(x)
  nx <- length(dx)
  graphics::plot.new()
  pin <- graphics::par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L] / pin[2L]) * xlim
  else ylim <- (pin[2L] / pin[1L]) * ylim
  grDevices::dev.hold()
  on.exit(grDevices::dev.flush())
  graphics::plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col))
    col <- if (is.null(density))
      c("white", "lightblue", "mistyrose", "lightcyan",
        "lavender", "cornsilk")
  else graphics::par("fg")
  if (!is.null(col))
    col <- rep_len(col, nx)
  if (!is.null(border))
    border <- rep_len(border, nx)
  if (!is.null(lty))
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density))
    density <- rep_len(density, nx)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi / 180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    graphics::polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i],
                      border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      graphics::lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      graphics::text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE,
                     adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  graphics::title(main = main, ...)
  invisible(NULL)
}


####################
##HELPER FUNCTIONS##
####################
.gdist_matrix <- function(xy){
  xy <- as.matrix(xy)
  return(stats::as.dist(geodist::geodist(xy)))
}



.Mi <- function(i,k,x,ts){
  MM <- NA
  if(k==1){
    if(ts[i]<=x && x<ts[i+1]){MM <- 1/(ts[i+1]-ts[i])}
    else{MM <- 0}
  }else{
    tss <-ts[i+k]
    if(i+k > length(ts)){tss <- 1}
    if(tss-ts[i]==0){
      MM <- 0
    }else{
      MM <- k*((x-ts[i])*.Mi(i,k-1,x,ts)+(tss-x)*.Mi(i+1,k-1,x,ts))/((k-1)*(tss-ts[i]))
    }
  }
  return(MM)
}


.Ii <- function(i,k,x,ts){
  j <- which(ts >= x)[1]-1
  II <- 0
  if((j-k+1)>i){
    II <- 1
  }else if(i <= j){
    for(m in i:j){
      if(m+k+1 > length(ts)){
        II <- II+(1-ts[m])*.Mi(m,k+1,x,ts)/(k+1)
      }else{
        II <- II+(ts[m+k+1]-ts[m])*.Mi(m,k+1,x,ts)/(k+1)
      }
    }
  }
  return(II)
}









####################
##DATA DESCRIPTION##
####################


#' South-East Australia Environmental Dataset at Coarse Scale
#'
#' Projected coordinates and environmental variables in 123, 100 x 100 km sites.
#'
#'
#' The data set contains the following variables:
#'
#' \itemize{
#' \item{x}: x-position in meters in UTM 53 South projection
#' \item{y}: y-position in meters in UTM 53 South projection
#' \item{Natural}: Proportion of area of conservation and natural environments
#' \item{Irrigated}: Proportion of area of production from irrigated agriculture and plantations
#' \item{Water}: Proportion of area of water features
#' \item{Elevation}: Elevation
#' \item{ApP}: Area per person
#' \item{Temp}: Temperature
#' \item{Precip}: Precipitation
#' }
#'
#' Location: Australia -- 51° 27' 2.27" S, 135° 21' 35.19" E
#'
#' Data owners: ABARES, Australian Bureau of Statistics,GEBCO, WorldClim
#' @name bird.env.coarse
#' @usage data(bird.env.coarse)
#' @docType data
#' @references http://data.daff.gov.au/anrdl/metadata_files/pa_luav4g9abl07811a00.xml
#' @references http://www.gebco.net/
#' @references http://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.0072011?
#' @references http://www.worldclim.org/
#' @references Hijmans, R.J., Cameron, S.E., Parra, J.L., Jones, P.G. & Jarvis, A. (2005) Very high resolution interpolated climate surfaces for global land areas. International journal of climatology, 25, 1965-1978.
#' @keywords data
#' @format A data frame with 123 rows (sites) and 9 columns (xy coordinates and environmental variables).
"bird.env.coarse"


#' South-East Australia Environmental Dataset at Fine Scale
#'
#' Projected coordinates and environmental variables in 604, 25 x 25 km contiguous sites.
#'
#'
#' The data set contains the following variables:
#'
#' \itemize{
#' \item{x}: x-position in meters in UTM 53 South projection
#' \item{y}: y-position in meters in UTM 53 South projection
#' \item{Natural}: Proportion of area of conservation and natural environments
#' \item{Irrigated}: Proportion of area of production from irrigated agriculture and plantations
#' \item{Water}: Proportion of area of water features
#' \item{Elevation}: Elevation
#' \item{ApP}: Area per person
#' \item{Temp}: Temperature
#' \item{Precip}: Precipitation
#' }
#'
#' Location: Australia -- 50° 33' 5.03" S, 135° 21' 10.40" E
#'
#' Data owners: ABARES, Australian Bureau of Statistics,GEBCO, WorldClim
#' @name bird.env.fine
#' @usage data(bird.env.fine)
#' @docType data
#' @references http://data.daff.gov.au/anrdl/metadata_files/pa_luav4g9abl07811a00.xml
#' @references http://www.gebco.net/
#' @references http://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.0072011?
#' @references http://www.worldclim.org/
#' @references Hijmans, R.J., Cameron, S.E., Parra, J.L., Jones, P.G. & Jarvis, A. (2005) Very high resolution interpolated climate surfaces for global land areas. International journal of climatology, 25, 1965-1978.
#' @keywords data
#' @format A data frame with 604 rows (sites) and 9 columns (xy coordinates and environmental variables).
"bird.env.fine"


#' Australia Bird Atlas Species Occurrence Dataset at Coarse Scale over South-East Australia
#'
#' Inventory of bird species occurrence in 123, 100 x 100 km sites.
#'
#'\itemize{
#' \item{x}: x-position in meters in UTM 53 South projection
#' \item{y}: y-position in meters in UTM 53 South projection
#' \item{columns 3-193}: bird species occurrence
#' }
#'
#' The original bird occurrence data were arranged into a continuous grid covering South-East Australia. Only cells whose richness was within 10 percents of real estimated richness are included here, so that the data corresponds to presence-absence data.
#'
#' Location: Australia -- 51° 27' 2.27" S, 135° 21' 35.19" E
#'
#' Data owner: BirdLife Australia
#' @name bird.spec.coarse
#' @usage data(bird.spec.coarse)
#' @docType data
#' @references Barrett, G., Silcocks, A., Barry, S., Cunningham, R. & Poulter, R. (2003) The new atlas of Australian birds. Royal Australasian Ornithologists Union, Melbourne, 1-824.
#' @keywords data
#' @format A data frame with 123 rows (sites) and 193 columns (xy coordinates and species).
"bird.spec.coarse"


#' Australia Bird Atlas Species Occurrence Dataset at Fine Scale over South-East Australia
#'
#' Inventory of bird species occurrence in 604, 25 x 25 km sites.
#'
#' \itemize{
#' \item{x}: x-position in meters in UTM 53 South projection
#' \item{y}: y-position in meters in UTM 53 South projection
#' \item{columns 3-192}: bird species occurrence
#' }
#'
#' Location: Australia -- 50° 33' 5.03" S, 135° 21' 10.40" E
#'
#' Data owner: BirdLife Australia
#'
#' The original bird occurrence data were arranged into a continuous grid covering South-East Australia. Only cells whose richness was within 10 percents of real estimated richness are included here, so that the data corresponds to presence-absence data.
#'
#' @name bird.spec.fine
#' @usage data(bird.spec.fine)
#' @docType data
#' @references Barrett, G., Silcocks, A., Barry, S., Cunningham, R. & Poulter, R. (2003) The new atlas of Australian birds. Royal Australasian Ornithologists Union, Melbourne, 1-824.
#' @keywords data
#' @format A data frame with 604 rows (sites) and 193 columns (xy coordinates and species).
"bird.spec.fine"



#' Marion Island Species Presence-Absence Dataset
#'
#' Inventory of springtails and mite species presence-absence in 12 plots (4 transects and 3 altitudes) on Marion Island.
#'
#' The data set contains the following variables:
#'
#' \itemize{
#' \item{x}: x-position in meters in UTM 37 South projection
#' \item{y}: y-position in meters in UTM 37 South projection
#' \item{columns 3-24}: mite species presence absence
#' \item{columns 25-33}: springtail species presence absence
#' }
#'
#' Location: Marion Island -- 46° 53' 34.2" S, 37 degrees 45' 02.3" E
#'
#' Data owner: Melodie A. McGeoch
#' @name Marion.species
#' @usage data(Marion.species)
#' @docType data
#' @references Nyakatya, M.J. & McGeoch, M.A. (2008). Temperature variation across Marion Island associated with a keystone plant species (\emph{Azorella selago} Hook. (Apiaceae)). Polar Biology, 31, 139-151.
#' @references McGeoch, M.A., Le Roux, P.C., Hugo, E.A. & Nyakatya, M.J. (2008). Spatial variation in the terrestrial biotic system. The Prince Edward Islands: Land-Sea Interactions in a Changing World (ed. by S.L. Chown and P.W. Froneman), pp. 245-276. African SunMedia, Stellenbosch.
#' @keywords data
#' @format A data frame with 12 rows (plots) and 33 columns (species).
"Marion.species"

#' Marion Island Environmental Dataset
#'
#' Geographic coordinates, altitude and island side (East, West) at 12 plots (4 transects and 3 altitudes) on Marion Island.
#'
#' The data set contains the following variables:
#'
#' \itemize{
#' \item{x}: x-position in meters in UTM 37 projection
#' \item{y}: y-position in meters in UTM 37 projection
#' \item{Altitude}: mean elevation
#' \item{Side}: cardinal (East or West) side of the island
#' }
#'
#' Location: Marion Island -- 46° 53' 34.2" S, 37 degrees 45' 02.3" E
#'
#' Data owner: Melodie A. McGeoch
#' @name Marion.env
#' @usage data(Marion.env)
#' @docType data
#' @references Nyakatya, M.J. & McGeoch, M.A. (2008). Temperature variation across Marion Island associated with a keystone plant species (\emph{Azorella selago} Hook. (Apiaceae)). Polar Biology, 31, 139-151.
#' @references McGeoch, M.A., Le Roux, P.C., Hugo, E.A. & Nyakatya, M.J. (2008). Spatial variation in the terrestrial biotic system. The Prince Edward Islands: Land-Sea Interactions in a Changing World (ed. by S.L. Chown and P.W. Froneman), pp. 245-276. African SunMedia, Stellenbosch.
#' @keywords data
#' @format A data frame with 12 rows (plots) and 4 columns (variables).
"Marion.env"









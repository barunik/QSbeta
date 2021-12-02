#' @include generics.R
#' @import quantspec
#' @import rugarch
#' @importMethodsFrom quantspec getPointwiseCIs
#' @useDynLib QSbeta
#' @importFrom Rcpp sourceCpp
NULL

################################################################################
#' Compute quantile spectral beta from a smoothed quantile periodogram.
#'
#' Returns quantile spectral beta defined as
#' \deqn{\frac{G^{j_1, j_2}(\omega; \tau_1, \tau_2)}{G^{j_1, j_1}(\omega; \tau_1, \tau_1) }}
#' where \eqn{G^{j_1, j_2}(\omega; \tau_1, \tau_2)} is the smoothed quantile
#' periodogram.
#' 
#' For the mechanism of selecting frequencies, dimensions and/or levels see,
#' for example, \code{\link{getValues-SmoothedPG}}.
#'
#' @name getBeta-SmoothedPG
#' @aliases getBeta,SmoothedPG-method
#'
#' @keywords Access-functions
#'
#' @param object \code{SmoothedPG} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#'
#' @return Returns data from the array \code{values} that's a slot of
#'          \code{object}.
#'
#' @seealso
#' An example on how to use this function is analogously to the example given in
#' \code{\link{getValues-QuantilePG}}.
#' 
#' @export
################################################################################
#' @importClassesFrom quantspec SmoothedPG
setMethod(f = "getBeta",
          signature = signature(
            "SmoothedPG"),
          definition = function(object,
                                frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
                                levels.1=getLevels(object,1),
                                levels.2=getLevels(object,2),
                                d1 = 1:(dim(object@values)[2]),
                                d2 = 1:(dim(object@values)[4])) {
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- 2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y)
            }
            if (!hasArg(levels.1)) {
              levels.1 <- object@levels[[1]]
            }
            if (!hasArg(levels.2)) {
              levels.2 <- object@levels[[2]]
            }
            if (!hasArg(d1)) {
              d1 <- 1:(dim(object@values)[2])
            }
            if (!hasArg(d2)) {
              d2 <- 1:(dim(object@values)[4])
            }
            # end: workaround
            
            J <- length(frequencies)
            #D <- dim(object@values)[2]
            D1 <- length(d1)
            D2 <- length(d2)
            K1 <- length(levels.1)
            K2 <- length(levels.2)
            res <- array(dim=c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1))
            
            
            if (class(object@weight) != "KernelWeight") {
              stop("QS beta can only be determined if weight is of type KernelWeight.")
            }
            d <- union(d1,d2)
            V <- array(getValues(object, d1 = d, d2 = d, frequencies = frequencies), dim=c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1))
            
            d1.pos <- closest.pos(d, d1)
            d2.pos <- closest.pos(d, d2)
            
            res <- .computeBeta(V, d1.pos, d2.pos)
            
            if (D1 == 1 && D2 == 1) {
              final.dim.res <- c(J, K1, K2, object@qPG@freqRep@B+1)
            } else {
              final.dim.res <- c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1)
            }
            
            res[is.nan(res)] <- NA
            res <- array(res, dim=final.dim.res)
            
            return(res)
          }
)

################################################################################
#' Get estimates for the standard deviation of the quantile-spectral beta computed
#' from smoothed quantile periodogram. This code is heavily based on the implementation
#' of \code{"getCoherencySdNaive"} from the original \code{quantspec} package.
#' It also uses slot sdCohNaive and sdCohSqNaive of the \code{SmoothedPG} class.
#'
#' Determines and returns an array of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)}, and
#' \code{K2=length(levels.2))}. Whether
#' available or not, boostrap repetitions are ignored by this procedure.
#' At position \code{(j,k1,k2)}
#' the returned value is the standard deviation estimated corresponding to
#' \code{frequencies[j]}, \code{levels.1[k1]} and \code{levels.2[k2]} that are
#' closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means.
#' 
#' If not only one, but multiple time series are under study, the dimension of
#' the returned vector is of dimension \code{[J,P,K1,P,K2]}, where \code{P}
#' denotes the dimension of the time series.
#'
#' Requires that the \code{\link{SmoothedPG}} is available at all Fourier
#' frequencies from \eqn{(0,\pi]}{(0,pi]}. If this is not the case the missing
#' values are imputed by taking one that is available and has a frequency
#' that is closest to the missing Fourier frequency; \code{closest.pos} is used
#' to determine which one this is.
#'
#' A precise definition on how the standard deviations of the smoothed quantile
#' periodogram are estimated is given in Barunik and Kley (2015). The estimate
#' returned is denoted by
#' \eqn{\sigma(\tau_1, \tau_2; \omega)}{sigma(tau1, tau2; omega)} on p. 26 of
#' the arXiv preprint.
#'
#' Note the ``standard deviation'' estimated here is not the square root of the
#' complex-valued variance. It's real part is the square root of the variance
#' of the real part of the estimator and the imaginary part is the square root
#' of the imaginary part of the variance of the estimator.
#'
#' @name getBetaSdNaive-SmoothedPG
#' @aliases getBetaSdNaive,SmoothedPG-method
#'
#' @keywords Access-functions
#'
#' @param object \code{\link{SmoothedPG}} of which to get the estimates for the
#'                standard deviation.
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#' @param type can be "1", where cov(Z, Conj(Z)) is subtracted, or "2", where
#' 					   it's not
#' @param impl choose "R" or "C" for one of the two implementations available
#'
#' @return Returns the estimate described above.
#'
#' @references
#' Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2016).
#' Quantile Spectral Processes: Asymptotic Analysis and Inference.
#' \emph{Bernoulli}, \bold{22}(3), 1770--1807.
#' [cf. \url{http://arxiv.org/abs/1401.8104}]
#' 
#' Barunik, J. & Kley, T. (2015).
#' Quantile Cross-Spectral Measures of Dependence between Economic Variables.
#' [preprint available from the authors]
#' 
#' Barunik, J. & Nevrla, M. (2021).
#' Quantile Spectral Beta: A Tale of Tail Risks, Investment Horizons, and Asset Prices.
#' [preprint available from the authors]
#' 
#' @export
################################################################################
setMethod(f = "getBetaSdNaive",
          signature = signature(
            object = "SmoothedPG"),
          definition = function(object,
                                frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
                                levels.1=getLevels(object,1),
                                levels.2=getLevels(object,2),
                                d1 = 1:(dim(object@values)[2]),
                                d2 = 1:(dim(object@values)[4]),
                                type = c("1", "2"),
                                impl=c("R","C")) {
            
            if (class(getWeight(object)) != "KernelWeight") {
              stop("getSdNaive currently only available for 'KernelWeight'.")
            }
            
            Y <- object@qPG@freqRep@Y
            objLevels1 <- object@levels[[1]]
            objLevels2 <- object@levels[[2]]
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- 2*pi*(0:(lenTS(Y)-1))/lenTS(Y)
            }
            if (!hasArg(levels.1)) {
              levels.1 <- objLevels1
            }
            if (!hasArg(levels.2)) {
              levels.2 <- objLevels2
            }
            if (!hasArg(d1)) {
              d1 <- 1:(dim(object@values)[2])
            }
            if (!hasArg(d2)) {
              d2 <- 1:(dim(object@values)[4])
            }
            if (!hasArg(type)) {
              type <- "1"
            }
            if (!hasArg(impl)) {
              impl <- "R"
            }
            # end: workaround
            
            N <- lenTS(Y)
            K1 <- length(objLevels1)
            K2 <- length(objLevels2)
            D1 <- length(d1)
            D2 <- length(d2)
            
            ## First determine which frequencies still have to be computed:
            #  i. e. for which j is 2*pi*j / N is to be computed.
            
            J.requested <- round(frequenciesValidator(frequencies, N = N) * N / (2*pi))  ## TODO: check if rounding required??
            J.toCompute <- J.requested ## Get SD for all frequencies setdiff(J.requested, object@env$sdCohNaive.freq)
            
            J <- length(J.toCompute)
            
            if ( J > 0 ) {
              
              # TODO: Make this work with type = 1 and type = 2
              #if (object@env$sdCohNaive.done == FALSE) {
              
              weight <- object@weight
              
              # Define a list which covariances to estimate
              # List shall be a matrix with rows
              # (ja ta jb tb jc tc jd td)
              #
              # by default all combinations shall be computed
              
              if (impl == "R") {
                
                lC <- matrix(ncol = 8, nrow = 7 * K1*K2*D1*D2 )
                
                i <- 1
                for (k1 in 1:K1) {
                  for (i1 in 1:D1) {
                    for (k2 in 1:K2) {
                      for (i2 in 1:D2) {
                        jt_1 <- c(i1,k1)
                        jt_2 <- c(i2,k2)
                        lC[i,] <- c(jt_1, jt_1, jt_1, jt_1)
                        i <- i + 1
                        lC[i,] <- c(jt_1, jt_1, jt_1, jt_2)
                        i <- i + 1
                        lC[i,] <- c(jt_1, jt_1, jt_2, jt_2)
                        i <- i + 1
                        lC[i,] <- c(jt_1, jt_2, jt_1, jt_2)
                        i <- i + 1
                        lC[i,] <- c(jt_1, jt_2, jt_2, jt_1)
                        i <- i + 1
                        lC[i,] <- c(jt_1, jt_2, jt_2, jt_2)
                        i <- i + 1
                        lC[i,] <- c(jt_2, jt_2, jt_2, jt_2)
                        i <- i + 1
                      }
                    }
                  }
                }
                
                lC <- matrix(lC[1:(i-1),], ncol = 8, nrow = i-1)
                lC <- unique(lC)
                
                #####
                ## Variant 1: more or less vectorized... 
                #####
                WW <- getValues(weight, N = N)[c(2:N,1)] # WW[j] corresponds to W_n(2 pi j/n)
                WW3 <- rep(WW,4) 
                
                # TODO: check whether it works for all d1, d2!!
                V <- array(getValues(object, frequencies = 2*pi*(1:(N-1))/N, d1=d1, d2=d2), dim=c(N-1,D1,K1,D2,K2))     
                res_coh <- array(0,dim=c(J,D1,K1,D2,K2))
                res_cohSq <- array(0,dim=c(J,D1,K1,D2,K2))
                auxRes <- array(0,dim=c(nrow(lC), J))
                
                M1 <- matrix(0, ncol=N-1, nrow=N)
                M2 <- matrix(0, ncol=N-1, nrow=N)
                #M3 <- matrix(0, ncol=N-1, nrow=N)
                #M4 <- matrix(0, ncol=N-1, nrow=N)
                for (j in 0:(N-1)) { # Test 1:N
                  M1[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N+j-(1:(N-1))]
                  M2[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N+j+(1:(N-1))]
                  #M3[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N-j-(1:(N-1))]
                  #M4[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N-j+(1:(N-1))]
                }
                
                M1 <- M1[J.toCompute+1,]
                M2 <- M2[J.toCompute+1,]
                
                for (r in 1:nrow(lC)) {
                  jt <- lC[r,]
                  V1 <- matrix(V[,jt[1],jt[2],jt[5],jt[6]] * Conj(V[,jt[3],jt[4],jt[7],jt[8]]), ncol=1)
                  V2 <- matrix(V[,jt[1],jt[2],jt[7],jt[8]] * Conj(V[,jt[3],jt[4],jt[5],jt[6]]), ncol=1)
                  
                  auxRes[r,] <- rowSums(M1 %*% V1) + rowSums(M2 %*% V2)
                }
                
                V <- array(getValues(object, frequencies = 2*pi*(J.toCompute)/N, d1=d1, d2=d2), dim=c(J,D1,K1,D2,K2))
                Coh <- array(getBeta(object, frequencies = 2*pi*(J.toCompute)/N, d1=d1, d2=d2), dim=c(J,D1,K1,D2,K2))
                
                
                for (i1 in 1:D1) {
                  for (k1 in 1:K1) {
                    for (i2 in 1:D2) {
                      for (k2 in 1:K2) {
                        #r <- which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i2 & lC[,8] == k2)
                        if (i1 == i2 && k1 == k2) {
                          S_coh <- 0
                          S_cohSq <- 0
                        } else {
                          H1111 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i1 & lC[,8] == k1),]
                          H1112 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i2 & lC[,8] == k2),]
                          H1122 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i2 & lC[,8] == k2),]
                          H1212 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i2 & lC[,8] == k2),]
                          H1221 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i1 & lC[,8] == k1),]
                          H1222 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i2 & lC[,8] == k2),]
                          H2222 <- auxRes[which(lC[,1] == i2 & lC[,2] == k2 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i2 & lC[,8] == k2),]
                          
                          #H1121 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i1 & lC[,8] == k1),]
                          H2111 <- auxRes[which(lC[,1] == i2 & lC[,2] == k2 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i1 & lC[,8] == k1),]
                          
                          f12 <- V[,i1,k1,i2,k2]
                          f11 <- V[,i1,k1,i1,k1]
                          f22 <- V[,i2,k2,i2,k2]
                          
                          A <- 1/(f11^2)*H1111 - Conj(f12)/(f11^3)*H1112 - f12/(f11^3)*Conj(H1112) + abs(f12)^2/(f11^4)*H1212
                          
                          #B <- 1/(f11*f22)*H1122 - f12/(f11*f22^2)*H1121 - f12/(f11^2*f22)*H1222 + f12^2/(f11^2*f22^2)*H1221 # Correct v1
                          B <- 1/(f11*f22)*H1122 - f12/(f11*f22^2)*Conj(H2111) - f12/(f11^2*f22)*H1222 + f12^2/(f11^2*f22^2)*H1221 # Correct v2
                          
                          #if (type == "1") {
                          S_coh <- (1/2) * complex(real = Re(A + B), imaginary = Re(A - B))
                          #} else {
                          # CORRECT??
                          # Recall, A == L1212 and B == L1221 ??
                          #R12 <- f12 / sqrt(abs(f11) * abs(f22)) # note that f11, f22 can be negative if weights are negative
                          # Get beta
                          #R12 <- f12 / f11
                          R12 <- Coh[,i1,k1,i2,k2]
                          S_cohSq <- 2 * abs(R12)^2 * A
                          S_cohSq <- S_cohSq + 2 * (Re(R12)^2 - Im(R12)^2) * Re(B)
                          S_cohSq <- S_cohSq + 4 * Re(R12) * Im(R12) * Im(B)
                          #}
                          
                        }
                        ## TODO: Comment on the next line!!
                        ## Is this because S is always a linear combination of the
                        ## Cov(Lab, Lcd) terms??
                        
                        res_coh[,i1,k1,i2,k2] <- (2*pi/N)^2 * S_coh / ((weight@env$Wnj[c(N,1:(N-1))])[J.toCompute+1])^2
                        res_coh[,i2,k2,i1,k1] <- res_coh[,i1,k1,i2,k2]
                        
                        res_cohSq[,i1,k1,i2,k2] <- (2*pi/N)^2 * S_cohSq / ((weight@env$Wnj[c(N,1:(N-1))])[J.toCompute+1])^2
                        res_cohSq[,i2,k2,i1,k1] <- res_cohSq[,i1,k1,i2,k2]
                      }
                    }
                  }
                }
                
                sqrt.cw <- function(z) {
                  return(complex(real=sqrt(max(Re(z),1e-9)), imaginary=sqrt(max(Im(z),1e-9))))
                }
                #if (type == "1") {
                res_coh <- array(apply(res_coh,c(1,2,3,4,5),sqrt.cw), dim = c(J, D1, K1, D2, K2))  
                #} else {
                res_cohSq <- array(apply(res_cohSq,c(1,2,3,4,5),sqrt), dim = c(J, D1, K1, D2, K2))
                #}
                
                
                #####
                ## END Variant 1: more or less vectorized... 
                #####
                
              }
              
              
              object@env$sdCohNaive.freq <- union(object@env$sdCohNaive.freq, J.toCompute)
              object@env$sdCohNaive[J.toCompute+1,,,,] <- res_coh
              object@env$sdCohSqNaive[J.toCompute+1,,,,] <- res_cohSq
              
            } # End of if (J > 0)
            #      } # End of 'if (object@env$sdNaive.done == FALSE) {'
            
            if (type == "1") {
              resObj <- object@env$sdCohNaive
            } else {
              resObj <- object@env$sdCohSqNaive
            }
            
            ##############################
            ## (Similar) Code also in Class-FreqRep!!!
            ##############################
            
            # Transform all frequencies to [0,2pi)
            frequencies <- frequencies %% (2*pi)
            
            # Create an aux vector with all available frequencies
            oF <- object@frequencies
            f <- frequencies
            
            # returns TRUE if x c y
            subsetequal.approx <- function(x,y) {
              X <- round(x, .Machine$double.exponent-2)
              Y <- round(y, .Machine$double.exponent-2)
              return(setequal(X,intersect(X,Y)))
            }
            
            C1 <- subsetequal.approx(f[f <= pi], oF)
            C2 <- subsetequal.approx(f[f > pi], 2*pi - oF[which(oF != 0 & oF != pi)])
            
            # Ist dann der Fall wenn die frequencies nicht geeignete FF sind!!
            if (!(C1 & C2)) {
              warning("Not all 'values' for 'frequencies' requested were available. 'values' for the next available Fourier frequencies are returned.")
            }
            
            # Select columns
            c.1.pos <- closest.pos(objLevels1,levels.1)
            c.2.pos <- closest.pos(objLevels2,levels.2)
            
            # Select rows
            r1.pos <- closest.pos(oF,f[f <= pi])
            r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])
            
            J <- length(frequencies)
            K1 <- length(levels.1)
            K2 <- length(levels.2)
            res <- array(dim=c(J, D1, K1, D2, K2))
            
            if (length(r1.pos) > 0) {
              res[which(f <= pi),,,,] <- resObj[r1.pos, , c.1.pos, , c.2.pos]
            }
            if (length(r2.pos) > 0) {
              res[which(f > pi),,,,] <- Conj(resObj[r2.pos, , c.1.pos, , c.2.pos])
            }
            
            return(res)
          }
)

################################################################################
#' Get pointwise confidence intervals for the quantile spectral density kernel,
#' quantile coherency, quantile coherence or quantile-spectral beta. This is an updated
#' version of the original \code{getPointwiseCIs} function from the \code{quantspec}
#' package which includes CIs for the quantile-spectral beta measure.
#'
#' Returns a list of two arrays \code{lowerCIs} and \code{upperCIs} that contain
#' the upper and lower limits for a level \code{1-alpha} confidence interval of
#' the quantity of interest. Each array is of dimension \code{[J,K1,K2]} if a
#' univariate time series is being analysed or of dimension \code{[J,D1,K1,D2,K2]},
#' where \code{J=length(frequencies)}, \code{D1=length(d1)}, \code{D2=length(d2)},
#' \code{K1=length(levels.1)}, and \code{K2=length(levels.2))}.
#' At position \code{(j,k1,k2)} or \code{(j,i1,k1,i2,k2)} the real (imaginary)
#' part of the returned values are the bounds of the confidence interval for the
#' the real (imaginary) part of the quantity under anlysis, which corresponds to
#' \code{frequencies[j]}, \code{d1[i1]}, \code{d2[i2]}, \code{levels.1[k1]} and
#' \code{levels.2[k2]} closest to the Fourier frequencies, \code{levels.1} and
#' \code{levels.2} available in \code{object}; \code{\link{closest.pos}} is used
#' to determine what closest to means.
#'
#' Currently, pointwise confidence bands for three different \code{quantity}
#' are implemented:
#' \itemize{
#'   \item \code{"spectral density"}: confidence intervals for the quantile spectral
#' 					 density as described in Kley et. al (2016) for the univariate case and
#' 					 in Barunik and Kley (2015) for the multivariate case.
#'   \item \code{"coherency"}: confidence intervals for the quantile coherency as
#' 					 described in Barunik and Kley (2015).
#'   \item \code{"beta"}: confidence intervals for the quantile-spectral betas
#'           as described in Barunik and Nevrla (2021).  
#' }
#' 
#' Currently, three different \code{type}s of confidence intervals are
#' available:
#' \itemize{
#'   \item \code{"naive.sd"}: confidence intervals based on the asymptotic
#'           normality of the smoothed quantile periodogram; standard deviations
#'           are estimated using \code{\link{getSdNaive}}.
#'   \item \code{"boot.sd"}: confidence intervals based on the asymptotic
#'           normality of the smoothed quantile periodogram; standard deviations
#'           are estimated using \code{\link{getSdBoot}}.
#'   \item \code{"boot.full"}: confidence intervals determined by estimating the
#'           quantiles of he distribution of the smoothed quantile periodogram,
#'           by the empirical quantiles of the sample of bootstrapped
#'           replications.
#' }
#'
#' @name getPointwiseCIs-SmoothedPG
#' @aliases getPointwiseCIs,SmoothedPG-method
#' 
#' @importFrom stats qnorm
#' @importFrom stats quantile
#'
#' @keywords Access-functions
#'
#' @param object \code{SmoothedPG} of which to get the confidence intervals
#' @param quantity a flag indicating for which the pointwise confidence bands
#' 								 will be determined. Can take one of the possible values
#' 								 discussed above. 
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#' @param alpha the level of the confidence interval; must be from \eqn{(0,1)}
#' @param type a flag indicating which type of confidence interval should be
#'         returned; can take one of the three values discussed above.
#'
#' @return Returns a named list of two arrays \code{lowerCIS} and \code{upperCIs}
#'          containing the lower and upper bounds for the confidence intervals.
#'
#' @examples
#' sPG <- smoothedPG(rnorm(2^10), levels.1=0.5)
#' CI.upper <- Re(getPointwiseCIs(sPG)$upperCIs[,1,1])
#' CI.lower <- Re(getPointwiseCIs(sPG)$lowerCIs[,1,1])
#' freq = 2*pi*(0:1023)/1024
#' plot(x = freq, y = rep(0.25/(2*pi),1024),
#'    ylim=c(min(CI.lower), max(CI.upper)),
#'    type="l", col="red") # true spectrum
#' lines(x = freq, y = CI.upper)
#' lines(x = freq, y = CI.lower)
#' 
#' @export
################################################################################
setMethod(f = "getPointwiseCIs",
          signature = signature(
            object = "SmoothedPG"),
          definition = function(object,
                                quantity = c("spectral density", "coherency", "coherence", "beta"),
                                frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
                                levels.1=getLevels(object,1),
                                levels.2=getLevels(object,2),
                                d1 = 1:(dim(object@values)[2]),
                                d2 = 1:(dim(object@values)[4]),
                                alpha=.1, type=c("naive.sd", "boot.sd", "boot.full")) {
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- 2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y)
            }
            if (!hasArg(levels.1)) {
              levels.1 <- object@levels[[1]]
            }
            if (!hasArg(levels.2)) {
              levels.2 <- object@levels[[2]]
            }
            if (!hasArg(alpha)) {
              alpha <- 0.1
            }
            if (!hasArg(d1)) {
              d1 <- 1:(dim(object@values)[2])
            }
            if (!hasArg(d2)) {
              d2 <- 1:(dim(object@values)[4])
            }
            if (!hasArg(type)) {
              type <- "naive.sd"
            }
            # end: workaround
            
            type <- match.arg(type)[1]
            quantity <- match.arg(quantity)[1]
            switch(type,
                   "naive.sd" = {
                     switch(quantity,
                            "spectral density" = {
                              sdEstim <- getSdNaive(object,
                                                    frequencies = frequencies,
                                                    levels.1 = levels.1,
                                                    levels.2 = levels.2,
                                                    d1 = d1, d2 = d2)
                            },
                            "coherency" = {
                              sdEstim <- getCoherencySdNaive(object,
                                                             frequencies = frequencies,
                                                             levels.1 = levels.1,
                                                             levels.2 = levels.2,
                                                             d1 = d1, d2 = d2, type="1")
                            },
                            "coherence" = {
                              sdEstim <- getCoherencySdNaive(object,
                                                             frequencies = frequencies,
                                                             levels.1 = levels.1,
                                                             levels.2 = levels.2,
                                                             d1 = d1, d2 = d2, type="2")
                            },
                            "beta" = {
                              sdEstim <- getBetaSdNaive(object,
                                                        frequencies = frequencies,
                                                        levels.1 = levels.1,
                                                        levels.2 = levels.2,
                                                        d1 = d1, d2 = d2, type="1")
                            })
                   },
                   "boot.sd" = {
                     if (quantity == "spectral density") {
                       sdEstim <- getSdBoot(object,
                                            frequencies = frequencies,
                                            levels.1 = levels.1,
                                            levels.2 = levels.2)
                     } else {
                       stop("boot.sd is so far only implemented for quantity = 'spectral density'.")
                     }
                   }
            )
            
            J <- length(frequencies)
            K1 <- length(levels.1)
            K2 <- length(levels.2)
            D1 <- length(d1)
            D2 <- length(d2)
            
            
            if (type == "naive.sd" || type == "boot.sd") {
              
              switch(quantity,
                     "spectral density" = {
                       v <- array(getValues(object,
                                            frequencies = frequencies,
                                            levels.1 = levels.1,
                                            levels.2 = levels.2, d1=d1, d2=d2), dim = c(J, D1, K1, D2, K2))
                       sdEstim <- array(sdEstim, dim = c(J, D1, K1, D2, K2))
                       upperCIs <- array(v + sdEstim * qnorm(1-alpha/2), dim = c(J, D1, K1, D2, K2))
                       lowerCIs <- array(v + sdEstim * qnorm(alpha/2), dim = c(J, D1, K1, D2, K2))
                     },
                     "coherency" = {
                       v <- array(getCoherency(object,
                                               frequencies = frequencies,
                                               levels.1 = levels.1,
                                               levels.2 = levels.2, d1=d1, d2=d2), dim = c(J, D1, K1, D2, K2))
                       upperCIs <- array(v + sdEstim * qnorm(1-alpha/2), dim = c(J, D1, K1, D2, K2))
                       lowerCIs <- array(v + sdEstim * qnorm(alpha/2), dim = c(J, D1, K1, D2, K2))
                     },
                     "coherence" = {
                       v <- array(getCoherency(object,
                                               frequencies = frequencies,
                                               levels.1 = levels.1,
                                               levels.2 = levels.2, d1=d1, d2=d2), dim = c(J, D1, K1, D2, K2))
                       
                       upperCIs <- array(abs(v)^2 + Re(sdEstim) * qnorm(1-alpha/2), dim = c(J, D1, K1, D2, K2))
                       lowerCIs <- array(abs(v)^2 + Re(sdEstim) * qnorm(alpha/2), dim = c(J, D1, K1, D2, K2))
                     },
                     "beta" = {
                       v <- array(getBeta(object,
                                          frequencies = frequencies,
                                          levels.1 = levels.1,
                                          levels.2 = levels.2, d1=d1, d2=d2), dim = c(J, D1, K1, D2, K2))
                       upperCIs <- array(v + sdEstim * qnorm(1-alpha/2), dim = c(J, D1, K1, D2, K2))
                       lowerCIs <- array(v + sdEstim * qnorm(alpha/2), dim = c(J, D1, K1, D2, K2))
                     })
              
            } else if (type == "boot.full") {
              
              if (quantity == "spectral density") {
                
                # TODO: Error Msg ausgeben falls B == 0
                B <- object@qPG@freqRep@B
                # TODO: fix...
                
                switch(quantity,
                       "spectral density" = {
                         v <- getValues(object,
                                        frequencies = frequencies,
                                        levels.1 = levels.1,
                                        levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
                         
                       },
                       "coherency" = {
                         v <- getCoherency(object,
                                           frequencies = frequencies,
                                           levels.1 = levels.1,
                                           levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
                       },
                       "coherence" = {
                         v <- getCoherency(object,
                                           frequencies = frequencies,
                                           levels.1 = levels.1,
                                           levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
                         v <- abs(v)^2
                       },
                       "beta" = {
                         v <- getBeta(object,
                                      frequencies = frequencies,
                                      levels.1 = levels.1,
                                      levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
                       })
                
                v <- getValues(object,
                               frequencies = frequencies,
                               levels.1 = levels.1,
                               levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
                uQuantile <- function(x) {complex(real = quantile(Re(x),1-alpha/2),
                                                  imaginary = quantile(Im(x),1-alpha/2))}
                lQuantile <- function(x) {complex(real = quantile(Re(x),alpha/2),
                                                  imaginary = quantile(Im(x),alpha/2))}
                
                upperCIs <- apply(v, c(1,2,3), uQuantile)
                lowerCIs <- apply(v, c(1,2,3), lQuantile)
                
              }
            }
            
            if (D1 == 1 && D2 == 1) {
              final.dim.res <- c(J, K1, K2)
            } else {
              final.dim.res <- c(J, D1, K1, D2, K2)
            }
            
            lowerCIs <- array(lowerCIs, dim=final.dim.res)
            upperCIs <- array(upperCIs, dim=final.dim.res)
            
            res <- list(lowerCIs = lowerCIs, upperCIs = upperCIs)
            return(res)
            
          }
)
##' Fit a smoothing spline with constraints on derivatives
##'
##' For further details see the references.
##' @title constrained.smoothspline
##' @param y a vector of response variables.
##' @param times a vector of time points at which y was measured, same
##'   length as y.
##' @param pen.degree desired degree of the derivative in smoothing
##'   penalty.
##' @param constraint one of 'none', 'fixed' or 'bounded' depending on
##'   whether the derivatives should not be constrained, fixed to a
##'   constant or bounded in an interval, respectively.
##' @param derivative.values either a vector with the same length as y
##'   if contraint=='fixed' or a matrix with 2 columns conatining the
##'   lower and upper bounds on the derivatives if
##'   contraint=='bounded'.
##' @param initial.value optional paramter that specifies an initial
##'   value of the spline incase it should be fixed.
##' @param times.new optional vector of new time points at which the
##'   spline should be evaluated.
##' @param num.folds either a numeric value of the number of folds or
##'   the string "leave-one-out" for a leave one out type
##'   cross-validation in determining the penalty parameter.
##' @param lambda either a numeric value if the penalty parameter is
##'   fixed explicitely or one of the values 'optim' or 'grid.search'
##'   depending on the desired optimization procedure.
##' 
##' @return a list consisting of the following
##'   elements
##'
##' \item{smooth.vals}{predicted values at points times}
##' \item{residual}{residuals}
##' \item{smooth.vals.new}{predicted values at time points times.new}
##' \item{smooth.deriv}{predicted derivative values at points times}
##' \item{pen.par}{penality parameter used for smoothing}
##' 
##' @export
##'
##' @import fda cvTools quadprog
##'
##' @author Niklas Pfister, Stefan Bauer and Jonas Peters
##'
##' @references
##' Pfister, N., S. Bauer, J. Peters (2018).
##' Identifying Causal Structure in Large-Scale Kinetic Systems
##' ArXiv e-prints (arXiv:1810.11776).
##'
##' @seealso 
##'
##' @examples
##' ## Example
##' x <- seq(0,4,length.out=200)
##' x.long <- seq(0,4,length.out=200)
##' y <- x^2+rnorm(200,0,2)
##' dy <- 2*x
##' dybdd <- cbind(dy-0.5,dy+0.5)
##' plot(x,y)
##'
##' ptm <- proc.time()
##' fit <- constrained.smoothspline(y=y,
##'                                 times=x,
##'                                 pen.degree=2,
##'                                 constraint="none",
##'                                 derivative.values=NA,
##'                                 times.new=x.long,
##'                                 num.folds=5,
##'                                 lambda="optim")
##' print(proc.time()-ptm)
##' fit2 <- smooth.spline(x,y)
##' lines(x.long,fit[[3]],col="blue")
##' lines(fit2, col="green")
##' lines(x.long, x.long^2, col="black")


constrained.smoothspline <- function(y,
                                     times,
                                     pen.degree,
                                     constraint="fixed",
                                     derivative.values=NA,
                                     initial.value=NA,
                                     times.new,
                                     num.folds="leave-one-out",
                                     lambda="optim"){

  ## Reorder data and deal with repetitions
  ux <- sort(unique(times))
  w <- rep(1, length(y))
  tmp <- matrix(unlist(tapply(seq(along = y), match(times, ux), 
                              function(x, y, w) c(mean(y[x]), sum(w[x])), y = y, w = w)), 
              , 2, byrow = TRUE)
  y <- tmp[, 1]
  times <- ux
    
  ## Initialize some variables
  order_splines <- pen.degree+2
  
  ## Construct folds for CV

  if(num.folds=="leave-one-out"){
    folds <- vector("list", length(times))
    for(i in 1:length(times)){
      folds[[i]] <- i
    }
    num.folds <- length(folds)
  }
  else if(is.numeric(num.folds)){
    if(num.folds>1){
      folds_tmp <- cvFolds(length(times), K=num.folds, type="interleaved")
      folds <- vector("list", num.folds)
      for(i in 1:num.folds){
        folds[[i]] <- folds_tmp$subsets[folds_tmp$which == i]
      }
    }
    else{
      stop("num.folds should be at least 2")
    }
  }
  else{
    stop("num.folds was specified incorrectly")
  }

  ## Function to ensure matrix is positive definite
  make.posdef <- function(A, mineig=10^(-15)){
    aa <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    if(min(aa) < mineig){
      warning("Spline matrix is not (numerically) positive definite and has been adjusted.")
      if(min(aa) < 0){
        A <- A + diag(-min(aa) + mineig, dim(A)[1])
      }
      else{
        A <- A + diag(mineig, dim(A)[1])
      }
    }
    return(A)
  }


  ##############################
  # 
  # Step 1: Set up spline basis for CV
  #
  ##############################

  # intialize basis
  basis <- create.bspline.basis(norder=order_splines, breaks=times)
  penmat <- getbasispenalty(basis,Lfdobj=pen.degree)

  # CV variables
  if(!is.numeric(lambda)){
    meq <- vector("numeric", num.folds)
    dvec <- vector("list", num.folds)
    Amat <- vector("list", num.folds)
    bvec <- vector("list", num.folds)
    Dmat_1 <- vector("list", num.folds)
    Bmat.val <- vector("list", num.folds)
    validation <- vector("list", num.folds)
    # compute the basis variables for each fold
    for(i in 1:num.folds){
      train <- (1:length(times))[-folds[[i]]]
      validation[[i]] <- folds[[i]]
      train_y <- matrix(y[train], nrow=1)
      Bmat <- getbasismatrix(times[train], basis, nderiv=0, returnMatrix=FALSE)
      Bmat.val[[i]] <- getbasismatrix(times[validation[[i]]], basis, nderiv=0, returnMatrix=FALSE)
      Bmat.deriv <- getbasismatrix(times[train], basis, nderiv=1, returnMatrix=FALSE)
      dvec[[i]] <- train_y %*% Bmat
      Dmat_1[[i]] <- t(Bmat) %*% Bmat
      # define QP-variables according to constraint
      if(constraint=="fixed"){
        meq[i] <- length(train)
        Amat[[i]] <- t(Bmat.deriv)
        bvec[[i]] <- derivative.values[train]
      }
      else if(constraint=="bounded"){
        meq[i] <- 0
        Amat[[i]] <- t(rbind(Bmat.deriv,-Bmat.deriv))
        bvec[[i]] <- c(derivative.values[train,1],-derivative.values[train,2])
      }
      else if(constraint=="none"){
        meq[i] <- 0
        Amat[[i]] <- matrix(0, ncol(Bmat.deriv), nrow(Bmat.deriv))
        bvec[[i]] <- matrix(0, length(train))
      }
      else{
        stop("Specified constraint does not exist. Use either fixed of bounded.")
      }
    }
  }

  if(constraint=="none"){
    sol_index <- 3
  }
  else{
    sol_index <- 1
  }


  ##############################
  # 
  # Step 2: Define cost function
  #
  ##############################

  # Initialize upper and lower bound for penalty (see 5.4.1 in Functional Data Analysis)
  Bmat <- getbasismatrix(times, basis, nderiv=0, returnMatrix=FALSE)
  BBnorm <- sum(diag(t(Bmat)%*%Bmat))
  Rnorm <- sum(diag(penmat))
  r <- BBnorm/Rnorm
  lower.spar <- 0
  upper.spar <- 5
  lower.lambda <- r*256^(3*lower.spar-1)
  upper.lambda <- r*256^(3*upper.spar-1)


  cost_function <- function(spar){
    lambda <- r*256^(3*spar-1)
    Dmat_2 <- lambda * penmat
    rss <- 0
    for(i in 1:num.folds){
      Dmat <- Dmat_1[[i]] + Dmat_2
      sc <- norm(Dmat, "F")
      Dmatsc <- Dmat/sc
      # make sure that Dmatsc is pos.def
      Dmatsc <- make.posdef(Dmatsc)
      # solve quadratic program with equality constraints
      solutions_qp <- solve.QP(Dmatsc,
                               dvec[[i]]/sc,
                               Amat[[i]],
                               bvec[[i]],
                               meq[i])
      csol <- matrix(solutions_qp[[sol_index]], ncol=1)
      # compute validation error
      rss <- rss + sum((y[validation[[i]]]-Bmat.val[[i]]%*%csol)^2)
    }
    return(rss/num.folds)
  }
  
    
  ##############################
  # 
  # Step 3: Compute lambda
  #
  ##############################

  if(lambda=="optim"){
    solutions.optim <- optimize(cost_function, c(lower.spar, upper.spar))
    spar <- as.numeric(solutions.optim$minimum)
    lambda <- r*256^(3*spar-1)
  }
  else if(lambda=="grid.search"){
    #lambda.vec <- 10^(seq(log(lower.lambda)/log(10),log(upper.lambda)/log(10), by = 0.05))
    spar.vec <- seq(lower.spar, upper.spar, length.out=100)
    RSS.vec <- rep(NA, length(spar.vec))
    for(kkk in 1:length(spar.vec)){
      RSS.vec[kkk] <- cost_function(spar.vec[kkk])
    }
    index.best.spar <- length(RSS.vec)
    best.RSS <- RSS.vec[index.best.spar]
    for(j in (index.best.spar-1):1){
      if(RSS.vec[j] < 0.9*best.RSS){
        index.best.spar <- j
        best.RSS <- RSS.vec[index.best.spar]
      }
    }
    spar <- spar.vec[index.best.spar]
    lambda <- r*256^(3*spar-1)
  }
  else if(is.numeric(lambda)){
    lambda <- lambda
  }
  else{
    stop("specified lambda is not a valid parameter or method")
  }
  
  # Check whether lambda is attained at boundaries
  if(abs(lambda - lower.lambda) < 10^-16){
    warning("There was at least one case in which CV yields a lambda at the lower boundary.")
  }
  if(abs(lambda - upper.lambda) < 10^-16){
    warning("There was at least one case in which CV yields a lambda at the upper boundary.")
  }

    
  ##############################
  # 
  # Step 4: Smoothing based on lambda
  #
  ##############################

  # set up quadratic program
  Bmat <- getbasismatrix(times, basis, nderiv=0, returnMatrix=FALSE)
  Bmat.deriv <- getbasismatrix(times, basis, nderiv=1, returnMatrix=FALSE)
  dvec <- t(y %*% Bmat)
  Dmat_1 <- t(Bmat) %*% Bmat
  Dmat_2 <- lambda * penmat
  Dmat <- Dmat_1 + Dmat_2
  sc <- norm(Dmat, "F")
  Dmatsc <- Dmat/sc
  if(constraint=="fixed"){
    meq <- length(times)
    Amat <- t(Bmat.deriv)
    bvec <- derivative.values
  }
  else if(constraint=="bounded"){
    meq <- 0
    Amat <- t(rbind(Bmat.deriv,-Bmat.deriv))
    bvec <- c(derivative.values[,1],-derivative.values[,2])
  }
  else if(constraint=="none"){
    meq <- 0
    Amat <- matrix(0, ncol(Bmat.deriv), nrow(Bmat.deriv))
    bvec <- matrix(0, length(times))
  }

  # add initial value as constraint
  if(!is.na(initial.value)){
    meq <- meq+1
    Amat <- cbind(Amat, c(1, rep(0, nrow(Amat)-1)))
    bvec <- c(bvec, initial.value)
  }

  # make sure Dmatsc is pos.def.
  Dmatsc <- make.posdef(Dmatsc)
  # solve quadratic program depending on constraint
  solutions_qp <- solve.QP(Dmatsc,
                           dvec/sc,
                           Amat,
                           bvec,
                           meq)
  csol <- matrix(solutions_qp[[sol_index]], ncol=1)
  # smoothed values at times.new
  Bmat.new <- getbasismatrix(times.new, basis, nderiv=0, returnMatrix=FALSE)
  predict_values <- Bmat.new %*% csol
  # smoothed values at times
  smoothed_spline_values <- Bmat %*% csol
  # smoothed derivative values
  smooth_deriv <- Bmat.deriv %*% csol
  # residuals
  residuals <- (y-smoothed_spline_values)
  
  return(list(smooth.vals=smoothed_spline_values,
              residuals=residuals,
              smooth.vals.new=predict_values,
              smooth.deriv=smooth_deriv,
              pen.par=lambda))
}

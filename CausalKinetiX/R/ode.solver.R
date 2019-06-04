##' Solves a mass-action ODE for a target variable Y by using smooth
##' approximations of the predictor variables X.
##'
##' For further details see the references.
##' @title ode.solver
##' @param time_vec numeric vector. Specifices the points at which to
##'   evaluate the target trajectory.
##' @param initial_value numeric value. Specifies the value of the
##'   target at initial time point.
##' @param times numeric vector of length L. Specifices the time
##'   points at which the predictors where observed.
##' @param X predictor matrix of dimension L x d. Each column
##'   corresponds to a different predictor observed at the time points
##'   \code{times}.
##' @param model list of mass-action terms. Each element in the list
##'   consists of a vector of predictor variables which are multiplied
##'   to a single term in the mass-action equation.
##' @param target integer specifing which variable is the target.
##' @param coefs numeric vector. Specifies the parameter values for
##'   each term in \code{model}.
##' @param smooth.type string. Specifies which type of smoothing to
##'   use. The following options exist: "smoothing.spline", "loess",
##'   "linear", "constant".
##' @param reltol numeric. Relative tolarance used in CVODE.
##' @param abstol numeric. Absolute tolarance used in CVODE.
##' 
##' @return object returned by CVODE.
##' 
##' @export
##'
##' @import sundialr
##' 
##' @author Niklas Pfister, Stefan Bauer and Jonas Peters
##'
##' @references
##' Pfister, N., S. Bauer, J. Peters (2018).
##' Identifying Causal Structure in Large-Scale Kinetic Systems
##' ArXiv e-prints (arXiv:1810.11776).
##'
##' @examples
##'
##' ## Generate data from Maillard reaction
##' simulation.obj <- generate.data.maillard(target=4,
##'                                          env=rep(1:5, each=5),
##'                                          L=20)
##'
##' D <- simulation.obj$simulated.data
##' time <- simulation.obj$time
##' env <- simulation.obj$env
##' target <- simulation.obj$target
##'
##' ## Solve for Melanoidin
##' X <- do.call(cbind, split(as.vector(t(D[1:5,])), rep(1:11, each=length(unique(time)))))
##' times <- rep(unique(time), 5)
##' odefit <- ode.solver(time, 0, times, X, list(c(8)), 11, 0.12514)
##' plot(odefit[,1], odefit[,2], type="l")
##' points(times, X[,11])


ode.solver <- function(time_vec, initial_value, times, X, model,
                       target, coefs, smooth.type="smoothing.spline",
                       reltol=10^(-10), abstol=10^(-16)){

  ## Remove all variables not contained in model
  included.vars <- sort(unique(unlist(model)))

  ## Fit spline on each predictor
  if(smooth.type == "smoothing.spline"){
    splinefun <- lapply(1:length(included.vars), function(j) smooth.spline(times, X[,j]))
    splinefun <- lapply(splinefun, function(fit){function(t) predict(fit, t)$y})
  }
  else if(smooth.type == "loess"){
    splinefun <- lapply(1:length(included.vars), function(j) loess(X[,j] ~ times, span=0.50))
    splinefun <- lapply(splinefun, function(fit){function(t) predict(fit, t)})
  }
  else if(smooth.type == "linear"){
    splinefun <- lapply(1:length(included.vars),
                        function(j) approxfun(times, X[,j], method="linear"))
  }
  else if(smooth.type == "constant"){
    splinefun <- lapply(1:length(included.vars),
                        function(j) approxfun(times, X[,j], method="constant"))
  }

  ## Construct RHS
  ## odefun <- function(t, y, par){
  ##   deriv <- 0
  ##   for(term in 1:length(model)){
  ##     tmp <- 1
  ##     for(var in model[[term]]){
  ##       if(var == target){
  ##         tmp <- tmp*y
  ##       }
  ##       else{
  ##         tmp <- tmp*splinefun[[which(var == included.vars)]](t)
  ##       }
  ##     }
  ##     deriv <- deriv + par[term]*tmp
  ##   }
  ##   return(deriv)
  ## }
  par <- coefs
  odefun <- function(t, y){
    deriv <- 0
    for(term in 1:length(model)){
      tmp <- 1
      for(var in model[[term]]){
        if(var == target){
          tmp <- tmp*y
        }
        else{
          tmp <- tmp*splinefun[[which(var == included.vars)]](t)
        }
      }
      deriv <- deriv + par[term]*tmp
    }
    return(deriv)
  }
  

  ## Solve ODE
  ## odefit <- cvode(time_vec, initial_value, odefun, coefs, reltol, abstol)
  odefit <- cvode(time_vec, initial_value, odefun, reltol, abstol)

  
  return(odefit)
}

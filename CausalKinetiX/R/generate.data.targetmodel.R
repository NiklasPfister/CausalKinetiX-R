##' Generate sample data from the target model based on predictor trajectories.
##'
##' For further details see the references.
##' @title Target model based on predictor trajectories
##' @param env integer vector of length n encoding to which experiment
##'   each repetition belongs.
##' @param noise.sd numerical value specifying the standard deviation
##'   of the noise.
##' @param L number of time points for evaluation.
##' @param d number of total variables (d-1 preditor variables).
##' @param seed random seed. Does not work if a "Detected blow-up"
##'   warning shows up.
##'
##' @return list consisting of the following elements
##' 
##' \item{simulated.data}{D-matrix of noisy data.}
##' \item{time}{vector containing time points}
##' \item{env}{vector specifying the experimental environment.}
##' \item{true.model}{vector specifying the target equation model.}
##' \item{target}{target variable.}
##' 
##' @export
##'
##' @import deSolve
##'
##' @author Niklas Pfister, Stefan Bauer and Jonas Peters
##'
##' @references
##' Pfister, N., S. Bauer, J. Peters (2018).
##' Identifying Causal Structure in Large-Scale Kinetic Systems
##' ArXiv e-prints (arXiv:1810.11776).
##'
##' @seealso The functions \code{\link{generate.data.maillard}} and
##'   \code{\link{generate.data.hidden}} allow to simulate ODE data
##'   from two additional models.
##'
##' @examples
##'
##' simulation.obj <- generate.data.targetmodel(env=rep(1:5, 3),
##'                                             L=15,
##'                                             d=5)
##'
##' D <- simulation.obj$simulated.data
##' fulldata <- simulation.obj$simulated.model
##' time <- simulation.obj$time
##' plot(time, D[1,1:length(time)], col="red", pch=19)
##' legend("topright", c("observations"),
##'        col=c("red"), pch=c(19))


generate.data.targetmodel <- function(env=rep(1,10),
                                      noise.sd=0.01,
                                      L=15,
                                      d=7,
                                      seed=NA){

  # set seed
  if(is.numeric(seed)){
    set.seed(seed)
  }

  # read out parameters
  n <- 10000
  n.env <- length(unique(env))
  reps <- length(env)
  time <- seq(0, 10, length.out=L)
  time_index <- floor(seq(1, n, length.out=L))
  simulated.data <- matrix(NA, reps, d*L)
  target <- d
  env.size <- as.vector(table(env))

  # Define function to generate random smooth functions
  smooth_fun <- function(par, n){
    fun <- function(t){
      return(par[1]/(1+exp(par[2]*t))+par[3]/(1+exp(par[4]*t)))
    }
    tvec <- seq(-3,3, length.out=n)
    fvec <- fun(tvec)
    return(fvec)
  }

  # Generate data
  for(a in 1:n.env){
    current_env <- unique(env)[a]
    # Generate random smooth predictor functions
    Xmat <- matrix(NA, n, d-1)
    for(i in 1:(d-1)){
      Xmat[, i] <- smooth_fun(rnorm(4, 0, 1), n)
    }
    
    # read out data for predictors
    noise_var <- apply(Xmat, 2, function(x) noise.sd*diff(range(x)))
    noiseterm <- matrix(rnorm(L*(d-1)*env.size[a], 0, rep(rep(noise_var, each=L), env.size[a])), env.size[a], L*(d-1), byrow=TRUE)
    simulated.data[(1:reps)[env==current_env], 1:((d-1)*L)] <- matrix(
      rep(as.vector(Xmat[time_index,]), env.size[a]), env.size[a], (d-1)*L, byrow=TRUE) + noiseterm

    # use numerical integration to generate Y
    Y1 <- c(0,cumsum(0.5*(Xmat[1:(n-1), 1] + Xmat[2:n, 1])))
    Y2 <- c(0,cumsum(0.5*(Xmat[1:(n-1), 2] + Xmat[2:n, 2])))
    Y <- (0.1*Y1+0.2*Y2)*0.001
    noise_var <- noise.sd*diff(range(Y))
    noiseterm <- matrix(rnorm(L*env.size[a], 0, noise_var), env.size, L)
    simulated.data[env==current_env, ((d-1)*L+1):(d*L)] <- matrix(
      rep(Y[time_index], env.size[a]), env.size[a], L, byrow=TRUE) + noiseterm
 
  }
  
  return(list(simulated.data=simulated.data,
              time=time,
              env=env,
              target=target,
              true.model=c(list(1), list(2), list(c(1, 2)))))
}

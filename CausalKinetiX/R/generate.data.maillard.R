##' Generate sample data from the Maillard reaction as specified by
##' Bio52 in the BioModel data base.
##'
##' For further details see the references.
##' @title Maillard reaction
##' @param target specifies which species is used as a target, needs
##'   to be an integer between 1 and 11.
##' @param env integer vector of length n encoding to which experiment
##'   each repetition belongs.
##' @param L number of time points for evaluation.
##' @param par.noise list of parameters that specify the added
##'   noise. \code{noise.sd} specifies the standard deviation of
##'   noise, \code{only.target.noise} specifies whether to only add
##'   noise to target and \code{relative} specifies if the size of the noise
##'   should be relative to size of variable (if TRUE standard
##'   deviation is given by par.noise$noise.sd*(x(t)-x(t-1))).
##' @param intervention string specifying type of
##'   intervention. Currently three type of interventions are
##'   implemented "initial" (only intervene on intial values),
##'   "blockreactions" (intervene by blocking random reactions) or
##'   "intial_blockreactions" (intervene on both initial values and
##'   blockreactions").
##' @param ode.solver specifies which ODE solver to use when solving
##'   ODE. Should be one of the methods from the \code{deSolve}
##'   package ("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
##'   "euler", "rk4", "ode23", "ode45", "radau", "bdf", "bdf_d",
##'   "adams", "impAdams", "impAdams_d", "iteration").
##' @param seed random seed. Does not work if a "Detected blow-up"
##'   warning shows up.
##' @param silent set to TRUE if no status output should be produced.
##' 
##' @return list consisting of the following elements
##' 
##' \item{simulated.data}{D-matrix of noisy data.}
##' \item{time}{vector containing time points}
##' \item{env}{vector specifying the experimental environment.}
##' \item{simulated.model}{object returned by ODE solver.}
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
##' Brands C. and van Boekel M. (2002).
##' Kinetic modeling of reactions in heated monosaccharide-casein systems.
##' Journal of agricultural and food chemistry, 50(23):6725–6739.
##'
##' @seealso The functions \code{\link{generate.data.hidden}} and
##'   \code{\link{generate.data.targetmodel}} allow to simulate ODE data
##'   from two additional models.
##'
##' @examples
##'
##' simulation.obj <- generate.data.maillard(target=1,
##'                                          env=rep(1:5, 3),
##'                                          L=15)
##'
##' D <- simulation.obj$simulated.data
##' fulldata <- simulation.obj$simulated.model
##' time <- simulation.obj$time
##' plot(fulldata[[1]][,1], fulldata[[1]][,2], type= "l", lty=2,
##'      xlab="time", ylab="concentration")
##' points(time, D[1,1:length(time)], col="red", pch=19)
##' legend("topright", c("true trajectory", "observations"),
##'        col=c("black", "red"), lty=c(2, NA), pch=c(NA, 19))


generate.data.maillard <- function(target,
                                   env=rep(1,10),
                                   L=15,
                                   par.noise=list(noise.sd=0.01,
                                                  only.target.noise=FALSE,
                                                  relativ=FALSE),
                                   intervention="initial_blockreactions",
                                   ode.solver="lsoda",
                                   seed=NA,
                                   silent=FALSE){
  # set seed
  if(is.numeric(seed)){
    set.seed(seed)
  }
    
  ###
  # Setting default parameters
  ###
  
  if(!exists("noise.sd",par.noise)){
    par.noise$noise.sd <- 0.01
  }
  if(!exists("only.target.noise",par.noise)){
    par.noise$only.target.noise <- TRUE
  }
  if(!exists("relativ",par.noise)){
    par.noise$relativ <- FALSE
  }

  
  ######################################
  #
  # Simulate data
  #
  ######################################
  
  
  ###
  # Initialize RHS
  ###

  reactions <- function(t, x, theta){
    
    dx1 <- -(theta[1]+theta[3])*x[1] + theta[2]*x[2] - theta[7]*x[1]*x[10]
    dx2 <- -(theta[2]+theta[4]+theta[5])*x[2] + theta[1]*x[1] - theta[10]*x[2]*x[10]
    dx3 <- theta[3]*x[1] + theta[4]*x[2]
    dx4 <- 2*theta[5]*x[2] - theta[6]*x[4]
    dx5 <- theta[6]*x[4] + theta[8]*x[7]
    dx6 <- theta[6]*x[4]
    dx7 <- -(theta[8]+theta[9])*x[7] + theta[7]*x[1]*x[10]
    dx8 <- theta[9]*x[7] - theta[11]*x[8] + theta[10]*x[2]*x[10]
    dx9 <- theta[3]*x[1] + theta[4]*x[2]
    dx10 <- theta[8]*x[7] - theta[7]*x[1]*x[10] - theta[10]*x[2]*x[10]
    dx11 <- theta[11]*x[8]
    
    return(list(c(dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx10,dx11)))
  }

  
  ###
  # Set parameters
  ###
  
  d <- 11
  time.grid <- seq(0, 100, by=0.005)
  time.index <- floor(seq(1, sqrt(length(time.grid)), length.out=L)^2)
  initial_obs <- c(160, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0)
  theta_obs <- c(0.01, 0.00509, 0.00047, 0.0011, 0.00712, 0.00439,
                 0.00018, 0.11134, 0.14359, 0.00015, 0.12514)
  n.env <- length(unique(env))
  
  included_reactions <- list(c(1, 2, 3, 7),
                             c(1, 2, 4, 5 ,10),
                             c(3, 4),
                             c(5, 6),
                             c(6, 8),
                             c(6),
                             c(7, 8, 9),
                             c(9, 10, 11),
                             c(3, 4),
                             c(7, 8, 10),
                             c(11))
  included_reactions <- included_reactions[[target]]
  
  true_set <- list(c(list(1), list(2), list(c(1, 10))),
                   c(list(1), list(2), list(c(2, 10))),
                   c(list(1), list(2)),
                   c(list(2), list(4)),
                   c(list(4), list(7)),
                   c(list(4)),
                   c(list(7), list(c(1, 10))),
                   c(list(7), list(8), list(c(2, 10))),
                   c(list(1), list(2)),
                   c(list(7), list(c(1, 10)), list(c(2, 10))),
                   c(list(8)))
  true_set <- true_set[[target]]
  
  ###
  # Define interventions
  ###  
  
  if(intervention == "initial"){
    intervention_fun <- function(){
      initial_int <- c(runif(1, 0, 360), runif(1, 0, 360), 0, 0, 0, 0, 0, 0, 0,
                       runif(1, 0, 30), 0)
      return(list(initial=initial_int,
                  theta=theta_obs))
    }
  }
  else if(intervention == "blockreactions"){
    intervention_fun <- function(){
      num_reactions <- d-length(included_reactions)
      r_vec <- (1:d)[-included_reactions]
      theta_int <- theta_obs
      theta_int[r_vec] <- theta_int[r_vec]*rbinom(num_reactions, 1, 1-2/num_reactions)
      return(list(initial=initial_obs,
                  theta=theta_int))
    }
  }
  else if(intervention == "initial_blockreactions"){
    intervention_fun <- function(){
      num_reactions <- d-length(included_reactions)
      r_vec <- (1:d)[-included_reactions]
      theta_int <- theta_obs
      theta_int[r_vec] <- theta_int[r_vec]*rbinom(num_reactions, 1, 1-3/num_reactions)
      initial_int <- c(runif(1, 0, 800), runif(1, 0, 800), 0, 0, 0, 0, 0, 0, 0,
                       runif(1, 0, 75), 0)
      return(list(initial=initial_int,
                  theta=theta_int))
    }
  }

  ###
  # Generate data from exact ODE and generate observations
  ###
  
  # set up required variables
  simulated.data <- matrix(NA, length(env), L*d)
  simulated.model <- vector("list", n.env)
  time <- time.grid[time.index]
  
  # iterate over environments
  for(i in 1:n.env){
    env.size <- sum(env == i)
    # perform intervention: initial condition
    if(i == 1){
      initial <- initial_obs
      theta <- theta_obs
    }
    else{
      int_data <- intervention_fun()
      initial <- int_data$initial
      theta <- int_data$theta
    }
    
    # solve ODE numerically
    if(!silent){
      print(paste("Currently solving ODE-system on environment", i))
    }
    
    # catch warnings from ODE solver --> not a good model
    simulated.model[[i]] <- ode(initial, time.grid, reactions, theta, method=ode.solver)
    if(!is.numeric(simulated.model[[i]]) | nrow(simulated.model[[i]])<length(time.grid)){
      if(!silent){
        print("Problem in ode-solver")
      }
      return(NA)
    }
    
    # select time points and add noise
    if(par.noise$only.target.noise){
      target.ind <- ((target-1)*L+1):(target*L)
      tmp <- simulated.model[[i]][time.index, -1, drop=FALSE]
      if(par.noise$relativ){
        noise_var <- par.noise$noise.sd*diff(range(tmp[, target]))+0.0000001
        noiseterm <- matrix(rnorm(L*env.size, 0, noise_var), env.size, L)
      }
      else{
        noiseterm <- matrix(rnorm(L*env.size, 0, target), env.size, L)
      }
      simulated.data[env==i, ] <- matrix(rep(tmp, env.size), env.size, L*d, byrow=TRUE)
      simulated.data[env==i, target.ind] <- simulated.data[env==i, target.ind] + noiseterm
    }
    else{
      tmp <- simulated.model[[i]][time.index, -1]
      if(par.noise$relativ){
        noise_var <- apply(tmp, 2, function(x) par.noise$noise.sd*diff(range(x)))+0.0000001
        noiseterm <- matrix(rnorm(L*d*env.size, 0, rep(rep(noise_var, each=L), env.size)), env.size, L*d, byrow=TRUE)
      }
      else{
        noiseterm <- matrix(rnorm(L*d*env.size, 0, par.noise$noise.sd), env.size, L*d)
      }
      simulated.data[env==i,] <-  matrix(rep(tmp, env.size), env.size, L*d, byrow=TRUE) + noiseterm
    }
  }

  ###
  # Check if a blow-up occured
  ###

  blowup <- vector("logical", n.env)
  for(i in 1:n.env){
    blowup[i] <- sum(abs(simulated.model[[i]][length(time.grid),-1])>10^8 |
                       is.na(abs(simulated.model[[i]][length(time.grid),-1])))>0
  }
  if(sum(blowup)>0){
    if(!silent){
      print("Detected blow-up")
    }
    return(NA)
  }
  

  return(list(simulated.data=simulated.data,
              time=time,
              env=env,
              simulated.model=simulated.model,
              true.model=true_set,
              target=target))
}

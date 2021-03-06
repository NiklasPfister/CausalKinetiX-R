##' Applies CausalKinetiX framework to rank a list models according to
##' their stability.
##'
##' This function scores a specified list of models and does not
##' include a variable ranking.
##' @title CausalKinetix.modelranking
##' @param D data matrix. Should have dimension n x (L*d), where n is
##'   the number of repetitions (over all experiments), L is the
##'   number of time points and d is the number of predictor
##'   variables.
##' @param times vector of length L specifying the time points at
##'   which data was observed.
##' @param env integer vector of length n encoding to which experiment
##'   each repetition belongs.
##' @param target integer specifing which variable is the target.
##' @param models list of models. Each model is specified by a list of
##'   vectors specifiying the variables included in the interactions
##'   of each term.
##' @param pars list of the following parameters: \code{pen.degree}
##'   (default 2) specifies the penalization degree in the smoothing
##'   spline, \code{num.folds} (default 2) number of folds used in
##'   cross-validation of smoothing spline, \code{include.vars}
##'   (default NA) specifies variables that should be included in each
##'   model, \code{include.intercept} (default FALSE) specifies
##'   whether to include a intercept in models, \code{average.reps}
##'   (default FALSE) specifies whether to average repetitions in each
##'   environment, \code{smooth.X} (default FALSE) specifies whether
##'   to smooth predictor observations before fitting, \code{smooth.Y}
##'   (default FALSE) specifies whether to smooth target observations
##'   before fitting, \code{regression.class} (default OLS) other
##'   options are signed.OLS, optim, random.forest,
##'   \code{sample.splitting} (default "loo") either leave-one-out
##'   (loo) or no splitting (none), \code{score.type} (default
##'   "mean_absolute") specifies the type of score funtion to use
##'   (note that "mean" and "max" are proportional scores which should
##'   be used if the noise variance is expected to change across
##'   experiments, if this is not the case "mean_absolute" and
##'   "max_absolute" are prefered as they also work for perfect
##'   unconstrained spline fits), \code{integrated.model} (default
##'   TRUE) specifies whether to fit the integrated or the derived
##'   model, \code{splitting.env} (default NA) an additonal
##'   environment vector used for scoring, \code{weight.vec} (default
##'   rep(1, length(env)) a weight vector used when
##'   scoring.type=="weighted.mean", \code{set.initial} (default
##'   FALSE) specifies whether to fix the initial value, \code{silent}
##'   (default TRUE) turn of additional output, \code{show.plot}
##'   (default FALSE) show diagnostic plots.
##' 
##' @return returns a list with the entries "scores" and "parameter"
##'
##' @export
##' 
##' @import quadprog randomForest pspline
##' @importFrom graphics plot lines
##'
##' @author Niklas Pfister, Stefan Bauer and Jonas Peters
##'
##' @references
##' Pfister, N., S. Bauer, J. Peters (2018).
##' Identifying Causal Structure in Large-Scale Kinetic Systems
##' ArXiv e-prints (arXiv:1810.11776).
##'
##' @seealso The function \code{\link{CausalKinetiX}} is a wrapper for
##'   this function that also computes the variable ranking and
##'   generates sensible classes of models.
##'
##' @examples
##' ## Generate data from Maillard reaction
##' simulation.obj <- generate.data.maillard(target=1,
##'                                          env=rep(1:5, 3),
##'                                          L=20,
##'                                          par.noise=list(noise.sd=1))
##' D <- simulation.obj$simulated.data
##' time <- simulation.obj$time
##' env <- simulation.obj$env
##' target <- simulation.obj$target
##'
##' ## Fit data to the following two models using CausalKinetiX:
##' ## 1: dy = theta_1*x_1 + theta_2*x_2 + theta_3*x_1*x_10 (true model)
##' ## 2: dy = theta_1*x_2 + theta_2*x_4 + theta_3*x_3*x_10 (wrong model)
##' ck.fit <- CausalKinetiX.modelranking(D, time, env, target,
##'                                      list(list(1, 2, c(1, 10)), list(2, 4, c(3, 10))))
##' print(ck.fit$scores)



CausalKinetiX.modelranking <- function(D,
                                       times,
                                       env,
                                       target,
                                       models,
                                       pars=list()){

  ############################
  #
  # setting the default values
  #
  ############################
  

  if(!exists("pen.degree", pars)){
    pars$pen.degree <- 2
  }
  if(!exists("num.folds", pars)){
    pars$num.folds <- 2
  }
  if(!exists("include.vars",pars)){
    pars$include.vars <- NA
  }
  if(!exists("include.intercept",pars)){
    pars$include.intercept <- FALSE
  }
  if(!exists("average.reps",pars)){
    pars$average.reps <- FALSE
  }
  if(!exists("smooth.X",pars)){
    pars$smooth.X <- FALSE
  }
  if(!exists("smooth.Y",pars)){
    pars$smooth.Y <- FALSE
  }
  if(!exists("regression.class", pars)){
    pars$regression.class <- "OLS"
  }
  if(!exists("sample.splitting",pars)){
    if(sum(length(unique(env))) == 1){
      pars$sample.splitting <- "none"
    }
    else{
      pars$sample.splitting <- "loo"
    }
  }
  if(!exists("score.type",pars)){
    pars$score.type <- "mean_absolute"
  }
  if(!exists("integrated.model",pars)){
    pars$integrated.model <- TRUE
  }
  if(!exists("splitting.env",pars)){
    pars$splitting.env <- NA
  }
  if(!exists("weight.vec",pars)){
    pars$weight.vec <- rep(1, length(env))
  }
  if(!exists("set.initial",pars)){
    pars$set.initial <- FALSE
  }
  if(!exists("silent",pars)){
    pars$silent <- TRUE
  }
  if(!exists("show.plot",pars)){
    pars$show.plot <- FALSE
  }

  ## Parameter consistency checks
  if(pars$smooth.Y & !is.na(pars$splitting.env)){
    stop("If smooth.Y is TRUE, splitting.env needs to be NA.")
  }
  if(pars$sample.splitting != "none" & length(unique(env)) == 1){
    warning(paste("pars$sample.splitting == ", pars$sample.splitting,
                  " does not make sense since only one environment",
                  " has been given. It has been set to 'none'.",
                  sep=""))
    pars$sample.splitting <- "none"
  }
  
  ############################
  #
  # initialize
  #
  ############################

  # remove all variables not involved in any models
  unique_vars <- unique(c(unlist(lapply(models, unlist)), target, pars$include.vars))
  unique_vars <- sort(unique_vars[!is.na(unique_vars)])
  L <- length(times)
  matching_perm <- order(c(unique_vars, setdiff(1:(ncol(D)/L), unique_vars)))
  target <- matching_perm[target]
  models <- lapply(models, function(x) lapply(x, function(y) matching_perm[y]))
  keep_ind <- do.call(c, lapply(unique_vars, function(i) ((i-1)*L+1):(i*L)))
  D <- D[,keep_ind]
  if(!is.na(pars$include.vars)[1]){
    pars$include.vars <- matching_perm[pars$include.vars]
  }
  
  # read out parameters (as in paper)
  L <- length(times)
  n <- nrow(D)
  d <- ncol(D)/L
  num.folds <- pars$num.folds
  pen.degree <- pars$pen.degree
  if(length(pen.degree)==1){
    pen.degree <- rep(pen.degree, 2)
  }
  include.vars <- pars$include.vars
  include.intercept <- pars$include.intercept
  average.reps <- pars$average.reps
  smooth.X <- pars$smooth.X
  smooth.Y <- pars$smooth.Y
  score.type <- pars$score.type
  silent <- pars$silent
  sample.splitting <- pars$sample.splitting
  splitting.env <- pars$splitting.env
  weight.vec <- pars$weight.vec
  regression.class <- pars$regression.class
  show.plot <- pars$show.plot

  # check whether to include products and interactions
  non_zero_mod <- sapply(models, function(x) length(x)>0)
  products <- sum(sapply(models[non_zero_mod], function(model)
    sum(sapply(model, function(term) length(term) > length(unique(term)))) > 0)) > 0
  interactions <- sum(sapply(models[non_zero_mod], function(model)
    sum(sapply(model, function(term) length(unique(term)) > 1)) > 0)) > 0

  # sort environments to increasing order
  if(is.na(splitting.env)[1]){
    splitting.env <- env
  }
  env_order <- order(env)
  D <- D[env_order,]
  splitting.env <- splitting.env[env_order]
  env <- env[env_order]
  if(smooth.Y){
    weight.vec <- weight.vec[order(unique(env))]
  }
  else{
    weight.vec <- weight.vec[env_order]
  }
  
  # construct DmatY
  target_ind <- ((target-1)*L+1):(target*L)
  DmatY <- D[,target_ind]

  # add interactions to D-matrix if interactions==TRUE
  if(interactions | products | !is.na(include.vars)[1]){
    if(!is.na(include.vars)[1]){
      interactions <- TRUE
      products <- TRUE
    }
    include.obj <- extend_Dmat(D, L, d, n,
                               products=products,
                               interactions=interactions,
                               include.vars=include.vars)
    D <- include.obj$Dnew
    ordering <- include.obj$ordering
    dtot <- ncol(D)/L
  }
  else{
    dtot <- d
    ordering <- as.list(1:d)
  }
  
  ##################################################
  #
  # Step 1: Average Repetitions & Smooth predictors
  #
  ##################################################

  ## Averaging
  if(average.reps){
    D <- apply(D, 2, function(col) sapply(split(col, env),
                                          function(x) mean(x, na.rm=TRUE)))
    DmatY <- apply(DmatY, 2, function(col) sapply(split(col, env),
                                                  function(x) mean(x, na.rm=TRUE)))
    splitting.env <- splitting.env[!duplicated(env)]
    env <- unique(env)
    n <- length(env)
  }

  ## Smoothing
  Dlist <- vector("list", n)
  if(smooth.X){
    for(i in 1:n){
      Dlist[[i]] <- matrix(D[i,], nrow=L, ncol=dtot, byrow=FALSE)
      # smooth X-values
      for(j in 1:dtot){
        na_ind <- is.na(Dlist[[i]][,j])
        fit <- sm.spline(times[!na_ind], Dlist[[i]][!na_ind,j], norder=2, cv=TRUE)
        Dlist[[i]][,j] <- predict(fit, times, nderiv=0)
      }
    }
  }
  else{
    for(i in 1:n){
      Dlist[[i]] <- matrix(D[i,], nrow=L, ncol=dtot, byrow=FALSE)
    }
  }

  ######################################
  #
  # Step 2: Fit reference model
  #
  #####################################
  
  if(smooth.Y){
    # initialize variables
    unique_env <- unique(env)
    num_env <- length(unique_env)
    Ylist <- vector("list", num_env)
    envtimes <- vector("list", num_env)
    dYlist <- vector("list", num_env)
    RSS_A <- vector("numeric", num_env)
    UpDown_A <- vector("numeric", num_env)
    RSS3_A <- vector("numeric", num_env)
    lambda <- vector("numeric", num_env)
    initial_values <- vector("numeric", num_env)
    times_new <- min(times)
    if(!silent | show.plot){
      Ya <- vector("list", num_env)
      Yb <- vector("list", num_env)
      times_new <- seq(min(times), max(times), length.out=100)
    }
    else{
      Ya <- NA
      Yb <- NA
    }
    for(i in 1:length(unique_env)){
      # fit higher order model to compute derivatives (based on sm.smooth)
      Ylist[[i]] <- as.vector(DmatY[env==unique_env[i],])
      len_env <- sum(env==unique_env[i])
      envtimes[[i]] <- rep(times, each=len_env)
      fit <- constrained.smoothspline(Ylist[[i]],
                                      envtimes[[i]],
                                      pen.degree[2],
                                      constraint="none",
                                      times.new=times_new,
                                      num.folds=num.folds,
                                      lambda="optim")
      if(!silent | show.plot){
        Yb[[i]] <- fit$smooth.vals.new
      }
      dYlist[[i]] <- rep(fit$smooth.deriv, len_env)
      # compute differences for intergrated model fit
      if(pars$integrated.model){
        dYlist[[i]] <- as.vector(apply(DmatY[env==unique_env[i],,drop=FALSE], 1, diff))
      }
      if(!silent | show.plot){
        Ya[[i]] <- fit$smooth.vals.new
      }
      lambda[i] <- fit$pen.par
      if(pars$set.initial){
        initial_values[i] <- fit$smooth.vals[1,1]
      }
      else{
        initial_values[i] <- NA
      }
      RSS_A[i] <- sum((rep(fit$smooth.vals, each=len_env)-Ylist[[i]])^2)
      UpDown_A[i] <- fit$smooth.vals[length(fit$smooth.vals)]
      RSS3_A[i] <- NA
    }
  }
  else{
    # initialize variables
    Ylist <- vector("list", n)
    dYlist <- vector("list", n)
    RSS_A <- vector("numeric", n)
    UpDown_A <- vector("numeric", n)
    RSS3_A <- vector("numeric", n)
    lambda <- vector("numeric", n)
    initial_values <- vector("numeric", n)
    times_new <- min(times)
    if(!silent | show.plot){
      Ya <- vector("list", n)
      Yb <- vector("list", n)
      times_new <- seq(min(times), max(times), length.out=100)
    }
    else{
      Ya <- NA
      Yb <- NA
    }
    for(i in 1:n){
      # fit higher order model to compute derivatives (based on sm.smooth)
      Ylist[[i]] <- as.vector(DmatY[i,])
      na_ind <- is.na(Ylist[[i]])
      fit <- constrained.smoothspline(Ylist[[i]][!na_ind],
                                      times[!na_ind],
                                      pen.degree[1],
                                      constraint="none",
                                      times.new=times_new,
                                      num.folds=num.folds,
                                      lambda="optim")
      if(!silent | show.plot){
        Yb[[i]] <- fit$smooth.vals.new
      }
      dYlist[[i]] <- fit$smooth.deriv
      # compute differences for intergrated model fit
      if(pars$integrated.model){
        dYlist[[i]] <- diff(Ylist[[i]])
      }
      
      # fit the reference model and compute penalty par and RSS
      if(pen.degree[1] != pen.degree[2]){
        fit <- constrained.smoothspline(Ylist[[i]][!na_ind],
                                        times[!na_ind],
                                        pen.degree[2],
                                        constraint="none",
                                        times.new=times_new,
                                        num.folds=num.folds,
                                        lambda="optim")
      }
      if(!silent | show.plot){
        Ya[[i]] <- fit$smooth.vals.new
      }
      lambda[i] <- fit$pen.par
      if(pars$set.initial){
        initial_values[i] <- fit$smooth.vals[1,1]
      }
      else{
        initial_values[i] <- NA
      }
      RSS_A[i] <- sum(fit$residuals^2)
      UpDown_A[i] <- fit$smooth.vals[length(fit$smooth.vals)]
      RSS3_A[i] <- sum(fit$residuals[c(1, floor(L/2), L)]^2)
    }
  }
  
  ######################################
  #
  # Step 3: Pre-compute models 
  #
  ######################################

  # convert models to modelstot
  if(is.list(models[[1]])){
    modelstot <- models
  }
  else{
    modelstot <- list()
    for(k in 1:length(models)){
      modelstot <- append(modelstot, list(models[k]))
    }
  }
  
  # intialize variables
  num.models <- length(modelstot)
  data_list <- vector("list", num.models)
  varlist <- vector("list", num.models)
  data_list2 <- vector("list", num.models)
  
  # iteration over all potential models
  for(model in 1:num.models){
    if(length(unlist(modelstot[[model]]))==0){
      model.index <- numeric()
    }
    else{
      model.index <- match(modelstot[[model]], ordering)
    }
    
    ## collect predictors X
    Xlist <- vector("list", n)
    num.pred <- length(unlist(modelstot[[model]]))+include.intercept
    if(num.pred==0){
      num.pred <- 1
      for(i in 1:n){
        Xlist[[i]] <- matrix(1, L, 1)
      }
    }
    else{
      for(i in 1:n){
        Xlist[[i]] <- Dlist[[i]][, model.index, drop=FALSE]
        if(include.intercept){
          Xlist[[i]] <- cbind(Xlist[[i]], matrix(1, L, 1))
        }
      }
    }
    data_list[[model]] <- Xlist
    
    # compute predictors for integrated model fit
    if(pars$integrated.model){
      Xlist2 <- vector("list", n)
      num.pred <- length(unlist(modelstot[[model]]))+include.intercept
      if(num.pred==0){
        num.pred <- 1
        for(i in 1:n){
          Xlist2[[i]] <- matrix(1, L-1, 1)
        }
      }
      else{
        for(i in 1:n){
          Xlist2[[i]] <- Dlist[[i]][, model.index, drop=FALSE]
          tmp <- (Xlist2[[i]][1:(L-1),,drop=FALSE]+Xlist2[[i]][2:L,,drop=FALSE])/2
          Xlist2[[i]] <- tmp*matrix(rep(diff(times), ncol(tmp)), L-1, ncol(tmp))
          if(include.intercept){
            Xlist2[[i]] <- cbind(Xlist2[[i]], matrix(1, L-1, 1))
          }
        }
      }
      data_list2[[model]] <- Xlist2
    }
  }  
  
  ######################################
  #
  # Step 4: Compute score
  #
  ######################################
  
  # Iterate over all models and compute score
  scores <- vector("numeric", num.models)
  parameters <- vector("list", num.models)
  for(model in 1:num.models){
    Xlist <- data_list[[model]]
    Xlist2 <- data_list2[[model]]
    # output
    if(!silent){
      print(paste("Scoring model ", toString(modelstot[model])))
    }
    ## Compute constraint and RSS on constrained model
    if(sample.splitting=="none"){
      ###
      # Without sample splitting
      ###
      unique_env <- unique(splitting.env)
      num.env <- length(unique_env)
      RSS_B <- vector("numeric", length(RSS_A))
      UpDown_B <- vector("numeric", length(RSS_A))
      RSS3_B <- vector("numeric", length(RSS_A))
      X <- do.call(rbind, Xlist)
      Xpred <- do.call(rbind, Xlist)
      dY <- unlist(dYlist)
      subenv.ind <- lapply(1:num.env, function(k) rep(1:sum(splitting.env==unique_env[k]), each=L))
      loo.ind <- rep(splitting.env, each=L)
      # adjust for missing obs in integrated model fit
      if(!is.null(Xlist2)){
        X <- do.call(rbind, Xlist2)
        loo.ind2 <- loo.ind
        loo.ind <- rep(splitting.env, each=L-1)
      }
      count <- 1
      ### PLOT
      if(show.plot){
        minyplot <- vector("numeric", length(RSS_A))
        maxyplot <- vector("numeric", length(RSS_A))
        constrained_fit <- vector("list", length(RSS_A))
      }
      # Fit model on all data
      # Classical OLS regression (main effects and interactions)
      if(regression.class=="OLS"){
        fit <- lm(dY ~ -1 + X)
        # remove coefficients resulting from singular fits (perfectly correlated predictors)
        coefs <- coefficients(fit)
        coefs[is.na(coefs)] <- 0
        fitted_dY <- Xpred %*% matrix(coefs, ncol(Xpred), 1)
        parameters[[model]] <- coefs
      }
      # OLS regression with sign constraints on parameters
      else if(regression.class=="signed.OLS"){
        len_model <- length(modelstot[[model]])
        tmp <- sapply(modelstot[[model]], function(x) x[1])
        ind <- rep(1, len_model)
        ind[tmp==target] <- -1
        # define quadProg parameters
        bvec <- rep(0, len_model)
        Amat <- diag(ind)
        dvec <- matrix(dY, 1, nrow(X)) %*% X
        Dmat <- (t(X) %*% X)
        sc <- norm(Dmat,"2")
        fit <- solve.QP(Dmat/sc, dvec/sc, Amat, bvec, meq=0)
        coefs <- fit$solution
        fitted_dY <- Xpred %*% matrix(coefs, ncol(Xpred), 1)
        parameters[[model]] <- coefs
      }
      # OLS regression with pruning based on score
      else if(regression.class=="optim"){
        # define loss function
        loss_fun <- function(beta, Y, X, env_vec, ind){
          coefs <- (10^beta)*ind
          tmp_vec <- (Y-X %*% matrix(coefs, ncol(X), 1))^2
          return(max(sapply(split(tmp_vec, env_vec), mean)))
        }
        # compute starting value using quadratic program
        len_model <- length(modelstot[[model]])
        tmp <- sapply(modelstot[[model]], function(x) x[1])
        ind <- rep(1, len_model)
        ind[tmp==target] <- -1
        bvec <- rep(10^(-10), len_model)
        Amat <- diag(ind)
        dvec <- matrix(dY, 1, nrow(X)) %*% X
        Dmat <- (t(X) %*% X)
        sc <- norm(Dmat,"2")
        fit <- solve.QP(Dmat/sc, dvec/sc, Amat, bvec, meq=0)
        coefs <- fit$solution
        # perform optimization
        opt.res <- optim(log(coefs*ind)/log(10), loss_fun,
                         Y=dY,
                         X=X,
                         env_vec=splitting.env,
                         ind=ind,
                         lower=-10, upper=5)
        coefs <- (10^opt.res$par)*ind
        fitted_dY <- Xpred %*% matrix(coefs, ncol(Xpred), 1)
        parameters[[model]] <- coefs
      }
      # Random forest regression
      else if(regression.class=="random.forest"){
        fit <- randomForest(y ~ ., data=cbind(as.data.frame(X), data.frame(y=dY)))
        fitted_dY <- predict(fit, newdata=as.data.frame(Xpred))
        parameters[[model]] <- fit
      }
      # Wrong regression.class
      else{
        stop("Specified regression.class does not exist. Use OLS, signed.OLS, optim or random.forest.")
      }
      # compute score using splitting environment
      if(!smooth.Y){
        for(i in 1:num.env){
          num.reps <- sum(splitting.env==unique_env[i])
          env_ind <- loo.ind2==unique_env[i]
          for(j in 1:num.reps){
            fitted_dY_tmp <- fitted_dY[env_ind][subenv.ind[[i]]==j]
            fit <- constrained.smoothspline(Ylist[[count]],
                                            times,
                                            pen.degree[2],
                                            constraint="fixed",
                                            derivative.values=fitted_dY_tmp,
                                            initial.value=initial_values[count],
                                            times.new=times_new,
                                            num.folds=num.folds,
                                            lambda=lambda[count])
            
            RSS_B[count] <- sum(fit$residuals^2)
            UpDown_B[count] <- fit$smooth.vals[length(fit$smooth.vals)]
            RSS3_B[count] <- sum(fit$residuals[c(1, floor(L/2), L)]^2)
            ### PLOT
            if(show.plot){
              constrained_fit[[count]] <- fit$smooth.vals.new
            }
            count <- count+1
          }
        }
      }
      else{
        for(i in 1:num.env){
          num.reps <- sum(splitting.env==unique_env[i])
          env_ind <- loo.ind2==unique_env[i]
          fitted_dY_tmp <- sapply(split(fitted_dY[env_ind], rep(1:L, num.reps)), mean)
          fit <- constrained.smoothspline(Ylist[[i]],
                                          rep(times, each=num.reps),
                                          pen.degree[2],
                                          constraint="fixed",
                                          derivative.values=fitted_dY_tmp,
                                          initial.value=initial_values[i],
                                          times.new=times_new,
                                          num.folds=num.folds,
                                          lambda=lambda[i])
          RSS_B[i] <- sum((rep(fit$smooth.vals, each=num.reps)-Ylist[[i]])^2)
          UpDown_B[i] <- fit$smooth.vals[length(fit$smooth.vals)]
          RSS3_B[i] <- NA
          if(show.plot){
            constrained_fit[[i]] <- fit$smooth.vals.new
          }
        }
      }
      ### PLOT
      if(show.plot){
        if(!smooth.Y){
          for(i in 1:num.env){
            env_ind <- splitting.env == unique_env[i]
            Y1plot <- unlist(Ylist[env_ind])
            times1 <- rep(times, sum(env_ind))
            miny <- min(c(Y1plot, unlist(Ya[env_ind]),
                          unlist(Yb[env_ind]), unlist(constrained_fit[env_ind])))
            maxy <- max(c(Y1plot, unlist(Ya[env_ind]),
                          unlist(Yb[env_ind]), unlist(constrained_fit[env_ind])))
            # plot
            plot(times1, Y1plot, xlab="times", ylab="concentration",
                 ylim = c(miny, maxy))
            which_ind <- which(env_ind)
            for(k in which_ind){
              lines(times_new, Ya[[k]], col="red")
              lines(times_new, Yb[[k]], col="blue")
              lines(times_new, constrained_fit[[k]], col="green")
            }
            readline("Press enter")
          }
        }
        else{
          for(i in 1:num.env){
            env_ind <- splitting.env == unique_env[i]
            times1 <- rep(times, each=sum(env_ind))
            L <- length(times)
            miny <- min(c(Ylist[[i]], Ya[[i]],
                          Yb[[i]], constrained_fit[[i]]))
            maxy <- max(c(Ylist[[i]], Ya[[i]],
                          Yb[[i]], constrained_fit[[i]]))
            # plot
            plot(times1, Ylist[[i]], xlab="times", ylab="concentration",
                 ylim = c(miny, maxy))
            for(k in 1:sum(env_ind)){
              lines(times_new, Ya[[i]], col="red")
              lines(times_new, Yb[[i]], col="blue")
              lines(times_new, constrained_fit[[i]], col="green")
            }
            readline("Press enter")
          }
        }
      }
    }
    else if(sample.splitting=="loo"){
      ###
      # With leave-one-out sample splitting
      ###
      unique_env <- unique(splitting.env)
      num.env <- length(unique_env)
      RSS_B <- vector("numeric", length(RSS_A))
      UpDown_B <- vector("numeric", length(RSS_A))
      RSS3_B <- vector("numeric", length(RSS_A))
      X <- do.call(rbind, Xlist)
      dY <- unlist(dYlist)
      subenv.ind <- lapply(1:num.env, function(k) rep(1:sum(splitting.env==unique_env[k]), each=L))
      loo.ind <- rep(splitting.env, each=L)
      # adjust for missing obs in integrated model fit
      if(!is.null(Xlist2)){
        X2 <- X                            
        X <- do.call(rbind, Xlist2)
        loo.ind2 <- loo.ind
        loo.ind <- rep(splitting.env, each=L-1)                               
      }
      count <- 1
      ### PLOT
      if(show.plot){
        minyplot <- vector("numeric", length(RSS_A))
        maxyplot <- vector("numeric", length(RSS_A))
        constrained_fit <- vector("list", length(RSS_A))
      }
      parameters[[model]] <- vector("list", num.env)
      for(i in 1:num.env){
        # compute derivative constraint using OLS
        dYout <- dY[loo.ind!=unique_env[i]]
        Xout <- X[loo.ind!=unique_env[i],,drop=FALSE]
        env_out <- loo.ind[loo.ind!=unique_env[i]]
        Xin <- X[loo.ind==unique_env[i],,drop=FALSE]
        # adjust for missing vlaue in integrated model fit
        if(!is.null(Xlist2)){
          Xin <- X2[loo.ind2==unique_env[i],,drop=FALSE]
        }
        # dealing with na
        Xout_nona <- Xout[rowSums(is.na(Xout))==0,,drop=FALSE]
        dYout_nona <- dYout[rowSums(is.na(Xout))==0]
        env_out <- env_out[rowSums(is.na(Xout))==0]
        # Classical OLS regression (main effects and interactions)
        if(regression.class=="OLS"){
          fit <- lm(dYout_nona ~ -1 + Xout_nona)
          # remove coefficients resulting from singular fits (perfectly correlated predictors)
          coefs <- coefficients(fit)
          coefs[is.na(coefs)] <- 0
          fitted_dY <- Xin %*% matrix(coefs, ncol(Xin), 1)
          parameters[[model]][[i]] <- coefs
        }
        # OLS regression with sign constraints on parameters
        else if(regression.class=="signed.OLS"){
          len_model <- length(modelstot[[model]])
          tmp <- sapply(modelstot[[model]], function(x) x[1])
          ind <- rep(1, len_model)
          ind[tmp==target] <- -1
          # define quadProg parameters
          bvec <- rep(0, len_model)
          if(length(ind)==1){
            Amat <- as.matrix(ind)
          }
          else{
            Amat <- diag(ind)
          }
          dvec <- matrix(dYout_nona, 1, nrow(Xout_nona)) %*% Xout_nona
          Dmat <- (t(Xout_nona) %*% Xout_nona)
          sc <- norm(Dmat,"2")
          fit <- solve.QP(Dmat/sc, dvec/sc, Amat, bvec, meq=0)
          coefs <- fit$solution
          fitted_dY <- Xin %*% matrix(coefs, ncol(Xin), 1)
          parameters[[model]][[i]] <- coefs
        }
        # OLS regression with pruning based on score
        else if(regression.class=="optim"){
          # define loss function
          loss_fun <- function(beta, Y, X, env_vec, ind){
            coefs <- (10^beta)*ind
            tmp_vec <- (Y-X %*% matrix(coefs, ncol(X), 1))^2
            return(mean(sapply(split(tmp_vec, env_vec), mean)))
          }
          # compute starting value using quadratic program
          len_model <- length(modelstot[[model]])
          tmp <- sapply(modelstot[[model]], function(x) x[1])
          ind <- rep(1, len_model)
          ind[tmp==target] <- -1
          bvec <- rep(10^(-10), len_model)
          Amat <- diag(ind)
          dvec <- matrix(dYout_nona, 1, nrow(Xout_nona)) %*% Xout_nona
          Dmat <- (t(Xout_nona) %*% Xout_nona)
          sc <- norm(Dmat,"2")
          fit <- solve.QP(Dmat/sc, dvec/sc, Amat, bvec, meq=0)
          coefs <- fit$solution
          # perform optimization
          opt.res <- optim(log(coefs*ind)/log(10), loss_fun,
                           Y=dYout_nona,
                           X=Xout_nona,
                           env_vec=env_out,
                           ind=ind,
                           lower=-10, upper=5)
          coefs <- (10^opt.res$par)*ind
          fitted_dY <- Xin %*% matrix(coefs, ncol(Xin), 1)
          parameters[[model]][[i]] <- coefs
        }
        # Random forest regression
        else if(regression.class=="random.forest"){
          fit <- randomForest(y ~ ., data=cbind(as.data.frame(Xout_nona), data.frame(y=dYout_nona)))
          fitted_dY <- predict(fit, newdata=as.data.frame(Xin))
          parameters[[model]][[i]] <- fit
        }
        # Wrong regression.class
        else{
          stop("Specified regression.class does not exist. Use OLS, signed.OLS, optim or random.forest.")
        }
        
        # Fit individual models with derivative constraint
        if(!smooth.Y){
          for(j in 1:sum(splitting.env==unique_env[i])){
            na_ind <- is.na(Ylist[[count]])|(rowSums(is.na(Xin[subenv.ind[[i]]==j,,drop=FALSE]))>0)
            fitted_dY_tmp <- fitted_dY[subenv.ind[[i]]==j]
            fit <- constrained.smoothspline(Ylist[[count]][!na_ind],
                                            times[!na_ind],
                                            pen.degree[2],
                                            constraint="fixed",
                                            derivative.values=fitted_dY_tmp[!na_ind],
                                            initial.value=initial_values[count],
                                            times.new=times_new,
                                            num.folds=num.folds,
                                            lambda=lambda[count])
            
            RSS_B[count] <- sum(fit$residuals^2)
            UpDown_B[count] <- fit$smooth.vals[length(fit$smooth.vals)]
            RSS3_B[count] <- sum(fit$residuals[c(1, floor(L/2), L)]^2)
            ### PLOT
            if(show.plot){
              constrained_fit[[count]] <- fit$smooth.vals.new
            }
            count <- count+1
          }
        }
        else{
          len_env <- sum(splitting.env==unique_env[i])
          fitted_dY_tmp <- sapply(split(fitted_dY, rep(1:L, len_env)), mean)
          fit <- constrained.smoothspline(Ylist[[i]],
                                          rep(times, each=len_env),
                                          pen.degree[2],
                                          constraint="fixed",
                                          derivative.values=fitted_dY_tmp,
                                          initial.value=initial_values[i],
                                          times.new=times_new,
                                          num.folds=num.folds,
                                          lambda=lambda[i])
          RSS_B[i] <- sum((rep(fit$smooth.vals, each=len_env)-Ylist[[i]])^2)
          UpDown_B[i] <- fit$smooth.vals[length(fit$smooth.vals)]
          RSS3_B[i] <- NA
          if(show.plot){
            constrained_fit[[i]] <- fit$smooth.vals.new
          }
        }
      }
      ### PLOT
      if(show.plot){
        if(!smooth.Y){
          for(i in 1:num.env){
            env_ind <- splitting.env == unique_env[i]
            Y1plot <- unlist(Ylist[env_ind])
            times1 <- rep(times, sum(env_ind))
            miny <- min(c(Y1plot, unlist(Ya[env_ind]),
                          unlist(Yb[env_ind]), unlist(constrained_fit[env_ind])))
            maxy <- max(c(Y1plot, unlist(Ya[env_ind]),
                          unlist(Yb[env_ind]), unlist(constrained_fit[env_ind])))
            # plot
            plot(times1, Y1plot, xlab="times", ylab="concentration",
                 ylim = c(miny, maxy))
            which_ind <- which(env_ind)
            for(k in which_ind){
              lines(times_new, Ya[[k]], col="red")
              lines(times_new, Yb[[k]], col="blue")
              lines(times_new, constrained_fit[[k]], col="green")
            }
            readline("Press enter")
          }
        }
        else{
          for(i in 1:num.env){
            env_ind <- splitting.env == unique_env[i]
            times1 <- rep(times, each=sum(env_ind))
            L <- length(times)
            miny <- min(c(Ylist[[i]], Ya[[i]],
                          Yb[[i]], constrained_fit[[i]]))
            maxy <- max(c(Ylist[[i]], Ya[[i]],
                          Yb[[i]], constrained_fit[[i]]))
            # plot
            plot(times1, Ylist[[i]], xlab="times", ylab="concentration",
                 ylim = c(miny, maxy))
            lines(times_new, Ya[[i]], col="red")
            lines(times_new, Yb[[i]], col="blue")
            lines(times_new, constrained_fit[[i]], col="green")
            readline("Press enter")
          }
        }
      }
    }
    else{
      stop("Specified sample.splitting does not exist. Use none or loo.")
    }
    ## compute score
    if(score.type=="max"){
      if(min(abs(RSS_A))<10^-10){
        warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.type max_absolut.")
      }
      score <- max((RSS_B-RSS_A)/RSS_A)
    }
    else if(score.type=="mean"){
      if(min(abs(RSS_A))<10^-10){
        warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.type mean_absolut.")
      }
      score <- mean((RSS_B-RSS_A)/RSS_A)
    }
    else if(score.type=="mean_absolute"){
      score <- mean(RSS_B)
    }
    else if(score.type=="max_absolute"){
      score <- max(RSS_B)
    }
    else if(score.type=="mean.weighted"){
      if(min(abs(RSS_A))<10^-10){
        warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.types mean_absolut or max_absolute.")
      }
      score <- mean(weight.vec*(RSS_B-RSS_A)/RSS_A)
    }
    else if(score.type=="stability"){
      pairs <- combn(1:length(RSS_B), 2)
      score <- max(apply(pairs, 2, function(x) abs(RSS_B[x[1]]-RSS_B[x[2]])))
    }
    else{
      stop("Specified score.type does not exist. Use max, mean or max-mean.")
    }
    # Output
    if(!silent){
      print(paste("Model has a score of", score))
    }
    scores[model] <- score
  }

  ## Return results
  return(list(scores=scores,
              parameters=parameters))

}

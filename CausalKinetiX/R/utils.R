####
# Small helper functions that are not exported
###

##' @import utils glmnet

## Lasso based on deltaY vs intX
ode_integratedlasso_rank_vars <- function(D,
                                          times,
                                          env=NULL,
                                          target,
                                          pars = list()){
  
  if(!exists("silent",pars)){
    pars$silent <- TRUE
  }
  if(!exists("interactions",pars)){
    pars$interactions <- TRUE
  }
  if(!exists("rm.target",pars)){
    pars$rm.target <- TRUE
  }
  if(!exists("signed",pars)){
    pars$signed <- FALSE
  }
  
  L <- length(times)
  d <- ncol(D)/L
  n <- nrow(D)

  Xint <- matrix(NA, (L-1)*n, d)
  deltaY <- vector("numeric", (L-1)*n)
  for(i in 1:n){
    deltaY[((i-1)*(L-1)+1):(i*(L-1))] <- diff(D[i, ((target-1)*L+1):(target*L)])
  }
  for(j in 1:d){
    for(i in 1:n){
      tmp <- D[i, ((j-1)*L+1):(j*L)]
      Xint[((i-1)*(L-1)+1):(i*(L-1)),j] <- (tmp[1:(L-1)]+tmp[2:L])/2*diff(times)
    }
  }

  # remove NAs
  na_ind <- is.na(deltaY) | (rowSums(is.na(Xint))>0)
  deltaY <- deltaY[!na_ind]
  Xint <- Xint[!na_ind,]
  
  # Perform lasso
  if(pars$interactions){
    var.names <- c(as.list(1:d), combn(d, 2, simplify = FALSE))
    Xint.interactions <- matrix(NA, n*(L-1), length(var.names))
    Xint.interactions[,1:d] <- Xint
    for(i in (d+1):length(var.names)){
      Xint.interactions[,i] <- Xint[,var.names[[i]][1]] * Xint[,var.names[[i]][2]]
    }
    if(pars$signed){
      fit <- nnlasso(Xint.interactions, deltaY, family="normal", path=TRUE, min.lambda = 1e-15)
      sel.matrix <- (t(fit$coef) > max(fit$coef[1,]))
    }
    else{
      fit <- glmnet(Xint.interactions, deltaY)
      sel.matrix <- (fit$beta != 0)
    }
    first.entrance <- apply(sel.matrix, MARGIN = 1, FUN = which.max)
    # find all rows without ones and set first entrance to Inf
    first.entrance[which(apply(sel.matrix, MARGIN = 1, FUN = sum) == 0)] <- Inf
    ranking <- order(first.entrance)
    ranking <- unique(unlist(var.names[ranking]))
  }
  else{
    if(pars$signed){
      fit <- nnlasso(Xint, deltaY, family="normal", path=TRUE, min.lambda = 1e-15)
      sel.matrix <- (t(fit$coef) > max(fit$coef[1,]))
    }
    else{
      fit <- glmnet(Xint, deltaY)
      sel.matrix <- (fit$beta != 0)
    }
    first.entrance <- apply(sel.matrix, MARGIN = 1, FUN = which.max)
    # find all rows without ones and set first entrance to Inf
    first.entrance[which(apply(sel.matrix, MARGIN = 1, FUN = sum) == 0)] <- Inf
    ranking <- order(first.entrance)
  }

  if(pars$rm.target){
    ranking <- ranking[ranking != target]
  }
  
  
  return(list(ranking = ranking))
}



construct_models <- function(D, L, d, n, target, times,
                             maineffects.models, screening,
                             interactions, products, include.vars,
                             max_preds, expsize, env=NULL, signed=FALSE){

  ## Main-Effect and Full-Effect models depends on maineffects.models
  if(!maineffects.models){
    # construct variable vector
    if(!is.na(include.vars)[1]){
      vv <- (1:d)[-include.vars[include.vars!=0]]
    }
    else{
      vv <- 1:d
    }
    ## Decide which terms to keep depending on screening, interactions and products
    if(is.numeric(screening)){
      tmp <- extend_Dmat(D, L, d, n, products, interactions, include.vars)
      Dfull <- tmp$Dnew
      ordering <- tmp$ordering
      res <- ode_integratedlasso_rank_vars(Dfull,
                                           times,
                                           env=env,
                                           target, list(interactions=FALSE,
                                                        rm.target=FALSE,
                                                        signed=signed))$ranking
      keep_terms <- ordering[res[1:screening]]
      ## print(keep_terms)
      num_terms <- screening
    }
    else{
      keep_terms <- as.list(vv)
      # add interactions
      if(interactions){
        keep_terms <- append(keep_terms, combn(vv, 2, simplify=FALSE))
      }
      # add products
      if(products){
        keep_terms <- append(keep_terms, lapply(vv, function(i) c(i,i)))
      }
      # include.vars
      if(!is.na(include.vars)[1]){
        keep_terms.new <- list()
        for(i in 1:length(keep_terms)){
          for(var in include.vars){
            if(var == 0){
              keep_terms.new <- append(keep_terms.new, keep_terms[[i]])
            }
            else{
              tmp_term <- sort(c(var, keep_terms[[i]]))
              keep_terms.new <- append(keep_terms.new, list(tmp_term))
            }
          }
        }
        keep_terms <- append(keep_terms.new, as.list(include.vars))
      }
      num_terms <- length(keep_terms)
    }
    ## Construct models
    if(max_preds){
      models <- list(list())
      for(k in 1:(expsize+1)){
        models <- append(models, lapply(combn(keep_terms, k, simplify=FALSE), as.list))
      }
    }
    else{
      models <- lapply(combn(keep_terms, expsize+1, simplify=FALSE), as.list)
    }
  }
  else{
    if(!is.na(include.vars)){
      warning("include.vars is not defined for maineffects.models==FALSE")
    }
    ## Construct models
    if(max_preds){
      models <- list(list())
      for(k in 1:(expsize+1)){
        models <- append(models, lapply(combn(as.list(1:d), k, simplify=FALSE), as.list))
      }
    }
    else{
      models <- lapply(combn(as.list(1:d), expsize+1, simplify=FALSE), as.list)
    }
    # add interactions and products
    for(i in 1:length(models)){
      if(length(models[[i]])>1){
        if(interactions){
          models[[i]] <- c(models[[i]],
                           lapply(combn(unlist(models[[i]]), 2, simplify=FALSE), c))
        }
        if(products){
          models[[i]] <- c(models[[i]],
                           lapply(lapply(unique(unlist(models[[i]])), function(i) c(i, i)), c))
        }
      }
    }
    num_terms <- d  
  }
    
  # return output
  result <- list(models=models,
                 num_terms=num_terms)
  return(result)
}


extend_Dmat <- function(D, L, d, n,
                        products,
                        interactions,
                        include.vars){

  # construct variable vector
  if(!is.na(include.vars)[1]){
    vv <- (1:d)[-include.vars[include.vars!=0]]
    vv_ind <- numeric()
    for(var in include.vars[include.vars!=0]){
      vv_ind <- c(vv_ind, ((var-1)*L+1):(var*L))
    }
    vv_ind <- (1:(L*d))[-vv_ind]
  }
  else{
    vv <- 1:d
    vv_ind <- (1:(L*d))
  }
  dnew <- length(vv)

  # initialize
  num_vars <- dnew + products*dnew + interactions*choose(dnew, 2)
  if(length(vv) == 1){
    num_vars <- dnew + products*dnew
  }
  Dnew <- matrix(NA, n, L*num_vars)
  Dnew[,1:(dnew*L)] <- D[,vv_ind]
  ordering <- vector("list", num_vars)
  ordering[1:dnew] <- as.list(vv)
  count <- dnew
  
  # add interactions
  if(interactions & length(vv) >1){
    for(i in 1:(dnew-1)){
      indi <- ((vv[i]-1)*L+1):(vv[i]*L)
      for(j in (i+1):dnew){
        ordering[[count+1]] <- c(vv[i], vv[j])
        indj <- ((vv[j]-1)*L+1):(vv[j]*L)
        Dnew[,(count*L+1):((count+1)*L)] <- D[, indi]*D[, indj]
        count <- count+1
      }
    }
  }
  
  # add products
  if(products){
    for(i in 1:dnew){
      indi <- ((vv[i]-1)*L+1):(vv[i]*L)
      ordering[[count+1]] <- c(vv[i], vv[i])
      Dnew[,(count*L+1):((count+1)*L)] <- D[, indi]^2
      count <- count + 1
    }
  }

  # include variables to every term
  if(!is.na(include.vars)[1]){
    Dfinal <- matrix(NA, nrow(Dnew), length(include.vars)*ncol(Dnew)+sum(include.vars!=0)*L)
    ordering_final <- vector("list", length(include.vars)*length(ordering)+sum(include.vars!=0))
    count <- 1
    for(var in include.vars){
      if(var != 0){
        Dfinal[,((count-1)*L+1):(count*L)] <- D[,((var-1)*L+1):(var*L)]
        ordering_final[[count]] <- var
        count <- count + 1
      }
    }
    for(j in 1:length(ordering)){
      for(var in include.vars){
        if(var == 0){
          ordering_final[[count]] <- ordering[[j]]
          Dfinal[,((count-1)*L+1):(count*L)] <- Dnew[,((j-1)*L+1):(j*L)]
        }
        else{
          ordering_final[[count]] <- sort(c(var, ordering[[j]]))
          Dfinal[,((count-1)*L+1):(count*L)] <- (D[,((var-1)*L+1):(var*L)]*
                                                   Dnew[,((j-1)*L+1):(j*L)])
        }
        count <- count + 1
      }
    }
    ordering <- ordering_final
    Dnew <- Dfinal
  }

  return(list(Dnew=Dnew,
              ordering=ordering))
}

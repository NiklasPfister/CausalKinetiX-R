##' Applies CausalKinetiX framework to rank variables and models according to their stability
##'
##' For further details see the references.
##' @title CausalKinetiX
##' @param D data matrix. Should have dimension n x (L*d), where n is
##'   the number of repetitions (over all experiments), L is the
##'   number of time points and d is the number of predictor
##'   variables.
##' @param times vector of length L specifying the time points at
##'   which data was observed.
##' @param env integer vector of length n encoding to which experiment
##'   each repetition belongs.
##' @param target integer specifing which variable is the target.
##' @param pars list of parameters.
##' @param models list of models. Each model is specified by a list of
##'   vectors specifiying the variables included in the interactions
##'   of each term.
##' 
##' @return object of class 'CausalKinetiX' consisting of the following
##'   elements
##'
##' \item{models}{list of the individually scored models.}
##' \item{model.scores}{vector containing the score for each model.}
##' \item{variable.scores}{vector containing the score of each variable.}
##' \item{ranking}{vector specifying the ranking of each variable.}
##' 
##' @export
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
##'
##' x <- 4

CausalKinetiX <- function(D,
                          times,
                          env,
                          target,
                          models=NA,
                          pars){

  # set defaults for pars
  if(!exists("max.preds", pars)){
    pars$max.preds <- TRUE
  }
  if(!exists("expsize", pars)){
    pars$expsize <- 2
  }
  if(!exists("interactions", pars)){
    pars$interactions <- TRUE
  }
  if(!exists("include.vars", pars)){
    pars$include.vars <- NA
  }
  if(!exists("products", pars)){
    pars$products <- FALSE
  }
  if(!exists("stability.cutoff", pars)){
    pars$stability.cutoff <- 3/4
  }
  if(!exists("stability.selection", pars)){
    pars$stability.selection <- TRUE
  }
  if(!exists("individual.models", pars)){
    pars$individual.models <- FALSE
  }
  if(!exists("rm.target",pars)){
    pars$rm.target <- TRUE
  }
  if(!exists("screening",pars)){
    pars$screening <- NA
  }
  if(!exists("K",pars)){
    pars$K <- NA
  }

  
  # for random forest regression remove interactions
  if(exists("regression.class", pars)){
    if(pars$regression.class == "random.forest"){
      pars$interactions <- FALSE
      pars$interactions.Y <- FALSE
      pars$include.intercept <- FALSE
    }
  }

  # read out variables
  n <- NROW(D)
  L <- length(times)
  d <- NCOL(D)/L  

  # Check whether a list of models was specified else generate models
  if(is.na(models)){
    constructed_mods <- construct_models(D, L, d, n, target, times,
                                         pars$individual.models,
                                         pars$screening,
                                         pars$interactions,
                                         pars$products,
                                         pars$include.vars,
                                         pars$max.preds,
                                         pars$expsize,
                                         env)
    models <- constructed_mods$models
    if(is.na(pars$K)){
      pars$K <- constructed_mods$num_terms-pars$expsize
    }
  }

  # check whether parameter K was specified
  if(is.na(pars$K)){
    warning("K was not specified and the default does not make sense for arbitrary lists of models. It was set 1, but this can invalidate the variable ranking.")
    pars$K <- 1
  }
      
  ###
  # Compute model scores
  ###
  
  model.scores <- CausalKinetiX.modelranking(D, times, env, target, models, pars)

  ###
  # Rank variables
  ###
  
  Mlen <- length(models)
  Mjlen <- sapply(1:d, function(j) sum(sapply(models, function(mod) j %in% unlist(mod))))
  # compute p-values based on hypergeometric distribution
  best_mods <- models[order(model.scores)[1:K]]
  var_scores <- sapply(1:d,
                       function(x) sum(x==unlist(
                         lapply(best_mods, function(y) unique(unlist(y))))))/K
  var_pvals <- sapply(1:d, function(j) phyper(var_scores[j]*K, Mjlen[j],
                                              Mlen-Mjlen[j],
                                              K, lower.tail=FALSE))
  var_pvals[Mjlen==0] <- Inf
  if(pars$rm.target){
    idx <- order(var_pvals[-target])
    ranking <- ((1:d)[-target])[idx]
    scores <- (var_pvals[-target])[idx]
  }
  else{
    idx <- order(var_pvals)
    ranking <- ((1:d))[idx]
    scores <- (var_pvals)[idx]
  }

  # output results
  return(list(models=models,
              model.scores=model.scores,
              variable.scores=scores,
              ranking=ranking))
  
}

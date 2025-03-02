
#' @keywords internal
get_likelis = function(mus, ws, X, y){
  XW = X %*% diag(ws)
  likelis = (t(y) %*% XW) %*% mus - sum(log(1 + exp(XW %*% mus )))
  return(likelis)
}

#' @keywords internal
expit = function(x){
  1/(1+exp(-1*x))
}

#' @keywords internal
my_norm = function(x){
  sqrt(sum(x^2))
}

#' @keywords internal
grad_likelis = function(mus, ws, X, y){
  XW = X %*% diag(ws)
  t(XW) %*% (y - expit(XW %*% mus))
}

#' @keywords internal
hes_likelis = function(mus, ws, X, y){
  XW = X %*% diag(ws)
  m = expit(XW %*% mus)
  w = as.vector(m * (1-m))
  hessian = -1*t(XW*w) %*% XW
  return(hessian)
}

#' @keywords internal
get_grad = function(mus, ws, sigma, X, y){
  as.vector(grad_likelis(mus, ws, X, y)) - 1/sigma^2 * mus
}

get_hes = function(mus, ws, sigma, X, y){
  hes = hes_likelis(mus, ws, X, y)
  hes - diag(1/sigma^2, nrow = nrow(hes), ncol = ncol(hes))
}

#' @keywords internal
iterate_mu = function(mus, ws, sigma, X, y)
{
  inv = solve(get_hes(mus, ws, sigma, X, y))
  return(mus - inv %*% get_grad(mus, ws, sigma, X, y))
}

#' @keywords internal
solve_mus = function(ws, sigma, X, y)
{
  MAXITER = 100
  mus_hats = rep(0, ncol(X))

  eps = 10e-5; diff = 1

  for (ITER in 1:MAXITER) {
    mus_hats_new = iterate_mu(mus_hats, ws, sigma, X, y)
    diff = my_norm(mus_hats_new - mus_hats)
    mus_hats = mus_hats_new

    if (diff < eps) {
      break;
    }
  }

  mus_hats[ws < 0.2] = 0
  return(as.vector(mus_hats))
}

#' @keywords internal
get_tau <- function(mus, w_jk, sigma, rho, X, y, k){
  tau = get_likelis(mus, w_jk, X, y) -1/(2*sigma^2)*my_norm(mus)^2 - 1/2 * log(det(-1*get_hes(mus, w_jk, sigma, X, y))) + k*log(rho/(1-rho))
  return(tau)
}

#' @keywords internal
solve_ws = function(mus, ws,sigma, rho, X, y){

  w_updated = ws

  for(j in 1:ncol(X)){
    w1 = ws; w1[j] = 1
    w0 = ws; w0[j] = 0
    tau1 = get_tau(mus, w1, sigma, rho, X, y, 1)
    tau0 = get_tau(mus, w0, sigma, rho, X, y, 0)
    w_updated[j] = 1/(1+exp(tau0- tau1))
  }
  return(w_updated)
}

#' @title Collapsed VB solver for coefficients
#'
#' @description Estimates the regression coefficients \eqn{ \boldsymbol \mu } and model selection mask \eqn{ \mathbf w } by taking the Maximum a priori estimates of the posterior density of the  random variables \eqn{ \boldsymbol \beta } and \eqn{ \mathbf w } with prior densities \eqn{\beta_j \sim \mathrm{Normal(0,\sigma^2)} } and \eqn{w_j \sim \mathrm{Bernoulli(\rho)} } using the reverse collapsed variational bayesian approximation.
#' @usage cvb_model(sigma, rho, X, y, add.intercept = TRUE, MAXITER = 100,
#'  EPS = 1e-5, verboseIter = FALSE)
#' @param sigma A hyperparameter tuning the regression coefficients' Normal(0,\code{sigma^2}) prior density.
#' @param rho A hyperparameter representing the sparsity of the Bernoulli(\code{rho}) prior density for the model selection mask.
#' @param X A matrix or vector containing the variables in the model.
#' @param y A vector containing the observed binary response for fitting the model.
#' @param add.intercept An optional logical for whether to add an intercept column of \code{1}s to the \code{X} input.
#' @param MAXITER The maximum number of iterations allowed before concluding that the estimated coefficients do not converge.
#' @param EPS The threshold for residuals under which the model coefficients are considered converged.
#' @param verboseIter An optional logical for whether to print information about each iteration.
#'
#' @return A list with the regression coefficients \code{mus}, model selection coefficients \code{ws} and logical \code{added.intercept}.
#' Can be used in the \code{predict_cvb} function to perform classification.
#'
#' @examples
#'
#' # Initial Data
#' X = c(66,70,69,68,67,72,73,70,57,63,70,78,67)
#' y = c(0,1,0,0,0,0,0,0,1,1,1,0,0)
#' new.X = c(66,70,68,59)
#'
#' # Fit the model
#' fitted_model <- cvb_model(sigma = 10, rho = 0.1, X, y,
#'                 add.intercept = TRUE, MAXITER = 100, EPS = 1e-5, verboseIter = FALSE)
#' fitted_model
#'
#' @export

cvb_model <- function(sigma, rho, X, y, add.intercept = TRUE, MAXITER = 100, EPS = 1e-5, verboseIter = FALSE){

  if(sigma <= 0 | rho<= 0 | rho >= 1){
    if(rho<= 0 | rho >= 1){
      cat("Sparsity parameter must be strictly in interval (0,1).\n")
    }
    if(sigma <= 0){
      cat("Regression prior density variance must be postive.\n")
    }
    stop("Please check your inputs.")
  }

  X = as.matrix(X)

  if(dim(X)[1] != length(y)){
    cat("Error: number of observations in predictor matrix X does not response vector y.\n")
    stop("Please check your training data.")
  }

  if(add.intercept){
    X = cbind(rep(1,nrow(X)), X)
  }

  mus = rep(0, ncol(X))
  ws = rep(1, ncol(X))

  diff = 1
  its = 1

  for(ITER in 1:MAXITER){
    mus_next = solve_mus(ws, sigma, X, y)
    ws_next = solve_ws(mus_next, ws, sigma, rho, X, y)

    diff_next = sum(c(mus_next - mus, ws_next - ws)^2)
    change = 100*(diff_next - diff)/diff
    diff = diff_next
    mus = mus_next
    ws = ws_next
    if(verboseIter){
      cat(paste("Finished iteration step", ITER, "with residual", ifelse(diff < 0.01, round(diff,6), round(diff,3)) ,". Change:", round(change,3), "%\n"))
    }
    its = its + 1

    if (diff < EPS)
    {
      s = (paste("Model converged in", its - 1, "iterations.\n"))
      message(s)
      break
    }
  }
  if (its - 1 == MAXITER){
    s = paste("Model was not able to converge in", its-1, " iterations. Results may not be accurate.")
    message(s)
  }

  res = list(mus = mus, ws = ws, added.intercept = add.intercept)
  return(res)
}

#' @title Classifier for fitted models
#'
#' @description Uses a fitted logistic model from \code{cvb_model} and new data to perform binary classification. The probability threshold is \eqn{0.5}.
#'
#' @param fitted_model A fitted model from \code{cvb_model}
#' @param data A matrix or vector containing new observations for classification.
#' @param add.intercept An optional logical for whether to add an intercept column of 1s to the X input.
#' If left blank, the function will assume the option used to train the model.
#' @param MAXITER The maximum number of iterations allowed before concluding that the estimated coefficients do not converge.
#' @param EPS The threshold for residuals under which the model coefficients are considered converged.
#'
#' @return A vector of predictions.
#' @examples
#' # Initial Data
#' X = c(66,70,69,68,67,72,73,70,57,63,70,78,67)
#' y = c(0,1,0,0,0,0,0,0,1,1,1,0,0)
#' new.X = c(66,70,68,59)
#'
#' # Fit the model
#' cvbdl::cvb_model(sigma = 10, rho = 0.1, X, y,
#'                 add.intercept = TRUE, MAXITER = 100, EPS = 1e-5, verboseIter = FALSE) -> m; m
#'
#' # Perform Prediction
#' cvbdl::predict_cvb(fitted_model = m, X = new.X)
#'
#' @export

predict_cvb <- function(fitted_model, X, add.intercept = NA){

  a.i = fitted_model$added.intercept
  if(is.na(add.intercept)){
    if(a.i == TRUE){
      X = cbind(rep(1,nrow(as.matrix(X))), X)
      }
  } else if(add.intercept == TRUE) {
      X = cbind(rep(1,nrow(as.matrix(X))), X)
  } else if(add.intercept == FALSE) {
    a.i = FALSE
  } else{
    cat("Error: add.intercept should be a logical or NA value \n")
    stop("Please check your inputs")
  }

  ws = ifelse(fitted_model$ws < 0.5, 0, 1)

  if(dim(as.matrix(X))[2] != length(ws)){
    cat("Error: dimension mismatch.\n")
    if(dim(as.matrix(X))[2] - length(ws) == 1){
      stop("Did you intend add.intercept = FALSE ?")
    } else if(dim(as.matrix(X))[2] - length(ws) == -1){
      stop("Did you intend add.intercept = TRUE ?")
    }
  }

  mus = fitted_model$mus
  probs = expit(X %*% diag(ws) %*% mus)
  classified = ifelse(probs > 0.5, 1, 0)
  return(as.vector(classified))
}

#' @title Cross-validation for fitted models
#'
#' @description Estimates the MSE of the model via repeated cross-validation.
#'
#' @param sigma A hyperparameter tuning the regression coefficients' Normal(0,\code{sigma^2}) prior density.
#' @param rho A hyperparameter representing the sparsity of the Bernoulli(\code{rho}) prior density for the model selection mask.
#' @param X A matrix or vector containing the variables in the model.
#' @param y A vector containing the observed binary response for fitting the model.
#' @param add.intercept An optional logical for whether to add an intercept column of \code{1}s to the \code{X} input.
#' @param reps The number of repetitions for cross-validation.
#' @param seed The seed for drawing the folds from the supplied data.
#' @param MAXITER The maximum number of iterations allowed before concluding that the estimated coefficients do not converge.
#' @param EPS The threshold for residuals under which the model coefficients are considered converged.
#' @param verboseIter An optional logical for whether to print information about each cross-validation test.
#'
#' @return A vector of cross-validation errors corresponding to each repetition.
#' @examples
#' # Initial Data
#' X = c(66,70,69,68,67,72,73,70,57,63,70,78,67)
#' y = c(0,1,0,0,0,0,0,0,1,1,1,0,0)
#'
#' # Perform cross-validation
#' cross_validate(folds = 5, sigma = 10, rho = 0.1, X = X, y = y,
#' reps = 3, seed = 0, add.intercept = TRUE, MAXITER = 10, EPS = 1e-5, verboseIter = TRUE)
#' @export

cross_validate <- function(folds, sigma, rho, X, y, add.intercept = TRUE,
                           reps = 50, seed = 0, MAXITER = 100, EPS = 1e-5, verboseIter = FALSE){

  if(sigma <= 0 | rho<= 0 | rho >= 1){
    if(rho<= 0 | rho >= 1){
      cat("Sparsity parameter must be strictly in interval (0,1).\n")
    }
    if(sigma <= 0){
      cat("Regression prior density variance must be postive.\n")
    }
    stop("Please check your inputs.")
  }

  X = as.matrix(X)

  if(dim(X)[1] != length(y)){
    cat("Error: number of observations in predictor matrix X does not response vector y.\n")
    stop("Please check your training data.")
  }

  if(add.intercept){
    X = cbind(rep(1,nrow(X)), X)
  }

  K = folds
  mi = MAXITER
  ep = EPS
  error_vector = vector(length = reps)
  shift = seed

  samp_num = nrow(X) %/% K
  remainder = nrow(X) %% K
  n = nrow(X)-remainder

  for (j in 1:reps){

    set.seed(j + shift)
    rand_samp = sample(c(1:n), n, replace = FALSE)
    rand_mat = matrix(rand_samp, nrow = K, ncol = samp_num, byrow = FALSE)

    err = vector(length = K)

    for (i in 1:K){
      test_set = X[rand_mat[i,],]
      test_outcomes = y[rand_mat[i,]]
      train_set = X[rand_mat[-i,],]
      train_outcomes = y[rand_mat[-i,]]

      fit = suppressMessages(cvb_model(sigma, rho, train_set, train_outcomes, MAXITER = mi, EPS = ep))
      results = predict_cvb(fit, test_set)
      err[i] =  sum(results != test_outcomes)
    }
    em = sum(err)/n
    if(verboseIter){
      cat(paste("Finished iteration ", j, " and obtained cross-validation error ", round(em,3),".\n"))
    }
    error_vector[j] = em
  }

  s = paste("Completed with mean", K, "- fold cross-validation error: ", round(mean(error_vector),3),".")
  message(s)
  return(error_vector)
}

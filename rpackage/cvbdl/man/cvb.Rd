\name{cvbdl}
\alias{cvbdl}
\title{Collapsed Variational Bayes Logistic Regression}
\description{
Given continuous predictor data and binary response, provides methods for performing simultaneous model selection and regression for logistic regression using Reverse Collapsed Variational Bayes methods.
}
\usage{
cvb_model(sigma, rho, X, y,
          add.intercept = FALSE, MAXITER = 100,
          EPS = 1e-5, verboseIter = FALSE)

predict_cvb(fitted_model, X, add.intercept = NA)

cross_validate(folds, sigma, rho, X, y,
               reps = 1, seed = 0, MAXITER = 100,
               EPS = 1e-5, verboseIter = FALSE)
}

\arguments{
  \item{sigma}{a hyperparameter tuning the regression coefficients' Normal(0,\code{sigma^2}) prior density.}
  \item{rho}{a hyperparameter representing the sparsity of the Bernoulli(\code{rho}) prior density for the model selection mask. }
  \item{X}{a matrix or vector containing the variables in the model.}
  \item{y}{a vector containing the observed binary response for fitting the model.}
  \item{add.intercept}{an optional logical for whether to add an intercept column of \code{1}s to the \code{X} input.}
  \item{MAXITER}{maximum number of iterations allowed before concluding that the estimated coefficients do not converge.}
  \item{EPS}{residuals threshold under which the model coefficients are considered converged.}
  \item{verboseIter}{an optional logical for whether to print information about each iteration.}
  \item{fitted_model}{a fitted model from the \code{cvb_model} funtion.}
  \item{folds}{the number of folds to be used in cross-validation.}
  \item{reps}{the number of repetitions for cross-validation.}
  \item{seed}{the seed for drawing the folds from the supplied data.}
}

\details{
For information encoded in the format
\deqn{
 \mathbf{X} = \begin{bmatrix}
  1 & \mathrm{predictor\_1}_1 & \cdots & \mathrm{predictor\_p}_1 \\
  \vdots &  \vdots &  \ddots & \vdots \\
  1 &  \mathrm{predictor\_1}_n&  \cdots & \mathrm{predictor\_p}_n  \\
\end{bmatrix} \text{ and } \mathbf{y} =\begin{bmatrix} y_1 \\ \vdots \\ y_n \,,\\
\end{bmatrix}
}
where \eqn{ y_i \in \{0,1\} } represents the binary response, we simultaneously regression coefficients \eqn{\boldsymbol \mu = (\mu_0,\dotsb, \mu_p)} and model selection coefficients \eqn{ \mathbf{w} = (w_0,\dotsb, w_p)} for the logistic model
\deqn{ P(\mathbf{Y} = \mathbf{1}) = \mathrm{expit}(\mathbf X \mathrm{diag}(\mathbf w) \boldsymbol \mu)  := \frac{e^{\mathbf X \mathrm{diag}(\mathbf w)  \boldsymbol \beta}}{ 1 -e^{ \mathbf X \mathrm{diag}(\mathbf w)  \boldsymbol \mu}} \,. }

We estimate the vectors \eqn{ \boldsymbol \mu } and \eqn{ \mathbf w } by taking the Maximum a priori estimates of the posterior density of the  random variables \eqn{ \boldsymbol \beta } and \eqn{ \mathbf w } with prior densities \eqn{\beta_j \sim \mathrm{Normal(0,\sigma^2)} } and \eqn{w_j \sim \mathrm{Bernoulli(\rho)} } using the reverse collapsed variational bayesian approximation.
}

\examples{
# Initial Data
X = c(66,70,69,68,67,72,73,70,57,63,70,78,67)
y = c(0,1,0,0,0,0,0,0,1,1,1,0,0)
new.X = c(66,70,68,59)

# Fit the model
cvbdl::cvb_model(sigma = 10, rho = 0.1, X, y,
  add.intercept = TRUE, MAXITER = 100, EPS = 1e-5, verboseIter = FALSE) -> m; m

# Perform Prediction
cvbdl::predict_cvb(fitted_model = m, X = new.X)
}


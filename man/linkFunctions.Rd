\name{Links}
\alias{Links}
\alias{links}
\alias{links}

\title{ Logit and probit links for generalised linear modelling }

\description{
  
Generalised linear models with binomial resonse variables ("logistic regression") use a link function to link the response on a (0,1) scale to a linear predictor on (-Inf, Inf). The cannonical link is the logistic ("logit") function, which has some nice theoretical properties and can be interpreted as the log of the odds of success. Other links are available, notably the cumulative standard normal ("probit") link, which allows for Gibbs sampling with truncated normal distributions. For that reason, several of the Bayesian estimation functions in \pkg{wiqid} use the probit link.

The form of the logit and probit links are shown in the figure below.

\figure{linkFunctions.jpg}

Both curves are symmetric, with probability = 0.5 when the linear predictor = 0. The pobit curve is steeper, so coefficients in a probit regression will be smaller than those in a logit regression (by a factor of about 1.7).

}


\section{Standardisation}{

Continuous variables are standardised to facilitate writing models, and optimisation. Standardisation also means that the size of regression coefficients directly reflect the importance of the corresponding variables.

Binary variables coded as 0/1 and dummy variables are not changed. To make continuous variables comparable with these, they are centred by subtracting the mean, and then divided by \emph{two times} their standard deviation (Gelman, 2008).

Note that all numerical inputs (ie, \code{is.numeric == TRUE}) that appear in the \code{data} argument will be standardised, including binary variables coded as 0/1. Variables coded as TRUE/FALSE or as factors are not affected.

The same standardisation strategy is used for both Bayesian and maximum likelihood functions.
}

\section{Priors for logistic regression coefficients}{

Following Gelman et al (2008), we use independent Cauchy priors with centre 0 and scale 10 for the intercept and scale 2.5 for all other coefficients.

}

\section{Priors for probabilities}{

We use independent Uniform(0, 1) priors for proababilities.

}

\references{
  
Gelman, A. (2008) Scaling regression inputs by dividing by two standard deviations. \emph{Statistics in Medicine}, 27, 2865-2873

Gelman, Jakulin, Pittau and Su (2008) A weakly informative default prior distribution for logistic and other regression models. \emph{Annals of Applied Statistics} 2, 1360-1383.

}

\keyword{ manip }
\keyword{ models }


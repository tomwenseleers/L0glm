\name{extended_families}
\alias{extended_families}
\title{Additional Likelihood Families}
\description{
Use \code{mgcv::nb()} for negative binomial fitting with fixed or estimated
shape. See \link{L0glm_negative_binomial} for mathematics and an example.
Use \code{mgcv::multinom(K)} for joint multinomial IRLS with CLR
coefficients; see \link{L0glm_multinomial}.
These family-specific loops are serial. Ordinary GLMs retain their existing paths.
}

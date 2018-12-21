\name{recrisk.mixture}
\alias{recrisk.mixture}
\title{
  Recurrence Risk based on a mixure of low and high risk families  
}
\description{
 Recurrence risk assuming a mixture of low and high risk families based on a mixture of truncated binomial densities
}
\usage{
 recrisk.mixture(s, a, k, max.iter=1e4, eps=1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{vector of sizes of families}
  \item{a}{vector of number of affected members per family; same length
    as vector s}
  \item{k}{number of affecteds to condition on; the ascertainment criterion,
    e.g., at least k affecteds}
  \item{max.iter}{maximum number of EM iterations}
  \item{eps}{EM stop criterion value}
}
\details{
Model heterogeneity  of risk across families by a mixture of truncated
binomial distributions, assuming two groups: high and low risk.
}
\value{
A list with the following values
  \item{recrisk.high}{disease recurrence risk for high-risk group}
  \item{recrisk.low}{disease recurrence risk for low-risk group}
  \item{prob.high.group}{the mixture probability that a familiy is in
    the high-risk group (a scalar value)}
  \item{prob.fam.high}{a vector of posterior probabilities that families
  are in the high-risk group (same length as vector s)}
  \item{iter}{Number of EM iterations}
}
\references{
Schaid DJ, McDonnell SK, Thibodeau SN. Familial recurrence risk with
varying amount of family history. Under review at Genetic Epidemiology, 2018.
}
\author{
Schaid DJ, Sinnwell JP.
}
\examples{
s <- c(4,  8,  6,  6,  10, 8,  7,  8,  5,  6,  10, 4,  4,  8,  6,  8,  4,  5,  9, 
       9,  5,  4,  4,  7,  3,  3,  9,  5,  3,  3,  4,  6,  8,  3,  5,  8,  6,  8, 9,  3)
a <- c(3, 4, 2, 5, 4, 5, 3, 4, 3, 3, 8, 2, 0, 3, 4, 5, 2, 3, 4, 2, 1, 1,
       1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 1, 0, 2, 3, 0)
recrisk.mixture(s,a,k=1)
}


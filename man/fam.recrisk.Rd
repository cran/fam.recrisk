\name{fam.recrisk}
\alias{fam.recrisk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Familial risk of recurrent disease
}
\description{
  Given vectors of family sizes and number of affecteds per family, calculates the risk
  of disease recurrence
}
\usage{
fam.recrisk(s, a, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{vector of sizes of families}
  \item{a}{vector of number of affected members per family; same length as vector s}
  \item{k}{number of affecteds to condition on; the ascertainment criterion,
    e.g., at least k affecteds}
}
\details{
Estimate risk of disease for a person with at least k affected relatives, where k can be specified in order to determine how risk
varies with k. Estimates are based on a truncated binomial distribution
where families with less than k affected relatives are truncated from
the sample.   If the sample of families in not a random sample from a population, the prevalence will
be biased, and merely represent the frequency of disease in the sample.
}
\value{
A list with the following values
  \item{phat}{disease recurrence risk}
  \item{prev}{estimate of disease prevalence}
  \item{var.phat}{variance of phat}
  \item{var.prev}{variance of prev}
  \item{covar}{covariance between phat and prev}
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
fam.recrisk(s,a,1)

}


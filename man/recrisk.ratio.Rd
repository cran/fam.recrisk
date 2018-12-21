\name{recrisk.ratio}
\alias{recrisk.ratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Recurrence risk ratio of familial recurrence risk
}
\description{
 Estimate the recurrence risk ratio of familial recurrence risk divided by population 
 prevalence, based on results returned from fam.recrisk
}
\usage{
recrisk.ratio(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{ list returned from fam.recrisk}
}
\details{
The recurrence risk is based on a truncated binomial density, and the
prevalence is the frequency of disease in the sample. If the sample of
families in not a random sample from a population, the prevalence will
be biased, and merely represent the frequency of disease in the sample.
}
\value{
A list with the following values
  \item{rr}{ratio of recurrence risk divided by population prevalence}
  \item{var.rr}{variance of rr}
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
  fit <- fam.recrisk(s,a,1)
  recrisk.ratio(fit)
}
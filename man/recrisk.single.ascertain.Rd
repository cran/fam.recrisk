\name{recrisk.single.ascertain}
\alias{recrisk.single.ascertain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Familial recurrence risk assuming single proband 
}
\description{
 Estimate familial recurrence risk, assuming there is a single
 identified proband in a family
}
\usage{
recrisk.single.ascertain(s, a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{vector of sizes of families}
  \item{a}{vector of number of affected members}
}
\details{
When affected subjects (i.e., probands) attend a clinic independent of
their family size and independent of their family history, and if the
family history provided by probands is complete (i.e., reporting of
family history does not depend on the number of affected members), then
we can view the probands and their families as sampled under complete
ascertainment. This function estimates the disease recurrence risk when
there is a single identified proband per family.
}
\value{
  A list with the following values
  \item{phat}{disease recurrence risk}
  \item{var.phat}{variance of phat}
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
recrisk.single.ascertain(s,a)

}


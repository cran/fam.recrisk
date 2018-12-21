
fam.recrisk <- function(s,a,k){
  ## Description: estimate familial recurrence risk, the probability a person will have disease, conditional
  ## on family size (s) and at least k affected relatives, when there are a affecteds in a family.
  
  
  prev <- sum(a)/sum(s)
  g <- a/prev - (s-a)/(1-prev)
  dg <- -a/prev^2 -(s-a)/(1-prev)^2
  
  asc  <- 1*(a >= k)
  
  phat <- sum(  asc*(a*(a-k)) ) / sum(  asc*(a*(s-1) - s*(k-1)) )
  f <- asc * ( a*(a-k)/phat - ( (s-a)*(a-k + 1) ) /(1-phat) )
  df  <-  asc*( -a*(a-k)/phat^2 - ( (s-a)*(a-k + 1) )  / (1-phat)^2 )
  
  var.phat  <- sum(f^2)/sum(df)^2
  var.prev <-  sum(g^2)/sum(dg)^2
  covar <- sum(f*g)/(sum(df)*sum(dg))
  
  
  return(list(phat=phat, prev=prev,
              var.phat=var.phat, var.prev=var.prev, covar=covar))
}

#print.fam.recrisk <- function(x, ...) {
#}

summary.fam.recrisk <- function(x, ...) {
  
}

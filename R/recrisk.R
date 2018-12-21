recrisk.single.ascertain <- function(s, a){

  ## Description:  estimate familial recurrence risk, assuming there is a single identified proband in a family
    
  asc  <- 1*(a >= 1)
  
  phat <- sum(  asc*(a-1) ) / sum(  asc*(s-1) )
  
  f <- asc * ( (a-1)/phat - (s - a) /(1-phat) )
  
  df  <-  asc*( -(a-1)/phat^2 - ( s-a ) / (1-phat)^2 )
  
  var.phat  <- sum(f^2)/sum(df)^2
  
  return(list(phat=phat, var.phat=var.phat))
}

recrisk.ratio <- function(fit){
  ## Descrription: estimate the recurrence risk ratio of familial recurrence risk divided by population 
  ## prevalence, based on results (fit) returned from fam.recrisk
  
  rr <- fit$phat/fit$prev
  correction <- 1 + fit$covar/(fit$phat*fit$prev) - fit$var.prev/fit$prev^2
  rr.corrected <- rr*correction
  var.rr <- rr^2*(fit$var.phat/fit$phat^2 -
                    2*fit$covar/(fit$phat*fit$prev) +
                    fit$var.prev/fit$prev^2)
  var.rr.corrected <- var.rr * correction^2
  return(list(rr=rr.corrected,  var.rr=var.rr.corrected))
}

recrisk.mixture <- function(s, a, k, max.iter=1e4, eps=1e-6){
  ## Description: fit a mixure of low and high risk families based on a mixture of 
  ## truncasted binomial densities
  
  ascert <-  (a >= k)
  if(sum(ascert) < 10)
  {
    stop("fewer than 10 families ascertained")
  }
  
  s <- s[ascert]
  a <- a[ascert]

  n <- length(a)
  post <- runif(n)
  prior <- mean(post)
  
  delta <- 10^6
  phat <- fam.recrisk(s,a,k)$phat
  
  ## initial lnlike assuming no mixure
  like <- dbinom(a, s, phat) / (1 - pbinom(k - 1, s, phat))
  
  lnlike.old <- sum(log(like))
  iter <- 0
  
  while(delta > eps & iter < max.iter){
    iter <- iter + 1
    num1 <-  post * a * (a - k)
    den1 <-  post * (s * (a - k + 1) - a)
    p1 <- sum(num1) / sum(den1)
    num2 <-  (1 - post) * a * (a - k)
    den2 <-  (1 - post) * (s * (a - k + 1) - a)
    p2 <- sum(num2) / sum(den2)
    
    prob1 <- dbinom(a, s, p1) / (1 - pbinom(k - 1, s, p1))
    prob2 <- dbinom(a, s, p2) / (1 - pbinom(k - 1, s, p2))
    
    like <- (prior * prob1 + (1 - prior) * prob2)
    lnlike <- sum(log(like))
    post.new <- prior * prob1 / like
    post <- post.new
    prior <- mean(post)
    delta <- abs(lnlike.old - lnlike)
    lnlike.old <- lnlike
  }
  if(delta > eps){warning("EM failed to converge")}
  
  if(p1 > p2){
    prob.high <- p1
    prob.low <- p2
    prob.fam.high <- post
    prob.high.group <- prior
  } else {
    prob.high <- p2
    prob.low <- p1
    prob.fam.high <- 1- post
    prob.high.group <- 1-prior
  }
  
  return(list(recrisk.high=prob.high, recrisk.low=prob.low, prob.high.group = prob.high.group,
              prob.fam.high=prob.fam.high, iter=iter))
  
}

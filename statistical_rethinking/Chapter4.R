## Chapter 4
## 11 May 2021

library(rethinking)
library(dplyr)
library(tidybayes)
library(tidybayes.rethinking)
library(modelr)

## 4.1 Normal by addition
  # coin flip example 4.1 (1000 different, independent random walks converge to a normal distribution)
  pos <- replicate(1000,sum(runif(n = 16,-1,1)))
  plot(density(pos),xlim=c(-6,6)) # can change n in line above; larger numbers increasingly resemble gaussian distribution

## 4.2 Normal by multiplication
  # alleles interaction example 4.2 (effects of loci interactions multiply)
  # sample 12 random numbers between 1 and 1.1, each representing a proportional increase in growth:
  prod(1+runif(12,0,0.1)) 
  
  # generate 10,000 random products:
  growth <- replicate(10000,prod(1+runif(12,0,0.1)))
  dens(growth,norm.comp=TRUE)
  # small effects that multiply together approx. additive (comparison below):
  big <- replicate(10000,prod(1+runif(12,0,0.5)))
  small <- replicate(10000,prod(1+runif(12,0,0.01)))  
  dens(big,norm.comp=TRUE)
  dens(small,norm.comp=TRUE)

## 4.3 Normal by log-multiplication
  # Large deviates multiplied together do produce gaussian dist. on log scale (4.5)
  # (adding logs equivalent to multiplying original numbers)
  log.big <- replicate(10000,log(prod(1+runif(12,0.05))))
  dens(log.big,norm.comp = TRUE)  
  
## 4.6 Model definition to Bayes theorem (globe tossing example)
  w <- 6; n <- 9
  p_grid <- seq(from=0,to=1,length.out=100)
  posterior <- dbinom(x=w,size=n,prob = p_grid)*dunif(x=p_grid,min=0,max=1)
  posterior <- posterior/sum(posterior)
  p_grid[which.max(posterior)]
  
## 4.7 A Gaussian model of height
  data("Howell1")
  d <- Howell1  
  str(d)
  
  # filter out adults (adult heights approximately normal):
  d2 <- d[d$age>=18,]
  dens(d2$height,norm.comp = TRUE,col="blue")
  
  # plot priors to see the assumption they build into the model:
  curve(dnorm(x,178,20),from=100,to=250)
  curve(dunif(x,0,50),from=0.10,to=60)  
  
  # what do these priors imply about distribution of individual heights?
  sample_mu <- rnorm(1e4,178,20)
  sample_sigma <- runif(1e4,0,50)
  prior_h <- rnorm(1e4,sample_mu,sample_sigma)    
  dens(prior_h,norm.comp = TRUE,col="blue")  
  
## 4.14 Find posterior using grid approximation
  mu.list <- seq(from=140,to=160,length.out=200)
  sigma.list <- seq(from=4,to=9,length.out=200)  
  post <- expand.grid(mu=mu.list,sigma=sigma.list)  
  post$LL <- sapply(1:nrow(post),function(i) sum(dnorm(x=d2$height,mean=post$mu[i],sd=post$sigma[i],log=TRUE)))
  post$prod <- post$LL + dnorm(x = post$mu,mean=178,sd=20,log=TRUE) + dunif(post$sigma,min=0,max=500,log=TRUE)
  post$prob <- exp(post$prod - max(post$prod))
  
  # inspect posterior:
  contour_xyz(post$mu,post$sigma,post$prob)
  image_xyz(post$mu,post$sigma,post$prob)
  
  # sample parameter combinations (and plot, showing most plausible values of mu and sigma):
  sample.rows <- sample(1:nrow(post),size=1e4,replace=TRUE,prob=post$prob)
  sample.mu <- post$mu[sample.rows]  
  sample.sigma <- post$sigma[sample.rows]  
  plot(sample.mu,sample.sigma,cex=0.5,pch=16,col=col.alpha(rangi2,0.1))
    
  # summarize widths of densities with highest posterior density intervals:
  HPDI(sample.mu)
  HPDI(sample.sigma)

## Find maximum a posteriori estimate (MAP):
  # R code 4.25
  flist <- alist(
    height ~ dnorm(mu,sigma),
    mu ~ dnorm(178,20),
    sigma ~ dunif(0,50)
  )

  # MAP quadratic approximation of the posterior:
  m4.1 <- map(flist,data=d2)
  
  # Gaussian approximations for each parameter's marginal distribution:
  precis(m4.1)
  
  # Give map a place to start searching the prior:
  start <- list(mu = mean(d2$height), sigma = sd(d2$height))
  
  # R code 4.29: Now use a more informative prior for mu:
  m4.2 <- map(alist(
    height ~ dnorm(mu,sigma),
    mu ~ dnorm(178,0.1),
    sigma ~ dunif(0,50)
  ),data=d2)
  precis(m4.2)
  
  # inspect multi-dimensional Gaussian distribution:
  vcov(m4.1)
  diag(vcov(m4.1))
  cov2cor(vcov(m4.1))  
  
  # Now sample vectors of parameter values from multidimensional Gaussian distribution:
  post <- extract.samples(m4.1,n=1e4)
  # each value is a sample from posterior:
  head(post)
  precis(post)  
  plot(post)

  # Estimating log(sigma) instead:
  m4.1_logsigma <- map(alist(
    height ~ dnorm(mu,exp(log_sigma)),
    mu ~ dnorm(178,20),
    log_sigma ~ dnorm(2,10)
  ),data=d2)
  
  post <- extract.samples(m4.1_logsigma)
  sigma <- exp(post$log_sigma) 
  
## Linear regression: adding a weight predictor to the linear model of height
  plot(d2$height ~ d2$weight)
  m4.3 <- map(alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b*weight,
    a ~ dnorm(156,100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)),data=d2)
  precis(m4.3)
  
  # return variance-covariance matrix:
  precis(m4.3,corr=TRUE) # alpha and beta almost perfectly negatively correlated

  # A few tricks for avoiding correlations among model parameters:
  # 4.42 Centering
  d2$weight.c <- d2$weight - mean(d2$weight)

  m4.4 <- map(alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b*weight.c,
    a ~ dnorm(178,100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)),data=d2)
  precis(m4.4,corr=TRUE)
  
  # Plotting the posterior distribution:
  plot(height~weight,data=d2)
  abline(a=coef(m4.3)["a"],b=coef(m4.3)["b"])   # maximum a posteriori line mean height at each weight
  
  # To better visualize uncertainty in the regression relationship:
  post <- extract.samples(m4.3)
  post[1:5,]  
  
  N <- 200
  dN <- d2[1:N,]
  mN <- map(alist(
    height ~ dnorm(mu,sigma),
    mu <- a +b*weight,
    a ~ dnorm(178,100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,50)),data=dN)
  
  # extract 20 samples from posterior:
  post <- extract.samples(mN,n=20)
  plot(dN$weight,dN$height,xlim=range(d2$weight),ylim=range(d2$height),
       col=rangi2,xlab="weightc",ylab="heightc")
  mtext(concat("N = ",N))  
  for(i in 1:20){
    abline(a=post$a[i],b=post$b[i],col=col.alpha("black",0.3))
  }
  
  # compute regression intervals:
  mu_at_50 <- post$a + post$b*50
  dens(mu_at_50,col=rangi2,lwd=2,xlab="mu|weight=50")
  HPDI(mu_at_50,prob=0.89)

  mu <- link(m4.3)
  
  # define sequence of weights to computer predictions for:
  weight.seq <- seq(from=25,to=70,by=1)
  
  # use link to compute mu for each sample from posterior & each weight:
  mu <- link(m4.3,data=data.frame(weight=weight.seq))
  
  plot(height~weight,d2,type="n")
  for(i in 1:100){
    points(weight.seq,mu[i,],pch=16,col=col.alpha(rangi2,0.1))
  }    
  
  # summarize the distribution of mu:
  mu.mean <- apply(mu,2,mean)
  mu.HPDI <- apply(mu,2,HPDI,prob=0.89)
  
  plot(height ~ weight, data=d2,col=col.alpha(rangi2,0.5))
  lines(weight.seq,mu.mean)
  shade(mu.HPDI,weight.seq)
  
  # prediction intervals
  sim.height <- sim(m4.3,data=list(weight=weight.seq),n=1e4)
  height.PI <- apply(sim.height,2,PI,prob=0.89)
  
  plot(height~weight,d2,col=col.alpha(rangi2,0.5))
  lines(weight.seq,mu.mean)   # MAP line
  shade(mu.HPDI,weight.seq)   # HDPI region for line
  shade(height.PI,weight.seq) # prediction interval for simulated heights
  
## Polynomial regression
  # standardize weight:
  d$weight.s <- (d$weight - mean(d$weight))/sd(d$weight)
  
  m4.5 <- map(alist(
    height ~ dnorm(mu,sigma),
    mu <- a +b1*weight.s + b2*weight.s^2,
    a ~ dnorm(140,100),
    b1 ~ dnorm(0,10),
    b2 ~ dnorm(0,10),
    sigma ~ dunif(0,50)),data=d)
  precis(m4.5)
  
  weight.seq <- seq(from=-2.2,to=2,length.out=30)
  pred_dat <- list(weight.s=weight.seq,weight.s2 = weight.seq^2)
  mu <- link(m4.5,data=pred_dat)
  mu.mean <- apply(mu,2,mean)
  mu.PI <- apply(mu,2,PI,prob=0.89)
  sim.height <- sim(m4.5,data=pred_dat)
  height.PI <- apply(sim.height,2,PI,prob=0.89)
  plot(height~weight.s,d,col=col.alpha(rangi2,0.5))
  lines(weight.seq,mu.mean)
  shade(mu.PI,weight.seq)
  shade(height.PI,weight.seq)
  
  # converting back to natural scale (from z scores)
  plot(height~weight.s,d,col=col.alpha(rangi2,0.5),xaxt="n")
  at <- c(-2,-1,0,1,2)
  labels <- at*sd(d$weight)+mean(d$weight) # convert zscores back to original scale
  axis(side=1,at=at,labels=round(labels,1))
  
  
## 4.7 Practice problems
  
  # 4M1. Simulate observed heights from the prior:
  mu <- 0; sd <- 10
  sim.heights <- data.frame(mu = rnorm(1e4,mean=0,sd=10),
                            sigma = runif(1e4,min = 0,max = 10))
  sim.heights$sim.y <- rnorm(1e4,mean=sim.heights$mu,sd=sim.heights$sigma)                          
  ggplot() + geom_density(data=sim.heights,aes(x=sim.y))+labs(x="y",y="Density")
  
  # 4M2. Translate model into map formula:
  mod.4M2 <- alist(
    y ~ dnorm(mu,sigma),
    mu ~ dnorm(0,10),
    sigma ~ dunif(0,10))
  
  # 4H1. 
  weights<- c(46.95,43.72,64.78,32.59,54.63)
  weights.c <- weights-mean(weights)
  pred_dat <- list(weight.c=weights.c)
  print(m4.3)
  mu <- link(m4.4,data=pred_dat)
  mu.mean <- apply(mu,2,mean)
  mu.PI <- apply(mu,2,PI,prob=0.89)
  sim.height <- sim(m4.4,data=pred_dat)
  height.PI <- apply(sim.height,2,PI,prob=0.89)
  
  pred.missing <- data.frame(individual = seq(from=1,to=5,by=1),
                             weight = weights,
                             expected_height = mu.mean,
                             low_int = NA,
                             high_int = NA)
  for(i in 1:length(pred.missing$individual)){
    pred.missing$low_int[i] <- as.numeric(height.PI[1,i])
    pred.missing$high_int[i] <- as.numeric(height.PI[2,i])
  }
  
  # 4H2. Select Howell rows with ages below 18 years
  d3 <- d[d$age<18,]
  d3$weight.c <- d3$weight-mean(d3$weight)
  
  # a. Fit linear regression using map
  m4H2 <- map(alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b*weight,
    a ~ dnorm(50,15),
    b ~ dnorm(0,10),
    sigma ~ dunif(0,20)),data=d3)
  precis(m4H2)
  
  # b. Plot raw data, with height on vertical axis and weight on horizontal axis.
  # Superimpose MAP regression line and 89% HPDI for the mean.
  
  weight.seq <- modelr::seq_range(d3$weight,100)
  
  mu <- link(m4H2,data=data.frame(weight=weight.seq))
  mu.mean <- apply(mu,2,mean)
  mu.HPDI <- apply(mu,2,HPDI,prob=0.89)
  
  # prediction intervals
  sim.height <- sim(m4H2,data=list(weight=weight.seq),n=1e4)
  height.PI <- apply(sim.height,2,PI,prob=0.89)
  
  plot(d3$height~d3$weight,col=col.alpha(rangi2,0.5))
  lines(weight.seq,mu.mean)  # MAP regression
  shade(mu.HPDI,weight.seq) # 89% HPDI for the mean
  shade(height.PI,weight.seq) # prediction interval for predicted heights
  
  # try using tidybayes:
  mod_fits <- d3 %>%
    data_grid(weight = seq_range(weight, 100)) %>%
    add_fitted_draws(m4H2) %>%
    group_by(weight) %>%
    mean_qi(.value, .width = 0.89)
  
  mod_preds <- d3 %>%
    data_grid(weight = seq_range(weight, 100)) %>%
    add_predicted_draws(m4H2) %>%
    group_by(weight) %>%
    mean_qi(.prediction, .width = 0.89)
  
  ggplot(d3, aes(x = weight)) +
    geom_point(aes(y = height), alpha = 0.4) +
    geom_ribbon(data = mod_preds, aes(ymin = .lower, ymax = .upper),
                alpha = 0.2) +
    geom_lineribbon(data = mod_fits,
                    aes(y = .value, ymin = .lower, ymax = .upper),
                    fill = "grey60", size = 1) +
    labs(x = "Weight", y = "Height")
  
  
# 4H3. Assume log of body weight scales with height
  
  # a. Fit model using Howell data set:
  d$weight.log <- log(d$weight)
  
  m4H3 <- map(alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b*weight.log,
    a ~ dnorm(178,100),
    b ~ dnorm(0,100),
    sigma ~ dunif(0,50)),data=d)
  precis(m4H3) 
  # beta est. ~ avg. expected increase in height assoc. w/ 1-unit increase in log-weight, conditional on the data and the model
  
  weight.seq <- modelr::seq_range(d$weight.log,100)
  mu <- link(m4H3,data=data.frame(weight.log=weight.seq))
  mu.mean <- apply(mu,2,mean)
  mu.HPDI <- apply(mu,2,HPDI,prob=0.97)
  
  # prediction intervals
  sim.height <- sim(m4H3,data=list(weight.log=weight.seq),n=1e4)
  height.PI <- apply(sim.height,2,PI,prob=0.97)
  
  plot(d$height~d$weight,col=col.alpha(rangi2,0.5))
  lines(exp(weight.seq),mu.mean)   # MAP regression
  shade(mu.HPDI,exp(weight.seq))   # 97% HPDI for the mean
  shade(height.PI,exp(weight.seq)) # prediction interval for predicted heights
  
  # now try with tidybayes
  mod_fits <- d %>%
    data_grid(weight.log = seq_range(weight.log, 100)) %>%
    add_fitted_draws(m4H3) %>%
    group_by(weight.log) %>%
    mean_qi(.value, .width = 0.97)
  
  mod_preds <- d %>%
    data_grid(weight.log = seq_range(weight.log, 100)) %>%
    add_predicted_draws(m4H3) %>%
    group_by(weight.log) %>%
    mean_qi(.prediction, .width = 0.97)
  
  ggplot(d, aes(x = exp(weight.log))) +
    geom_point(aes(y = height), alpha = 0.4) +
    geom_ribbon(data = mod_preds, aes(ymin = .lower, ymax = .upper),
                alpha = 0.2) +
    geom_lineribbon(data = mod_fits,
                    aes(y = .value, ymin = .lower, ymax = .upper),
                    fill = "grey60", size = 1) +
    labs(x = "Weight", y = "Height")
  
  

  
  
  
  
  
  
  
  
  
  
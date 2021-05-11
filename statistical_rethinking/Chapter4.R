## Chapter 4
## 11 May 2021

library(rethinking)

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
  
  # inpect posterior:
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
  
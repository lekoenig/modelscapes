## Chapter 3
## 11 May 2021

## Problem 3M1:
# Suppose globe tossing data had 8 water in 15 tosses. 
# Construct the posterior distribution, using grid approximation (reference R code section 3.2 in book)
# Use same flat prior as before.

p_grid <- seq(from=0,to=1,length.out = 1000)
prior <- rep(1,1000)
likelihood <- dbinom(8,size=15,prob=p_grid)
posterior <- likelihood * prior
posterior <- posterior/sum(posterior)

globe.df <- data.frame(p=p_grid,posterior=posterior)
ggplot() + geom_line(data=globe.df,aes(x=p,y=posterior))

## Problem 3M2:
# Draw 10000 samples from the grid approximation from above. 
# Then use the samples to calculate the 90% HPDI for p.

sample.post <- sample(x=p_grid,prob=posterior,size=10000,replace=TRUE)
HPDI(sample.post,prob=0.9)

## Problem 3M3:
# Construct a posterior predictive check for this model and data. 
# Simulate the distribution of samples, averaging over the posterior uncertainty in p.
# What is the probability of observing 8 water in 15 tosses?

w <- rbinom(n=10000,size=15,prob=sample.post)
hist(w)
length(which(w==8))/length(w)
mean(w==8) # should match above line






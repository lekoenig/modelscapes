## Chapter 5: Multivariate linear models
## 4 June 2021

library(rethinking)
library(dplyr)
library(tidybayes)
library(tidybayes.rethinking)
library(modelr)
library(brms)


## 5.1 Spurious association
data(WaffleDivorce)
d <- WaffleDivorce

# standardize predictor
d$MedianAgeMarriage.s <- (d$MedianAgeMarriage - mean(d$MedianAgeMarriage))/sd(d$MedianAgeMarriage)

# fit model
m5.1 <- map(
    alist(
      Divorce ~ dnorm(mu,sigma),
      mu <- a + bA * MedianAgeMarriage.s,
      a ~ dnorm(10,10),
      bA ~ dnorm(0,1),
      sigma ~ dunif(0,10)),
    data = d
    )

# compute percentile interval of mean 
MAM.seq <- seq(from=-3,to=3.5,length.out = 30)
mu <- link(m5.1,data=data.frame(MedianAgeMarriage.s = MAM.seq))
mu.PI <- apply(mu,2,PI)

# plot
plot(Divorce ~ MedianAgeMarriage.s,data=d,col=rangi2)
abline(m5.1)
shade(mu.PI,MAM.seq)

# Fit model for marriage rate.s (standardized)
d$Marriage.s <- (d$Marriage - mean(d$Marriage))/sd(d$Marriage)
m5.2 <- map(
  alist(
    Divorce ~ dnorm(mu,sigma),
    mu <- a + bR * Marriage.s,
    a ~ dnorm(10,10),
    bR ~ dnorm(0,1),
    sigma ~ dunif(0,10)),
  data = d
)
precis(m5.2)

# plot
plot(Divorce ~ Marriage.s,data=d,col=rangi2)
abline(m5.2)

# Text section 5.1.2: Fitting the multivariate linear model
m5.3 <- map(
  alist(
    Divorce ~ dnorm(mu,sigma),
    mu <- a + bR*Marriage.s + bA*MedianAgeMarriage.s,
    a ~ dnorm(10,10),
    bR ~ dnorm(0,1),
    bA ~ dnorm(0,1),
    sigma ~ dunif(0,10)),
  data = d
)
precis(m5.3)

# visualize posterior distribution estimates:
plot(precis(m5.3))
# Note: once we know median age at marriage for a state, there is little or 
# no additional predictive pwoer in also knowing rate of marriage in that state

# PREDICTOR RESIDUAL PLOTS: avg. pred. error when we use the other pred. var. to model predictor of interest
# (leave in var. not expected by the model of the mean as a fxn of the other predictors)
m5.4 <- map(
  alist(
    Marriage.s ~ dnorm(mu, sigma),
    mu <- a + b*MedianAgeMarriage.s,
    a ~ dnorm(0, 10),
    b ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)),
  data = d)

# then compute residuals by sustracting obs. marriage rate in each state from predicted rate, based on using age at marriage:
mu <- coef(m5.4)["a"] + coef(m5.4)["b"]*d$MedianAgeMarriage.s

# compute residual for each state:
m.resid <- d$Marriage.s - mu

plot(Marriage.s ~ MedianAgeMarriage.s, d, col=rangi2)
abline(m5.4)
# loop over States
for (i in 1:length(m.resid)) {
  x <- d$MedianAgeMarriage.s[i] # x location of line segment
  y <- d$Marriage.s[i] # observed endpoint of line segment
  # draw the line segment
  lines(c(x,x), c(mu[i],y), lwd=0.5, col=col.alpha("black",0.7) )
}
# residuals = variation in marriage rate that is left over after accounting for linear assoc. with median age at marriage

# COUNTERFACTUAL PLOTS
A.avg <- mean(d$MedianAgeMarriage.s)
R.seq <- seq(from=-3,to=3,length.out=30)
pred.data <- data.frame(
      Marriage.s = R.seq,
      MedianAgeMarriage.s = A.avg
)

# compute counterfacutal mean divorce (mu)
mu <- link(m5.3,data=pred.data)
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI)

# simulate counterfactual divorce outcomes:
R.sim <- sim(m5.3,data=pred.data,n=1e4)
R.PI <- apply(R.sim,2,PI)

# display predictions, hiding raw data with type="n"
plot(Divorce ~ Marriage.s,data=d,type="n")
mtext("MedianAgeMarriage.s = 0")
lines(R.seq,mu.mean)
shade(mu.PI,R.seq)
shade(R.PI,R.seq)

# same strategy, now with Marriage.s set to is avg. and MedianAgeMarriage.s allowed to vary
R.avg <- mean(d$Marriage.s)
A.seq <- seq(from=-3,to=3.5,length.out=30)
pred.data2 <- data.frame(
    Marriage.s = R.avg,
    MedianAgeMarriage.s = A.seq
)
mu <- link(m5.3,data=pred.data2)
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI)
A.sim <- sim(m5.3,data=pred.data2,n=1e4)
A.PI <- apply(A.sim,2,PI)

plot(Divorce ~ MedianAgeMarriage.s,data=d,type="n")
mtext("Marriage.s = 0")
lines(A.seq,mu.mean)
shade(mu.PI,A.seq)
shade(A.PI,A.seq)

# POSTERIOR PREDICTION PLOTS
# call link without specifying new data (so it uses original data):
mu <- link(m5.3)

# summarize samples across cases:
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu,2,PI)

# simulate observations (no new data, uses original data)
divorce.sim <- sim(m5.3,n=1e4)
divorce.PI <- apply(divorce.sim,2,PI)

# plot predictions against observed:
plot(mu.mean ~ d$Divorce,col=rangi2,ylim=range(mu.PI),xlab="Observed divorce",ylab="Predicted divorce")
abline(a=0,b=1,lty=2)
for(i in 1:nrow(d)){
  lines(rep(d$Divorce[i],2),c(mu.PI[1,i],mu.PI[2,i]),col=rangi2)
}

# compute residuals in order to show the mean prediction error for each row:
divorce.resid <- d$Divorce - mu.mean
# get ordering by divorce rate:
o <- order(divorce.resid)
# plot
dotchart(divorce.resid[o], labels = d$Loc[o], xlim = c(-6, 5), cex = 0.6)
abline(v = 0, col = col.alpha ("black",0.2))
for (i in 1:nrow(d)) {
  j <- o[i] # which State in order
  lines(d$Divorce[j]-c(mu.PI[1,j], mu.PI[2,j]), rep(i,2))
  points(d$Divorce[j]-c(divorce.PI[1,j], divorce.PI[2,j]), rep(i,2),
         pch = 3, cex = 0.6, col = "gray")
}

# Simulate spurious association:
N <- 100                            # number of cases
x_real <- rnorm(N)                  # x_real as Gaussian with mean 0 and stddev 1
x_spur <- rnorm(N,x_real)           # x_spur as Gaussian with mean=x_real
y <- rnorm(N,x_real)                # y as Gaussian with mean = x_real
d <- data.frame(y,x_real,x_spur)    # bind all together in data frame
pairs(d)
mspur <- map(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a + b*x_real + c*x_spur,
    a ~ dnorm(0, 10),
    b ~ dnorm(0, 1),
    c ~ dnorm(0,1),
    sigma ~ dunif(0, 10)),
  data = d)
precis(mspur)

# Text section 5.2: Masked relationship
data(milk)
d <- milk
str(d)

# simple bivariate regression between kilocalories and neocortext percent:
m5.5 <- map(alist(
  kcal.per.g ~ dnorm(mu, sigma),
  mu <- a + bn*neocortex.perc,
  a ~ dnorm(0, 100),
  bn ~ dnorm(0, 1),
  sigma ~ dunif(0, 1)),
  data=d)

# look at missing values:
d$neocortex.perc
dcc <- d[complete.cases(d), ]

m5.5 <- map(alist(
  kcal.per.g ~ dnorm(mu, sigma),
  mu <- a + bn*neocortex.perc,
  a ~ dnorm(0, 100),
  bn ~ dnorm(0, 1),
  sigma ~ dunif(0, 1)),
  data=dcc)
precis(m5.5,digits=3)

# change from smallest neocortext % in the data (55%) to the largest (76%) would result in an expected change of:
coef(m5.5)["bn"] * (76 - 55)

np.seq <- 0:100
pred.data <- data.frame(neocortex.perc=np.seq)
mu <- link(m5.5, data = pred.data, n = 1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(kcal.per.g ~ neocortex.perc, data = dcc, col = rangi2)
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty=2)
lines(np.seq, mu.PI[2,], lty=2)

dcc$log.mass <- log(dcc$mass)

m5.6 <- map(alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bm*log.mass,
    a ~ dnorm(0, 100),
    bm ~ dnorm(0, 1),
    sigma ~ dunif(0, 1)),
  data=dcc)
precis(m5.6)

np.seq <- -2:4
pred.data <- data.frame(log.mass=np.seq)
mu <- link(m5.6, data = pred.data, n = 1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(kcal.per.g ~ log.mass, data = dcc, col = rangi2)
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty=2)
lines(np.seq, mu.PI[2,], lty=2)

# Fit joint model:
m5.7 <- map(alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bn*neocortex.perc + bm*log.mass,
    a ~ dnorm(0, 100),
    bn ~ dnorm(0, 1),
    bm ~ dnorm(0, 1),
    sigma ~ dunif(0, 1)),
  data=dcc)
precis(m5.7)

# counterfactual plot: neocortex percent
mean.log.mass <- mean(log(dcc$mass) )
np.seq <- 0:100
pred.data <- data.frame(
  neocortex.perc=np.seq,
  log.mass=mean.log.mass
)
mu <- link(m5.7, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(kcal.per.g ~ neocortex.perc, data=dcc, type="n")
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty=2)
lines(np.seq, mu.PI[2,], lty=2)

# counterfactual plot: magnitude of female body mass
mean.neo.mass <- mean(dcc$neocortex.perc)
np.seq <- -10:100
pred.data <- data.frame(
  log.mass=np.seq,
  neocortex.perc=mean.neo.mass
)
mu <- link(m5.7, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(kcal.per.g ~ log.mass, data=dcc, type="n")
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty=2)
lines(np.seq, mu.PI[2,], lty=2)

# Simulate masked associations:
N <- 100                                   # number of cases
rho <- 0.7                                 # correlation between x_pos and x_neg
x_pos <- rnorm(N)                          # x_pos as Gaussian
x_neg <- rnorm(N,rho*x_pos,sqrt(1-rho^2))  # x_neg correlated with x_pos
y <- rnorm(N,x_pos - x_neg)                # y equally associated with x_pos, x_neg
d <- data.frame(y,x_pos,x_neg)             # bind all together in data frame
pairs(d)

mmask <- map(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a + b*x_pos + c*x_neg,
    a ~ dnorm(0, 10),
    b ~ dnorm(0, 1),
    c ~ dnorm(0,1),
    sigma ~ dunif(0, 10)),
  data = d)
precis(mmask)

# Text section 5.3.1: Multicollinearity
N <- 100                                      # number of individuals
height <- rnorm(N,10,2)                       # sim total height of each
leg_prop <- runif(N,0.4,0.5)                  # leg as proportion of height
leg_left <- leg_prop*height + rnorm(N,0,0.02)  # sim left leg as proportion + error
leg_right <- leg_prop*height + rnorm(N,0,0.02) # sim right leg as proportion + error
d <- data.frame(height,leg_left,leg_right)    # combine into data frame

m5.8 <- map(alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2, 10),
    br ~ dnorm(2, 10),
    sigma ~ dunif(0, 10)),
  data=d)
precis(m5.8)
plot(precis(m5.8))

# look at the bivariate posterior distribution for bl and br:
post <- extract.samples(m5.8)
plot(bl~br,post,col=col.alpha(rangi2,0.1),pch=16)

# posterior distribution of sum of bl and br:
sum_blbr <- post$bl + post$br
dens(sum_blbr, col = rangi2, lwd = 2, xlab = ("sum of bl and br"))
     
# fit regression with only one of leg length variables:
m5.9 <- map(alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left,
    a ~ dnorm(10, 100),
    bl ~ dnorm(2, 10),
    sigma ~ dunif(0, 10)),
  data=d)
precis(m5.9) 
# take-home: when two predictor variables are very strongly correlated, including both in a model may lead to confusion

# look at multicollinearity with real dataset (primate milk)
data(milk)
d <- milk

# model kcal as a function of perc.fat and perc.lactose, but in two bivariate regressions:
# kcal.per.g regressed on perc.fat
m5.10 <- map(alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bf*perc.fat,
    a ~ dnorm(0.6, 10),
    bf ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)),
  data=d)

# kcal.per.g regressed on perc.lactose
m5.11 <- map(alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bl*perc.lactose,
    a ~ dnorm(0.6, 10),
    bl ~ dnorm(0,1),
    sigma ~ dunif(0, 10)),
  data=d)
precis(m5.10, digits=3)
precis(m5.11, digits=3)
# posterior means are essentially mirror images of one another

# place both predictors in the same model:
m5.12 <- map(alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bf*perc.fat + bl*perc.lactose,
    a ~ dnorm(0.6, 10),
    bf ~ dnorm(0, 1),
    bl ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)),
  data=d)
precis(m5.12, digits=3)
# perc.fat and perc.lactose contain much of the same information:
pairs(~kcal.per.g + perc.fat + perc.lactose,data=d,col=rangi2)
cor(d$perc.fat,d$perc.lactose)

# simulating collinearity:
sim.coll <- function(r=0.9) {
  d$x <- rnorm(nrow(d), mean=r*d$perc.fat,sd=sqrt((1-r^2)*var(d$perc.fat)))
  m <- lm(kcal.per.g ~ perc.fat + x, data=d)
  sqrt(diag(vcov(m)))[2] # stddev of parameter
}
rep.sim.coll <- function(r=0.9, n=100) {
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}
r.seq <- seq(from=0,to=0.99,by=0.01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r=z,n=100) )
plot(stddev ~ r.seq, type = "l", col = rangi2, lwd=2, xlab = "correlation")

# Text section 5.3.2: Post-treatment bias

N <- 100                              # number of plants
h0 <- rnorm(N, 10, 2)                 # simulate initial heights
treatment <- rep(0:1, each=N/2)       # assign treatments and simulate fungus and growth
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5 - 3*fungus)
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus) 

m5.13 <- map(alist(
    h1 ~ dnorm(mu,sigma),
    mu <- a + bh*h0 + bt*treatment + bf*fungus,
    a ~ dnorm(0,100),
    c(bh,bt,bf) ~ dnorm(0,10),
    sigma ~ dunif(0,10)),
  data=d)
precis(m5.13)
# fungus is mostly a consequence of treatment (fungus is a post-treatment variable)
# implicitly asking: once we already know whether or not a plant developed fungus, does soil treatment matter?

# to properly measure the impact of soil treatment on growth, omit post-treatment variable fungus:
m5.14 <- map(alist(
    h1 ~ dnorm(mu,sigma),
    mu <- a + bh*h0 + bt*treatment,
    a ~ dnorm(0,100),
    c(bh,bt) ~ dnorm(0,10),
    sigma ~ dunif(0,10)),
  data=d)
precis(m5.14)

# Text section 5.4, Categorical variables
data(milk)
d <- milk
unique(d$clade)

# make dummy variables:
d$clade.NWM <- ifelse(d$clade == "New World Monkey", 1, 0)
d$clade.OWM <- ifelse(d$clade == "Old World Monkey",1,0)
d$clade.S <- ifelse(d$clade == "Strepsirrhine", 1, 0)

m5.16 <- map(alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + b.NWM*clade.NWM + b.OWM*clade.OWM + b.S*clade.S,
    a ~ dnorm(0.6, 10),
    b.NWM ~ dnorm(0,1),
    b.OWM ~ dnorm(0,1),
    b.S ~ dnorm(0,1),
    sigma ~ dunif(0, 10)),
  data=d)
precis(m5.16)

# estimate a = avg. milk energy for apes. To get posterior dist. of avg. milk energy in each category, use samples:
# sample posterior
post <- extract.samples(m5.16)

# compute averages for each category
mu.ape <- post$a
mu.NWM <- post$a + post$b.NWM
mu.OWM <- post$a + post$b.OWM
mu.S <- post$a + post$b.S

# summarize using precis
precis(data.frame(mu.ape,mu.NWM,mu.OWM,mu.S) )

# another approach: unique intercepts (creating an index variable that says which param. goes with each case)
d$clade_id <- coerce_index(d$clade)

m5.16_alt <- map(alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm(0.6, 10),
    sigma ~ dunif(0, 10 )),
  data=d)
precis(m5.16_alt, depth=2) 

# lm and OLS
# rethinking package has a fxn, glimmer, that translates design formulates into map-style model formulas:
data(cars)
glimmer(dist~speed,data=cars) # glimmer automatically adds default priors

## Practice Problems
# 5M1. Invent your own example of a spurious association.
n <- 100                          
x_moisture <- rnorm(n)                  
x_rain <- rnorm(n,x_moisture)          
Resp <- rnorm(n,x_moisture)                
df_ex <- tibble(Resp,x_moisture,x_rain) %>% mutate(across(everything(),standardize))

# Fit models:
m_moist <- map(alist(
    Resp ~ dnorm(mu, sigma),
    mu <- a + b*x_moisture,
    a ~ dnorm(0, 10),
    b ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)),
  data = df_ex)
m_rain <- map(alist(
  Resp ~ dnorm(mu, sigma),
  mu <- a + b*x_rain,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  sigma ~ dunif(0, 10)),
  data = df_ex)
m_all <- map(alist(
  Resp ~ dnorm(mu, sigma),
  mu <- a + b*x_moisture + c*x_rain,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  c ~ dnorm(0,1),
  sigma ~ dunif(0, 10)),
  data = df_ex)

precis(m_moist)
precis(m_rain)
precis(m_all)


# 5M2. Invent your own example of a masked association.
n <- 100                            
rho <- 0.7
x_temp <- rnorm(n)                  
x_predation <- rnorm(n,rho*x_temp,sqrt(1-rho^2))          
GPP <- rnorm(n,x_temp-x_predation)                
df_ex <- tibble(GPP,x_temp,x_predation) %>% mutate(across(everything(),standardize))
pairs(df_ex)

# Fit models:
m_temp <- map(alist(
  GPP ~ dnorm(mu, sigma),
  mu <- a + b*x_temp,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  sigma ~ dunif(0, 10)),
  data = df_ex)
m_predation <- map(alist(
  GPP ~ dnorm(mu, sigma),
  mu <- a + b*x_predation,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  sigma ~ dunif(0, 10)),
  data = df_ex)
m_all <- map(alist(
  GPP ~ dnorm(mu, sigma),
  mu <- a + b*x_temp + c*x_predation,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  c ~ dnorm(0,1),
  sigma ~ dunif(0, 10)),
  data = df_ex)

precis(m_temp)
precis(m_predation)
precis(m_all)

# 5H1.

data(foxes)
d <- foxes
str(d)
d_fox <- d %>% mutate(across(everything(),standardize))

# body weight as a linear function of territory size
m_fox_terr <- map(alist(
  weight ~ dnorm(mu, sigma),
  mu <- a + b*area,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  sigma ~ dunif(0, 10)),
  data = d_fox)

# body weight as a linear function of group size
m_fox_group <- map(alist(
  weight ~ dnorm(mu, sigma),
  mu <- a + b*groupsize,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  sigma ~ dunif(0, 10)),
  data = d_fox)

n.seq <- seq(from=-3,to=4,by=1/16.5)
pred.data <- data.frame(
  area=n.seq,
  groupsize = n.seq,
  weight=d_fox$weight
)
# territory size:
mu <- link(m_fox_terr, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(weight ~ area, data=d_fox, type="n")
lines(n.seq, mu.mean)
lines(n.seq, mu.PI[1,], lty=2)
lines(n.seq, mu.PI[2,], lty=2)

# group size:
mu <- link(m_fox_group, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
plot(weight ~ groupsize, data=d_fox, type="n")
lines(n.seq, mu.mean)
lines(n.seq, mu.PI[1,], lty=2)
lines(n.seq, mu.PI[2,], lty=2)

precis(m_fox_terr)
precis(m_fox_group) 

# 5H2. Now fit a multiple linear regression with weight as outcome and area + group size as predictor variables
m_fox_MR <- map(alist(
  weight ~ dnorm(mu, sigma),
  mu <- a + b*groupsize + c*area,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  c ~ dnorm(0,1),
  sigma ~ dunif(0, 10)),
  data = d_fox)
precis(m_fox_MR)

# prepare counterfactual plots, holding the other predictor constant at its mean:
# prepare new counterfactual data
A.avg <- mean(d_fox$area)
F.seq <- seq(from=-3, to=3, length.out=30)
pred.data <- data.frame(
  groupsize=F.seq,
  area=A.avg
)

# compute counterfactual 
mu <- link(m_fox_MR, data=pred.data)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# simulate counterfactual outcomes
F.sim <- sim(m_fox_MR, data=pred.data, n=1e4)
F.PI <- apply(F.sim, 2, PI)

# display predictions, hiding raw data with type = “n”
plot(weight ~ groupsize, data=d_fox, type = "n")
mtext("Mean territory size = 0")
lines(F.seq, mu.mean)
shade(mu.PI, F.seq)
shade(F.PI, F.seq)

# prepare new counterfactual data
G.avg <- mean(d_fox$area)
F.seq <- seq(from=-3, to=3, length.out=30)
pred.data <- data.frame(
  groupsize=G.avg,
  area=F.seq
)

# compute counterfactual mean divorce (mu)
mu <- link(m_fox_MR, data=pred.data)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# simulate counterfactual outcomes
F.sim <- sim(m_fox_MR, data=pred.data, n=1e4)
F.PI <- apply(F.sim, 2, PI)

# display predictions, hiding raw data with type = “n”
plot(weight ~ area, data=d_fox, type = "n")
mtext("Mean group size = 0")
lines(F.seq, mu.mean)
shade(mu.PI, F.seq)
shade(F.PI, F.seq)

pairs(d_fox)
# Territory size and group size are correlated and have opposing effects on body weight (masked association)

# 5H3. 
# body weight as an additive function of avgfood and group size:
m_fox_MR1 <- map(alist(
  weight ~ dnorm(mu, sigma),
  mu <- a + b*groupsize + c*avgfood,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  c ~ dnorm(0,1),
  sigma ~ dunif(0, 10)),
  data = d_fox)

# body weight as a function of avgfood, group size, and area
m_fox_MR2 <- map(alist(
  weight ~ dnorm(mu, sigma),
  mu <- a + b*groupsize + c*avgfood + d*area,
  a ~ dnorm(0, 10),
  b ~ dnorm(0, 1),
  c ~ dnorm(0,1),
  d ~ dnorm(0,1),
  sigma ~ dunif(0, 10)),
  data = d_fox)

precis(m_fox_MR1)
precis(m_fox_MR2)

# area (sd is lower for b and c compared to avgfood model)
# when both avgfood or area area in the same model, their effects are reduced (closer to zero) and their std. errors
# area larger than when they area included in separate models because they contain very similar information (they are substitutes for one anotehr)
# posterior distribution ends up describing long ridge of combinations of b and c that are equally plausible




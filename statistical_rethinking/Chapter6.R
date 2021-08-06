## Chapter 6 (2nd edition): The haunted DAG and causal terror
## 21 July 2021

library(rethinking)
library(dplyr)
library(tidybayes)
library(tidybayes.rethinking)
library(modelr)
library(brms)

# Rcode 6.2 - multicollinearity example
N <- 100                      # number of individuals
set.seed(909)
height <- rnorm(N,10,2)       # sim total height of each person
leg_prop <- runif(N,0.4,0.5)  # leg length as proportion of height
leg_left <- leg_prop*height + rnorm(N,0,0.02)   # add error
leg_right <- leg_prop*height + rnorm(N,0,0.02)  # add error
d <- data.frame(height,leg_left,leg_right)     # combine into a data frame

m6.1 <- quap(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    br ~ dnorm(2,10),
    sigma ~ dexp(1)),
  data = d)
precis(m6.1)
plot(precis(m6.1))

# Look at join posterior prediction for bl and br:
post <- extract.samples(m6.1)
plot(bl ~ br, post, col=col.alpha(rangi2,0.1), pch=16 )

# Compute posterior distribution of sum of bl and br, and plot:
sum_blbr <- post$bl + post$br
dens(sum_blbr, col=rangi2, lwd=2, xlab="sum of bl and br")

# Now fit a regression with only one of the leg length variables
m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    sigma ~ dexp(1)), data=d)
precis(m6.2)
# Take-home: When 2 predictor vars are very strongly correlated 
# (conditional on other variables in the model), including both 
# in a model may lead to confusion.


# Rcode 6.8: Multicollinearity example (primate milk data)
data(milk)
d <- milk
d$K <- standardize(d$kcal.per.g)
d$F <- standardize(d$perc.fat)
d$L <- standardize(d$perc.lactose)

# kcal.per.g regressed on perc.fat
m6.3 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bF*F,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    sigma ~ dexp(1)), data=d)

# kcal.per.g regressed on perc.lactose
m6.4 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bL*L,
    a ~ dnorm(0,0.2),
    bL ~ dnorm(0,0.5),
    sigma ~ dexp(1)), data=d)
precis(m6.3) # Posterior distributions for bF and bL essentially mirror images of one another
precis(m6.4)

# Now place both predictors in the same regression model:
m6.5 <- quap(
  alist(
    K ~ dnorm(mu,sigma),
    mu <- a + bF*F + bL*L,
    a ~ dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    bL ~ dnorm(0,0.5),
    sigma ~ dexp(1)),
  data=d)
precis(m6.5)
plot(precis(m6.5))
# Two variables essentially form a single axis of variation:
pairs(~kcal.per.g + perc.fat + perc.lactose, data=d, col=rangi2)


# Rcode 6.12 Simulating collinearity
data(milk)
d <- milk
sim.coll <- function(r=0.9) {
  d$x <- rnorm(nrow(d), mean=r*d$perc.fat,
                sd=sqrt( (1-r^2)*var(d$perc.fat)))
  m <- lm( kcal.per.g ~ perc.fat + x, data=d)
  sqrt(diag(vcov(m)))[2] # stddev of parameter
}

rep.sim.coll <- function(r=0.9, n=100){
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}
r.seq <- seq(from=0,to=0.99,by=0.01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r=z,n=100))
plot(stddev ~ r.seq, type="l", col=rangi2, lwd=2, xlab="correlation") # uses flat priors


# R code 6.13 - an example of post-treatment bias
set.seed(71)

# number of plants
N <- 100

# simulate initial heights
h0 <- rnorm(N,10,2)

# assign treatments and simulate fungus and growth
treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5 - 3*fungus)

# compose a clean data frame
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d)

# plant growth as a proportion of h1/h0:
sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))

m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p ~ dlnorm(0, 0.25),
    sigma ~ dexp(1)), data=d)
precis(m6.6)

# Now include treatment and fungus as predictors of proportion growth, p:
m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)), data=d)
precis(m6.7)
# fungus is a post-trtmt variable, so when we control for fungus, the model is 
# answering question: once we already know whether or not a plant developed
# fungus, does soil trtmt matter? (no, bc soil trtmt has its effects on growth
# through reducing fungus)

# omit post-treatment variable:
m6.8 <- quap(
          alist(
                h1 ~ dnorm(mu, sigma),
                mu <- h0 * p,
                p <- a + bt*treatment,
                a ~ dlnorm(0, 0.2),
                bt ~ dnorm(0, 0.5),
                sigma ~ dexp(1)), data=d)
precis(m6.8)

# Rcode 6.18 - Look at DAG for plant treatment/fungus example:
library(dagitty)
plant_dag <- dagitty("dag {
  H_0 -> H_1
  F -> H_1
  T -> F
}")
coordinates(plant_dag) <- list(x = c(H_0=0,T=2,F=1.5,H_1=1),
                               y = c(H_0=0,T=0,F=0,H_1=0))
drawdag(plant_dag)
impliedConditionalIndependencies(plant_dag)


# Rcode 6.20 - fungus has no influence on growth, but M influences both H1 and F:
set.seed(71)
N <- 1000
h0 <- rnorm(N,10,2)
treatment <- rep(0:1, each=N/2)
M <- rbern(N)
fungus <- rbinom(N, size=1, prob=0.5 - treatment*0.4 + 0.4*M)
h1 <- h0 + rnorm(N, 5 + 3*M)
d2 <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)

m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    bf ~ dnorm(0, 0.5),
    sigma ~ dexp(1)), data=d2)
precis(m6.7)

m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p <- a + bt*treatment,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm(0, 0.5),
    sigma ~ dexp(1)), data=d2)
precis(m6.8)


# Rcode 6.21 - Collider bias with happiness - marriage - age:
d <- sim_happiness(seed=1977, N_years=1000)
precis(d)

# rescale data to make sense of model coefficients
d2 <- d[d$age>17,] # only adults
d2$A <- (d2$age - 18)/(65 - 18)
d2$mid <- d2$married + 1

m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)), data=d2)
precis(m6.9,depth=2)

# Compare to model that omits marriage status:
m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)), data=d2)
precis(m6.10)
# marriage status is a collider (common consequence of age and happiness,
# so when we condition on it, we create a spurious association between
# the two causes)


# Rcode 6.25 - grandparents - parents - children education DAG example
N <- 200   # number of grandparent-parent-child triads
b_GP <- 1  # direct effect of G on P
b_GC <- 0  # direct effect of G on C
b_PC <- 1  # direct effect of P on C
b_U <- 2   # direct effect of U on P and C

# use defined slopes to draw random observations:
set.seed(1)
U <- 2*rbern(N, 0.5) - 1
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U*U)
C <- rnorm(N, b_PC*P + b_GC*G + b_U*U)
d <- data.frame(C=C, P=P, G=G, U=U)

m6.11 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G,
    a ~ dnorm(0, 1),
    c(b_PC,b_GC) ~ dnorm(0, 1),
    sigma ~ dexp(1)), data=d)
precis(m6.11)

# Regression that conditions on U (neighborhood effect):
m6.12 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm(0, 1),
    c(b_PC,b_GC,b_U) ~ dnorm(0, 1),
    sigma ~ dexp(1)), data=d)
precis(m6.12)


# Rcode 6.29 - Closing the backdoor (confounds):
dag_6.1 <- dagitty("dag {
  U [unobserved]
  X -> Y
  X <- U <- A -> C -> Y
  U -> B <- C
}")
adjustmentSets( dag_6.1 , exposure="X", outcome="Y")

# wafflehouse-divorce example:
dag_6.2 <- dagitty("dag{
  A -> D
  A -> M -> D
  A <- S -> M
  S -> W -> D
}")
adjustmentSets(dag_6.2, exposure="W",outcome="D")
impliedConditionalIndependencies(dag_6.2)
# first: median age of marriage should be independent of waffle houses,
# conditioning on a state being in the south


# 6.6 Practice problems
#6M1.
library(ggdag)
dag_coords <- tibble(name = c("X", "U", "A", "B", "C", "Y", "V"),
                     x = c(1, 1, 2, 2, 3, 3, 3.5),
                     y = c(1, 2, 2.5, 1.5, 2, 1, 1.5))
dagify(Y ~ X + C + V,
       X ~ U,
       U ~ A,
       B ~ U + C,
       C ~ A + V,
       coords = dag_coords) %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(data = . %>% filter(name %in% c("U", "V")),
                 shape = 1, stroke = 2, color = "black") +
  geom_dag_text(color = "black", size = 10) +
  geom_dag_edges(edge_color = "black", edge_width = 2,
                 arrow_directed = grid::arrow(length = grid::unit(15, "pt"),
                                              type = "closed")) +
  theme_void()

new_dag <- dagitty("dag { U [unobserved]
                          V [unobserved]
                          X -> Y
                          X <- U <- A -> C -> Y
                          U -> B <- C
                          C <- V -> Y }")
adjustmentSets(new_dag, exposure = "X", outcome = "Y")

#6M2

set.seed(1984)
n <- 100
dat <- tibble(x = rnorm(n)) %>%
  mutate(z = rnorm(n, mean = x, sd = 0.4),
         y = rnorm(n, mean = z),
         across(everything(), standardize))

sim_cor <- cor(dat$x, dat$z)
sim_cor

b6m2 <- brm(y ~ 1 + x + z, data = dat, family = gaussian,
            prior = c(prior(normal(0, 0.2), class = Intercept),
                      prior(normal(0, 0.5), class = b),
                      prior(exponential(1), class = sigma)),
            iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234,
            file = here("fits", "chp6", "b6m2"))

posterior_samples(b6m2) %>%
  as_tibble() %>%
  select(-lp__) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value, y = name)) +
  stat_halfeye(.width = c(0.67, 0.89, 0.97))



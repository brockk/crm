
# Introduction
skeleton = c(0.1, 0.15, 0.25, 0.35, 0.5)
target = 0.25
num_doses = length(skeleton)
dose_indices = 1:num_doses

# The skeleton raised to a power on (0, 1) increases the values
skeleton^0.7

# The skeleton raised to a power >1 decreases the values
skeleton^1.3

# Task: What happens to the skeleton when raised to a power <0?
# Is that OK?



# dfcrm

# We fit a model to 2NNN 3NTN
library(dfcrm)
fit = crm(prior = skeleton, target = target, 
          tox = c(0,0,0, 0,1,0), 
          level = c(2,2,2, 3,3,3))
fit

fit$ptox
fit$mtd


# Task. Fit a CRM model to outcomes 2NNN 3NNN 4TTN
# What is the posterior probability of toxicity?
# What does is recommended for the next cohort?

# The one-parameter logistic CRM can be fit by:
fit3 = crm(prior = skeleton, target = target, 
           tox = c(0,0,0, 0,1,0), 
           level = c(2,2,2, 3,3,3), model = 'logistic')
fit3

# Task: Fit a logistic CRM to 2NNN 3NNN 4TTN
# Is the dose recommended the same?
# Do the posterior Prob(DLT) curves differ?


# getprior
# Two different priors
prior1 = getprior(halfwidth = 0.05, target = target, nu = 2, nlevel = num_doses)
prior2 = getprior(halfwidth = 0.1, target = target, nu = 2, nlevel = num_doses)
# And an ugly plot
plot(1:num_doses, prior2, type = 'b', col = 'blue')
points(1:num_doses, prior1, type = 'b', col = 'green')
# A nicer plot
bind_rows(data.frame(Dose = 1:num_doses, ProbDLT = prior1, HalfWidth = 0.05),
          data.frame(Dose = 1:num_doses, ProbDLT = prior2, HalfWidth = 0.1)) %>% 
  ggplot(aes(x = Dose, y = ProbDLT, group = HalfWidth, col = HalfWidth)) + 
  geom_point() + geom_line()


# Priors on the slope parameter, beta
fit5 = crm(prior = skeleton, target = target, 
           tox = c(0,0,0, 1,0,1), 
           level = c(2,2,2, 3,3,3), scale = sqrt(0.75))
fit6 = crm(prior = skeleton, target = target, 
           tox = c(0,0,0, 1,0,1), 
           level = c(2,2,2, 3,3,3), scale = sqrt(1.34))
fit7 = crm(prior = skeleton, target = target, 
           tox = c(0,0,0, 1,0,1), 
           level = c(2,2,2, 3,3,3), scale = sqrt(2))

bind_rows(
  tibble(Dose = 1:num_doses, ProbDLT = fit5$ptox, Series = 'BetaVar = 0.75'),
  tibble(Dose = 1:num_doses, ProbDLT = fit6$ptox, Series = 'BetaVar = 1.34'),
  tibble(Dose = 1:num_doses, ProbDLT = fit7$ptox, Series = 'BetaVar = 2'),
  tibble(Dose = 1:num_doses, ProbDLT = skeleton, Series = 'Skeleton')
) %>% 
  ggplot(aes(x = Dose, y = ProbDLT, group = Series, col = Series)) + 
  geom_point() + geom_line()
# The tighter prior on beta hugs the posterior closer to the skeleton.


# Numerical solution to estimate beta_hat
empiric_lik <- function(beta, skeleton, doses, tox) {
  lik <- 1
  for(i in 1:length(doses)) {
    lik = lik * (skeleton[doses[i]]^exp(beta))^tox[i] * 
      (1 - skeleton[doses[i]]^exp(beta))^(1 - tox[i])
  }
  lik
}

# This is the likelihood function.
# Notice that the maximum likelihood estimate is at roughly -0.3,
# not far from the posterior mean estimate of fit7
f = function(x) empiric_lik(x, skeleton, doses = c(2,2,2, 3,3,3), 
                            tox = c(0,0,0, 1,0,1))
curve(f, from = -2, to = 1)

# Incorporating the prior shifts the estimate slightly
g = function(x) f(x) * dnorm(x, mean = 0, sd = sqrt(2))
curve(g, from = -2, to = 1)

# Monte Carlo integration
beta_samp <- runif(1000, -2, 2)
g_samp = g(beta_samp)
data.frame(beta_samp, g_samp)
mean(beta_samp * g_samp) / mean(g_samp)
fit7$estimate
# Close. Larger n will get closer.

beta_hat_df = tibble(x = beta_samp, y = g_samp)
beta_hat_head_df = beta_hat_df %>% head(30)

# How to solve an integral numerically, graphical guide:

# Function
beta_hat_df %>% 
  ggplot(aes(x, y)) + geom_line()
# Sample points on x-axis
beta_hat_df %>% 
  ggplot(aes(x, y)) + geom_line() + 
  geom_point(data = beta_hat_head_df, aes(y = 0))
# Evaluate them at y
beta_hat_df %>% 
  ggplot(aes(x, y)) + geom_line() + 
  geom_point(data = beta_hat_head_df, aes(y = 0)) + 
  geom_point(data = beta_hat_head_df)
# Calculate heights  
beta_hat_df %>% 
  ggplot(aes(x, y)) + geom_line() + 
  geom_point(data = beta_hat_head_df, aes(y = 0)) + 
  geom_point(data = beta_hat_head_df) +
  geom_linerange(data = beta_hat_df %>% head(30), aes(ymin = 0, ymax = y))
# Take their mean
beta_hat_df %>% 
  ggplot(aes(x, y)) + geom_line() + 
  geom_point(data = beta_hat_head_df, aes(y = 0)) + 
  geom_point(data = beta_hat_head_df) +
  geom_linerange(data = beta_hat_df %>% head(30), aes(ymin = 0, ymax = y)) + 
  geom_hline(yintercept = mean(beta_hat_df$y), col = 'blue', 
             linetype = 'dashed')

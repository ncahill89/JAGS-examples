
library(tidyverse)
library(splines)
library(R2jags)
library(rjags)

GetSplines <- function( # a function to get spline basis 
  x.i, # Vector of x-values (without NAs) for which splines need to be calculated
  x0 = NULL, # x-value which determines knot placement. By default, knot is placed half-interval before last observation
  I = 2.5, # Interval length between two knots during observation period
  degree = 3
) {
  if (is.null(x0)) {
    x0 <- max(x.i)-0.1*I
  }
  # get knots, given that one knot needs to be in year0
  knots <- seq(x0-1000*I, x0+1000*I, I)
  while (min(x.i) < knots[1]) knots <- c(seq(knots[1]-1000*I, knots[1]-I,I), knots)
  while (max(x.i) > knots[length(knots)]) knots <- c(knots, seq(knots[length(knots)]+I,
                                                                knots[length(knots)]+1000*I, I))
  Btemp.ik <- bs(x.i, knots = knots[-c(1, length(knots))],  degree = degree,
                 Boundary.knots = knots[c(1, length(knots))])
  indicesofcolswithoutzeroes <- which(apply(Btemp.ik, 2, sum) > 0)
  # only remove columns with zeroes at start and end
  startnonzerocol <- indicesofcolswithoutzeroes[1]
  endnonzerocol <- indicesofcolswithoutzeroes[length(indicesofcolswithoutzeroes)]
  B.ik <- Btemp.ik[,startnonzerocol:endnonzerocol]
  colnames(B.ik) <- paste0("spline", seq(1, dim(B.ik)[2]))
  knots.k <- knots[startnonzerocol:endnonzerocol]
  names(knots.k) <- paste0("spline", seq(1, dim(B.ik)[2]))
  ##value<< List of B-splines containing:
  return(list(B.ik = B.ik, ##<< Matrix, each row is one observation, each column is one B-spline.
              knots.k = knots.k ##<< Vector of knots.
  ))
}

# read in database
Common_Era_Database_2020 <- read_csv("Common Era Database 2020.csv")

# get regions
regions <- Common_Era_Database_2020$Region %>% unique()

# select a region
sl_dat <- Common_Era_Database_2020 %>% 
            filter(Region == regions[30]) %>% 
            arrange(`Age (CE)`)
  

# plot the data (just focusing on points w/o uncertainty for this example)
ggplot(sl_dat, aes(x = `Age (CE)`, y = `RSL (m)`)) +
  geom_point()


# splines -----------------------------------------------------------------

# the arguments I and x0 should be changed depending on the data
res <- GetSplines(sl_dat$`Age (CE)`, 
                  I = 300, # Knots every 300 years across the whole time seq 
                  x0 = 0.1 # knot placed 1/10-interval before last observation
                  ) 

K <- res$knots.k %>% length
B.ik <- res$B.ik

# plot basis functions
plot(sl_dat$`Age (CE)`,res$B.ik[,1], type= "n",
     xlab = "Year",
     ylim = c(0,1), ylab ="Splines",
     xlim = range(sl_dat$`Age (CE)`))
abline(v=res$knots.k, col = seq(1, K), lwd = 1)
for (k in 1:K){
  lines(sl_dat$`Age (CE)`,res$B.ik[,k], type= "l", col = k, lwd = 1)
}

## using a different implementation for the penalisation so that we can put a prior directly on the coefficient differences
## Eilers (1999) proposed Z = BD'(DD')^(-1) (where D is the differencing matrix)
## then instead of y = B*beta, we have y = alpha + Z*delta
## where delta are the first order differences (if delta is 0 then splines stay flat)
## this works much better for convergence

D.hk <- diff(diag(K), diff = 1) # first order difference matrix (h = k-1)
Q.kh <- t(D.hk)%*%solve(D.hk%*%t(D.hk))
Z.ih <- B.ik%*%Q.kh 
H <- dim(Z.ih)[2]

# get basis for prediction

year_pred <- pretty(sl_dat$`Age (CE)`,n = 50)
n_pred <- length(year_pred)

res_pred <- GetSplines(year_pred, 
                  I = 300, # Knots every 300 years across the whole time seq 
                  x0 = 0.1 # knot placed 1/10-interval before last observation
) 

Kpred <- res_pred$knots.k %>% length
Bpred.ik <- res_pred$B.ik
Dpred.hk <- diff(diag(Kpred), diff = 1) # first order difference matrix (h = k-1)
Qpred.kh <- t(Dpred.hk)%*%solve(Dpred.hk%*%t(Dpred.hk))
Zpred.ih <- Bpred.ik%*%Qpred.kh 


# Run JAGS ----------------------------------------------------------------

# required JAGS data
jags_data <- list(y.i = sl_dat$`RSL (m)`, 
                         Z.ih = Z.ih,
                         n_obs = nrow(sl_dat),
                         H = H,
                         Zpred.ih = Zpred.ih,
                         n_pred = n_pred)

# parameters to look at
jags_pars <- c("s_pred",
               "delta.h",
               "alpha",
               "sigma_delta",
               "sigma_err")

# specify the model
spline_mod = "model
{
  # Process Model (P-Spline)
  for(i in 1:n_obs)
  {
    s[i] <- alpha + inprod(Z.ih[i,], delta.h)
  }
  
  # Data Model
  for(i in 1:n_obs)
  {
    y.i[i]~dnorm(s[i],sigma_err^(-2))
  }
  
  # Priors
  for (h in 1:H)
  {
  delta.h[h] ~ dnorm(0,sigma_delta^-2) # penalising the 1st order differences
  }
  sigma_delta ~ dt(0,1,1)T(0,)

  alpha ~ dnorm(0, 0.01) 
  sigma_err ~ dt(0,1,1)T(0,)
  
  # Prediction
  for(j in 1:n_pred)
  {
      s_pred[j] <- alpha + inprod(Zpred.ih[j,], delta.h)
  }
  
}##End model
"
  
  
# run the model
mod <- jags(data=jags_data, 
                parameters.to.save=jags_pars, 
                model.file = textConnection(spline_mod),
                n.chains = 2, 
                n.iter = 10000,
                n.burnin = 1000,
                n.thin = 5)

# (quick) check for convergence
plot(mod)

# get posterior samples and format as a tibble
post_samps <-  mod$BUGSoutput$sims.list
s_samps <- as_tibble(post_samps$s_pred) %>%
  rename_at(vars(everything()),~ as.character(seq(1:n_pred))) %>%
  pivot_longer(everything(), names_to = "year_index") %>%
  mutate(year_index = as.numeric(year_index))

# get predictions and 95% UI
res_dat <- s_samps %>% 
  group_by(year_index) %>% 
  summarise(rsl = mean(value),
            l95 = quantile(value, probs = 0.025),
            u95 = quantile(value, probs = 0.975)) %>% 
  mutate(year = year_pred)

# plot the results 
ggplot(res_dat, aes(x = year, y = rsl))+
  geom_line()+
  geom_ribbon(aes(ymin = l95, ymax = u95, fill = "RSL"),alpha = 0.5)+
  geom_point(data = sl_dat, aes(x = `Age (CE)`, y = `RSL (m)`), alpha = 0.3)+
  ylab("RSL (m)")+
  xlab("Year CE") +
  theme_bw() +
  labs(fill = "95% UI")








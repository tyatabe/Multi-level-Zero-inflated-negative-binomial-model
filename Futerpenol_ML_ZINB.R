## Setting working directory
# PC
setwd("C:/Users/tyatabe/OneDrive/Docs/Projects/Futerpenol/Code&data")
# Laptop
#setwd("C:/Users/Tadaishi/SkyDrive/Docs/Projects/Futerpenol/Code&data")
# Read in data set

d <- read.csv("data.csv")

# Create dummy variables indicating cumulative weeks of treatment
d$fut0 <- ifelse(d$futc_prev < 4, 1,0)
d$fut4 <- ifelse(d$futc_prev >= 4 & d$futc_prev<8, 1,0)
d$fut8 <- ifelse(d$futc_prev >= 8 & d$futc_prev<12, 1,0)
d$fut12 <- ifelse(d$futc_prev >= 12, 1,0)
#d$fut16 <- ifelse(d$futc_prev >= 16, 1,0)

# Create an ordinal variable for cumulative exposure
d$futord <- ifelse(d$futc_prev < 4, 0,ifelse(d$futc_prev >=4 & d$futc_prev <8, 1,
                                             ifelse(d$futc_prev >=8 & d$futc_prev <12, 2, 
                                                    ifelse(d$futc_prev >=12 & d$futc_prev <16,3,4))))
d$fut.f <- as.factor(d$futord)

# Create an ordinal variable for time since the last futerpenol treatment
d$tsince <- ifelse(d$time_since <0, "NT", ifelse(d$time_since==0, "CT", 
                                                 ifelse(d$time_since >0 & d$time_since <=2, "2W", 
                                                        ifelse(d$time_since >2 & d$time_since <=4, "4W", 
                                                               ifelse(d$time_since >4 & d$time_since <=6, "6W", "6W+")))))
d$tsince.f <- as.factor(d$tsince)
# Changing the reference level to no treatment
d<- within(d, tsince.f <- relevel(tsince.f, ref = "CT"))
d<- within(d, tsince.f <- relevel(tsince.f, ref = "NT"))

# Dummy variables for time since last treatment
d$tsinceNT <- ifelse(d$time_since < 0, 1, 0)
d$tsinceCT <- ifelse(d$time_since == 0, 1, 0)
d$tsince2W <- ifelse(d$time_since >0 & d$time_since <=2, 1, 0)
d$tsince2PW <- ifelse(d$time_since >2, 1, 0)
#d$tsince4PW <- ifelse(d$time_since >4, 1, 0)
#d$tsince7WP <- ifelse(d$time_since >6, 1, 0)

# Last week antibiotic treatment binary variable
d$ab_lw <- rep(NA, nrow(d))
# index of when cages start
index <- which(d$time<2)
# For each unique cage
for(j in 1:length(unique(d$cage))){
  # For each unique time within cage j
  for(k in 1:length(unique(d[d$cage==unique(d$cage)[j],"time"]))){
    #The value of the index - 1 (this will be time 0) + time k equals 0 if k==1 (time 0), else the value of ab_cat of the week before (k - 1), staring from tim zero (index[j]-1)
    d$ab_lw[index[j]-1+k] <- ifelse(k==1,0,d$ab_cat[index[j]-1+k-1])
  }
}

# Last week antiparasitic treatment binary variable
d$em_lw <- rep(NA, nrow(d))
# index of when cages start
index <- which(d$time<2)

for(j in 1:length(unique(d$cage))){
  for(k in 1:length(unique(d[d$cage==unique(d$cage)[j],"time"]))){
    d$em_lw[index[j]-1+k] <- ifelse(k==1,0,d$em_cat[index[j]-1+k-1])
  }
}

# Last week stocking density (fish/m3)
d$dens_lw <- rep(NA, nrow(d))
# index of when cages start
index <- which(d$time<2)

for(j in 1:length(unique(d$cage))){
  for(k in 1:length(unique(d[d$cage==unique(d$cage)[j],"time"]))){
    d$dens_lw[index[j]-1+k] <- ifelse(k==1,d$density_fish[index[j]],d$density_fish[index[j]-1+k-1])
  }
}
# Season variable
d$season <- as.factor(ifelse(d$week>=26 & d$week <39, "Winter", ifelse(d$week>=39 & d$week <52, "Spring", ifelse(d$week>=12 & d$week <26, "Fall", "Summer"))))
# Dummy variable for season
d$Fall <- ifelse(d$season=="Fall", 1, 0)
d$Spring <- ifelse(d$season=="Spring", 1, 0)
d$Summer <- ifelse(d$season=="Summer", 1, 0)

# Creating an offset for poisson regression
d$lpop <- log(d$pop_mid)
# Creating obs id within cage
d$obs <- d$time

# Creatig a pairwise difference in time matrix (time distance b/w obs)
dist.u <- dist(unique(d$time), method = "maximum", diag = T, upper = T, p = 2)
dist.m.u <- as.matrix(dist.u)

## Extracting water temperature data from NASA Earth Observations (NEO)
library(raster)
library(rgdal)
library(ncdf4)

# Reading weekly temp data from NOAA (https://www.esrl.noaa.gov/psd/repository/entry/show?entryid=12159560-ab82-48a1-b3e4-88ace20475cd)
temp <- stack("sst.wkmean.1990-present.nc", varname="sst")# Temp data
world <- raster("lsmask.nc")# world mask
# Rotate to change log from 0 to 360 to -180 to 180
temp <- rotate(temp)
world <- rotate(world)
# Subsetting the stack for the required weeks
temp15_16 <- subset(temp, 1326:1374)
# Downloading Chile shp file
cl <- getData('GADM', country='CL', level=1)
# Reproject CL to the raster projection
cl <- spTransform(cl, crs(temp))
# Subsetting the regions of interest
cl <- cl[cl$NAME_1 %in% c('Los Lagos','Aisén del General Carlos Ibáñez del Campo'), ]
# Tipying in farm coordinates
site_coord <- data.frame(c("A1", "A2", "A3"), c(-72.4655, -72.6385, -74.116), c(-42.10, -42.055, -44.803))
colnames(site_coord) <- c("site", "x", "y")
# Creating spatialPoints object
site.pts <- SpatialPoints(data.frame(site_coord$x, site_coord$y))
site.pts <- SpatialPointsDataFrame(site.pts, site_coord)
projection(site.pts) <- "+proj=longlat +datum=WGS84"
site.pts <- spTransform(site.pts, crs(temp))

# Plotting to see if it makes sense
plot(cl, col="grey")
plot(temp15_16[[1]], add=T)
plot(cl, add=T, col="grey")
plot(site.pts, col="red", pch=20, add=T)

# Extract temp values for farms
sitetemp <- extract(temp15_16, site.pts)
# Transpose the data frame
sitetemp <- t(sitetemp)
# Change from matrix to dataframe
sitetemp <- as.data.frame(sitetemp)
# Renaming to match mort data set
colnames(sitetemp) <-as.factor(c("A1", "A2", "A3"))
rownames(sitetemp) <- NULL
sitetemp$time <- seq(1:49)
# Making long data frame to use merge
temp_long <- data.frame(c(rep("A1", 49), rep("A2", 49), rep("A3", 49)), 
                        c(sitetemp$A1, sitetemp$A2, sitetemp$A3), rep(seq(1:49), 3))
colnames(temp_long)<- c("site", "temp", "time")
# Adding temp data to mort data set
d <- merge(d, temp_long, by=c("site", "time"))

# Preparing data for Stan
data=list(logpop = d$lpop, mort=d$srs_mort, futc=d$fut.f, tsince=d$tsince.f, 
          ablw=d$ab_lw, emlw=d$em_lw, prevmort=d$prev.mort, dens=d$dens_lw, 
          site = as.numeric(d$site), cage=as.numeric(d$cage), time=d$time, 
          N_time=length(unique(d$time)), Dmat=dist.m.u, N = nrow(d), 
          N_site = length(unique(d$site)), N_cage = length(unique(d$cage)), 
          fut=d$fut_cum, fut2=(d$fut_cum)^2, fut3=(d$fut_cum)^3, timesince=d$time_since, 
          timesince2=(d$time_since)^2, fut4 = d$fut4, fut8 = d$fut8, fut12 = d$fut12, 
          CT = d$tsinceCT, W2 = d$tsince2W, W2P = d$tsince2PW, weight = d$weight_mid,
          fall = d$Fall, spring = d$Spring, summer=d$Summer, temp=d$temp)

# Centered (binary) and scaled (numerical) data
data=list(logpop = d$lpop, mort=d$srs_mort, futc=d$fut.f, tsince=d$tsince.f, 
          ablw=d$ab_lw - mean(d$ab_lw) , emlw=d$em_lw - mean(d$em_lw), 
          prevmort= as.vector(scale(log(d$prev.mort + 0.001), scale=sd(log(d$prev.mort + 0.001))*2)),
          dens=as.vector(scale(d$dens_lw, scale = sd(d$dens_lw)*2)), 
          site = as.numeric(d$site), cage=as.numeric(d$cage), time=d$time, 
          N_time=length(unique(d$time)), Dmat=dist.m.u, N = nrow(d), 
          N_site = length(unique(d$site)), N_cage = length(unique(d$cage)), 
          fut=as.vector(scale(d$fut_cum, scale=sd(d$fut_cum)*2)),
          fut2=as.vector(scale(d$fut_cum, scale=sd(d$fut_cum)*2)^2), 
          fut3=as.vector(scale(d$fut_cum, scale=sd(d$fut_cum)*2)^3), 
          timesince=as.vector(scale(d$time_since, scale = sd(d$time_since)*2)), 
          timesince2=as.vector(scale(d$time_since, scale = sd(d$time_since)*2)^2), 
          fut4 = d$fut4 - mean(d$fut4), fut8 = d$fut8 - mean(d$fut8),
          fut12 = d$fut12 - mean(d$fut12), CT = d$tsinceCT - mean(d$tsinceCT),
          W2 = d$tsince2W - mean(d$tsince2W), W2P = d$tsince2PW - mean(d$tsince2PW),
          weight = as.vector(scale(d$weight_mid, scale = sd(d$weight_mid)*2)),
          temp=as.vector(scale(d$temp, scale = sd(d$temp)*2)))

# Fitting a zero-inflated negative binomial model with main exposures (cumulative weeks of
# treatment and cum weeks since last treatment) as categorical predictors
stancode= "data {
int<lower=1> N;
int<lower=1> N_site;
int<lower=1> N_cage;
int<lower=1> N_time;
int<lower=0> mort[N];
real logpop[N];
real temp[N];
real ablw[N];
real emlw[N];
real prevmort[N];
int site[N];
int cage[N];
int time[N];
real fut4[N];
real fut8[N];
real fut12[N];
real CT[N];
real W2[N];
real W2P[N];
real weight[N];
}
parameters {
real a;
vector[N_site] a_site_raw;
vector[N_cage] a_cage_raw;
vector[N_time] a_week_raw;
real b;
vector[N_site] b_site_raw;
vector[N_cage] b_cage_raw;
vector[N_time] b_week_raw;
//vector[N_site] a_site;
//vector[N_cage] a_cage;
//vector[N_time] a_week;
real bfut4;
real bfut8;
real bfut12;
real pfut4;
real pfut8;
real pfut12;
real bct;
real b2w;
real b2Pw;
real bablw;
real bemlw;
real bprev;
real bweight;
real btemp;
real bfut4_temp;
real bfut8_temp;
real bfut12_temp;
real pct;
real p2w;
real p2Pw;
real pablw;
real pemlw;
real pprev;
real pweight;
real ptemp;
real pfut4_temp;
real pfut8_temp;
real pfut12_temp;
real<lower=0> sigma_site;
real<lower=0> sigma_cage;
real<lower=0> sigma_week;
real<lower=0> psigma_site;
real<lower=0> psigma_cage;
real<lower=0> psigma_week;
real<lower=0> theta;
//real<lower=0, upper=1> p;
//real<lower=0> lambda;
}
transformed parameters{
vector[N_site] a_site;
vector[N_cage] a_cage;
vector[N_time] a_week;
vector[N_site] b_site;
vector[N_cage] b_cage;
vector[N_time] b_week;
a_site = sigma_site*a_site_raw;
a_cage = sigma_cage*a_cage_raw;
a_week = sigma_week*a_week_raw;
b_site = psigma_site*b_site_raw;
b_cage = psigma_cage*b_cage_raw;
b_week = psigma_week*b_week_raw;

}
model {
vector[N] p;
vector[N] lambda;
a ~ normal( 0 , 1);
b ~ normal( 0 , 1);
//lambda ~ cauchy (0,2);
theta ~ exponential( 2 );
//p ~ cauchy(0, 5);
sigma_site ~ cauchy( 0 , 2);
sigma_cage ~ cauchy( 0 , 2);
sigma_week ~ cauchy( 0 , 2);
a_site_raw ~ normal( 0 , 1);
a_cage_raw ~ normal( 0 , 1);
a_week_raw ~ normal( 0 , 1);
psigma_site ~ cauchy( 0 , 2);
psigma_cage ~ cauchy( 0 , 2);
psigma_week ~ cauchy( 0 , 2);
b_site_raw ~ normal( 0 , 1);
b_cage_raw ~ normal( 0 , 1);
b_week_raw ~ normal( 0 , 1);
bfut4 ~ normal(0,0.5);
bfut8 ~ normal(0,0.5);
bfut12 ~ normal(0,0.5);
pfut4 ~ normal(0,0.5);
pfut8 ~ normal(0,0.5);
pfut12 ~ normal(0,0.5);
bct ~ normal(0,0.5);
b2w ~ normal(0,0.5);
b2Pw ~ normal(0,0.5);
bablw ~ normal(0,0.5);
bemlw ~ normal(0,0.5);
bprev ~ normal(0,0.5);
bweight ~ normal(0,0.5);
btemp ~ normal(0,0.5);
bfut4_temp ~ normal(0,0.5);
bfut8_temp ~ normal(0,0.5);
bfut12_temp ~ normal(0,0.5);
pct ~ normal(0,0.5);
p2w ~ normal(0,0.5);
p2Pw ~ normal(0,0.5);
pablw ~ normal(0,0.5);
pemlw ~ normal(0,0.5);
pprev ~ normal(0,0.5);
pweight ~ normal(0,0.5);
ptemp ~ normal(0,0.5);
pfut4_temp ~ normal(0,0.5);
pfut8_temp ~ normal(0,0.5);
pfut12_temp ~ normal(0,0.5);


for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]] + bfut4*fut4[i] + bfut8*fut8[i] + bfut12*fut12[i] + bct*CT[i] + b2w*W2[i] + b2Pw*W2P[i] + bablw*ablw[i] + bemlw*emlw[i] + bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i] + bfut4_temp*temp[i]*fut4[i] + bfut8_temp*temp[i]*fut8[i] + bfut12_temp*temp[i]*fut12[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_site[site[i]] + b_cage[cage[i]] + b_week[time[i]] + pfut4*fut4[i] + pfut8*fut8[i] + pfut12*fut12[i] + pct*CT[i] + p2w*W2[i] + p2Pw*W2P[i] + pablw*ablw[i] + pemlw*emlw[i] + pprev*prevmort[i] + pweight*weight[i] + ptemp*temp[i] + pfut4_temp*temp[i]*fut4[i] + pfut8_temp*temp[i]*fut8[i] + pfut12_temp*temp[i]*fut12[i]);
}
for (n in 1:N) {
if (mort[n] == 0){
target += log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
target += bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
}
}

generated quantities {
real lambda[N];
real p[N];
real pred[N];
real resid[N];
real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + 
a_week[time[i]] + bfut4*fut4[i] + bfut8*fut8[i] + bfut12*fut12[i] + 
bct*CT[i] + b2w*W2[i] + b2Pw*W2P[i] + bablw*ablw[i] + bemlw*emlw[i] + 
bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i] + 
bfut4_temp*temp[i]*fut4[i] + bfut8_temp*temp[i]*fut8[i] + 
bfut12_temp*temp[i]*fut12[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_site[site[i]] + b_cage[cage[i]] + 
b_week[time[i]] + pfut4*fut4[i] + pfut8*fut8[i] + pfut12*fut12[i] + 
pct*CT[i] + p2w*W2[i] + p2Pw*W2P[i] + pablw*ablw[i] + pemlw*emlw[i] + 
pprev*prevmort[i] + pweight*weight[i] + ptemp*temp[i] + 
pfut4_temp*temp[i]*fut4[i] + pfut8_temp*temp[i]*fut8[i] + 
pfut12_temp*temp[i]*fut12[i]);
}

for (n in 1:N) {
if (mort[n] == 0){
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
log_lik[n] = bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
pred[n] = (1-p[n]) * lambda[n];
resid[n] = (mort[n] - pred[n])/sqrt(pred[n] + theta * pred[n]^2);
}
}
"

# Setting up inits
n_chains <- 4
start <- list(a=0, b=0, a_site = rep(0, data$N_site), b_site = rep(0, data$N_site), a_cage = rep(0, data$N_cage), b_cage = rep(0, data$N_cage), a_week = rep(0, data$N_time), b_week = rep(0, data$N_time), bfut4=0, bfut8=0, bfut12=0, bct=0, b2w=0, b2Pw=0, bablw=0, bemlw=0, bprev=0, bweight=0, btemp=0, bfut4_temp=0, bfut8_temp=0, bfut12_temp=0, pfut4=0, pfut8=0, pfut12=0, pct=0, p2w=0, p2Pw=0, pablw=0, pemlw=0, pprev=0, pweight=0, ptemp=0, pfut4_temp=0, pfut8_temp=0, pfut12_temp=0)
init <- list()
for ( i in 1:n_chains ) init[[i]] <- start

#Running model
library(rstan)
library(Rcpp)
m1.4 <- stan(model_code = stancode, iter=4000, chains=4, cores=4, warmup=2000, control = list(adapt_delta = 0.95, stepsize=1), data=data, init=init, algorithm = "NUTS")

# Printing estimates of interest
#Population level effects
print(m1.4, probs=c(.025, .975), 
      pars=c("a", "bfut4", "bfut8", "bfut12", "bct", "b2w", "b2Pw", "bablw", "bemlw", 
             "bprev", "bweight", "btemp", "bfut4_temp", "bfut8_temp", "bfut12_temp", "b", 
             "pfut4", "pfut8", "pfut12", "pct", "p2w", "p2Pw", "pablw", "pemlw", "pprev", 
             "pweight", "ptemp", "pfut4_temp", "pfut8_temp", "pfut12_temp", "theta",
             "sigma_site", "sigma_cage", "sigma_week", "psigma_site", "psigma_cage", 
             "psigma_week"))
#Varying effects rate model: site, cage, week
print(m1.4, probs=c(.025, .975), pars=c("a_site", "a_cage", "a_week"))
#Varying effects logit model: site, cage, week
print(m1.4, probs=c(.025, .975), pars=c("b_site", "b_cage", "b_week"))


# Plotting population level effects
plot(m1.4, 
     pars=c("a", "bfut4", "bfut8", "bfut12", "bct", "b2w", "b2Pw", "bablw", "bemlw", 
            "bprev", "bweight", "btemp", "bfut4_temp", "bfut8_temp", "bfut12_temp", "b", 
            "pfut4", "pfut8", "pfut12", "pct", "p2w", "p2Pw", "pablw", "pemlw", "pprev", 
            "pweight", "ptemp", "pfut4_temp", "pfut8_temp", "pfut12_temp", "theta"))

# Traceplots
# Population level effects
traceplot(m1.4, pars=c("a", "bfut4", "bfut8", "bfut12", "bct", "b2w", "b2Pw", "bablw", 
                       "bemlw", "bprev", "bweight", "btemp", "bfut4_temp", "bfut8_temp", 
                       "bfut12_temp", "b", "pfut4", "pfut8", "pfut12", "pct", "p2w", 
                       "p2Pw", "pablw", "pemlw", "pprev", "pweight", "ptemp", "pfut4_temp",
                       "pfut8_temp", "pfut12_temp", "theta", "sigma_site", "sigma_cage", 
                       "sigma_week", "psigma_site", "psigma_cage", "psigma_week"))
# Varying effects rate model
traceplot(m1.4, pars=c("a_site", "a_cage"))
traceplot(m1.4, pars=c("a_week"))
# Varying effects logit model
traceplot(m1.4, pars=c("b_site", "b_cage"))
traceplot(m1.4, pars=c("b_week"))

# Checking looic
library(loo)
log_lik <- extract_log_lik(m1.4)
(Cat_loo <- loo(log_lik))




# Fitting a zero-inflated negative binomial model with main exposures (cumulative weeks of
# treatment and cum weeks since last treatment) as numerical predictors
stancode2= "data {
int<lower=1> N;
int<lower=1> N_site;
int<lower=1> N_cage;
int<lower=1> N_time;
int<lower=0> mort[N];
real logpop[N];
real temp[N];
real ablw[N];
real emlw[N];
real prevmort[N];
int site[N];
int cage[N];
int time[N];
real fut[N];
real fut2[N];
real timesince[N];
real timesince2[N];
real weight[N];
}
parameters {
real a;
vector[N_site] a_site_raw;
vector[N_cage] a_cage_raw;
vector[N_time] a_week_raw;
real b;
vector[N_site] b_site_raw;
vector[N_cage] b_cage_raw;
vector[N_time] b_week_raw;
//vector[N_site] a_site;
//vector[N_cage] a_cage;
//vector[N_time] a_week;
real bfut;
real bfut2;
real pfut;
real pfut2;
real btsince;
real btsince2;
real bablw;
real bemlw;
real bprev;
real bweight;
real btemp;
real bfut_temp;
real bfut2_temp;
real ptsince;
real ptsince2;
real pablw;
real pemlw;
real pprev;
real pweight;
real ptemp;
real pfut_temp;
real pfut2_temp;
real<lower=0> sigma_site;
real<lower=0> sigma_cage;
real<lower=0> sigma_week;
real<lower=0> psigma_site;
real<lower=0> psigma_cage;
real<lower=0> psigma_week;
real<lower=0> theta;
//real<lower=0, upper=1> p;
//real<lower=0> lambda;
}
transformed parameters{
vector[N_site] a_site;
vector[N_cage] a_cage;
vector[N_time] a_week;
vector[N_site] b_site;
vector[N_cage] b_cage;
vector[N_time] b_week;
a_site = sigma_site*a_site_raw;
a_cage = sigma_cage*a_cage_raw;
a_week = sigma_week*a_week_raw;
b_site = psigma_site*b_site_raw;
b_cage = psigma_cage*b_cage_raw;
b_week = psigma_week*b_week_raw;

}
model {
vector[N] p;
vector[N] lambda;
a ~ normal( 0 , 1);
b ~ normal( 0 , 1);
//lambda ~ cauchy (0,2);
theta ~ exponential( 2 );
//p ~ cauchy(0, 5);
sigma_site ~ cauchy( 0 , 2);
sigma_cage ~ cauchy( 0 , 2);
sigma_week ~ cauchy( 0 , 2);
a_site_raw ~ normal( 0 , 1);
a_cage_raw ~ normal( 0 , 1);
a_week_raw ~ normal( 0 , 1);
psigma_site ~ cauchy( 0 , 2);
psigma_cage ~ cauchy( 0 , 2);
psigma_week ~ cauchy( 0 , 2);
b_site_raw ~ normal( 0 , 1);
b_cage_raw ~ normal( 0 , 1);
b_week_raw ~ normal( 0 , 1);
bfut ~ normal(0,0.5);
bfut2 ~ normal(0,0.5);
pfut ~ normal(0,0.5);
pfut2 ~ normal(0,0.5);
btsince ~ normal(0,0.5);
btsince2 ~ normal(0,0.5);
bablw ~ normal(0,0.5);
bemlw ~ normal(0,0.5);
bprev ~ normal(0,0.5);
bweight ~ normal(0,0.5);
btemp ~ normal(0,0.5);
bfut_temp ~ normal(0,0.5);
bfut2_temp ~ normal(0,0.5);
ptsince ~ normal(0,0.5);
ptsince2 ~ normal(0,0.5);
pablw ~ normal(0,0.5);
pemlw ~ normal(0,0.5);
pprev ~ normal(0,0.5);
pweight ~ normal(0,0.5);
ptemp ~ normal(0,0.5);
pfut_temp ~ normal(0,0.5);
pfut2_temp ~ normal(0,0.5);


for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ bfut*fut[i] + bfut2*fut2[i] + btsince*timesince[i] + btsince2*timesince2[i] 
+ bablw*ablw[i] + bemlw*emlw[i] + bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i] 
+ bfut_temp*temp[i]*fut[i] + bfut2_temp*temp[i]*fut2[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_site[site[i]] + b_cage[cage[i]] + b_week[time[i]] + pfut*fut[i] 
+ pfut2*fut2[i] + ptsince*timesince[i] + ptsince2*timesince2[i] + pablw*ablw[i] 
+ pemlw*emlw[i] + pprev*prevmort[i] + pweight*weight[i] + ptemp*temp[i] 
+ pfut_temp*temp[i]*fut[i] + pfut2_temp*temp[i]*fut2[i]);
}
for (n in 1:N) {
if (mort[n] == 0){
target += log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
target += bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
}
}

generated quantities {
  real lambda[N];
  real p[N];
  real pred[N];
  real resid[N];
  real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ bfut*fut[i] + bfut2*fut2[i] + btsince*timesince[i] + btsince2*timesince2[i] 
+ bablw*ablw[i] + bemlw*emlw[i] + bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i] 
+ bfut_temp*temp[i]*fut[i] + bfut2_temp*temp[i]*fut2[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_site[site[i]] + b_cage[cage[i]] + b_week[time[i]] + pfut*fut[i] 
+ pfut2*fut2[i] + ptsince*timesince[i] + ptsince2*timesince2[i] + pablw*ablw[i] 
+ pemlw*emlw[i] + pprev*prevmort[i] + pweight*weight[i] + ptemp*temp[i] 
+ pfut_temp*temp[i]*fut[i] + pfut2_temp*temp[i]*fut2[i]);
}

  for (n in 1:N) {
    if (mort[n] == 0){
      log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p[n]),
      bernoulli_lpmf(0 | p[n])
      + neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
      }else{
      log_lik[n] = bernoulli_lpmf(0 | p[n])
      + neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
      }
   pred[n] = (1-p[n]) * lambda[n];
   resid[n] = (mort[n] - pred[n])/sqrt(pred[n] + theta * pred[n]^2);
  }
}
"

# Setting up inits
n_chains <- 4
start <- list(a=0, b=0, a_site = rep(0, data$N_site), b_site = rep(0, data$N_site), 
              a_cage = rep(0, data$N_cage), b_cage = rep(0, data$N_cage), 
              a_week = rep(0, data$N_time), b_week = rep(0, data$N_time), bfut=0, bfut2=0,
              btsince=0, btsince2=0, bablw=0, bemlw=0, bprev=0, bweight=0, btemp=0, 
              bfut_temp=0, bfut2_temp=0, pfut=0, pfut2=0, ptsince=0, ptsince2=0, 
              pablw=0, pemlw=0, pprev=0, pweight=0, ptemp=0, pfut_temp=0, pfut2_temp=0)
init <- list()
for ( i in 1:n_chains ) init[[i]] <- start

#Running model
library(rstan)
library(Rcpp)
m1.5 <- stan(model_code = stancode2, iter=4000, chains=4, cores=4, warmup=2000, control = list(adapt_delta = 0.95, stepsize=1), data=data, init=init, algorithm = "NUTS")

# Printing estimates of interest
#Population level effects
print(m1.5, probs=c(.025, .975), 
      pars=c("a", "bfut", "bfut2", "btsince", "btsince2", "bablw", "bemlw", 
             "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp", "b", 
             "pfut", "pfut2", "ptsince", "ptsince2", "pablw", "pemlw", "pprev", 
             "pweight", "ptemp", "pfut_temp", "pfut2_temp", "theta", "sigma_site", 
             "sigma_cage", "sigma_week", "psigma_site", "psigma_cage", "psigma_week"))
#Varying effects rate model: site, cage, week
print(m1.5, probs=c(.025, .975), pars=c("a_site", "a_cage", "a_week"))
#Varying effects logit model: site, cage, week
print(m1.5, probs=c(.025, .975), pars=c("b_site", "b_cage", "b_week"))


# Plotting population level effects
plot(m1.5, 
     pars=c("a", "bfut", "bfut2", "btsince", "btsince2", "bablw", "bemlw", 
            "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp", "b", 
            "pfut", "pfut2", "ptsince", "ptsince2", "pablw", "pemlw", "pprev", 
            "pweight", "ptemp", "pfut_temp", "pfut2_temp", "theta"))

# Traceplots
# Population level effects
traceplot(m1.5, pars=c("a", "bfut", "bfut2", "btsince", "btsince2", "bablw", "bemlw", 
                       "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp", "b", 
                       "pfut", "pfut2", "ptsince", "ptsince2", "pablw", "pemlw", "pprev", 
                       "pweight", "ptemp", "pfut_temp", "pfut2_temp", "theta", "sigma_site", 
                       "sigma_cage", "sigma_week", "psigma_site", "psigma_cage", "psigma_week"))
# Varying effects rate model
traceplot(m1.5, pars=c("a_site", "a_cage"))
traceplot(m1.5, pars=c("a_week"))
# Varying effects logit model
traceplot(m1.5, pars=c("b_site", "b_cage"))
traceplot(m1.5, pars=c("b_week"))

# Checking looic
library(loo)
log_lik1 <- extract_log_lik(m1.5)
(Cont_loo <- loo(log_lik1))

# Obs Vs predicted...need to plot on the log-scale
post.cont <- extract(m1.5)
pred <- post.cont$pred
pred.med <- apply(pred, 2, median)

# ZNIB parameter estimates
p <- post.cont$p
lambda <- post.cont$lambda
theta <- post.cont$theta


# Save image of workspace
save.image()
# Save list of objects
save(m1.4, m1.5, d,file="models.RData")

# Simulating observations from a ZINB
library(ZIM)
kk <- matrix(rep(NA,8000*2480), ncol=2480)

for (i in 1:2480){   
  kk[,i] <- rzinb(1, theta, lambda[,i], omega=p[,i])
}

# Mean and variance observed Vs predicted
# Observed
mean(d$srs_mort)
median(d$srs_mort)
var(d$srs_mort)
# Predicted
mean(as.vector(kk))
median(as.vector(kk))
var(as.vector(kk))
# Proportion of zeroes
length(d$srs_mort[d$srs_mort<1])/length(d$srs_mort)
length(kk[kk<1])/length(kk)

# Comparing histograms
par(mfrow=c(1,2))
simplehist(d$srs_mort)
simplehist(kk)


# q-q plot
qqplot(d$srs_mort, as.vector(kk), ylim=c(0,600), xlim=c(0,600))
abline(a=0, b=1)
qqplot(d$srs_mort[d$srs_mort<35], as.vector(kk[kk<35]), ylim=c(0,35), xlim=c(0,35))
abline(a=0, b=1)

# Permutation test to compare distributions? They look different

# Experiment with model implications
a <- as.vector(post.cont$a)
a_site1 <- post.cont$a_site[,1]
a_site2 <- post.cont$a_site[,2]
a_site3 <- post.cont$a_site[,3]
a_cage <- apply(post.cont$a_cage, 1, FUN = "mean")
a_week <- apply(post.cont$a_week, 1, FUN = "mean")
b <- as.vector(post.cont$b)
b_site1 <- post.cont$b_site[,1]
b_site2 <- post.cont$b_site[,2]
b_site3 <- post.cont$b_site[,3]
b_cage <- apply(post.cont$b_cage, 1, FUN = "min")
b_week <- apply(post.cont$b_week, 1, FUN = "min")
bfut <- as.vector(post.cont$bfut)
bfut2 <- post.cont$bfut2
pfut <- post.cont$pfut
pfut2 <- post.cont$pfut2
btsince <- post.cont$btsince
btsince2 <- post.cont$btsince2
bablw <- post.cont$bablw
bemlw <- post.cont$bablw
bprev <- post.cont$bablw
bweight <- post.cont$bweight
btemp <- post.cont$btemp
bfut_temp <- post.cont$bfut_temp
bfut2_temp <- post.cont$bfut2_temp
ptsince <- post.cont$ptsince
ptsince2 <- post.cont$ptsince2
pablw <- post.cont$pablw
pemlw <- post.cont$pemlw
pprev <- post.cont$pprev
pweight <- post.cont$pweight
ptemp <- post.cont$ptemp
pfut_temp <- post.cont$pfut_temp
pfut2_temp <- post.cont$pfut2_temp

# Effect of cumulative exposure to fut for differnt levels of temp
# Values of predictors to be evaluated
fut.seq <- seq(from=min(data$fut), to=max(data$fut), length.out=40)
temp.low <- min(data$temp)
temp.mid <- median(data$temp)
temp.max <- max(data$temp)
tsince.med <- median(data$timesince)
ablw <- mean(data$ablw)
emlw <- mean(data$emlw)
prevmort <- median(data$prevmort)
weight <- median(data$weight)
logpop <- median(data$logpop)


# Estimating mortality rate predicted by the model for a range of cum exposure and temps
library(boot)
# Low temp
lambda.low.temp <- matrix(rep(NA, 40*8000),ncol=40)
p.low.temp <- matrix(rep(NA, 40*8000),ncol=40)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.low.temp[,i] <- 
  logpop + a + a_site1 + a_cage + a_week + bfut*fut.seq[i] + bfut2*fut.seq[i]^2 + btsince*tsince.med 
    + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort + bweight*weight 
    + btemp*temp.low + bfut_temp*temp.low*fut.seq[i] + bfut2_temp*temp.low*fut.seq[i]^2
}  

# (logit)Prob of zero inflation
for (i in 1:length(fut.seq)){  
  p.low.temp[,i] <- 
  b + b_site1 + b_cage + b_week + pfut*fut.seq[i] + pfut2*fut.seq[i]^2 + ptsince*tsince.med 
    + ptsince2*tsince.med^2 + pablw*ablw + pemlw*emlw + pprev*prevmort + pweight*weight 
    + ptemp*temp.low + pfut_temp*temp.low*fut.seq[i] + pfut2_temp*temp.low*fut.seq[i]^2
}

# Transforming to rate and probability
lambda.low.temp <- exp(lambda.low.temp)
p.low.temp <- inv.logit(p.low.temp)

# Estimating the expectation
rate.low.temp <- lambda.low.temp*(1 - p.low.temp)


# Median, 66%, and 95% PI
median.rate.low.temp <- apply( rate.low.temp, 2 , median )
PI66.rate.low.temp <- apply( rate.low.temp, 2 , PI, prob=0.66 )
PI95.rate.low.temp <- apply( rate.low.temp, 2 , PI, prob=0.95 )


# Ploting

plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,200),
     main="Low temp", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
    length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.low.temp)
shade( PI95.rate.low.temp , fut.seq , col=col.alpha("pink",0.2) )
shade( PI66.rate.low.temp , fut.seq , col=col.alpha(rangi2,0.2) )


# Medium temp
lambda.med.temp <- matrix(rep(NA, 40*8000),ncol=40)
p.med.temp <- matrix(rep(NA, 40*8000),ncol=40)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.med.temp[,i] <- 
    logpop + a + a_site1 + a_cage + a_week + bfut*fut.seq[i] + bfut2*fut.seq[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort + bweight*weight 
  + btemp*temp.mid + bfut_temp*temp.mid*fut.seq[i] + bfut2_temp*temp.mid*fut.seq[i]^2
}  

# (logit)Prob of zero inflation
for (i in 1:length(fut.seq)){  
  p.med.temp[,i] <- 
    b + b_site1 + b_cage + b_week + pfut*fut.seq[i] + pfut2*fut.seq[i]^2 + ptsince*tsince.med 
  + ptsince2*tsince.med^2 + pablw*ablw + pemlw*emlw + pprev*prevmort + pweight*weight 
  + ptemp*temp.mid + pfut_temp*temp.mid*fut.seq[i] + pfut2_temp*temp.mid*fut.seq[i]^2
}

# Transforming to rate and probability
lambda.med.temp <- exp(lambda.med.temp)
p.med.temp <- inv.logit(p.med.temp)

# Estimating the expectation
rate.mid.temp <- lambda.med.temp*(1 - p.med.temp)

# Median, 66%, and 95% PI
median.rate.mid.temp <- apply( rate.mid.temp, 2 , median )
PI66.rate.mid.temp <- apply( rate.mid.temp, 2 , PI, prob=0.66 )
PI95.rate.mid.temp <- apply( rate.mid.temp, 2 , PI, prob=0.95 )


# Ploting

plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,200),
     main="Medium temp", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.mid.temp)
shade( PI95.rate.mid.temp , fut.seq , col=col.alpha("pink",0.2) )
shade( PI66.rate.mid.temp , fut.seq , col=col.alpha(rangi2,0.2) )


# High temp
lambda.high.temp <- matrix(rep(NA, 40*8000),ncol=40)
p.high.temp <- matrix(rep(NA, 40*8000),ncol=40)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.high.temp[,i] <- 
    logpop + a + a_site1 + a_cage + a_week + bfut*fut.seq[i] + bfut2*fut.seq[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort + bweight*weight 
  + btemp*temp.max + bfut_temp*temp.max*fut.seq[i] + bfut2_temp*temp.max*fut.seq[i]^2
}  

# (logit)Prob of zero inflation
for (i in 1:length(fut.seq)){  
  p.high.temp[,i] <- 
    b + b_site1 + b_cage + b_week + pfut*fut.seq[i] + pfut2*fut.seq[i]^2 + ptsince*tsince.med 
  + ptsince2*tsince.med^2 + pablw*ablw + pemlw*emlw + pprev*prevmort + pweight*weight 
  + ptemp*temp.max + pfut_temp*temp.max*fut.seq[i] + pfut2_temp*temp.max*fut.seq[i]^2
}

# Transforming to rate and probability
lambda.high.temp <- exp(lambda.high.temp)
p.high.temp <- inv.logit(p.high.temp)

# Estimating the expectation
rate.high.temp <- lambda.high.temp*(1 - p.high.temp)

# Median, 66%, and 95% PI
median.rate.high.temp <- apply( rate.high.temp, 2 , median )
PI66.rate.high.temp <- apply( rate.high.temp, 2 , PI, prob=0.66 )
PI95.rate.high.temp <- apply( rate.high.temp, 2 , PI, prob=0.95 )


# Ploting

plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,200),
     main="High temp", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.high.temp)
shade( PI95.rate.high.temp , fut.seq , col=col.alpha("pink",0.2) )
shade( PI66.rate.high.temp , fut.seq , col=col.alpha(rangi2,0.2) )



# Exploring other models
data2=data.frame(logpop = d$lpop, mort=d$srs_mort, futc=d$fut.f, tsince=d$tsince.f, 
          ablw=d$ab_lw - mean(d$ab_lw) , emlw=d$em_lw - mean(d$em_lw), 
          prevmort= as.vector(scale(log(d$prev.mort + 0.001), scale=sd(log(d$prev.mort + 0.001))*2)),
          dens=as.vector(scale(d$dens_lw, scale = sd(d$dens_lw)*2)), 
          site = as.numeric(d$site), cage=as.numeric(d$cage), time=d$time, 
          fut=as.vector(scale(d$fut_cum, scale=sd(d$fut_cum)*2)),
          fut2=as.vector(scale(d$fut_cum, scale=sd(d$fut_cum)*2)^2), 
          fut3=as.vector(scale(d$fut_cum, scale=sd(d$fut_cum)*2)^3), 
          timesince=as.vector(scale(d$time_since, scale = sd(d$time_since)*2)), 
          timesince2=as.vector(scale(d$time_since, scale = sd(d$time_since)*2)^2), 
          fut4 = d$fut4 - mean(d$fut4), fut8 = d$fut8 - mean(d$fut8),
          fut12 = d$fut12 - mean(d$fut12), CT = d$tsinceCT - mean(d$tsinceCT),
          W2 = d$tsince2W - mean(d$tsince2W), W2P = d$tsince2PW - mean(d$tsince2PW),
          weight = as.vector(scale(d$weight_mid, scale = sd(d$weight_mid)*2)),
          temp=as.vector(scale(d$temp, scale = sd(d$temp)*2)))

kk <- lm(log(mort+0.001) ~ fut + fut2 + timesince + timesince2 + fut:timesince + fut2:timesince2 + temp + fut:temp 
         + timesince:temp + fut2:temp + timesince2:temp + prevmort, data=data2)


# Adding interaction between fut and timesince
stancode3= "data {
int<lower=1> N;
int<lower=1> N_site;
int<lower=1> N_cage;
int<lower=1> N_time;
int<lower=0> mort[N];
real logpop[N];
real temp[N];
real ablw[N];
real emlw[N];
real prevmort[N];
int site[N];
int cage[N];
int time[N];
real fut[N];
real fut2[N];
real timesince[N];
real timesince2[N];
real weight[N];
}
parameters {
real a;
vector[N_site] a_site_raw;
vector[N_cage] a_cage_raw;
vector[N_time] a_week_raw;
real b;
vector[N_site] b_site_raw;
vector[N_cage] b_cage_raw;
vector[N_time] b_week_raw;
//vector[N_site] a_site;
//vector[N_cage] a_cage;
//vector[N_time] a_week;
real bfut;
real bfut2;
real pfut;
real pfut2;
real btsince;
real btsince2;
real bfut_tsince;
real bfut2_tsince2;
real bablw;
real bemlw;
real bprev;
real bweight;
real btemp;
real bfut_temp;
real bfut2_temp;
real ptsince;
real ptsince2;
real pfut_tsince;
real pfut2_tsince2;
real pablw;
real pemlw;
real pprev;
real pweight;
real ptemp;
real pfut_temp;
real pfut2_temp;
real<lower=0> sigma_site;
real<lower=0> sigma_cage;
real<lower=0> sigma_week;
real<lower=0> psigma_site;
real<lower=0> psigma_cage;
real<lower=0> psigma_week;
real<lower=0> theta;
//real<lower=0, upper=1> p;
//real<lower=0> lambda;
}
transformed parameters{
vector[N_site] a_site;
vector[N_cage] a_cage;
vector[N_time] a_week;
vector[N_site] b_site;
vector[N_cage] b_cage;
vector[N_time] b_week;
a_site = sigma_site*a_site_raw;
a_cage = sigma_cage*a_cage_raw;
a_week = sigma_week*a_week_raw;
b_site = psigma_site*b_site_raw;
b_cage = psigma_cage*b_cage_raw;
b_week = psigma_week*b_week_raw;

}
model {
vector[N] p;
vector[N] lambda;
a ~ normal( 0 , 1);
b ~ normal( 0 , 1);
//lambda ~ cauchy (0,2);
theta ~ exponential( 2 );
//p ~ cauchy(0, 5);
sigma_site ~ cauchy( 0 , 2);
sigma_cage ~ cauchy( 0 , 2);
sigma_week ~ cauchy( 0 , 2);
a_site_raw ~ normal( 0 , 1);
a_cage_raw ~ normal( 0 , 1);
a_week_raw ~ normal( 0 , 1);
psigma_site ~ cauchy( 0 , 2);
psigma_cage ~ cauchy( 0 , 2);
psigma_week ~ cauchy( 0 , 2);
b_site_raw ~ normal( 0 , 1);
b_cage_raw ~ normal( 0 , 1);
b_week_raw ~ normal( 0 , 1);
bfut ~ normal(0,0.5);
bfut2 ~ normal(0,0.5);
pfut ~ normal(0,0.5);
pfut2 ~ normal(0,0.5);
btsince ~ normal(0,0.5);
btsince2 ~ normal(0,0.5);
bfut_tsince ~ normal(0,0.5);
bfut2_tsince2 ~ normal(0,0.5);
bablw ~ normal(0,0.5);
bemlw ~ normal(0,0.5);
bprev ~ normal(0,0.5);
bweight ~ normal(0,0.5);
btemp ~ normal(0,0.5);
bfut_temp ~ normal(0,0.5);
bfut2_temp ~ normal(0,0.5);
ptsince ~ normal(0,0.5);
ptsince2 ~ normal(0,0.5);
pfut_tsince ~ normal(0,0.5);
pfut2_tsince2 ~ normal(0,0.5);
pablw ~ normal(0,0.5);
pemlw ~ normal(0,0.5);
pprev ~ normal(0,0.5);
pweight ~ normal(0,0.5);
ptemp ~ normal(0,0.5);
pfut_temp ~ normal(0,0.5);
pfut2_temp ~ normal(0,0.5);


for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ bfut*fut[i] + bfut2*fut2[i] + btsince*timesince[i] + btsince2*timesince2[i] 
+ bfut_tsince*fut[i]*timesince[i] + bfut2_tsince2*fut2[i]*timesince2[i] + bablw*ablw[i] + bemlw*emlw[i] + bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i] 
+ bfut_temp*temp[i]*fut[i] + bfut2_temp*temp[i]*fut2[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_site[site[i]] + b_cage[cage[i]] + b_week[time[i]] + pfut*fut[i] 
+ pfut2*fut2[i] + ptsince*timesince[i] + ptsince2*timesince2[i] + pfut_tsince*fut[i]*timesince[i] + pfut2_tsince2*fut2[i]*timesince2[i] + pablw*ablw[i] 
+ pemlw*emlw[i] + pprev*prevmort[i] + pweight*weight[i] + ptemp*temp[i] 
+ pfut_temp*temp[i]*fut[i] + pfut2_temp*temp[i]*fut2[i]);
}
for (n in 1:N) {
if (mort[n] == 0){
target += log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
target += bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
}
}

generated quantities {
real lambda[N];
real p[N];
real pred[N];
real resid[N];
real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ bfut*fut[i] + bfut2*fut2[i] + btsince*timesince[i] + btsince2*timesince2[i] 
+ bfut_tsince*fut[i]*timesince[i] + bfut2_tsince2*fut2[i]*timesince2[i] + bablw*ablw[i] + bemlw*emlw[i] + bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i] 
+ bfut_temp*temp[i]*fut[i] + bfut2_temp*temp[i]*fut2[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b + b_site[site[i]] + b_cage[cage[i]] + b_week[time[i]] + pfut*fut[i] 
+ pfut2*fut2[i] + ptsince*timesince[i] + ptsince2*timesince2[i] + pfut_tsince*fut[i]*timesince[i] + pfut2_tsince2*fut2[i]*timesince2[i] + pablw*ablw[i] 
+ pemlw*emlw[i] + pprev*prevmort[i] + pweight*weight[i] + ptemp*temp[i] 
+ pfut_temp*temp[i]*fut[i] + pfut2_temp*temp[i]*fut2[i]);
}

for (n in 1:N) {
if (mort[n] == 0){
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
log_lik[n] = bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
pred[n] = (1-p[n]) * lambda[n];
resid[n] = (mort[n] - pred[n])/sqrt(pred[n] + theta * pred[n]^2);
}
}
"

# Setting up inits
n_chains <- 4
start <- list(a=0, b=0, a_site = rep(0, data$N_site), b_site = rep(0, data$N_site), 
              a_cage = rep(0, data$N_cage), b_cage = rep(0, data$N_cage), 
              a_week = rep(0, data$N_time), b_week = rep(0, data$N_time), bfut=0, bfut2=0,
              btsince=0, btsince2=0, bfut_tsince=0, bfut2_tsince2=0, bablw=0, bemlw=0, bprev=0, bweight=0, btemp=0, 
              bfut_temp=0, bfut2_temp=0, pfut=0, pfut2=0, ptsince=0, ptsince2=0, 
              pfut_tsince=0, pfut2_tsince2=0, pablw=0, pemlw=0, pprev=0, pweight=0, ptemp=0, pfut_temp=0, pfut2_temp=0)
init <- list()
for ( i in 1:n_chains ) init[[i]] <- start

#Running model
m1.6 <- stan(model_code = stancode3, iter=4000, chains=4, cores=4, warmup=2000, control = list(adapt_delta = 0.95, stepsize=1), data=data, init=init, algorithm = "NUTS")

# Save Model
save(m1.6, file="m1.6.RData")

# Printing estimates of interest
#Population level effects
print(m1.6, probs=c(.025, .975), 
      pars=c("a", "bfut", "bfut2", "btsince", "btsince2", "bfut_tsince", "bfut2_tsince2", "bablw", "bemlw", 
             "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp", "b", 
             "pfut", "pfut2", "ptsince", "ptsince2", "pfut_tsince", "pfut2_tsince2", "pablw", "pemlw", "pprev", 
             "pweight", "ptemp", "pfut_temp", "pfut2_temp", "theta", "sigma_site", 
             "sigma_cage", "sigma_week", "psigma_site", "psigma_cage", "psigma_week"))
#Varying effects rate model: site, cage, week
print(m1.6, probs=c(.025, .975), pars=c("a_site", "a_cage", "a_week"))
#Varying effects logit model: site, cage, week
print(m1.6, probs=c(.025, .975), pars=c("b_site", "b_cage", "b_week"))


# Plotting population level effects
plot(m1.6, 
     pars=c("a", "bfut", "bfut2", "btsince", "btsince2", "bfut_tsince", "bfut2_tsince2", "bablw", "bemlw", 
            "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp", "b", 
            "pfut", "pfut2", "ptsince", "ptsince2", "pfut_tsince", "pfut2_tsince2", "pablw", "pemlw", "pprev", 
            "pweight", "ptemp", "pfut_temp", "pfut2_temp", "theta"))

# Traceplots
# Population level effects
traceplot(m1.5, pars=c("a", "bfut", "bfut2", "btsince", "btsince2", "bfut_tsince", "bfut2_tsince2", "bablw", "bemlw", 
                       "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp", "b", 
                       "pfut", "pfut2", "ptsince", "ptsince2", "pfut_tsince", "pfut2_tsince2", "pablw", "pemlw", "pprev", 
                       "pweight", "ptemp", "pfut_temp", "pfut2_temp", "theta", "sigma_site", 
                       "sigma_cage", "sigma_week", "psigma_site", "psigma_cage", "psigma_week"))
# Varying effects rate model
traceplot(m1.5, pars=c("a_site", "a_cage"))
traceplot(m1.5, pars=c("a_week"))
# Varying effects logit model
traceplot(m1.5, pars=c("b_site", "b_cage"))
traceplot(m1.5, pars=c("b_week"))

# Checking looic

log_lik2 <- extract_log_lik(m1.6)
(Cont_loo2 <- loo(log_lik2))


# Varying effects. Only treatment effects for negbin part (as model fit is bad)
stancode4= "data {
int<lower=1> N;
int<lower=1> N_site;
int<lower=1> N_cage;
int<lower=1> N_time;
int<lower=0> mort[N];
real logpop[N];
real temp[N];
real ablw[N];
real emlw[N];
real prevmort[N];
int site[N];
int cage[N];
int time[N];
real fut[N];
real fut2[N];
real timesince[N];
real timesince2[N];
real weight[N];
}
parameters {
real a;
//vector[N_site] a_site_raw;
vector[N_cage] a_cage_raw;
vector[N_time] a_week_raw;
real b;
vector[N_site] a_site;
real bfut;
real bfut2;
vector[N_site] bfut_site;
vector[N_site] bfut2_site;
real btsince;
real btsince2;
real bablw;
real bemlw;
real bprev;
real bweight;
real btemp;
real bfut_temp;
real bfut2_temp;
vector[N_site] bfut_temp_site;
vector[N_site] bfut2_temp_site;
vector<lower=0>[5] sigma_site;
real<lower=0> sigma_cage;
real<lower=0> sigma_week;
real<lower=0> theta;
corr_matrix[5] Rho_a;
}
transformed parameters{
vector[N_cage] a_cage;
vector[N_time] a_week;

// mortality rate model
vector[5] v_a_site_bfut_site_bfut2_site[N_site];
cov_matrix[5] SRS_sigma_siteRho;
for ( j in 1:N_site ) {
        v_a_site_bfut_site_bfut2_site[j,1] = a_site[j];
        v_a_site_bfut_site_bfut2_site[j,2] = bfut_site[j];
        v_a_site_bfut_site_bfut2_site[j,3] = bfut2_site[j];
        v_a_site_bfut_site_bfut2_site[j,4] = bfut_temp_site[j];
        v_a_site_bfut_site_bfut2_site[j,5] = bfut2_temp_site[j];
    }
    SRS_sigma_siteRho = quad_form_diag(Rho_a,sigma_site);

a_cage = sigma_cage*a_cage_raw;
a_week = sigma_week*a_week_raw;

}
model {
vector[N] p;
vector[N] lambda;
a ~ normal( 0 , 1);
b ~ normal( 0 , 1);
Rho_a ~ lkj_corr( 4 );
theta ~ exponential( 2 );
sigma_site ~ cauchy( 0 , 2);
sigma_cage ~ cauchy( 0 , 2);
sigma_week ~ cauchy( 0 , 2);
a_cage_raw ~ normal( 0 , 1);
a_week_raw ~ normal( 0 , 1);
bfut ~ normal(0,0.5);
bfut2 ~ normal(0,0.5);
bfut_site ~ normal(0,0.5);
bfut2_site ~ normal(0,0.5);
btsince ~ normal(0,0.5);
btsince2 ~ normal(0,0.5);
bablw ~ normal(0,0.5);
bemlw ~ normal(0,0.5);
bprev ~ normal(0,0.5);
bweight ~ normal(0,0.5);
btemp ~ normal(0,0.5);
bfut_temp ~ normal(0,0.5);
bfut2_temp ~ normal(0,0.5);
bfut_temp_site ~ normal(0,0.5);
bfut2_temp_site ~ normal(0,0.5);
v_a_site_bfut_site_bfut2_site ~ multi_normal( rep_vector(0,5) , SRS_sigma_siteRho );


for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ (bfut + bfut_site[site[i]])*fut[i] + (bfut2 + bfut2_site[site[i]])*fut2[i] 
+ btsince*timesince[i] + btsince2*timesince2[i] + bablw*ablw[i] + bemlw*emlw[i]
+ bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i]
+ (bfut_temp + bfut_temp_site[site[i]])*fut[i]*temp[i]
+ (bfut2_temp + bfut2_temp_site[site[i]])*fut2[i]*temp[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b);
}
for (n in 1:N) {
if (mort[n] == 0){
target += log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
target += bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
}
}

generated quantities {
real lambda[N];
real p[N];
real pred[N];
real resid[N];
real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ (bfut + bfut_site[site[i]])*fut[i] + (bfut2 + bfut2_site[site[i]])*fut2[i] 
+ btsince*timesince[i] + btsince2*timesince2[i] + bablw*ablw[i] + bemlw*emlw[i]
+ bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i]
+ (bfut_temp + bfut_temp_site[site[i]])*fut[i]*temp[i]
+ (bfut2_temp + bfut2_temp_site[site[i]])*fut2[i]*temp[i];
lambda[i] = exp(lambda[i]);
}

for ( i in 1:N ) {
p[i] = inv_logit(b);
}

for (n in 1:N) {
if (mort[n] == 0){
log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p[n]),
bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
}else{
log_lik[n] = bernoulli_lpmf(0 | p[n])
+ neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
}
pred[n] = (1-p[n]) * lambda[n];
resid[n] = (mort[n] - pred[n])/sqrt(pred[n] + theta * pred[n]^2);
}
}
"

n_chains <- 4
start <- list(a=0, b=0, a_site = rep(0, data$N_site), a_cage_raw = rep(0, data$N_cage), 
              a_week_raw = rep(0, data$N_time), bfut=0, bfut2=0, bfut_site = rep(0, data$N_site), 
              bfut2_site = rep(0, data$N_site), btsince=0, btsince2=0, bablw=0,
              bemlw=0, bprev=0, bweight=0, btemp=0, bfut_temp = 0, bfut2_temp=0,
              bfut_temp_site=rep(0, data$N_site), bfut2_temp_site=rep(0, data$N_site)
              )
init <- list()
for ( i in 1:n_chains ) init[[i]] <- start

#Running model
m1.7 <- stan(model_code = stancode4, iter=4000, chains=4, cores=4, warmup=2000, 
             control = list(adapt_delta = 0.95, stepsize=1, max_treedepth=15), data=data, init=init, algorithm = "NUTS")


# Save Model
save(m1.7, file="m1.7.RData")

# Printing estimates of interest
#Population level effects
print(m1.7, probs=c(.025, .975), 
      pars=c("a", "bfut", "bfut_site", "bfut2", "bfut2_site", "btsince", 
             "btsince2", "bablw", "bemlw", "bprev", "bweight", "btemp",
             "bfut_temp", "bfut_temp_site", "bfut2_temp", "bfut2_temp_site",
             "theta", "sigma_site", "sigma_cage", "sigma_week", "b"))
#Varying effects rate model: site, cage, week
print(m1.7, probs=c(.025, .975), pars=c("a_site", "a_cage", "a_week"))
#Varying effects logit model: site, cage, week
print(m1.7, probs=c(.025, .975), pars=c("b", "b_site", "b_cage", "b_week"))


# Plotting population level effects
plot(m1.7, 
     pars=c("a", "bfut", "bfut_site", "bfut2", "bfut2_site", "btsince", 
            "btsince2", "bablw", "bemlw", "bprev", "bweight", "btemp",
            "bfut_temp", "bfut_temp_site", "bfut2_temp", "bfut2_temp_site",
            "theta", "b"))
# Traceplots
# Population level effects
traceplot(m1.7, pars=c("a", "bfut", "bfut_site", "bfut2", "bfut2_site", "btsince", 
                       "btsince2", "bablw", "bemlw", "bprev", "bweight", "btemp",
                       "bfut_temp", "bfut_temp_site", "bfut2_temp", "bfut2_temp_site",
                       "theta", "sigma_site", "sigma_cage", "sigma_week", "b"))


# Varying effects rate model
traceplot(m1.7, pars=c("a_site", "a_cage"))
traceplot(m1.7, pars=c("a_week"))

log_lik3 <- extract_log_lik(m1.7)
(Cont_loo3 <- loo(log_lik3))

# Experiment with model implications
post.cont <- extract(m1.7)
a <- as.vector(post.cont$a)
a_site1 <- post.cont$a_site[,1]
a_site2 <- post.cont$a_site[,2]
a_site3 <- post.cont$a_site[,3]
a_cage <- apply(post.cont$a_cage, 1, FUN = "mean")
a_week <- apply(post.cont$a_week, 1, FUN = "mean")
b <- as.vector(post.cont$b)
bfut <- as.vector(post.cont$bfut)
bfut_site1 <- as.vector(post.cont$bfut_site[,1])
bfut_site2 <- as.vector(post.cont$bfut_site[,2])
bfut_site3 <- as.vector(post.cont$bfut_site[,3])
bfut2 <- post.cont$bfut2
bfut2_site1 <- as.vector(post.cont$bfut2_site[,1])
bfut2_site2 <- as.vector(post.cont$bfut2_site[,2])
bfut2_site3 <- as.vector(post.cont$bfut2_site[,3])
btsince <- post.cont$btsince
btsince2 <- post.cont$btsince2
bablw <- post.cont$bablw
bemlw <- post.cont$bablw
bprev <- post.cont$bablw
bweight <- post.cont$bweight
btemp <- post.cont$btemp
bfut_temp <- post.cont$bfut_temp
bfut_temp_site1 <- as.vector(post.cont$bfut_temp_site[,1])
bfut_temp_site2 <- as.vector(post.cont$bfut_temp_site[,2])
bfut_temp_site3 <- as.vector(post.cont$bfut_temp_site[,3])
bfut2_temp <- post.cont$bfut2_temp
bfut2_temp_site1 <- as.vector(post.cont$bfut2_temp_site[,1])
bfut2_temp_site2 <- as.vector(post.cont$bfut2_temp_site[,2])
bfut2_temp_site3 <- as.vector(post.cont$bfut2_temp_site[,3])

# Effect of cumulative exposure to fut for differnt levels of temp
# Values of predictors to be evaluated
fut.seq <- seq(from=min(data$fut), to=max(data$fut), length.out=80)
temp.low <- min(data$temp)
temp.mid <- median(data$temp)
temp.max <- max(data$temp)
tsince.med <- median(data$timesince)
ablw <- mean(data$ablw)
emlw <- mean(data$emlw)
prevmort <- median(data$prevmort)
weight <- median(data$weight)
logpop <- median(data$logpop)


# Estimating mortality rate predicted by the model for a range of cum exposure and temps
library(boot)
# Colors based on temp
data$tempcol <- ifelse(data$temp < -0.30744812, "green",
                  ifelse(data$temp >=-0.30744812 & data$temp< 0.09201039,"yellow", "red"))


#Site 1
# Low temp
lambda.low.temp <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.low.temp[,i] <- 
    (logpop + a + a_site1 + a_cage + a_week + (bfut + bfut_site1)*fut.seq[i] 
     + (bfut2 + bfut2_site1)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.low 
     + (bfut_temp + bfut_temp_site1)*temp.low*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site1)*temp.low*fut.seq[i]^2)
}  

# Transforming to rate and probability
lambda.low.temp <- exp(lambda.low.temp)
p.low.temp <- inv.logit(b)

# Estimating the expectation
rate.low.temp <- lambda.low.temp*(1 - p.low.temp)


# Median, 66%, and 95% PI
median.rate.low.temp <- apply( rate.low.temp, 2 , median )
PI66.rate.low.temp <- apply( rate.low.temp, 2 , PI, prob=0.66 )
PI95.rate.low.temp <- apply( rate.low.temp, 2 , PI, prob=0.95 )


# Ploting
# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,50),
     main="Low temp site 1", ylab="weekly mortality", 
     xlab="Cumulative weeks of treatment", xaxt = "n",  
     col=data$tempcol[data$site==1])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.low.temp)
shade( PI95.rate.low.temp , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.low.temp , fut.seq , col=col.alpha(rangi2,0.2) )



# Medium temp
lambda.med.temp <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.med.temp[,i] <- 
    (logpop + a + a_site1 + a_cage + a_week + (bfut + bfut_site1)*fut.seq[i] 
     + (bfut2 + bfut2_site1)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.mid 
     + (bfut_temp + bfut_temp_site1)*temp.mid*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site1)*temp.mid*fut.seq[i]^2)
}    


# Transforming to rate and probability
lambda.med.temp <- exp(lambda.med.temp)
p.med.temp <- inv.logit(b)

# Estimating the expectation
rate.mid.temp <- lambda.med.temp*(1 - p.med.temp)

# Median, 66%, and 95% PI
median.rate.mid.temp <- apply( rate.mid.temp, 2 , median )
PI66.rate.mid.temp <- apply( rate.mid.temp, 2 , PI, prob=0.66 )
PI95.rate.mid.temp <- apply( rate.mid.temp, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,50),
     main="Medium temp site 1", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n",
     col=data$tempcol[data$site==1])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.mid.temp)
shade( PI95.rate.mid.temp , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.mid.temp , fut.seq , col=col.alpha(rangi2,0.2) )


# High temp
lambda.high.temp <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.high.temp[,i] <- 
    (logpop + a + a_site1 + a_cage + a_week + (bfut + bfut_site1)*fut.seq[i] 
     + (bfut2 + bfut2_site1)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.max 
     + (bfut_temp + bfut_temp_site1)*temp.max*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site1)*temp.max*fut.seq[i]^2)
}    


# Transforming to rate and probability
lambda.high.temp <- exp(lambda.high.temp)
p.high.temp <- inv.logit(b)

# Estimating the expectation
rate.high.temp <- lambda.high.temp*(1 - p.high.temp)

# Median, 66%, and 95% PI
median.rate.high.temp <- apply( rate.high.temp, 2 , median )
PI66.rate.high.temp <- apply( rate.high.temp, 2 , PI, prob=0.66 )
PI95.rate.high.temp <- apply( rate.high.temp, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,50),
     main="High temp site 1", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n",
     col=data$tempcol[data$site==1])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.high.temp)
shade( PI95.rate.high.temp , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.high.temp , fut.seq , col=col.alpha(rangi2,0.2) )


#Site 2
# Low temp
lambda.low.temp2 <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.low.temp2[,i] <- 
    (logpop + a + a_site2 + a_cage + a_week + (bfut + bfut_site2)*fut.seq[i] 
     + (bfut2 + bfut2_site2)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.low 
     + (bfut_temp + bfut_temp_site2)*temp.low*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site2)*temp.low*fut.seq[i]^2)
}  

# Transforming to rate and probability
lambda.low.temp2 <- exp(lambda.low.temp2)
p.low.temp <- inv.logit(b)

# Estimating the expectation
rate.low.temp2 <- lambda.low.temp2*(1 - p.low.temp)


# Median, 66%, and 95% PI
median.rate.low.temp2 <- apply( rate.low.temp2, 2 , median )
PI66.rate.low.temp2 <- apply( rate.low.temp2, 2 , PI, prob=0.66 )
PI95.rate.low.temp2 <- apply( rate.low.temp2, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==2]), jitter(data$mort[data$site==2]), ylim=c(0,2),
     main="Low temp site 2", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n",
     col=data$tempcol[data$site==2])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.low.temp2)
shade( PI95.rate.low.temp2 , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.low.temp2 , fut.seq , col=col.alpha(rangi2,0.2) )



# Medium temp
lambda.mid.temp2 <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.mid.temp2[,i] <- 
    (logpop + a + a_site2 + a_cage + a_week + (bfut + bfut_site2)*fut.seq[i] 
     + (bfut2 + bfut2_site2)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.mid 
     + (bfut_temp + bfut_temp_site2)*temp.mid*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site2)*temp.mid*fut.seq[i]^2)
}  

# Transforming to rate and probability
lambda.mid.temp2 <- exp(lambda.mid.temp2)
p.mid.temp <- inv.logit(b)

# Estimating the expectation
rate.mid.temp2 <- lambda.mid.temp2*(1 - p.low.temp)


# Median, 66%, and 95% PI
median.rate.mid.temp2 <- apply( rate.mid.temp2, 2 , median )
PI66.rate.mid.temp2 <- apply( rate.mid.temp2, 2 , PI, prob=0.66 )
PI95.rate.mid.temp2 <- apply( rate.mid.temp2, 2 , PI, prob=0.95 )


# Ploting
# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==2]), jitter(data$mort[data$site==2]), ylim=c(0,2),
     main="Medium temp site 2", ylab="weekly mortality", 
     xlab="Cumulative weeks of treatment", xaxt = "n",  
     col=data$tempcol[data$site==2])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.mid.temp2)
shade( PI95.rate.mid.temp2 , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.mid.temp2 , fut.seq , col=col.alpha(rangi2,0.2) )


# High temp
lambda.high.temp2 <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.high.temp2[,i] <- 
    (logpop + a + a_site2 + a_cage + a_week + (bfut + bfut_site2)*fut.seq[i] 
     + (bfut2 + bfut2_site2)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.max 
     + (bfut_temp + bfut_temp_site2)*temp.max*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site2)*temp.max*fut.seq[i]^2)
}    


# Transforming to rate and probability
lambda.high.temp2 <- exp(lambda.high.temp2)
p.high.temp <- inv.logit(b)

# Estimating the expectation
rate.high.temp2 <- lambda.high.temp2*(1 - p.high.temp)

# Median, 66%, and 95% PI
median.rate.high.temp2 <- apply( rate.high.temp2, 2 , median )
PI66.rate.high.temp2 <- apply( rate.high.temp2, 2 , PI, prob=0.66 )
PI95.rate.high.temp2 <- apply( rate.high.temp2, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==2]), jitter(data$mort[data$site==2]), ylim=c(0,10),
     main="High temp site 2", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n",
     col=data$tempcol[data$site==2])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.high.temp2)
shade( PI95.rate.high.temp2 , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.high.temp2 , fut.seq , col=col.alpha(rangi2,0.2) )

#Site 3
# Low temp
lambda.low.temp3 <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.low.temp3[,i] <- 
    (logpop + a + a_site3 + a_cage + a_week + (bfut + bfut_site3)*fut.seq[i] 
     + (bfut2 + bfut2_site3)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.low 
     + (bfut_temp + bfut_temp_site3)*temp.low*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site3)*temp.low*fut.seq[i]^2)
}  

# Transforming to rate and probability
lambda.low.temp3 <- exp(lambda.low.temp3)
p.low.temp <- inv.logit(b)

# Estimating the expectation
rate.low.temp3 <- lambda.low.temp3*(1 - p.low.temp)


# Median, 66%, and 95% PI
median.rate.low.temp3 <- apply( rate.low.temp3, 2 , median )
PI66.rate.low.temp3 <- apply( rate.low.temp3, 2 , PI, prob=0.66 )
PI95.rate.low.temp3 <- apply( rate.low.temp3, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==3]), jitter(data$mort[data$site==3]), ylim=c(0,10),
     main="Low temp site 3", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n",
     col=data$tempcol[data$site==3])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.low.temp3)
shade( PI95.rate.low.temp3 , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.low.temp3 , fut.seq , col=col.alpha(rangi2,0.2) )



# Medium temp
lambda.mid.temp3 <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.mid.temp3[,i] <- 
    (logpop + a + a_site3 + a_cage + a_week + (bfut + bfut_site3)*fut.seq[i] 
     + (bfut2 + bfut2_site3)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.mid 
     + (bfut_temp + bfut_temp_site3)*temp.mid*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site3)*temp.mid*fut.seq[i]^2)
}  

# Transforming to rate and probability
lambda.mid.temp3 <- exp(lambda.mid.temp3)
p.mid.temp <- inv.logit(b)

# Estimating the expectation
rate.mid.temp3 <- lambda.mid.temp3*(1 - p.low.temp)


# Median, 66%, and 95% PI
median.rate.mid.temp3 <- apply( rate.mid.temp3, 2 , median )
PI66.rate.mid.temp3 <- apply( rate.mid.temp3, 2 , PI, prob=0.66 )
PI95.rate.mid.temp3 <- apply( rate.mid.temp3, 2 , PI, prob=0.95 )


# Ploting
# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==3]), jitter(data$mort[data$site==3]), ylim=c(0,10),
     main="Medium temp site 3", ylab="weekly mortality", 
     xlab="Cumulative weeks of treatment", xaxt = "n",  
     col=data$tempcol[data$site==3])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.mid.temp3)
shade( PI95.rate.mid.temp3 , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.mid.temp3 , fut.seq , col=col.alpha(rangi2,0.2) )


# High temp
lambda.high.temp3 <- matrix(rep(NA, 80*8000),ncol=80)

# (log)Lambda
for (i in 1:length(fut.seq)){
  lambda.high.temp3[,i] <- 
    (logpop + a + a_site1 + a_cage + a_week + (bfut + bfut_site3)*fut.seq[i] 
     + (bfut2 + bfut2_site3)*fut.seq[i]^2 + btsince*tsince.med 
     + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
     + bweight*weight + btemp*temp.max 
     + (bfut_temp + bfut_temp_site3)*temp.max*fut.seq[i] 
     + (bfut2_temp + bfut2_temp_site3)*temp.max*fut.seq[i]^2)
}    


# Transforming to rate and probability
lambda.high.temp3 <- exp(lambda.high.temp3)
p.high.temp <- inv.logit(b)

# Estimating the expectation
rate.high.temp3 <- lambda.high.temp3*(1 - p.high.temp)

# Median, 66%, and 95% PI
median.rate.high.temp3 <- apply( rate.high.temp3, 2 , median )
PI66.rate.high.temp3 <- apply( rate.high.temp3, 2 , PI, prob=0.66 )
PI95.rate.high.temp3 <- apply( rate.high.temp3, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==3]), jitter(data$mort[data$site==3]), ylim=c(0,10),
     main="High temp site 3", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n",
     col=data$tempcol[data$site==3])
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq, median.rate.high.temp3)
shade( PI95.rate.high.temp3 , fut.seq , col=col.alpha("black",0.2) )
shade( PI66.rate.high.temp3 , fut.seq , col=col.alpha(rangi2,0.2) )


# Negative binomial with Varying effects
stancode5= "data {
int<lower=1> N;
int<lower=1> N_site;
int<lower=1> N_cage;
int<lower=1> N_time;
int<lower=0> mort[N];
real logpop[N];
real temp[N];
real ablw[N];
real emlw[N];
real prevmort[N];
int site[N];
int cage[N];
int time[N];
real fut[N];
real fut2[N];
real timesince[N];
real timesince2[N];
real weight[N];
}
parameters {
real a;
vector[N_site] a_site;
vector[N_cage] a_cage_raw;
vector[N_time] a_week_raw;
real bfut;
real bfut2;
vector[N_site] bfut_site;
vector[N_site] bfut2_site;
real btsince;
real btsince2;
real bablw;
real bemlw;
real bprev;
real bweight;
real btemp;
real bfut_temp;
real bfut2_temp;
vector[N_site] bfut_temp_site;
vector[N_site] bfut2_temp_site;
vector<lower=0>[5] sigma_site;
real<lower=0> sigma_cage;
real<lower=0> sigma_week;
real<lower=0> theta;
corr_matrix[5] Rho_a;
}
transformed parameters{
vector[N_cage] a_cage;
vector[N_time] a_week;
// mortality rate model
vector[5] v_a_site_bfut_site_bfut2_site[N_site];
cov_matrix[5] SRS_sigma_siteRho;
for ( j in 1:N_site ) {
v_a_site_bfut_site_bfut2_site[j,1] = a_site[j];
v_a_site_bfut_site_bfut2_site[j,2] = bfut_site[j];
v_a_site_bfut_site_bfut2_site[j,3] = bfut2_site[j];
v_a_site_bfut_site_bfut2_site[j,4] = bfut_temp_site[j];
v_a_site_bfut_site_bfut2_site[j,5] = bfut2_temp_site[j];
}
SRS_sigma_siteRho = quad_form_diag(Rho_a,sigma_site);

a_cage = sigma_cage*a_cage_raw;
a_week = sigma_week*a_week_raw;
}

model {
vector[N] lambda;
a ~ normal( 0 , 1);
Rho_a ~ lkj_corr( 4 );
theta ~ exponential( 2 );
sigma_site ~ cauchy( 0 , 2);
sigma_cage ~ cauchy( 0 , 2);
sigma_week ~ cauchy( 0 , 2);
a_cage_raw ~ normal( 0 , 1);
a_week_raw ~ normal( 0 , 1);
bfut ~ normal(0,0.5);
bfut2 ~ normal(0,0.5);
btsince ~ normal(0,0.5);
btsince2 ~ normal(0,0.5);
bablw ~ normal(0,0.5);
bemlw ~ normal(0,0.5);
bprev ~ normal(0,0.5);
bweight ~ normal(0,0.5);
btemp ~ normal(0,0.5);
bfut_temp ~ normal(0,0.5);
bfut2_temp ~ normal(0,0.5);
v_a_site_bfut_site_bfut2_site ~ multi_normal( rep_vector(0,5) , SRS_sigma_siteRho );

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ (bfut + bfut_site[site[i]])*fut[i] + (bfut2 + bfut2_site[site[i]])*fut2[i] + btsince*timesince[i] + btsince2*timesince2[i] 
+ bablw*ablw[i] + bemlw*emlw[i] + bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i]
+ (bfut_temp + bfut_temp_site[site[i]])*fut[i]*temp[i] + (bfut2_temp + bfut2_temp_site[site[i]])*fut2[i]*temp[i];
lambda[i] = exp(lambda[i]);
}


for (n in 1:N) {
mort[n] ~ neg_binomial_2( lambda[n] , theta );
}
}

generated quantities {
real lambda[N];
real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_week[time[i]]
+ (bfut + bfut_site[site[i]])*fut[i] + (bfut2 + bfut2_site[site[i]])*fut2[i] + btsince*timesince[i] + btsince2*timesince2[i] 
+ bablw*ablw[i] + bemlw*emlw[i] + bprev*prevmort[i] + bweight*weight[i] + btemp*temp[i]
+ (bfut_temp + bfut_temp_site[site[i]])*fut[i]*temp[i] + (bfut2_temp + bfut2_temp_site[site[i]])*fut2[i]*temp[i];
lambda[i] = exp(lambda[i]);
}

for (n in 1:N) {
log_lik[n] = neg_binomial_2_lpmf( mort[n] | lambda[n] , theta );
}
}
"

n_chains <- 4
start <- list(a=0, a_site = rep(0, data$N_site), a_cage = rep(0, data$N_cage), 
              a_week = rep(0, data$N_time), bfut=0, bfut2=0,
              bfut_site = rep(0, data$N_site), bfut2_site = rep(0, data$N_site),
              btsince=0, btsince2=0, bablw=0, bemlw=0, bprev=0, bweight=0, btemp=0,
              bfut_temp = 0, bfut2_temp = 0, bfut_temp_site = rep(0, data$N_site), 
              bfut2_temp_site = rep(0, data$N_site))
init <- list()
for ( i in 1:n_chains ) init[[i]] <- start

#Running model
m1.8 <- stan(model_code = stancode5, iter=4000, chains=4, cores=4, warmup=2000, 
             control = list(adapt_delta = 0.95, stepsize=1, max_treedepth=15), data=data, init=init, algorithm = "NUTS")

save(m1.8, file = "m1.8.RData")

print(m1.8, probs=c(.025, .975), 
      pars=c("a", "bfut", "bfut2", "bfut_site", "bfut2_site", "btsince", "btsince2", 
             "bablw", "bemlw", "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp",
             "bfut_temp_site", "bfut2_temp_site", "theta", "sigma_site", "sigma_cage",
             "sigma_week"))

plot(m1.8, probs=c(.025, .975), 
      pars=c("a", "bfut", "bfut2", "bfut_site", "bfut2_site", "btsince", "btsince2", 
             "bablw", "bemlw", "bprev", "bweight", "btemp", "bfut_temp", "bfut2_temp",
             "bfut_temp_site", "bfut2_temp_site", "theta"))

log_lik4 <- extract_log_lik(m1.8)
(Cont_loo4 <- loo(log_lik4))


# Experiment with model implications
post.cont <- extract(m1.8)
a <- as.vector(post.cont$a)
a_site1 <- post.cont$a_site[,1]
a_site2 <- post.cont$a_site[,2]
a_site3 <- post.cont$a_site[,3]
a_cage <- apply(post.cont$a_cage, 1, FUN = "mean")
a_week <- apply(post.cont$a_week, 1, FUN = "mean")
bfut <- as.vector(post.cont$bfut)
bfut2 <- post.cont$bfut2
bfut_site1 <- post.cont$bfut_site[,1]
bfut_site2 <- post.cont$bfut_site[,2]
bfut_site3 <- post.cont$bfut_site[,3]
bfut2_site1 <- post.cont$bfut2_site[,1]
bfut2_site2 <- post.cont$bfut2_site[,2]
bfut2_site3 <- post.cont$bfut2_site[,3]
btsince <- post.cont$btsince
btsince2 <- post.cont$btsince2
bablw <- post.cont$bablw
bemlw <- post.cont$bablw
bprev <- post.cont$bablw
bweight <- post.cont$bweight
btemp <- post.cont$btemp
bfut_temp <- post.cont$bfut_temp
bfut2_temp <- post.cont$bfut2_temp
bfut_temp_site1 <- post.cont$bfut_temp_site[,1]
bfut_temp_site2 <- post.cont$bfut_temp_site[,2]
bfut_temp_site3 <- post.cont$bfut_temp_site[,3]
bfut2_temp_site1 <- post.cont$bfut2_temp_site[,1]
bfut2_temp_site2 <- post.cont$bfut2_temp_site[,2]
bfut2_temp_site3 <- post.cont$bfut2_temp_site[,3]

# Effect of cumulative exposure to fut for differnt levels of temp
# Values of predictors to be evaluated
fut.seq1 <- seq(from=min(data$fut[data$site==1]), to=max(data$fut[data$site==1]), length.out=80)
fut.seq2 <- seq(from=min(data$fut[data$site==2]), to=max(data$fut[data$site==2]), length.out=80)
fut.seq3 <- seq(from=min(data$fut[data$site==3]), to=max(data$fut[data$site==3]), length.out=80)
temp.low1 <- min(data$temp[data$site==1])
temp.mid1 <- median(data$temp[data$site==1])
temp.max1 <- max(data$temp[data$site==1])
temp.low2 <- min(data$temp[data$site==2])
temp.mid2 <- median(data$temp[data$site==2])
temp.max2 <- max(data$temp[data$site==2])
temp.low3 <- min(data$temp[data$site==3])
temp.mid3 <- median(data$temp[data$site==3])
temp.max3 <- max(data$temp[data$site==3])

tsince.med <- median(data$timesince)
ablw <- mean(data$ablw)
emlw <- mean(data$emlw)
prevmort <- median(data$prevmort)
weight <- median(data$weight)
logpop <- median(data$logpop)



# Estimating mortality rate predicted by the model for a range of cum exposure and temps
library(boot)

# Site 1
# Low temp
lambda.low.temp <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq1)){
  lambda.low.temp[,i] <- 
    (logpop + a + a_site1 + a_cage + a_week + (bfut + bfut_site1)*fut.seq1[i] 
  + (bfut2 + bfut2_site1)*fut.seq1[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.low1 
  + (bfut_temp + bfut_temp_site1)*temp.low1*fut.seq1[i] 
  + (bfut2_temp + bfut2_temp_site1)*temp.low1*fut.seq1[i]^2)
}  


# Transforming to rate
lambda.low.temp <- exp(lambda.low.temp)


# Estimating the expectation
rate.low.temp <- lambda.low.temp


# Median, 66%, and 95% PI
median.rate.low.temp <- apply( rate.low.temp, 2 , median )
PI66.rate.low.temp <- apply( rate.low.temp, 2 , PI, prob=0.66 )
PI95.rate.low.temp <- apply( rate.low.temp, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,50),
     main="Low temp site 1", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq1, median.rate.low.temp)
shade( PI95.rate.low.temp , fut.seq1 , col=col.alpha("black",0.2) )
shade( PI66.rate.low.temp , fut.seq1 , col=col.alpha(rangi2,0.2) )


# Medium temp
lambda.med.temp <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq1)){
  lambda.med.temp[,i] <- 
    (logpop + a + a_site1 + a_cage + a_week + (bfut + bfut_site1)*fut.seq1[i] 
  + (bfut2 + bfut2_site1)*fut.seq1[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.mid1 
  + (bfut_temp + bfut_temp_site1)*temp.mid1*fut.seq1[i] 
  + (bfut2_temp + bfut2_temp_site1)*temp.mid1*fut.seq1[i]^2)
}  

# Transforming to rate and probability
lambda.med.temp <- exp(lambda.med.temp)


# Estimating the expectation
rate.mid.temp <- lambda.med.temp

# Median, 66%, and 95% PI
median.rate.mid.temp <- apply( rate.mid.temp, 2 , median )
PI66.rate.mid.temp <- apply( rate.mid.temp, 2 , PI, prob=0.66 )
PI95.rate.mid.temp <- apply( rate.mid.temp, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,50),
     main="Medium temp site 1", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq1, median.rate.mid.temp)
shade( PI95.rate.mid.temp , fut.seq1 , col=col.alpha("black",0.2) )
shade( PI66.rate.mid.temp , fut.seq1, col=col.alpha(rangi2,0.2) )


# High temp
lambda.high.temp <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq1)){
  lambda.high.temp[,i] <- 
    (logpop + a + a_site1 + a_cage + a_week + (bfut + bfut_site1)*fut.seq1[i] 
  + (bfut2 + bfut2_site1)*fut.seq1[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.max1 
  + (bfut_temp + bfut_temp_site1)*temp.max1*fut.seq1[i] 
  + (bfut2_temp + bfut2_temp_site1)*temp.max1*fut.seq1[i]^2)
}  

# Transforming to rate and probability
lambda.high.temp <- exp(lambda.high.temp)

# Estimating the expectation
rate.high.temp <- lambda.high.temp

# Median, 66%, and 95% PI
median.rate.high.temp <- apply( rate.high.temp, 2 , median )
PI66.rate.high.temp <- apply( rate.high.temp, 2 , PI, prob=0.66 )
PI95.rate.high.temp <- apply( rate.high.temp, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==1]), jitter(data$mort[data$site==1]), ylim=c(0,50),
     main="High temp site 1", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq1, median.rate.high.temp)
shade( PI95.rate.high.temp , fut.seq1 , col=col.alpha("black",0.2) )
shade( PI66.rate.high.temp , fut.seq1 , col=col.alpha(rangi2,0.2) )




# Site 2
# Low temp
lambda.low.temp2 <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq2)){
  lambda.low.temp2[,i] <- 
    (logpop + a + a_site2 + a_cage + a_week + (bfut + bfut_site2)*fut.seq2[i] 
  + (bfut2 + bfut2_site2)*fut.seq2[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.low2 
  + (bfut_temp + bfut_temp_site2)*temp.low2*fut.seq2[i] 
  + (bfut2_temp + bfut2_temp_site2)*temp.low2*fut.seq2[i]^2)
}  


# Transforming to rate
lambda.low.temp2 <- exp(lambda.low.temp2)


# Estimating the expectation
rate.low.temp2 <- lambda.low.temp2


# Median, 66%, and 95% PI
median.rate.low.temp2 <- apply( rate.low.temp2, 2 , median )
PI66.rate.low.temp2 <- apply( rate.low.temp2, 2 , PI, prob=0.66 )
PI95.rate.low.temp2 <- apply( rate.low.temp2, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==2]), jitter(data$mort[data$site==2]), ylim=c(0,2),
     main="Low temp site 2", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq2, median.rate.low.temp2)
shade( PI95.rate.low.temp2 , fut.seq2 , col=col.alpha("black",0.2) )
shade( PI66.rate.low.temp2 , fut.seq2 , col=col.alpha(rangi2,0.2) )


# Medium temp
lambda.med.temp2 <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq2)){
  lambda.med.temp2[,i] <- 
    (logpop + a + a_site2 + a_cage + a_week + (bfut + bfut_site2)*fut.seq2[i] 
  + (bfut2 + bfut2_site2)*fut.seq2[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.mid2 
  + (bfut_temp + bfut_temp_site2)*temp.mid2*fut.seq2[i] 
  + (bfut2_temp + bfut2_temp_site2)*temp.mid2*fut.seq2[i]^2)
}  

# Transforming to rate and probability
lambda.med.temp2 <- exp(lambda.med.temp2)


# Estimating the expectation
rate.mid.temp2 <- lambda.med.temp2

# Median, 66%, and 95% PI
median.rate.mid.temp2 <- apply( rate.mid.temp2, 2 , median )
PI66.rate.mid.temp2 <- apply( rate.mid.temp2, 2 , PI, prob=0.66 )
PI95.rate.mid.temp2 <- apply( rate.mid.temp2, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==2]), jitter(data$mort[data$site==2]), ylim=c(0,2),
     main="Medium temp site 2", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq2, median.rate.mid.temp2)
shade( PI95.rate.mid.temp2 , fut.seq2 , col=col.alpha("black",0.2) )
shade( PI66.rate.mid.temp2 , fut.seq2 , col=col.alpha(rangi2,0.2) )


# High temp
lambda.high.temp2 <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq2)){
  lambda.high.temp2[,i] <- 
    (logpop + a + a_site2 + a_cage + a_week + (bfut + bfut_site2)*fut.seq2[i] 
  + (bfut2 + bfut2_site2)*fut.seq2[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.max2
  + (bfut_temp + bfut_temp_site2)*temp.max2*fut.seq2[i] 
  + (bfut2_temp + bfut2_temp_site2)*temp.max2*fut.seq2[i]^2)
}  

# Transforming to rate and probability
lambda.high.temp2 <- exp(lambda.high.temp2)

# Estimating the expectation
rate.high.temp2 <- lambda.high.temp2

# Median, 66%, and 95% PI
median.rate.high.temp2 <- apply( rate.high.temp2, 2 , median )
PI66.rate.high.temp2 <- apply( rate.high.temp2, 2 , PI, prob=0.66 )
PI95.rate.high.temp2 <- apply( rate.high.temp2, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==2]), jitter(data$mort[data$site==2]), 
     main="High temp site 2", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq2, median.rate.high.temp2)
shade( PI95.rate.high.temp2 , fut.seq2 , col=col.alpha("black",0.2) )
shade( PI66.rate.high.temp2 , fut.seq2 , col=col.alpha(rangi2,0.2) )


# Site 3
# Low temp
lambda.low.temp3 <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq3)){
  lambda.low.temp3[,i] <- 
    (logpop + a + a_site3 + a_cage + a_week + (bfut + bfut_site3)*fut.seq3[i] 
  + (bfut2 + bfut2_site3)*fut.seq3[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.low3 
  + (bfut_temp + bfut_temp_site3)*temp.low3*fut.seq3[i] 
  + (bfut2_temp + bfut2_temp_site3)*temp.low3*fut.seq3[i]^2)
}  


# Transforming to rate
lambda.low.temp3 <- exp(lambda.low.temp3)


# Estimating the expectation
rate.low.temp3 <- lambda.low.temp3


# Median, 66%, and 95% PI
median.rate.low.temp3 <- apply( rate.low.temp3, 2 , median )
PI66.rate.low.temp3 <- apply( rate.low.temp3, 2 , PI, prob=0.66 )
PI95.rate.low.temp3 <- apply( rate.low.temp3, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==3]), jitter(data$mort[data$site==3]), ylim=c(0,30),
     main="Low temp site 3", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq3, median.rate.low.temp3)
shade( PI95.rate.low.temp3 , fut.seq3 , col=col.alpha("black",0.2) )
shade( PI66.rate.low.temp3 , fut.seq3 , col=col.alpha(rangi2,0.2) )


# Medium temp
lambda.med.temp3 <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq3)){
  lambda.med.temp3[,i] <- 
    (logpop + a + a_site3 + a_cage + a_week + (bfut + bfut_site3)*fut.seq3[i] 
  + (bfut2 + bfut2_site3)*fut.seq3[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.mid3 
  + (bfut_temp + bfut_temp_site3)*temp.mid3*fut.seq3[i] 
  + (bfut2_temp + bfut2_temp_site3)*temp.mid3*fut.seq3[i]^2)
}  

# Transforming to rate and probability
lambda.med.temp3 <- exp(lambda.med.temp3)


# Estimating the expectation
rate.mid.temp3 <- lambda.med.temp3

# Median, 66%, and 95% PI
median.rate.mid.temp3 <- apply( rate.mid.temp3, 2 , median )
PI66.rate.mid.temp3 <- apply( rate.mid.temp3, 2 , PI, prob=0.66 )
PI95.rate.mid.temp3 <- apply( rate.mid.temp3, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==3]), jitter(data$mort[data$site==3]), ylim=c(0,30),
     main="Medium temp site 3", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq3, median.rate.mid.temp3)
shade( PI95.rate.mid.temp3 , fut.seq3 , col=col.alpha("black",0.2) )
shade( PI66.rate.mid.temp3 , fut.seq3 , col=col.alpha(rangi2,0.2) )


# High temp
lambda.high.temp3 <- matrix(rep(NA, 80*8000),ncol=80)


# (log)Lambda
for (i in 1:length(fut.seq3)){
  lambda.high.temp3[,i] <- 
    (logpop + a + a_site3 + a_cage + a_week + (bfut + bfut_site3)*fut.seq3[i] 
  + (bfut2 + bfut2_site3)*fut.seq3[i]^2 + btsince*tsince.med 
  + btsince2*tsince.med^2 + bablw*ablw + bemlw*emlw + bprev*prevmort 
  + bweight*weight + btemp*temp.max3 
  + (bfut_temp + bfut_temp_site3)*temp.max3*fut.seq3[i] 
  + (bfut2_temp + bfut2_temp_site3)*temp.max3*fut.seq3[i]^2)
}  

# Transforming to rate and probability
lambda.high.temp3 <- exp(lambda.high.temp3)

# Estimating the expectation
rate.high.temp3 <- lambda.high.temp3

# Median, 66%, and 95% PI
median.rate.high.temp3 <- apply( rate.high.temp3, 2 , median )
PI66.rate.high.temp3 <- apply( rate.high.temp3, 2 , PI, prob=0.66 )
PI95.rate.high.temp3 <- apply( rate.high.temp3, 2 , PI, prob=0.95 )


# Ploting
set.seed(0)
plot(jitter(data$fut[data$site==3]), jitter(data$mort[data$site==3]), ylim=c(0,100),
     main="High temp site 3", ylab="weekly mortality", xlab="Cumulative weeks of treatment", xaxt = "n")
axis(side = 1, cex.axis=8/12, at = seq(from=min(data$fut), to=max(data$fut), 
                                       length.out=20), lwd=0.5, labels =round(seq(from=0, to=19, length.out=20), digits=0))
lines(fut.seq3, median.rate.high.temp3)
shade( PI95.rate.high.temp3 , fut.seq3 , col=col.alpha("black",0.2) )
shade( PI66.rate.high.temp3 , fut.seq3 , col=col.alpha(rangi2,0.2) )


# Last chance: 2 groups treatment and control
data2 <- list(mort=d$srs_mort, site=as.numeric(d$site), cage=as.numeric(d$cage),
              time=d$time, time2=d$time^2, logpop=log(d$pop_mid), 
              group=as.numeric(d$group)-1, N=nrow(d), 
              N_site=length(unique(d$site)), N_cage=length(unique(d$cage)))

stancode6= "data {
int<lower=1> N;
int<lower=1> N_site;
int<lower=1> N_cage;
int<lower=0> mort[N];
real logpop[N];
real group[N];
int site[N];
int cage[N];
int time[N];
int time2[N];
}
parameters {
real a;
vector[N_site] a_site_raw;
real bgroup;
real btime;
real btime2;
real bfut;
real bfut2;
real<lower=0> sigma_site;
real<lower=0> theta;
}
transformed parameters{
vector[N_site] a_site;
a_site = a_site_raw*sigma_site;
}

model {
vector[N] lambda;
a ~ normal( 0 , 1);
a_site_raw ~ normal(0,1);
theta ~ exponential( 2 );
sigma_site ~ cauchy( 0 , 2);
bfut ~ normal(0,0.5);
bfut2 ~ normal(0,0.5);
bgroup ~ normal(0,0.5);
btime ~ normal(0,0.5);
btime2 ~ normal(0,0.5);


for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] 
+ bfut*group[i]*time[i] + bfut2*group[i]*time2[i] 
+ bgroup*group[i] + btime*time[i] + btime2*time2[i];
lambda[i] = exp(lambda[i]);
}


for (n in 1:N) {
mort[n] ~ neg_binomial_2( lambda[n] , theta );
}
}

generated quantities {
real lambda[N];
real log_lik[N];

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] 
+ bfut*group[i]*time[i] + bfut2*group[i]*time2[i] 
+ bgroup*group[i] + btime*time[i] + btime2*time2[i];
lambda[i] = exp(lambda[i]);
}

for (n in 1:N) {
log_lik[n] = neg_binomial_2_lpmf( mort[n] | lambda[n] , theta );
}
}
"

n_chains <- 4
start <- list(a=0, a_site_raw = rep(0, data2$N_site), bfut=0, bfut2=0,
              bgroup=0, btime=0, btime2=0)
init <- list()
for ( i in 1:n_chains ) init[[i]] <- start

#Running model
m1.9 <- stan(model_code = stancode6, iter=4000, chains=4, cores=4, warmup=2000, 
             control = list(adapt_delta = 0.95, stepsize=1, max_treedepth=15), 
             data=data2, init=init, algorithm = "NUTS")

save(m1.9, file = "m1.9.RData")
# It does not work by this measure either
print(m1.9, probs=c(0.025, 0.975), pars=c("a", "a_site", "bfut", "bfut2", "bgroup", "btime",
                                          "btime2", "sigma_site"))

# Bad fit (1000+ more looic)
log_lik5 <- extract_log_lik(m1.9)
(Cont_loo5 <- loo(log_lik5))


m6.1 <- map2stan(
  alist(mort ~ dgampois( lambda , theta),
        log(lambda) <- logpop + a + a_site[site] + a_cage[cage] + a_week[time] 
        + bfut*group + bfut_time*group*time + bfut_time2*group*time2,
        c(a, bfut, bfut_time, bfut_time2) ~ dnorm(0,1),
        theta ~ dexp(2),
        a_site[site] ~ dnorm(0,1),
        a_cage[cage] ~ dnorm(0,1),
        a_week[time] ~ dnorm(0,1)
        
  ), data=data,
  iter=3000 , warmup=1000 , chains=4 ,cores=4,
  control = list(adapt_delta = 0.8, stepsize=1))

precis(m6.1, depth=2)
plot(m6.1)
save(m6.1, file="m6.1.RData")

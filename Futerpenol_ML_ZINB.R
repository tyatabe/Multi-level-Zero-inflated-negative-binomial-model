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
             "pweight", "ptemp", "pfut4_temp", "pfut8_temp", "pfut12_temp", "theta"))
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
                       "pfut8_temp", "pfut12_temp", "theta"))
# Varying effects rate model
traceplot(m1.4, pars=c("a_site", "a_cage"))
traceplot(m1.4, pars=c("a_week"))
# Varying effects logit model
traceplot(m1.4, pars=c("b_site", "b_cage"))
traceplot(m1.4, pars=c("b_week"))




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
# treatment and cum weeks since last treatment) as numerical predictors
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
  real pred[N];
  real resid[N];
  real log_lik[N];
  for (n in 1:N) {
    if (mort[n] == 0){
      log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | p[n]),
      bernoulli_lpmf(0 | p[n])
      + neg_binomial_2_lpmf(mort[n] | lambda[n], theta));
      }else{
      log_lik[n] = bernoulli_lpmf(0 | p[n])
      + neg_binomial_2_lpmf(mort[n] | lambda[n], theta);
      }
   pred[n] <- (1-p[n]) * lambda[n];
   resid[n] <- (mort[n] - pred[n])/sqrt(pred[n] + theta * pred[n]^2);
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
             "pweight", "ptemp", "pfut4_temp", "pfut8_temp", "pfut12_temp", "theta"))
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
                       "pfut8_temp", "pfut12_temp", "theta"))
# Varying effects rate model
traceplot(m1.4, pars=c("a_site", "a_cage"))
traceplot(m1.4, pars=c("a_week"))
# Varying effects logit model
traceplot(m1.4, pars=c("b_site", "b_cage"))
traceplot(m1.4, pars=c("b_week"))



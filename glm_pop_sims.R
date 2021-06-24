
# Running Generalized Linear Models

# Code written by Alissa Brown.
# This code reads and restructures simulated data for inclusion in GLMs.
# We run 4 GLMs on 4 different sets of data: 
# 1) generic, 2) generic with bottlenecks, 3) oaks (case study), 4) oaks with bottlenecks

# Prepare work directory and packages
# setwd('C:/Users/abrow/Documents/pop_sims')
require(lsmeans)
require(tidyr)
require(dplyr)


#### GENERIC GLMS ####
# read in numbers of alleles (total and captured)
load('total_alleles.RData')
load('alleles_capt_highMig_highSamp.Rdata')
load('alleles_capt_highMig_lowSamp.Rdata')
load('alleles_capt_lowMig_lowSamp.Rdata')
load('alleles_capt_lowMig_highSamp.Rdata')

# get total number of alleles in populations
total_list <- list(total_alleles_highMig, total_alleles_lowMig)
for(i in 1:length(total_list)){
  total_list[[i]] <- as.data.frame(t(total_list[[i]]))
  total_list[[i]] <- pivot_longer(total_list[[i]], cols = colnames(total_list[[i]]),
                                  values_to = 'total_num', names_to = 'scenario')
}

# make total # alleles vectors for each of the 8 combinations of treatments
total_alleles <- c(rep(list(total_list[[1]]$total_num), 4),
                   rep(list(total_list[[2]]$total_num), 4))

# create list of 8 dataframes, one for each combination of treatments,
# add total number of alleles vectors, then combine all into one dataframe
captured_list <- list(samp_all_highMig_highSamp_equal, samp_all_highMig_highSamp_prop, 
                      samp_all_highMig_lowSamp_equal, samp_all_highMig_lowSamp_prop, 
                      samp_all_lowMig_highSamp_equal, samp_all_lowMig_highSamp_prop, 
                      samp_all_lowMig_lowSamp_equal, samp_all_lowMig_lowSamp_prop)

migration <- c(rep('high',4), rep('low',4))
intensity <- rep(c('high','high','low','low'), 2)
strategy <- rep(c('equal','prop'), 4)

for(i in 1:length(captured_list)){
  captured_list[[i]] <- as.data.frame(t(captured_list[[i]]))
  captured_list[[i]] <- pivot_longer(captured_list[[i]], cols = colnames(captured_list[[i]]), 
                                     values_to = 'success', names_to = 'scenario')
  captured_list[[i]]$scenario <- substr(captured_list[[i]]$scenario, start = 2, 
                                        stop = nchar(captured_list[[i]]$scenario))
  captured_list[[i]]$tot_alleles <- total_alleles[[i]]
  captured_list[[i]]$migration <- migration[i]
  captured_list[[i]]$intensity <- intensity[i]
  captured_list[[i]]$strategy <- strategy[i]
}
dat <- bind_rows(captured_list)
cols <- c('scenario', 'strategy','migration','intensity')
dat[,cols] <- lapply(dat[,cols], factor)

# calculate coefficient of variation for each scenario
cv1 <- sd(c(30, 100, 100, 100, 1170)) / mean(c(30, 100, 100, 100, 1170))
cv2 <- sd(c(40, 150, 150, 150, 1010)) / mean(c(40, 150, 150, 150, 1010))
cv3 <- sd(c(50, 200, 200, 200, 850)) / mean(c(50, 200, 200, 200, 850))
cv4 <- sd(c(100, 200, 200, 200, 800)) / mean(c(100, 200, 200, 200, 800))
cv5 <- sd(c(150, 200, 200, 200, 750)) / mean(c(150, 200, 200, 200, 750))
cv6 <- sd(c(200, 250, 250, 250, 550)) / mean(c(200, 250, 250, 250, 550))
cv7 <- sd(c(200, 300, 300, 300, 400)) / mean(c(200, 300, 300, 300, 400))
cv8 <- sd(c(290, 300, 300, 300, 310)) / mean(c(290, 300, 300, 300, 310))
cv9 <- sd(c(300, 300, 300, 300, 300)) / mean(c(300, 300, 300, 300, 300))
cv <- c(cv1,cv2,cv3,cv4,cv5,cv6,cv7,cv8,cv9)
dat$cv <- as.numeric(NA)
dat$scenario <- as.integer(dat$scenario)
for(i in 1:length(cv)){
  dat[dat$scenario == i,]$cv <- cv[i]
}
dat$fail <- dat$tot_alleles - dat$success
saveRDS(dat, 'data/generic/df_for_GLM.RDS')
dat <- readRDS('data/generic/df_for_GLM.RDS')

# prepare response variable
y <- as.matrix(dat[,c('success','fail')])

# run GLM
mod <- glm(y ~ scenario + migration + intensity + strategy + 
             scenario:strategy + 
             scenario:migration + 
             scenario:strategy:intensity:migration,
           family = binomial(link = 'logit'), data = dat)
summary(mod)

# construct pairwise contrasts
lsmeans(mod, poly ~ scenario)
# capture more alleles as scenario # increases

lsmeans(mod, pairwise ~ strategy|migration|scenario)

lsmeans(mod, pairwise ~ intensity|scenario)

lsmeans(mod, pairwise ~ strategy|scenario)
# diff'ce between equal/prop strategy significant for scenarios 1-6

lsmeans(mod, pairwise ~ migration)
# when migration rates low, capture more alleles (-0.156, p<0.0001)

lsmeans(mod, pairwise ~ intensity)
# when intensity high, capture more alleles
lsmeans(mod, pairwise ~ intensity|strategy|scenario)
# for all combinations of strategy and scenario, high intensity always captures more alleles
lsmeans(mod, pairwise ~ strategy|intensity|scenario)
# for scenarios 1-5, proportional captures more alleles (for both intensity levels, 
# but more so for high intensity levels)

lsmeans(mod, pairwise ~ migration|scenario)
# for almost all scenarios (except 8), high migration means capturing fewer alleles

lsmeans(mod, pairwise ~ strategy|intensity|migration|scenario)
# proportional better to use for all combinations of intensity and migration for sc 1-5
# proportional only better in scenario 6 when intensity high and migration low
# strategy doesn't matter for any intensity/migration combo for scenarios 7-9

### MODEL VALIDATION
anova(mod, test = 'Chisq') # goodness of fit measures

# check for overdispersion
deviance(mod) / mod$df.residual
# 0.81

# plot residuals
res <- residuals(mod)
qqline(res)
hist(res)
plot(y[,1], res)
abline(h = 0, col = 'red')

dat$res <- res
with(dat[dat$migration == 'high',], plot(success, res))
abline(h = 0, col = 'red')

with(dat[dat$migration == 'low',], plot(success, res))
abline(h = 0, col = 'red')

# plot predicted data against observations
pred <- predict.glm(mod6, type = 'response', se.fit = TRUE)
dat$prop <- dat$success / dat$tot_alleles
dat$fit <- pred$fit
plot(dat$prop, dat$fit, xlim = c(0.75, 1), ylim = c(0.75, 1))
abline(a = 0, b = 1, col = 'red')




#### GLM FOR OAKS ####
# read and restructure data
# QUERCUS ACERIFOLIA
load('alleles_capt_q_acer.RData')
load('combined_results_q_acerifolia.RData')
acer_h <- combined_q_acerifolia_high
acer_l <- combined_q_acerifolia_low
acer_capt_h <- c(alleles_cap_q_acerifolia_equal_high[1,],
                 alleles_cap_q_acerifolia_prop_high[1,])
acer_h$alleles_capt <- acer_capt_h
acer_h$intensity <- 'high'
acer_capt_l <- c(alleles_cap_q_acerifolia_equal_low[1,],
                 alleles_cap_q_acerifolia_prop_low[1,])
acer_l$alleles_capt <- acer_capt_l
acer_l$intensity <- 'low'
acer <- rbind(acer_l, acer_h)
acer$species <- 'q_acer'

# QUERCUS ENGELMANNII
load('alleles_capt_q_engel.RData')
load('combined_results_q_engelmannii.RData')
engel_h <- combined_q_engelmannii_high
engel_l <- combined_q_engelmannii_low
engel_capt_h <- c(alleles_cap_q_engel_equal_high[1,],
                  alleles_cap_q_engel_prop_high[1,])
engel_h$alleles_capt <- engel_capt_h
engel_h$intensity <- 'high'
engel_capt_l <- c(alleles_cap_q_engel_equal_low[1,],
                  alleles_cap_q_engel_prop_low[1,])
engel_l$alleles_capt <- engel_capt_l
engel_l$intensity <- 'low'
engel <- rbind(engel_l, engel_h)
engel$species <- 'q_engel'

# QUERCUS OGLETHORPENSIS
load('alleles_capt_q_ogle.RData')
load('combined_results_q_oglethorpensis.RData')
ogle_h <- combined_q_oglethorpensis_high
ogle_l <- combined_q_oglethorpensis_low
ogle_capt_h <- c(alleles_cap_q_ogle_equal_high[1,],
                 alleles_cap_q_ogle_prop_high[1,])
ogle_h$alleles_capt <- ogle_capt_h
ogle_h$intensity <- 'high'
ogle_capt_l <- c(alleles_cap_q_ogle_equal_low[1,],
                 alleles_cap_q_ogle_prop_low[1,])
ogle_l$alleles_capt <- ogle_capt_l
ogle_l$intensity <- 'low'
ogle <- rbind(ogle_l, ogle_h)
ogle$species <- 'q_ogle'

oak <- rbind(acer, engel, ogle)
oak$tot_alleles <- oak$alleles_capt / oak$prop_all
oak$fail <- oak$tot_alleles - oak$alleles_capt
cols <- c('strategy','intensity','species')
oak[,cols] <- lapply(oak[,cols], factor)
saveRDS(oak, 'data/oak/oak_df_for_GLM.RDS')
oak <- readRDS('data/oak/oak_df_for_GLM.RDS')

# prepare response variable
y <- as.matrix(oak[,c('alleles_capt','fail')])

# run GLM
mod_oak <- glm(y ~ species + strategy + intensity + 
                 species:strategy + 
                 species:intensity + 
                 species:strategy:intensity,
               family = binomial(link = 'logit'), data = oak)
anova(mod_oak, test = 'Chisq')
summary(mod_oak)

lsmeans(mod_oak, pairwise ~ strategy|species)
# for all species, proportional sampling works better, especially for engelmannii

lsmeans(mod_oak, pairwise ~ intensity|species)
# for all species, high intensity works better

lsm <- lsmeans(mod_oak, pairwise ~ strategy|intensity|species)
# proportional sampling works better, regardless of sampling intensity, for all species
# somewhat more effective when coupled with high intensity sampling for oglethorpensis

# summarize contrast results for supplemental file
lsm_df <- summary(lsm$contrasts)
lsm_df <- lsm_df %>% select(species, intensity, estimate, SE, z.ratio, p.value)
lsm_df[,3:6] <- signif(lsm_df[,3:6], 4)




#### GENERIC BOTTLENECKS ####
# BOTTLENECK 1 DATA
# read in numbers of alleles (total and captured)
load('data/generic/bottleneck1/total_alleles.RData')
load('data/generic/bottleneck1/alleles_capt_highMig_highSamp.Rdata')
load('data/generic/bottleneck1/alleles_capt_highMig_lowSamp.Rdata')
load('data/generic/bottleneck1/alleles_capt_lowMig_lowSamp.Rdata')
load('data/generic/bottleneck1/alleles_capt_lowMig_highSamp.Rdata')

# get total number of alleles in populations
total_list <- list(total_alleles_highMig, total_alleles_lowMig)
for(i in 1:length(total_list)){
  total_list[[i]] <- as.data.frame(t(total_list[[i]]))
  total_list[[i]] <- pivot_longer(total_list[[i]], cols = colnames(total_list[[i]]),
                                  values_to = 'total_num', names_to = 'scenario')
}

# make total # alleles vectors for each of the 8 combinations of treatments
total_alleles <- c(rep(list(total_list[[1]]$total_num), 4),
                   rep(list(total_list[[2]]$total_num), 4))

# create list of 8 dataframes, one for each combination of treatments,
# add total number of alleles vectors, then combine all into one dataframe
captured_list <- list(samp_all_highMig_highSamp_equal, samp_all_highMig_highSamp_prop, 
                      samp_all_highMig_lowSamp_equal, samp_all_highMig_lowSamp_prop, 
                      samp_all_lowMig_highSamp_equal, samp_all_lowMig_highSamp_prop, 
                      samp_all_lowMig_lowSamp_equal, samp_all_lowMig_lowSamp_prop)

migration <- c(rep('high',4), rep('low',4))
intensity <- rep(c('high','high','low','low'), 2)
strategy <- rep(c('equal','prop'), 4)

for(i in 1:length(captured_list)){
  captured_list[[i]] <- as.data.frame(t(captured_list[[i]]))
  captured_list[[i]] <- pivot_longer(captured_list[[i]], cols = colnames(captured_list[[i]]), 
                                     values_to = 'success', names_to = 'scenario')
  captured_list[[i]]$scenario <- substr(captured_list[[i]]$scenario, start = 2, 
                                        stop = nchar(captured_list[[i]]$scenario))
  captured_list[[i]]$tot_alleles <- total_alleles[[i]]
  captured_list[[i]]$migration <- migration[i]
  captured_list[[i]]$intensity <- intensity[i]
  captured_list[[i]]$strategy <- strategy[i]
}
bottle1 <- bind_rows(captured_list)

# BOTTLENECK 2 DATA
# read in numbers of alleles (total and captured)
load('data/generic/bottleneck2/total_alleles.RData')
load('data/generic/bottleneck2/alleles_capt_highMig_highSamp.Rdata')
load('data/generic/bottleneck2/alleles_capt_highMig_lowSamp.Rdata')
load('data/generic/bottleneck2/alleles_capt_lowMig_lowSamp.Rdata')
load('data/generic/bottleneck2/alleles_capt_lowMig_highSamp.Rdata')

# get total number of alleles in populations
total_list <- list(total_alleles_highMig, total_alleles_lowMig)
for(i in 1:length(total_list)){
  total_list[[i]] <- as.data.frame(t(total_list[[i]]))
  total_list[[i]] <- pivot_longer(total_list[[i]], cols = colnames(total_list[[i]]),
                                  values_to = 'total_num', names_to = 'scenario')
}

# make total # alleles vectors for each of the 8 combinations of treatments
total_alleles <- c(rep(list(total_list[[1]]$total_num), 4),
                   rep(list(total_list[[2]]$total_num), 4))

# create list of 8 dataframes, one for each combination of treatments,
# add total number of alleles vectors, then combine all into one dataframe
captured_list <- list(samp_all_highMig_highSamp_equal, samp_all_highMig_highSamp_prop, 
                      samp_all_highMig_lowSamp_equal, samp_all_highMig_lowSamp_prop, 
                      samp_all_lowMig_highSamp_equal, samp_all_lowMig_highSamp_prop, 
                      samp_all_lowMig_lowSamp_equal, samp_all_lowMig_lowSamp_prop)

migration <- c(rep('high',4), rep('low',4))
intensity <- rep(c('high','high','low','low'), 2)
strategy <- rep(c('equal','prop'), 4)

for(i in 1:length(captured_list)){
  captured_list[[i]] <- as.data.frame(t(captured_list[[i]]))
  captured_list[[i]] <- pivot_longer(captured_list[[i]], cols = colnames(captured_list[[i]]), 
                                     values_to = 'success', names_to = 'scenario')
  captured_list[[i]]$scenario <- substr(captured_list[[i]]$scenario, start = 2, 
                                        stop = nchar(captured_list[[i]]$scenario))
  captured_list[[i]]$tot_alleles <- total_alleles[[i]]
  captured_list[[i]]$migration <- migration[i]
  captured_list[[i]]$intensity <- intensity[i]
  captured_list[[i]]$strategy <- strategy[i]
}
bottle2 <- bind_rows(captured_list)

# combine bottleneck dataframes and prepare data for GLM
bottle1$bottleneck <- 1
bottle2$bottleneck <- 2
dat <- rbind(bottle1, bottle2)
cols <- c('scenario', 'strategy','migration','intensity', 'bottleneck')
dat[,cols] <- lapply(dat[,cols], factor)
saveRDS('data/generic_bottlenecks/bottleneck_df_for_GLM.RDS')
dat <- readRDS('data/generic_bottlenecks/bottleneck_df_for_GLM.RDS')

# prepare response variable
y <- as.matrix(dat[,c('success','tot_alleles')])

# run GLM
mod <- glm(y ~ scenario + migration + intensity + strategy + bottleneck +
              scenario:strategy + 
              scenario:migration + 
              scenario:strategy:intensity:migration + 
              scenario:strategy:migration:bottleneck,
            family = binomial(link = 'logit'), data = dat)
summary(mod)
anova(mod, test = 'Chisq')

# Use lsmeans to construct pairwise contrasts
lsmeans(mod, pairwise ~ strategy|scenario)
# diff'ce between equal/prop strategy significant for scenarios 1-6

lsmeans(mod, pairwise ~ migration)
# when migration rates low, capture more alleles (-0.156, p<0.0001)

lsmeans(mod, pairwise ~ intensity)
# when intensity high, capture more alleles
lsmeans(mod, pairwise ~ intensity|strategy|scenario)
# for all combinations of strategy and scenario, high intensity always captures more alleles
lsmeans(mod, pairwise ~ strategy|intensity|scenario)
# for scenarios 1-5, proportional captures more alleles (for both intensity levels, 
# but more so for high intensity levels)

lsmeans(mod, pairwise ~ migration|scenario)
# for almost all scenarios (except 8), high migration means capturing fewer alleles

lsm <- lsmeans(mod, pairwise ~ strategy|intensity|migration|scenario)
# proportional better to use for all combinations of intensity and migration for sc 1-5
# proportional only better in scenario 6 when intensity high and migration low
# strategy doesn't matter for any intensity/migration combo for scenarios 7-9

# summarize contrast results for supplemental file
lsm_df <- summary(lsm$contrasts)
lsm_df <- lsm_df %>% select(intensity, migration, scenario, p.value)
lsm_df$type <- paste0(lsm_df$migration, '_', lsm_df$intensity)
lsm_df <- lsm_df %>% select(-migration, -intensity)
lsm_df <- pivot_wider(lsm_df, names_from = 'type', values_from = 'p.value')
lsm_df[,2:5] <- signif(lsm_df[,2:5], 3)




#### OAK BOTTLENECKS ####
# Read and restructure data
# BOTTLENECK 1
# QUERCUS ACERIFOLIA
load('data/oak_bottlenecks/bottleneck1/alleles_capt_q_acer.RData')
load('data/oak_bottlenecks/bottleneck1/combined_results_q_acerifolia.RData')
acer_h <- combined_q_acerifolia_high
acer_l <- combined_q_acerifolia_low
acer_capt_h <- c(alleles_cap_q_acerifolia_equal_high[1,],
                 alleles_cap_q_acerifolia_prop_high[1,])
acer_h$alleles_capt <- acer_capt_h
acer_h$intensity <- 'high'
acer_capt_l <- c(alleles_cap_q_acerifolia_equal_low[1,],
                 alleles_cap_q_acerifolia_prop_low[1,])
acer_l$alleles_capt <- acer_capt_l
acer_l$intensity <- 'low'
acer <- rbind(acer_l, acer_h)
acer$species <- 'q_acer'

# QUERCUS ENGELMANNII
load('data/oak_bottlenecks/bottleneck1/alleles_capt_q_engel.RData')
load('data/oak_bottlenecks/bottleneck1/combined_results_q_engelmannii.RData')
engel_h <- combined_q_engelmannii_high
engel_l <- combined_q_engelmannii_low
engel_capt_h <- c(alleles_cap_q_engel_equal_high[1,],
                  alleles_cap_q_engel_prop_high[1,])
engel_h$alleles_capt <- engel_capt_h
engel_h$intensity <- 'high'
engel_capt_l <- c(alleles_cap_q_engel_equal_low[1,],
                  alleles_cap_q_engel_prop_low[1,])
engel_l$alleles_capt <- engel_capt_l
engel_l$intensity <- 'low'
engel <- rbind(engel_l, engel_h)
engel$species <- 'q_engel'

# QUERCUS OGLETHORPENSIS
load('data/oak_bottlenecks/bottleneck1/alleles_capt_q_ogle.RData')
load('data/oak_bottlenecks/bottleneck1/combined_results_q_oglethorpensis.RData')
ogle_h <- combined_q_oglethorpensis_high
ogle_l <- combined_q_oglethorpensis_low
ogle_capt_h <- c(alleles_cap_q_ogle_equal_high[1,],
                 alleles_cap_q_ogle_prop_high[1,])
ogle_h$alleles_capt <- ogle_capt_h
ogle_h$intensity <- 'high'
ogle_capt_l <- c(alleles_cap_q_ogle_equal_low[1,],
                 alleles_cap_q_ogle_prop_low[1,])
ogle_l$alleles_capt <- ogle_capt_l
ogle_l$intensity <- 'low'
ogle <- rbind(ogle_l, ogle_h)
ogle$species <- 'q_ogle'

# combine species for bottleneck 1
oak1 <- rbind(acer, engel, ogle)
oak1$bottleneck <- 1

# BOTTLENECK 2
# QUERCUS ACERIFOLIA
load('data/oak_bottlenecks/bottleneck2/alleles_capt_q_acer.RData')
load('data/oak_bottlenecks/bottleneck2/combined_results_q_acerifolia.RData')
acer_h <- combined_q_acerifolia_high
acer_l <- combined_q_acerifolia_low
acer_capt_h <- c(alleles_cap_q_acerifolia_equal_high[1,],
                 alleles_cap_q_acerifolia_prop_high[1,])
acer_h$alleles_capt <- acer_capt_h
acer_h$intensity <- 'high'
acer_capt_l <- c(alleles_cap_q_acerifolia_equal_low[1,],
                 alleles_cap_q_acerifolia_prop_low[1,])
acer_l$alleles_capt <- acer_capt_l
acer_l$intensity <- 'low'
acer <- rbind(acer_l, acer_h)
acer$species <- 'q_acer'

# QUERCUS ENGELMANNII
load('data/oak_bottlenecks/bottleneck2/alleles_capt_q_engel.RData')
load('data/oak_bottlenecks/bottleneck2/combined_results_q_engelmannii.RData')
engel_h <- combined_q_engelmannii_high
engel_l <- combined_q_engelmannii_low
engel_capt_h <- c(alleles_cap_q_engel_equal_high[1,],
                  alleles_cap_q_engel_prop_high[1,])
engel_h$alleles_capt <- engel_capt_h
engel_h$intensity <- 'high'
engel_capt_l <- c(alleles_cap_q_engel_equal_low[1,],
                  alleles_cap_q_engel_prop_low[1,])
engel_l$alleles_capt <- engel_capt_l
engel_l$intensity <- 'low'
engel <- rbind(engel_l, engel_h)
engel$species <- 'q_engel'

# QUERCUS OGLETHORPENSIS
load('data/oak_bottlenecks/bottleneck2/alleles_capt_q_ogle.RData')
load('data/oak_bottlenecks/bottleneck2/combined_results_q_oglethorpensis.RData')
ogle_h <- combined_q_oglethorpensis_high
ogle_l <- combined_q_oglethorpensis_low
ogle_capt_h <- c(alleles_cap_q_ogle_equal_high[1,],
                 alleles_cap_q_ogle_prop_high[1,])
ogle_h$alleles_capt <- ogle_capt_h
ogle_h$intensity <- 'high'
ogle_capt_l <- c(alleles_cap_q_ogle_equal_low[1,],
                 alleles_cap_q_ogle_prop_low[1,])
ogle_l$alleles_capt <- ogle_capt_l
ogle_l$intensity <- 'low'
ogle <- rbind(ogle_l, ogle_h)
ogle$species <- 'q_ogle'

# combine species for bottleneck 2
oak2 <- rbind(acer, engel, ogle)
oak2$bottleneck <- 2

# combine bottlenecks into one dataframe, calculate total # alleles and 'failures'
oak <- rbind(oak1,oak2)
oak$tot_alleles <- oak$alleles_capt / oak$prop_all
oak$fail <- oak$tot_alleles - oak$alleles_capt
cols <- c('strategy','intensity','species', 'bottleneck')
oak[,cols] <- lapply(oak[,cols], factor)
saveRDS(oak, 'data/oak_bottlenecks/oak_bottleneck_df_for_GLM.RDS')
oak <- readRDS('data/oak_bottlenecks/oak_bottleneck_df_for_GLM.RDS')

# run GLM
y <- as.matrix(oak[,c('alleles_capt','fail')])

mod_oak <- glm(y ~ species + strategy + intensity + bottleneck +
                 species:strategy + 
                 species:intensity + 
                 # species:strategy:intensity +
                 species:strategy:intensity:bottleneck,
               family = binomial(link = 'logit'), data = oak)
anova(mod_oak, test = 'Chisq')
# doesn't add meaning to include 3-way interaction, so remove it
summary(mod_oak)

# Use lsmeans to construct pairwise contrasts
lsmeans(mod_oak, pairwise ~ strategy|species)
# for all species, proportional sampling works better, especially for engelmannii

lsmeans(mod_oak, pairwise ~ intensity|species)
# for all species, high intensity works better

lsm <- lsmeans(mod_oak, pairwise ~ strategy|intensity|bottleneck|species)
# proportional sampling works better, regardless of sampling intensity, for all species
# somewhat more effective when coupled with high intensity sampling for oglethorpensis

# summarize contrast results for supplemental file
lsm_df <- summary(lsm$contrasts)
lsm_df <- lsm_df %>% select(species, bottleneck, intensity, estimate, SE, z.ratio, p.value)
lsm_df[,4:7] <- signif(lsm_df[,4:7], 4)

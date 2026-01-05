## Code associated with:
# Designing monitoring to determine efficacy of control efforts for small populations of invasive species 
# Katherine M. McClure*, Richard J. Camp*, Kristina M. McIntire^, Lucas Berio Fortini*, Dennis A. LaPointe*, Carter T. Atkinson*, Helen R. Sofaer* 

# *U.S. Geological Survey, Pacific Island Ecosystems Research Center, Hawaiʻi National Park, HI 96718, USA 
# ^Hawaiʻi Cooperative Studies Unit, University of Hawaiʻi at Hilo, Hilo, HI 96720, USA 

##  USGS disclaimer:
# This software has been approved for release by the U.S. Geological Survey (USGS). 
# Although the software has been subjected to rigorous review, 
# the USGS reserves the right to update the software as needed pursuant to further analysis and review. 
# No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality 
# of the software and related material nor shall the fact of release constitute any such warranty. 
# Furthermore, the software is released on condition that neither the USGS nor the U.S. Government 
# shall be held liable for any damages resulting from its authorized or unauthorized use.

## Script 1:
# Power analysis for estimating mosquito control efficacy - data simulation and model fitting
# This version includes models of abundance and logistic regression to analyze occurrence
# Includes realistic abund-occur relationship based on field data

# Code by Helen Sofaer
# 24 April 2025; ran in R v 4.5.0

library(MASS)
library(tidyverse)
library(glmmTMB)
library(boot)
library(DHARMa)
library(ggeffects)
library(broom.mixed)
library(cowplot)

##  Estimate parameters linking prevalence to abundance based on mosq data: -------------------------
# biocomplexity and kipahulu data summarized monthly (from Kat)
mosqData <- read_csv("C:/Users/hsofaer/OneDrive - DOI/ResearchProjects/Mosquito/PowerAnalysis/EcoApps/mosquito counts summary_Biocomplexity and Kipahulu.csv")
glimpse(mosqData)

mosqData <- mosqData %>%
  # drop observations that are all zero:
  filter(MeanObserved > 0) %>%
  mutate(logMeanObserved = log(MeanObserved),
         logMeanNonZero = log(MeanNonZero))

# estimate an intercept and slope, based on relationship between mean abundance where present and prevalence
# logit link is default for beta family
mNon0 <- glmmTMB(PropZero ~ logMeanNonZero,
                 data = mosqData,
                 family = beta_family())

# check fit:
simNon0 <- simulateResiduals(mNon0)
plot(simNon0)

predNon0 <- ggpredict(mNon0,
                      terms = "logMeanNonZero [0:3, by = .1]")
plot(predNon0)

p.estPvsA <- predNon0 %>%
  as_tibble() %>%
  mutate(MeanNonZero = exp(x)) %>%
  rename(PropZero = predicted) %>%
  ggplot() +
  geom_point(data = mosqData,
             aes(MeanNonZero, PropZero),
             color = "#5a2328") +
  geom_line(aes(MeanNonZero, PropZero),
            linewidth = 1) +
  scale_x_continuous(limits = c(1, round(max(mosqData$MeanNonZero))),
                     breaks = seq(1, round(max(mosqData$MeanNonZero)), by = 1)) +
  theme_bw() +
  ylab("Proportion of traps with zero captures") +
  xlab("Mean count in traps with captures") +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(size = 12))

###  Simulate based on these parameters: ------------------------------------------------
zi.int <- tidy(mNon0) %>%
  filter(term == "(Intercept)") %>%
  select(estimate) %>%
  as.numeric()
zi.slope <- tidy(mNon0) %>%
  filter(term == "logMeanNonZero") %>%
  select(estimate) %>%
  as.numeric()


# parameter setup in a tibble:
nSims <- 50
simData <- expand_grid(nbin_mu = c(2, 3, 5, 7, 10),
                       nbin_theta = 20,
                       nTraps = c(10, 25, 45),
                       nNights = c(4, 7, 10),
                       controlEfficacy = c(.6, .9),
                       detection = c(.5, .1),
                       sim_rep = seq(1, nSims)) %>%
  # add zero-inf prob:
  mutate(zeroinf_pi_pre = inv.logit(zi.int + zi.slope*log(nbin_mu)),
         zeroinf_pi_post = inv.logit(zi.int + zi.slope*log(nbin_mu * (1 - controlEfficacy))))


# functions for simulating ZINB count and imperfect obs-------------------------
sim_pre <- function(nbin_mu, nbin_theta, zeroinf_pi_pre, nTraps, nNights,
                    detection) {
  
  simData <- tibble(Trap = rep(1:nTraps, each = nNights),
                    Night = rep(1:nNights, times = nTraps),
                    MosqNum = rep(rnegbin(nTraps, 
                                          nbin_mu, 
                                          nbin_theta) * 
                                    rbinom(nTraps, 
                                           prob = 1 - zeroinf_pi_pre, 
                                           size = 1),
                                  each = nNights),
                    ObsCount = rbinom(n = nTraps*nNights,
                                      prob = detection,
                                      size = MosqNum),
                    ObsOcc = if_else(ObsCount > 0, 1, 0),
                    Treatment = "1pre")
}

sim_post <- function(nbin_mu, nbin_theta, zeroinf_pi_post, nTraps, nNights,
                     detection, controlEfficacy) {
  
  simData <- tibble(Trap = rep(1:nTraps, each = nNights),
                    Night = rep(1:nNights, times = nTraps),
                    MosqNum = rep(rnegbin(nTraps, 
                                          # control acts on the negbin mean
                                          nbin_mu * (1 - controlEfficacy), 
                                          nbin_theta) * 
                                    rbinom(nTraps, 
                                           prob = 1 - zeroinf_pi_post, 
                                           size = 1),
                                  each = nNights),
                    ObsCount = rbinom(n = nTraps*nNights,
                                      prob = detection,
                                      size = MosqNum),
                    ObsOcc = if_else(ObsCount > 0, 1, 0),
                    Treatment = "2post")
}



##  Simulate data: -------------------------------------------------------------

# Code was run 6 times with nSims == 50, using seeds 1234, 2345, ..., 6789
set.seed(6789)

simData <- simData %>%
  mutate(preControl = pmap(list(nbin_mu, 
                                nbin_theta,
                                zeroinf_pi_pre,
                                nTraps,
                                nNights,
                                detection), 
                           .f = sim_pre),
         postControl = pmap(list(nbin_mu, 
                                 nbin_theta,
                                 zeroinf_pi_post,
                                 nTraps,
                                 nNights,
                                 detection,
                                 controlEfficacy), 
                            .f = sim_post),
         # combine pre and post into one tibble
         PreAndPost = map2(preControl, postControl, bind_rows),
         # to confirm that initial obs'd counts approximate reasonable values:
         pre_meanObs = map_dbl(preControl, ~ mean(.x$ObsCount)),
         pre_meanObsNon0 = map_dbl(preControl, ~ mean(.x$ObsCount[.x$ObsCount != 0])),
         pre_prop0Obs = map_dbl(preControl, ~ sum(.x$ObsCount == 0)/nrow(.x)),
         pre_totalObs = map_dbl(preControl, ~ sum(.x$ObsCount)),
         post_meanObs = map_dbl(postControl, ~ mean(.x$ObsCount)),
         post_prop0Obs = map_dbl(postControl, ~ sum(.x$ObsCount == 0)/nrow(.x)),
         post_meanObsNon0 = map_dbl(postControl, ~ mean(.x$ObsCount[.x$ObsCount != 0])),
         post_totalObs = map_dbl(postControl, ~ sum(.x$ObsCount)),
         # also pull out what 'truth' was for occupancy in each grid:
         pre_prop0MosqNum = map_dbl(preControl, ~ sum(.x$MosqNum == 0)/nrow(.x)),
         post_prop0MosqNum = map_dbl(postControl, ~ sum(.x$MosqNum == 0)/nrow(.x))
  ) %>%
  dplyr::select(-preControl, -postControl)

##  Plot simulated and obs'd prevalence vs abundance relationships: ------------

p.simObs <- simData %>%
  select(nbin_mu, zeroinf_pi_pre, zeroinf_pi_post, sim_rep, 
         nTraps, nNights, detection, controlEfficacy,
         contains("meanObs"), contains("prop0Obs"),
         -contains("meanObsNon0")) %>%
  pivot_longer(c(contains("meanObs"), contains("prop0Obs")), 
               names_sep = "_",
               names_to = c("Dataset", "variable")) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(Dataset = case_when(
    Dataset == "pre" ~ "Simulated pre-treatment",
    Dataset == "post" ~ "Simulated post-treatment")) %>%
  bind_rows(mosqData %>%
              select(meanObs = MeanObserved, 
                     prop0Obs = PropZero) %>%
              mutate(Dataset = "Observed")) %>%
  mutate(Dataset = factor(Dataset,
                             levels = c("Observed",
                                        "Simulated pre-treatment",
                                        "Simulated post-treatment"))) %>%
  ggplot() +
  geom_point(aes(meanObs, prop0Obs,
                 color = Dataset)) +
  scale_color_manual(values = c("#5a2328",
                                "#d38c36",
                                "#e4bb25")) +
  theme_bw() +
  ylab("Proportion of traps with zero captures") +
  xlab("Mean count across all traps") +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(.6, .8))

png("C:/Users/hsofaer/OneDrive - DOI/ResearchProjects/Mosquito/PowerAnalysis/EcoApps/output/AbundPrev_6789.png",
    width = 8, height = 4, units = "in", res = 400)
plot_grid(p.estPvsA, p.simObs,
          labels = c("a", "b"))
dev.off()

##  Modeling functions: --------------------------------------------------------

# ZINB abundance model:
trt_model <- function(df) {
  
  m1 <- glmmTMB(ObsCount ~ Treatment,
                data = df,
                family = nbinom2,
                zi = ~ 1)
}


# due to convergence warnings, check if I get the same estimates (see glmmTMB troubleshooting vignette):
trt_model_bfgs <- function(df) {
  
  m1 <- glmmTMB(ObsCount ~ Treatment,
                data = df,
                family = nbinom2,
                zi = ~ 1,
                control = glmmTMBControl(optimizer = optim, 
                                         optArgs = list(method="BFGS"))
  )
}

# logistic regression: 
logistic_model <- function(df) {
  
  m1 <- glmmTMB(ObsOcc ~ Treatment,
                data = df,
                family = binomial)
}

logistic_model_bfgs <- function(df) {
  
  m1 <- glmmTMB(ObsOcc ~ Treatment,
                data = df,
                family = binomial,
                control = glmmTMBControl(optimizer = optim, 
                                         optArgs = list(method="BFGS"))
  )
}

##  Filter to the simulations with observed 1s in initial period ---------------

simData <- simData %>%
  filter(pre_prop0Obs != 1)


##  Model outcomes: -----------------------------------------------------------

simData <- simData %>%
  mutate(
    # estimate a ZINB model for each row:
    modZINB = map(PreAndPost, trt_model),
    
    # re-fit for convergence check:
    modZINB_BFGS = map(PreAndPost, trt_model_bfgs),
    
    # logistic model:
    modLogistic = map(PreAndPost, logistic_model),
    
    modLogistic_BFGS = map(PreAndPost, logistic_model_bfgs)
  ) 



##  Extract estimates: ---------------------------------------------------------
simData <- simData %>%
  mutate(
    ##  ZINB model: 
    # list column with estimates
    abund_tidyEst = map(modZINB, tidy),
    # pull out key pieces:
    abund_estTrtEff = map_dbl(abund_tidyEst, 
                              ~ filter(., term == "Treatment2post") %>%
                                select(., estimate) %>% 
                                as.numeric(.)),
    abund_estControlEfficacy = 1 - exp(abund_estTrtEff),
    # convergence checks:
    abund_estTrtP = map_dbl(abund_tidyEst, 
                            ~ filter(., term == "Treatment2post") %>%
                              select(., p.value) %>% 
                              as.numeric(.)),
    abund_bigBeta = map_dbl(abund_tidyEst,
                            ~ sum(abs(.$estimate) > 10)),
    abund_pdHess = map_lgl(modZINB, ~ .$sdr$pdHess),
    abund_diagnoseExplain = map(modZINB, \(x) capture.output(diagnose(x))),
    abund_diagnose1ok = map_dbl(modZINB, \(x) diagnose(x) %>%
                                  sum(.)),
    
    # ZINB model with different optimization method:
    abundBFGS_tidyEst = map(modZINB_BFGS, tidy),
    abundBFGS_estTrtEff = map_dbl(abundBFGS_tidyEst, 
                                  ~ filter(., term == "Treatment2post") %>%
                                    select(., estimate) %>% 
                                    as.numeric(.)),
    
    # Summarize my convergence checks:
    abund_Converge1ok = case_when(
      is.na(abund_estTrtP) ~ 0,
      isFALSE(abund_pdHess) ~ 0,
      abund_bigBeta > 0 ~ 0,
      abs(abund_estTrtEff - abundBFGS_estTrtEff) > .05 ~ 0,
      TRUE ~ 1),
    
    ##  Logistic model:
    logistic_tidyEst = map(modLogistic, tidy),
    logistic_estTrtEff = map_dbl(logistic_tidyEst,
                                 ~ filter(., term == "Treatment2post") %>%
                                   select(., estimate) %>%
                                   as.numeric(.)),
    logistic_estTrtP = map_dbl(logistic_tidyEst, 
                               ~ filter(., term == "Treatment2post") %>%
                                 select(., p.value) %>% 
                                 as.numeric(.)),
    logistic_bigBeta = map_dbl(logistic_tidyEst,
                               ~ sum(abs(.$estimate) > 10)),
    logistic_pdHess = map_lgl(modLogistic, ~ .$sdr$pdHess),
    logistic_diagnoseExplain = map(modLogistic, \(x) capture.output(diagnose(x))),
    logistic_diagnose1ok = map_dbl(modLogistic, \(x) diagnose(x) %>%
                                     sum(.)),
    
    #  Logistic model with different optimization method:
    logisticBFGS_tidyEst = map(modLogistic_BFGS, tidy),
    logisticBFGS_estTrtEff = map_dbl(logisticBFGS_tidyEst, 
                                     ~ filter(., term == "Treatment2post") %>%
                                       select(., estimate) %>% 
                                       as.numeric(.)),
    
    #  Summarize my convergence checks:
    logistic_Converge1ok = case_when(
      is.na(logistic_estTrtP) ~ 0,
      isFALSE(logistic_pdHess) ~ 0,
      logistic_bigBeta > 0 ~ 0,
      abs(logistic_estTrtEff - logisticBFGS_estTrtEff) > .05 ~ 0,
      TRUE ~ 1),
  ) %>%
  # drop extra models to keep object smaller:
  select(-modZINB_BFGS, -modLogistic_BFGS)



write_rds(simData, file = "C:/Users/hsofaer/Documents/MosqPower/linkedZeros/simAbundLogistic_LinkedZeros_seed6789.rds")

simData %>%
  select(-starts_with("mod")) %>%
  write_rds("C:/Users/hsofaer/Documents/MosqPower/linkedZeros/simAbundLogistic_LinkedZeros_seed6789noMod.rds")

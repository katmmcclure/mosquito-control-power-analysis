# Designing monitoring to determine efficacy of control efforts for small populations of invasive species
# Power analysis for estimating mosquito control efficacy - data simulation and model fitting

library(MASS)
library(tidyverse)
library(glmmTMB)
library(broom.mixed)
library(unmarked)


# parameter setup in a tibble:
nSims <- 50
simData <- expand_grid(nbin_mu = c(2, 3, 5, 7, 10),
                       nbin_theta = 20,
                       zeroinf_pi = c(.6),
                       # updated to match Rick's lower values:
                       nTraps = c(10, 25, 45),
                       nNights = c(4, 7, 10),
                       controlEfficacy = c(.6, .9),
                       detection = c(.5, .1),
                       sim_rep = seq(1, nSims)) 

# functions for simulating ZINB count and imperfect obs-------------------------
sim_pre <- function(nbin_mu, nbin_theta, zeroinf_pi, nTraps, nNights,
                    detection) {
  
  simData <- tibble(Trap = rep(1:nTraps, each = nNights),
                    Night = rep(1:nNights, times = nTraps),
                    MosqNum = rep(rnegbin(nTraps, 
                                          nbin_mu, 
                                          nbin_theta) * 
                                    rbinom(nTraps, 
                                           prob = 1 - zeroinf_pi, 
                                           size = 1),
                                  each = nNights),
                    ObsCount = rbinom(n = nTraps*nNights,
                                      prob = detection,
                                      size = MosqNum),
                    ObsOcc = if_else(ObsCount > 0, 1, 0),
                    Treatment = "1pre")
}

sim_post <- function(nbin_mu, nbin_theta, zeroinf_pi, nTraps, nNights,
                     detection, controlEfficacy) {
  
  simData <- tibble(Trap = rep(1:nTraps, each = nNights),
                    Night = rep(1:nNights, times = nTraps),
                    MosqNum = rep(rnegbin(nTraps, 
                                          # control acts on the negbin mean
                                          nbin_mu * (1 - controlEfficacy), 
                                          nbin_theta) * 
                                    rbinom(nTraps, 
                                           prob = 1 - zeroinf_pi, 
                                           size = 1),
                                  each = nNights),
                    ObsCount = rbinom(n = nTraps*nNights,
                                      prob = detection,
                                      size = MosqNum),
                    ObsOcc = if_else(ObsCount > 0, 1, 0),
                    Treatment = "2post")
}



##  Simulate data: -------------------------------------------------------------

# #  Set up for parallel processing: https://furrr.futureverse.org/
# plan(multisession, workers = 16) # this was slower than sequential

set.seed(4321)

simData <- simData %>%
  mutate(preControl = pmap(list(nbin_mu, 
                                nbin_theta,
                                zeroinf_pi,
                                nTraps,
                                nNights,
                                detection), 
                           .f = sim_pre),
         postControl = pmap(list(nbin_mu, 
                                 nbin_theta,
                                 zeroinf_pi,
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
         post_prop0MosqNum = map_dbl(postControl, ~ sum(.x$MosqNum == 0)/nrow(.x)),
         # wide data for occupancy modeling:
         PrePostWide = map(PreAndPost, ~ .x %>% 
                             dplyr::select(Trap, Treatment, Night, ObsOcc) %>%
                             pivot_wider(id_cols = c("Trap", "Treatment"),
                                         names_from = Night,
                                         names_prefix = "Night",
                                         values_from = ObsOcc) %>%
                             mutate(Treatment = factor(Treatment))),
         PrePostOcc = map(PrePostWide, ~ unmarkedFrameOccu(
           y = .x %>%
             dplyr::select(starts_with("Night")),
           siteCovs = .x %>%
             dplyr::select(Treatment)))
  ) %>%
  dplyr::select(-preControl, -postControl, -PrePostWide)



##  Modeling functions: --------------------------------------------------------

# ZINB abundance model:
trt_model <- function(df) {
  
  m1 <- glmmTMB(ObsCount ~ Treatment,
                data = df,
                family = nbinom2,
                zi = ~ 1)
}


# due to convergence warnings, check if I get the same estimates (see troubleshooting vignette):
trt_model_bfgs <- function(df) {
  
  m1 <- glmmTMB(ObsCount ~ Treatment,
                data = df,
                family = nbinom2,
                zi = ~ 1,
                control = glmmTMBControl(optimizer = optim, 
                                         optArgs = list(method="BFGS"))
  )
}


# occupancy model: unmarked formula is p first, psi second
occ_model <- function(unmarkedObj) {
  
  mOcc <- occu(~ 1 ~ Treatment,
               data = unmarkedObj)
  
}

# function for bootstrap occ model convergence check from Fisk and Chandler:
chisq <- function(fm) {
  umf <- fm@data
  y <- umf@y
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop = FALSE]
  fv <- fitted(fm, na.rm = TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm = TRUE)
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
    
    # occupancy model:
    modOcc = map(PrePostOcc, occ_model)
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
    abund_estTrtP = map_dbl(abund_tidyEst, 
                            ~ filter(., term == "Treatment2post") %>%
                              select(., p.value) %>% 
                              as.numeric(.)),
    abund_bigBeta = map_dbl(abund_tidyEst,
                            ~ sum(abs(.$estimate) > 10)),
    abund_pdHess = map_lgl(modZINB, ~ .$sdr$pdHess),
    ##  ZINB model with different optimization method:
    abundBFGS_tidyEst = map(modZINB_BFGS, tidy),
    abundBFGS_estTrtEff = map_dbl(abundBFGS_tidyEst, 
                                  ~ filter(., term == "Treatment2post") %>%
                                    select(., estimate) %>% 
                                    as.numeric(.)),
    ##  Occupancy model:
    occ_estTrtP = map_dbl(modOcc, ~ summary(.)$state[[2, 4]]),
    occ_pred = map(modOcc, ~ predict(., 
                                     newdata = data.frame(Treatment = c("1pre", "2post")), 
                                     type = "state", 
                                     appendData = TRUE)),
    occ_Pre_EstState = map_dbl(occ_pred, . %>% 
                                 filter(Treatment == "1pre") %>%
                                 select(Predicted) %>%
                                 as.numeric()),
    occ_Pre_EstLCI = map_dbl(occ_pred, . %>% 
                               filter(Treatment == "1pre") %>%
                               select(lower) %>%
                               as.numeric()),
    occ_Pre_EstUCI = map_dbl(occ_pred, . %>% 
                               filter(Treatment == "1pre") %>%
                               select(upper) %>%
                               as.numeric()),
    occ_Post_EstState = map_dbl(occ_pred, . %>% 
                                  filter(Treatment == "2post") %>%
                                  select(Predicted) %>%
                                  as.numeric()),
    occ_Post_EstLCI = map_dbl(occ_pred, . %>% 
                                filter(Treatment == "2post") %>%
                                select(lower) %>%
                                as.numeric()),
    occ_Post_EstUCI = map_dbl(occ_pred, . %>% 
                                filter(Treatment == "2post") %>%
                                select(upper) %>%
                                as.numeric()),
    occ_EstDet = map_dbl(modOcc, ~ backTransform(.x, 'det')@estimate),
    
    ##  Add convergence checks:
    
    abund_failConverge = case_when(
      is.na(abund_estTrtP) ~ 1,
      isFALSE(abund_pdHess) ~ 1,
      abund_bigBeta > 0 ~ 1,
      abs(abund_estTrtEff - abundBFGS_estTrtEff) > .05 ~ 1,
      TRUE ~ 0),
    occ_FailCheck = map_dbl(modOcc, ~ .x@opt$convergence),
    # bootstraping check:
    occ_convergeBoot = map(modOcc, ~ parboot(.x,
                                             statistic = chisq,
                                             nsim = 100,
                                             parallel = FALSE)),
    occ_convergeBootP = map_dbl(occ_convergeBoot, 
                                ~ sum(.x@t.star > .x@t0)/nrow(.x@t.star)),
    occ_failConverge = case_when(
      occ_FailCheck == 1 ~ 1,
      occ_convergeBootP < 0.05 ~ 1,
      TRUE ~ 0)
    
  ) %>%
  # drop extra model to keep it smaller:
  select(-modZINB_BFGS)


write_rds(simData, file = "./simAbundOccur_seed4321.rds")

simData %>%
  select(-starts_with("mod")) %>%
  write_rds("./simAbundOccur_sim4321noMod.rds")


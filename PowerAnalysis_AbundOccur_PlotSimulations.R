# Designing monitoring to determine efficacy of control efforts for small populations of invasive species
# Figures 1-2 & S Figures 1-5

# Power analysis for estimating mosquito control efficacy - plotting and summarizing results
# Note this script uses outputs produced by PowerAnalysis_AbundOccur_RunSimulations.R

library(tidyverse)
library(MASS)

theme_set(theme_bw() +
            theme(axis.text = element_text(color = "black")))
setwd("./Mosquito/PowerAnalysis/openPop")

simOut <- list.files(pattern = "noMod.rds$", full.names = TRUE)

sims <- simOut %>%
  map_dfr(read_rds)

nSims <- length(simOut) * 50


##  What was control efficacy (prop decline) incl zeros: -----------------------

controlEff <- sims %>%
  mutate(obsControlEff = 1 - post_meanObs/pre_meanObs,
         pre_meanMosqNum = map_dbl(PreAndPost,
                                   ~ mean(.x$MosqNum[.x$Treatment == "1pre"])),
         post_meanMosqNum = map_dbl(PreAndPost,
                                    ~ mean(.x$MosqNum[.x$Treatment == "2post"])),
         RatioMosqNum = post_meanMosqNum/pre_meanMosqNum,
         mosqNumControlEff = 1 - post_meanMosqNum/pre_meanMosqNum) 


controlEff %>%
  ggplot() +
  geom_violin(aes(factor(controlEfficacy), RatioMosqNum)) +
  #  ylim(0, 5) +
  facet_grid(nTraps ~ nNights)

controlEff %>%
  ggplot() +
  geom_violin(aes(factor(controlEfficacy), obsControlEff)) +
  #  ylim(0, 5) +
  facet_grid(nTraps ~ nNights)

sims %>%
  mutate(propObsPopChange = (pre_totalObs - post_totalObs)/pre_totalObs) %>%
  ggplot() +
  geom_violin(aes(factor(controlEfficacy), propObsPopChange))+
  facet_grid(~ nbin_theta) +
  geom_hline(yintercept = .6, linetype = "dashed") +
  geom_hline(yintercept = .9, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  ylim(c(-2, 1.2)) +
  ggtitle("Open population: resampling from rnegbin; cut ylim")

controlEff %>%
  ggplot() +
  geom_hex(aes(mosqNumControlEff, obsControlEff)) +
  scale_fill_gradient(trans = "log10") +
  ylim(c(-2, 1.2)) +
  xlim(c(-2, 1.2)) +  facet_wrap(~ controlEfficacy)

controlEff %>%
  ggplot() +
  geom_hex(aes(mosqNumControlEff, abund_estControlEfficacy)) +
  scale_fill_gradient(trans = "log10") +
  ylim(c(-2, 1.2)) +
  xlim(c(-2, 1.2)) +
  facet_wrap(~ controlEfficacy)

controlEff %>%
  mutate(controlEfficacy = factor(controlEfficacy, 
                                  levels = levels(factor(controlEfficacy)),
                                  labels = paste0("Control efficacy of ", levels(factor(controlEfficacy)))),
         detection = factor(detection, 
                            levels = levels(factor(detection)),
                            labels = paste0("Detection probability of ", levels(factor(detection))))
  ) %>%
  # just significant ones:
  filter(abund_estTrtP < .05,
         abund_failConverge == 0) %>%
  ggplot() +
  geom_hex(aes(obsControlEff, abund_estControlEfficacy)) +
  scale_fill_gradient(trans = "log10") +
  ylim(c(-5, 1.2)) +
  xlim(c(-5, 1.2)) +
  facet_grid(controlEfficacy ~ detection)

##  Calculate power: -----------------------------------------------------------

sims %>%
  filter(post_prop0Obs == 1,
         abund_failConverge == 0) %>%
  nrow()

sims %>%
  filter(post_prop0Obs == 1) %>%
  nrow()

abund_simPower <- sims %>%
  mutate(nTraps = factor(nTraps, 
                         levels = levels(factor(nTraps)),
                         labels = paste(levels(factor(nTraps)), "traps")),
         nNights = factor(nNights, 
                          levels = levels(factor(nNights)),
                          labels = paste(levels(factor(nNights)), "nights"))
  ) %>%
  mutate(sig = ifelse(abund_estTrtP < .05 & abund_failConverge == 0, 1, 0)) %>%
  group_by(nbin_mu, zeroinf_pi, nTraps, nNights, controlEfficacy, detection) %>%
  summarize(n_failConverge = sum(abund_failConverge == 1 & post_prop0Obs < 1),
            n_all0Post = sum(post_prop0Obs == 1),
            n_sig = sum(sig),
            propSig = sum(sig)/nSims,
            mean_pre_meanObs = mean(pre_meanObs),
            n = n()) %>%
  mutate(n_all0Pre = nSims - n,
         n_nonSig = nSims - n_all0Pre - n_all0Post - n_failConverge - n_sig)
head(abund_simPower)


## Stacked bar plots:  ----------------------------------------------------

# make data long:
abund_simPower.long <- abund_simPower %>%
  select(nbin_mu, zeroinf_pi, nTraps, nNights, controlEfficacy, detection,
         starts_with("n_"), mean_pre_meanObs) %>%
  mutate(mean_pre_meanObs = round(mean_pre_meanObs, 2)) %>%
  pivot_longer(starts_with("n_"), names_to = "outcome", values_to = "numSims") %>%
  # manually calculate proportion since it's easier to add line:
  mutate(propSims = numSims / nSims,
         outcome = factor(outcome,
                          levels = c("n_all0Pre",
                                     "n_all0Post",
                                     "n_failConverge",
                                     "n_nonSig",
                                     "n_sig"),
                          labels = c("Zero captures pre-control",
                                     "Zero captures post-control",
                                     "Model failed to converge",
                                     "Non-significant treatment effect",
                                     "Significant treatment effect"),
                          ordered = TRUE)) 

# ---------------------------------------------
# S Figure 5 ----------------------------------
# ---------------------------------------------
p.p5c9 <- abund_simPower.long %>%  
  filter(detection == .5,
         controlEfficacy == .9) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .05,
           position = "stack") +
  scale_fill_viridis_d(name = "Outcome when\np = 0.5\nControl efficacy = 0.9") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("power_stackBar_p5c9.png",
       width = 9, height = 6, dpi = 450, units = "in")  

# ---------------------------------------------
# S Figure 3 ----------------------------------
# ---------------------------------------------
p.p1c9 <- abund_simPower.long %>%  
  filter(detection == .1,
         controlEfficacy == .9) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .03,
           position = "stack") +
  scale_fill_viridis_d(name = "Outcome when\np = 0.1\nControl efficacy = 0.9") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("power_stackBar_p1c9.png",
       width = 9, height = 6, dpi = 450, units = "in") 

# ---------------------------------------------
# S Figure 2 ----------------------------------
# ---------------------------------------------
p.p1c6 <- abund_simPower.long %>%  
  filter(detection == .1,
         controlEfficacy == .6) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .02,
           position = "stack") +
  scale_fill_viridis_d(name = "Outcome when\np = 0.1\nControl efficacy = 0.6") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("power_stackBar_p1c6.png",
       width = 9, height = 6, dpi = 450, units = "in") 

# ---------------------------------------------
# S Figure 4 ----------------------------------
# ---------------------------------------------
p.p5c6 <- abund_simPower.long %>%  
  filter(detection == .5,
         controlEfficacy == .6) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .05,
           position = "stack") +
  scale_fill_viridis_d(name = "Outcome when\np = 0.5\nControl efficacy = 0.6") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("power_stackBar_p5c6.png",
       width = 9, height = 6, dpi = 450, units = "in")

# -------------------------------------------
# Figure 2 ----------------------------------
# -------------------------------------------
# prob = .1, efficacy = .9 on scale that best fits just that one:
p.p1c9.focus <- p.p1c9 +
  scale_x_continuous(limits = c(0, .45))
ggsave("power_stackBar_p1c9_focus.png",
       width = 9, height = 6, dpi = 450, units = "in")

##  comparison between scenarios for 7 nights 25 traps:
abund_simPower %>%
  mutate(controlEfficacy = factor(controlEfficacy, 
                                  levels = levels(factor(controlEfficacy)),
                                  labels = paste0("Control efficacy of ", levels(factor(controlEfficacy)))),
         detection = factor(detection, 
                            levels = levels(factor(detection)),
                            labels = paste0("Detection probability of ", levels(factor(detection))))
  ) %>%
  filter(nbin_mu <= 10) %>%
  select(nbin_mu, zeroinf_pi, nTraps, nNights, controlEfficacy, detection,
         starts_with("n_"), mean_pre_meanObs) %>%
  mutate(mean_pre_meanObs = round(mean_pre_meanObs, 2)) %>%
  pivot_longer(starts_with("n_"), names_to = "outcome", values_to = "numSims") %>%
  filter(nTraps == "25 traps",
         nNights == "7 nights") %>%
  ggplot(aes(mean_pre_meanObs, numSims)) +
  # geom_point(aes(mean_pre_meanObs, numSims)) +
  geom_bar(aes(fill = outcome),
           stat = "identity", width = .04) +
  scale_fill_viridis_d() +
  # TO DO: change to proportions
  #  geom_hline(yintercept = 80) +
  facet_grid(detection ~ controlEfficacy) +
  ggtitle("25 traps for 7 nights") +
  geom_line(data = . %>%
              filter(outcome == "n_sig")) +
  #  scale_x_continuous(limits = c(0, 2.5)) +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  #  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  ylab("Number of simulations")

##  Reason for failure across model runs: -------------------------------

#   For different control efficacies:
abund_simPower.long %>%
  pivot_wider(names_from = controlEfficacy, names_prefix = "controlEff_",
              values_from = propSims) %>%
  ggplot(aes(controlEff_0.6, controlEff_0.9)) +
  geom_point() +
  facet_wrap(~ outcome, scales = "free")

p.outcome <- abund_simPower.long %>%
  mutate(controlEfficacy = factor(controlEfficacy, 
                                  levels = levels(factor(controlEfficacy)),
                                  labels = paste0("Control efficacy of ", levels(factor(controlEfficacy)))),
         detection = factor(detection, 
                            levels = levels(factor(detection)),
                            labels = paste0("Detection probability of ", levels(factor(detection))))
  ) %>%
  # # recreate numeric trap nights:
  # mutate(nTraps.num = as.numeric(str_extract(nTraps, "\\d+")),
  #        nNights.num = as.numeric(str_extract(nNights, "\\d+")),
  #        nTrapNights = nTraps.num*nNights.num) %>%
  ggplot(aes(propSims, fill = nTraps)) +
  facet_grid(detection*controlEfficacy ~ outcome, scales = "free_x") +
  geom_histogram() +
  scale_fill_brewer(name = "Number of traps", palette = "Dark2") +
  xlab("Proportion of simulations with this outcome") +
  ylab("Number of simulation parameter combinations") +
  theme(strip.text = element_text(size = 8))
ggsave("outcomeByDetectControl.png",
       width = 12, height = 9, dpi = 450, units = "in")

# -------------------------------------------
# Figure 1 ----------------------------------
# -------------------------------------------
# panels showing significant outcomes:
p.sig <- abund_simPower.long %>%
  filter(outcome == "Significant treatment effect") %>%
  mutate(controlEfficacy = factor(controlEfficacy, 
                                  levels = levels(factor(controlEfficacy)),
                                  labels = paste0("Control efficacy of ", levels(factor(controlEfficacy)))),
         detection = factor(detection, 
                            levels = levels(factor(detection)),
                            labels = paste0("Detection probability of ", levels(factor(detection))))
  ) %>%
  ggplot(aes(propSims, fill = nTraps)) +
  facet_grid(detection ~ controlEfficacy, scales = "free_x") +
  geom_histogram() +
  scale_fill_brewer(name = "Number of traps", palette = "Dark2") +
  xlab("Power (proportion of simulations with significant treatment effect)") +
  ylab("Number of simulation parameter combinations") +
  theme(strip.text = element_text(size = 10))
ggsave("sigByDetectControl.png",
       width = 9, height = 6, dpi = 450, units = "in")


abund_simPower.long %>%
  group_by(detection, controlEfficacy, outcome) %>%
  count() 

abund_simPower.long %>%
  filter(nTraps == "25 traps",
         nNights == "7 nights",
         detection == .1,
         controlEfficacy == .9,
         outcome == "Significant treatment effect") 


#  Plot in terms of observed means: --------------------------------------------
sims %>%
  ggplot(aes(factor(nbin_mu), pre_meanObs,
             color = detection)) +
  geom_violin() +
  facet_grid(~ detection)

sims %>%
  ggplot(aes(pre_meanObs)) +
  geom_histogram() +
  facet_grid(controlEfficacy ~ detection)

# abund_simPower %>%
#   ggplot() +
#   geom_line(aes(mean_pre_meanObs, propSig,
#                  color = factor(controlEfficacy))) +
#   facet_grid(nTraps ~ nNights) +
#   geom_hline(yintercept = .8)
# 
# abund_simPower %>%
#   ggplot() +
#   geom_line(aes(mean_pre_meanObs, propSig,
#                 color = factor(controlEfficacy))) +
#   facet_grid(nTraps ~ nNights) +
#   geom_hline(yintercept = .8) +
#   geom_vline(xintercept = 0.48, lty = 2) +
#   ylab("Power\n(proportion of simulations with significant treatment effect)") +
#   xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
#   scale_color_manual(name = "True control efficacy",
#                      values = c("orange3", "navy"))

# p.abundPower <- abund_simPower %>%
#   ggplot() +
#   geom_line(aes(mean_pre_meanObs, propSig,
#                 color = factor(controlEfficacy),
#                 linetype = factor(detection)
#   )) +
#   facet_grid(nTraps ~ nNights) +
#   geom_hline(yintercept = .8) +
#   geom_vline(xintercept = 0.48, lty = 2) +
#   ylab("Power\n(proportion of simulations with significant treatment effect)") +
#   xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
#   scale_color_manual(name = "True control efficacy",
#                      values = c("orange3", "navy")) +
#   scale_linetype_discrete(name = "Detection probability")
# p.abundPower
# ggsave("./openPop/abundPower_openPop_fixpost1001_19Oct2023.png",
#        width = 11, height = 7)

#  look at control efficacy vs obs'd data: -------------------------------------
sims %>%
  ggplot(aes(factor(1 - controlEfficacy), post_meanObs/pre_meanObs)) +
  geom_violin() +
  facet_grid(detection ~ nbin_mu) +
#  ylim(0, 2) +
  geom_hline(yintercept = .1) +
  geom_hline(yintercept = .4, lty = 2)


#  group according to breaks in observed mosq: ---------------------------------
abund_simPower_obsdBreaks <- sims %>%
  filter(abund_failConverge == 0) %>%
  mutate(sig = ifelse(abund_estTrtP < .05, 1, 0),
         meanObs_cut = cut(pre_meanObs, breaks = c(seq(0, 2, by = .25), 3, 4, 8, 12, 25))) %>%
  group_by(meanObs_cut, nTraps, nNights, controlEfficacy) %>%
  summarize(propSig = sum(sig)/n(),
            nSims = n(),
            mean_pre_meanObs = mean(pre_meanObs)) 

abund_simPower_obsdBreaks %>%
  ggplot(aes(meanObs_cut, nSims)) +
  geom_point() +
  facet_grid(nTraps ~ nNights)

abund_simPower_obsdBreaks %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(meanObs_cut, propSig,
                 color = factor(controlEfficacy))) +
  facet_grid(nTraps ~ nNights) +
  geom_hline(yintercept = .8) +
  # geom_vline(xintercept = 0.48, lty = 2) +
  ylab("Power\n(proportion of simulations with significant treatment effect)") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)\n(Grouped by observed bins)") +
  scale_color_manual(name = "True control efficacy",
                     values = c("orange3", "navy")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))


#  Comparing estimated treatment effect to truth: ------------------------------

# show exp of estimated trt effect is 1-controlEfficacy, using one run as an example:
predict(simModel$mod[[6]], 
        newdata = tibble(Treatment = c("1pre", "2post")),
        type = "conditional") # this gives change in mean (not incl zi)
simModel$tidyEst[[6]]
simModel$tidyEst[[6]]$estimate[1] %>% exp(.) # est mean for pre (cond intercept is first param est)

(simModel$tidyEst[[6]]$estimate[1] + simModel$tidyEst[[6]]$estimate[2]) %>% exp(.) # est mean for post

simModel$tidyEst[[6]]$estimate[2] %>% exp(.) # exponentiating Trt effect gives EstPostTrt/EstPreTrt
((simModel$tidyEst[[6]]$estimate[1] + simModel$tidyEst[[6]]$estimate[2]) %>% exp(.)) / (simModel$tidyEst[[6]]$estimate[1] %>% exp(.))

# TrueMean_Post = TrueMean_Pre * (1-controlEfficacy), so ratio Post/Pre = (1-controlEfficacy) = exp(TrtEffect)

all.equal(sims$abund_estControlEfficacy, 1 - exp(sims$abund_estTrtEff))

# maybe make separate plots for each control efficacy or detection level
sims %>%
  filter(abund_failConverge == 0) %>%
  ggplot(aes(factor(nbin_mu), abund_estControlEfficacy,
             fill = factor(interaction(detection, controlEfficacy)))) +
  geom_violin() +
  facet_grid(nTraps ~ nNights) +
  ylim(-1, 1)

##  How tightly was control efficacy estimated, just when trt effect was significant:
sims %>%
  filter(abund_failConverge == 0,
         abund_estTrtP < .05) %>%
  ggplot(aes(factor(controlEfficacy), abund_estControlEfficacy)) +
  geom_violin()


# Justifying choice of parameters ----------------------------------------------------

# Justify parameters: confirm counts approx match distribution of obs'd data: -----
biocomplexity <- tibble(site = c("MAL",
                                 "NAN",
                                 "BRY",
                                 "WAI",
                                 "COO",
                                 "CRA",
                                 "PUU",
                                 "CJR",
                                 "SOL"
                                 # "DEL",
                                 # "PAL"
),
meanObs = c(0.61,
            1.74,
            0.91,
            4.52,
            3.54,
            0.16,
            7.61,
            0.03,
            0.00
            # 0.48,
            # 0.03
),
prop0Obs = c(89.7,
             76.9,
             86.9,
             66.9,
             65.9,
             97.5,
             57.4,
             99.2,
             99.9
             # 84.8,
             # 98.6
),
meanNon0 = c(1.5, 1.9, 1.7, 3.3, 2.6, 1.6, 4.5, 1.1, 1
             # NA, NA
)) %>%
  mutate(prop0Obs = prop0Obs/100)

##  Import simulations: --------------------------------------------------------
# note that this drops sims without any non-zeros, but those get dropped from plotting anyway

simOut <- list.files(pattern = "noMod.rds$", full.names = TRUE)

sims <- simOut %>%
  map_dfr(read_rds)

nSims <- length(simOut) * 50

# -------------------------------------------
# S Figure 1 --------------------------------
# -------------------------------------------
p.param <- sims %>%
  select(nbin_mu, zeroinf_pi, nTraps, nNights, controlEfficacy, detection,
         pre_meanObsNon0, pre_prop0Obs) %>%
  mutate(detection = factor(detection, 
                            levels = c(0.1, 0.5),
                            labels = c("Detection probability of 0.1",
                                       "Detection probability of 0.5"))) %>%
  ggplot() +
  geom_point(aes(pre_meanObsNon0, pre_prop0Obs,
                 color = factor(nbin_mu),
                 #                 shape = factor(detection)
  )) +
  scale_x_continuous(breaks = seq(1, 15),
                     labels = seq(1, 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  # geom_text(data = biocomplexity,
  #           aes(meanNon0, prop0Obs,
  #               label = site), color = "black", size = 4) +
  geom_point(data = biocomplexity,
             aes(meanNon0, prop0Obs), 
             color = "black", 
             size = 4, pch = 18) +
  facet_wrap(~ detection) +
  #  scale_color_brewer(name = "Simulated mean") +
  scale_color_viridis_d(name = "Simulated mean",
                        option = "plasma") +
  xlab("Mean of non-zero observations") +
  ylab("Proportion of zeros in observations") +
  theme(axis.title = element_text(size = 14),
        strip.text = element_text(size = 10))
p.param
ggsave("power_ParamJust_SuppFig.png", width = 8, height = 4, dpi = 450)




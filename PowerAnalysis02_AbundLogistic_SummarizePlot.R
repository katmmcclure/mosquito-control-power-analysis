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

## Script 2:
# Power analysis for estimating mosquito control efficacy - plotting and summarizing results

# Code by Helen Sofaer

library(tidyverse)
theme_set(theme_bw() +
            theme(axis.text = element_text(color = "black")))
setwd("C:/Users/hsofaer/OneDrive - DOI/ResearchProjects/Mosquito/PowerAnalysis/EcoApps/output/LinkedZeros")

simOut <- list.files(pattern = "noMod.rds$", 
                     full.names = TRUE)

sims <- simOut %>%
  map_dfr(read_rds)
nSims <- length(simOut) * 50


##  Level of zero inflation for each level of abundance: -----------------------
sims %>%
  select(nbin_mu, controlEfficacy, zeroinf_pi_pre, zeroinf_pi_post) %>%
  distinct()


##  Calculate power: -----------------------------------------------------------

abund_simPower <- sims %>%
  mutate(nTraps = factor(nTraps, 
                         levels = levels(factor(nTraps)),
                         labels = paste(levels(factor(nTraps)), "traps")),
         nNights = factor(nNights, 
                          levels = levels(factor(nNights)),
                          labels = paste(levels(factor(nNights)), "nights"))
  ) %>%
  mutate(sig = ifelse(abund_estTrtP < .05 & abund_Converge1ok == 1, 1, 0)) %>%
  group_by(nbin_mu, zeroinf_pi_pre, zeroinf_pi_post, nTraps, nNights, controlEfficacy, detection) %>%
  summarize(n_failConverge = sum(abund_Converge1ok == 0 & post_prop0Obs < 1),
            n_all0Post = sum(post_prop0Obs == 1),
            n_sig = sum(sig),
            propSig = sum(sig)/nSims,
            mean_pre_meanObs = mean(pre_meanObs),
            n = n()) %>%
  mutate(n_all0Pre = nSims - n, # I didn't try to fit models to all 0s pre-treatment
         n_nonSig = nSims - n_all0Pre - n_all0Post - n_failConverge - n_sig)
head(abund_simPower)


##  stacked bar plot:

# make data long:
abund_simPower.long <- abund_simPower %>%
  select(nbin_mu, zeroinf_pi_pre, zeroinf_pi_post, nTraps, nNights, controlEfficacy, detection,
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

p.p5c9 <- abund_simPower.long %>%  
  filter(detection == .5,
         controlEfficacy == .9) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .3,
           position = "stack") +
  scale_fill_viridis_d(name = "Abundance model\noutcome when\np = 0.5\nControl efficacy = 0.9") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  #  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/power_stackBar_p5c9.png",
       width = 9, height = 6, dpi = 450, units = "in")  

p.p1c9 <- abund_simPower.long %>%  
  filter(detection == .1,
         controlEfficacy == .9) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .05,
           position = "stack") +
  scale_fill_viridis_d(name = "Abundance model\noutcome when\np = 0.1\nControl efficacy = 0.9") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/power_stackBar_p1c9.png",
       width = 9, height = 6, dpi = 450, units = "in") 

p.p1c6 <- abund_simPower.long %>%  
  filter(detection == .1,
         controlEfficacy == .6) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .05,
           position = "stack") +
  scale_fill_viridis_d(name = "Abundance model\noutcome when\np = 0.1\nControl efficacy = 0.6") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  # scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/power_stackBar_p1c6.png",
       width = 9, height = 6, dpi = 450, units = "in") 

p.p5c6 <- abund_simPower.long %>%  
  filter(detection == .5,
         controlEfficacy == .6) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .2,
           position = "stack") +
  scale_fill_viridis_d(name = "Abundance model\noutcome when\np = 0.5\nControl efficacy = 0.6") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  #  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/power_stackBar_p5c6.png",
       width = 9, height = 6, dpi = 450, units = "in")

# prob = .1, efficacy = .9 on scale that best fits just that one:
p.p1c9.focus <- p.p1c9 +
  scale_x_continuous(limits = c(0, 1))
ggsave("./figs/power_stackBar_p1c9_focus.png",
       width = 9, height = 6, dpi = 450, units = "in")

##  Reason for failure across model runs: --------------------------------------

p.outcome <- abund_simPower.long %>%
  mutate(controlEfficacy = factor(controlEfficacy, 
                                  levels = levels(factor(controlEfficacy)),
                                  labels = paste0("Control efficacy of ", levels(factor(controlEfficacy)))),
         detection = factor(detection, 
                            levels = levels(factor(detection)),
                            labels = paste0("Detection probability of ", levels(factor(detection))))
  ) %>%
  ggplot(aes(propSims, fill = nTraps)) +
  facet_grid(detection*controlEfficacy ~ outcome, scales = "free_x") +
  geom_histogram() +
  scale_fill_brewer(name = "Number of traps", palette = "Dark2") +
  xlab("Proportion of simulations with this outcome") +
  ylab("Number of simulation parameter combinations") +
  theme(strip.text = element_text(size = 8))
ggsave("./figs/outcomeByDetectControl.png",
       width = 12, height = 9, dpi = 450, units = "in")

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
ggsave("./figs/sigByDetectControl.png",
       width = 9, height = 6, dpi = 450, units = "in")



##  logistic power: -----------------------------------------------------------

logistic_simPower <- sims %>%
  mutate(nTraps = factor(nTraps, 
                         levels = levels(factor(nTraps)),
                         labels = paste(levels(factor(nTraps)), "traps")),
         nNights = factor(nNights, 
                          levels = levels(factor(nNights)),
                          labels = paste(levels(factor(nNights)), "nights"))
  ) %>%
  mutate(sig = ifelse(!is.na(logistic_estTrtP) & logistic_estTrtP < .05 & logistic_Converge1ok == 1, 
                      1, 0)) %>%
  group_by(nbin_mu, zeroinf_pi_pre, zeroinf_pi_post, nTraps, nNights, controlEfficacy, detection) %>%
  summarize(n_failConverge = sum(logistic_Converge1ok == 0 & post_prop0Obs < 1),
            n_all0Post = sum(post_prop0Obs == 1),
            n_sig = sum(sig),
            propSig = sum(sig)/nSims,
            mean_pre_meanObs = round(mean(pre_meanObs), 2),
            mean_pre_Prevalence = round(mean(1 - pre_prop0Obs), 2),
            mean_PreMinusPost_meanObs = round(mean(pre_meanObs - post_meanObs), 2),
            mean_PreMinusPost_Prevalence = round(mean((1 - pre_prop0Obs) - (1 - post_prop0Obs)), 2),
            n = n()) %>%
  mutate(n_all0Pre = nSims - n,
         n_nonSig = nSims - n_all0Pre - n_all0Post - n_failConverge - n_sig)

with(logistic_simPower, plot(mean_pre_meanObs, mean_pre_Prevalence))

##  stacked bar plot:

# make data long:
logistic_simPower.long <- logistic_simPower %>%
  select(nbin_mu, zeroinf_pi_pre, zeroinf_pi_post, nTraps, nNights, controlEfficacy, detection,
         starts_with("n_"), 
         mean_pre_meanObs, mean_pre_Prevalence) %>%
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


p.logistic.p5c9 <- logistic_simPower.long %>%  
  filter(detection == .5,
         controlEfficacy == .9) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .2,
           position = "stack") +
  scale_fill_viridis_d(name = "Logistic model\noutcome when\np = 0.5\nControl efficacy = 0.9") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  #  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/logisticPower_stackBar_p5c9.png",
       width = 9, height = 6, dpi = 450, units = "in")

p.logistic.p5c6 <- logistic_simPower.long %>%  
  filter(detection == .5,
         controlEfficacy == .6) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .2,
           position = "stack") +
  scale_fill_viridis_d(name = "Logistic model\noutcome when\np = 0.5\nControl efficacy = 0.6") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  #  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/logisticPower_stackBar_p5c6.png",
       width = 9, height = 6, dpi = 450, units = "in")

p.logistic.p1c6 <- logistic_simPower.long %>%  
  filter(detection == .1,
         controlEfficacy == .6) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .05,
           position = "stack") +
  scale_fill_viridis_d(name = "Logistic model\noutcome when\np = 0.1\nControl efficacy = 0.6") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  #  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/logisticPower_stackBar_p1c6.png",
       width = 9, height = 6, dpi = 450, units = "in")

p.logistic.p1c9 <- logistic_simPower.long %>%  
  filter(detection == .1,
         controlEfficacy == .9) %>%
  ggplot(aes(x = mean_pre_meanObs, y = propSims,
             fill = outcome)) +
  geom_col(width = .05,
           position = "stack") +
  scale_fill_viridis_d(name = "Logistic model\noutcome when\np = 0.1\nControl efficacy = 0.9") +
  geom_hline(yintercept = .8, linetype = "dashed") +
  facet_grid(nTraps ~ nNights) +
  geom_line(data = . %>%
              filter(outcome == "Significant treatment effect")) +
  #  scale_x_continuous(limits = c(0, 2.5)) +
  ylab("Proportion of simulations") +
  xlab("Mean observed pre-treatment mosquito count\n(number per night per trap)") +
  theme(axis.title = element_text(size = 12))
ggsave("./figs/logisticPower_stackBar_p1c9.png",
       width = 9, height = 6, dpi = 450, units = "in")


##  Plot focusing on significant outcomes:

p.logistic.sig <- logistic_simPower.long %>%
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
ggsave("./figs/Logistic_sigByDetectControl.png",
       width = 9, height = 6, dpi = 450, units = "in")


###  Abundance vs logistic power: ----------------------------------------------

comparePower <- abund_simPower.long %>%
  select(-numSims) %>%
  rename(abund_propSims = propSims) %>%
  inner_join(logistic_simPower.long %>%
               select(-numSims) %>%
               rename(logistic_propSims = propSims)) 

head(comparePower)

p.comparePower <- comparePower %>%
  filter(!grepl("^Zero", outcome)) %>%
  ggplot() +
  geom_point(aes(abund_propSims, logistic_propSims,
                 color = outcome)) +
  # match colors of others:
  scale_color_manual(values = c("#21908CFF",
                                "#5DC863FF",
                                "#FDE725FF"),
                     name = "Outcome") +
  facet_wrap(~ outcome) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Abundance model (proportion of simulations)") +
  ylab("Logistic model (proportion of simulations)") +
  theme(axis.text = element_text(size = 8),
        strip.text = element_text(size = 8))
ggsave("./figs/comparePower_abundvslogistic.png",
       width = 8, height = 4)  

# Add post-control prevalence:
postProp0 <- sims %>%
  mutate(nTraps = factor(nTraps,
                         levels = levels(factor(nTraps)),
                         labels = paste(levels(factor(nTraps)), "traps")),
         nNights = factor(nNights,
                          levels = levels(factor(nNights)),
                          labels = paste(levels(factor(nNights)), "nights"))) %>%
  group_by(nbin_mu, zeroinf_pi_pre, zeroinf_pi_post, nTraps, nNights, controlEfficacy, detection) %>%
  summarize(mean_postProp0Obs = mean(post_prop0Obs))

comparePower <- comparePower %>%
  inner_join(postProp0)
# not too useful since so many are .99 or 1. Includes all sims, not just sig ones

# just the significant ones:
p.compareSig <- comparePower %>%
  filter(outcome == "Significant treatment effect") %>%
  ggplot() +
  geom_point(aes(abund_propSims, logistic_propSims,
                 # shape = factor(detection),
                 # color = factor(controlEfficacy)
                 ), fill = "#FDE725FF", shape = 21) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Abundance model:\nproportion of significant simulations") +
  ylab("Logistic model:\nproportion of significant simulations") +
  theme(axis.text = element_text(size = 8),
        strip.text = element_text(size = 8))
ggsave("./figs/comparePower_sigOnly.png",
       width = 6, height = 5)

comparePower %>%
  group_by(outcome) %>%
  summarize(cor = cor(logistic_propSims, abund_propSims))


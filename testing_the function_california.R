#example with using the code

#load in packages
library(tidyverse)
#tidysynth for the smoking dataset
library(tidysynth)
library(glmnet)
library(meboot)

data("smoking")


#load in my functions
source("Synthetic_control_functions/synthetic_control_functions_claude.R")


#check na
missing <- smoking |> 
  mutate(NA_cigsale = is.na(cigsale),
         NA_lnincome = is.na(lnincome),
         NA_beer = is.na(beer),
         NA_age = is.na(age15to24),
         NA_retprice = is.na(retprice)) |> 
  select(NA_cigsale:NA_retprice) |> 
  mutate(NA_cigsale = sum(NA_cigsale == "TRUE"),
         NA_lnincome = sum(NA_lnincome == "TRUE"),
         NA_beer = sum(NA_beer == "TRUE"),
         NA_age = sum(NA_age == "TRUE"),
         NA_retprice = sum(NA_retprice == "TRUE")) |> 
  distinct(NA_cigsale,.keep_all = TRUE)


#min and max 
min(smoking$cigsale)
max(smoking$cigsale)

#lets try the synthetic control

m1  <- synth_lasso_boot(
  data             = smoking,
  unit_col         = "state",
  time_cols        = "year",
  outcome_col      = "cigsale",
  treated_unit     = "California",
  treatment_period = c(year = 1989),
  boot_by_unit     = TRUE,
  seed             = 09272026,
  M                = 1,
  jack_conform     = TRUE,
  conform_level    = 0.95,
  reps             = 1000,
  trim             = list(trim = 0.10, xmin = 40),
  conf_int         = 0.95)



#lets plot this
ggplot(m1$effects, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = pred_jack_lo, ymax = pred_jack_hi),
              fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = treated_predicted, color = "Synthetic")) +
  geom_line(aes(y = treated_observed, color = "Observed")) +
  geom_vline(xintercept = "1989", linetype = 2, color = "gray40") +
  scale_color_manual(values = c(Observed = "black", Synthetic = "blue")) +
  labs(title = "California Cigarette Sales: observed vs synthetic",
       subtitle = "Shaded: 95% jackknife conformal prediction interval (far too narrow in this example)",
       x = "Period", y = "Cigarette Sales per capita", color = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#lets plot the sensitivity bounds
ggplot(m1$effects, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = bounds_lo, ymax = bounds_hi), fill = "tomato", alpha = 0.2) +
  geom_line(aes(y = att), color = "blue") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = "1989", linetype = 2, color = "gray40") +
  labs(title = "ATT with sensitivity bounds (M = 1)",
       subtitle = "Band = ATT ± the largest absolute pre-treatment gap",
       x = "Period", y = "ATT Cigarette Sales per capita") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#lets plot the att 
ggplot(m1$effects, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = eff_low, ymax = eff_high), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = att), color = "blue") +
  geom_point(aes(y = att), color = "blue", size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = "1989", linetype = 2, color = "gray40") +
  labs(title = "ATT: observed minus synthetic California (95% bootstrap intervals)",
       x = "Period", y = "ATT (Cigarette Sales per Capita)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#lets try with controls
min_retprice <- min(smoking$retprice)
min_cigsale <- min(smoking$cigsale)


#lets add in population controls and see what happens
m2 <- synth_lasso_boot(
  data             = smoking,
  unit_col         = "state",
  time_cols        = "year",
  outcome_col      = "cigsale",
  treated_unit     = "California",
  treatment_period = c(year = 1989),
  controls         = "retprice",
  control_units    = "donors",
  scale_controls   = TRUE,
  seed             = 09272026,
  boot_by_unit     = FALSE, #ensure that the bootstrap is done by unit
  M                = 1,
  jack_conform     = TRUE,
  conform_level    = 0.95,
  reps             = 1000,
  trim             = list(cigsale = list(xmin = min_cigsale), 
                          retprice = list(xmin = min_retprice)),
  conf_int         = 0.95)


view(m2$weights)


#lets plot this
ggplot(m2$effects, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = pred_jack_lo, ymax = pred_jack_hi),
              fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = treated_predicted, color = "Synthetic")) +
  geom_line(aes(y = treated_observed, color = "Observed")) +
  geom_vline(xintercept = "1989", linetype = 2, color = "gray40") +
  scale_color_manual(values = c(Observed = "black", Synthetic = "blue")) +
  labs(title = "California Cigarette Sales: observed vs synthetic",
       subtitle = "Shaded: 95% jackknife conformal prediction interval (far too narrow in this example)",
       x = "Period", y = "Cigarette Sales per capita", color = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#lets plot the sensitivity bounds
ggplot(m2$effects, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = bounds_lo, ymax = bounds_hi), fill = "tomato", alpha = 0.2) +
  geom_line(aes(y = att), color = "blue") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = "1989", linetype = 2, color = "gray40") +
  labs(title = "ATT with sensitivity bounds (M = 1)",
       subtitle = "Band = ATT ± the largest absolute pre-treatment gap",
       x = "Period", y = "ATT Cigarette Sales per capita") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#lets plot the att 
ggplot(m2$effects, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = eff_low, ymax = eff_high), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = att), color = "blue") +
  geom_point(aes(y = att), color = "blue", size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = "1989", linetype = 2, color = "gray40") +
  labs(title = "ATT: observed minus synthetic California (95% bootstrap intervals)",
       x = "Period", y = "ATT (Cigarette Sales per Capita)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))















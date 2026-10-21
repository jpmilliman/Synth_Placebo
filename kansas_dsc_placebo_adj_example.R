library(tidyverse)
library(augsynth)
library(fect)
library(CVXR)
library(augsynth)
library(forcats)
library(patchwork)
library(panelView)


source("synthetic_control_functions/synthetic_control_functions_claude.R")
kansas$emplvlcapita

kansas <- augsynth::kansas |>
  select(state, year, qtr, lngdp, avg_wkly_wage,
         estabscapita, emplvlcapita) |>
  as.data.frame()


#covarieates
# (i) log state and local revenue per capita, 
# (ii) log average weekly wages,
# (iii) number of establishments per capita, 
# (iv) the employment level, and (v) log GSP per capita

#unadjusted
fit <- dsc_placebo(
  data             = kansas,
  unit_col         = "state",
  time_cols        = c("year", "qtr"),
  outcome_col      = "lngdp",
  treated_unit     = "Kansas",
  treatment_period = c(year = 2012, qtr = 2)
)


fit$diagnostics


plot_full_path(fit, "Kansas vs synthetic Kansas, 1990-2016 (dsc_placebo defaults)", y_label = "Log(GDP) per capita ($)")
#adjusted

#lets see if we find an effect
plot_post_effect(fit = fit, measure = "cumulative")

plot_post_effect(fit = fit, measure = "effect")

#lets try this adjusted
#fitting with covariates can take a really long time -its essentially not working
#with certain datasets
fit_adj_ife <- dsc_placebo(
  data             = kansas,
  unit_col         = "state",
  time_cols        = c("year", "qtr"),
  outcome_col      = "lngdp",
  treated_unit     = "Kansas",
  adjust_covariates = c("estabscapita", "emplvlcapita", "avg_wkly_wage"),
  adjust_model = "ife",
    treatment_period = c(year = 2012, qtr = 2))



plot_full_path(fit_adj_ife, "Kansas vs synthetic Kansas, 1990-2016 (dsc_placebo defaults)", y_label = "Log(GDP) per capita ($)")



view(fit_adj_ife$test)


#lets see if we find an effect
plot_post_effect(fit = fit_adj_ife, 
                 measure = "cumulative")

plot_post_effect(fit = fit_adj_ife, 
                 measure = "effect")


#adjusted
pred_adj <- plot_prediction(fit=fit_adj_ife, title = "Adjusted")

#unadjusted
pred <- plot_prediction(fit = fit, title = "Unadjusted")

#plot adjusted and unadjusted treatment effects
plots <- pred_adj/pred

plots


#lets look at some plots
kansas <- kansas |> 
  group_by(state) |> 
  mutate(unit_id = seq(from =1, to = 105, by = 1)) |> 
  ungroup() |> 
  mutate(treated = case_when(state == "Kansas" & unit_id >= 90~ 1,
                             .default = 0))

treatment <- panelview(lngdp ~ treated,
            data = kansas, index = c("state", "unit_id"),
            type = "outcome", main = "",
            xlab = "QTR_Year", ylab = "Log (GDP per Capita)",
            theme = "red")


estab <- panelview(estabscapita ~ treated,
                       data = kansas, index = c("state", "unit_id"),
                       type = "outcome", main = "",
                       xlab = "QTR_Year", ylab = "Establishments per Capita",
                   theme = "red")



employ <- panelview(emplvlcapita ~ treated,
                   data = kansas, index = c("state", "unit_id"),
                   type = "outcome", main = "",
                   xlab = "QTR_Year", ylab = "Employment per Capita",
                   theme = "red")


week_wages <- panelview(avg_wkly_wage ~ treated,
                    data = kansas, index = c("state", "unit_id"),
                    type = "outcome", main = "",
                    xlab = "QTR_Year", ylab = "Employment per Capita",
                    theme = "red")



#regional plots
treat_region <- panelview(lngdp ~ treated,
                       data = kansas, index = c("state", "unit_id"),
                       type = "outcome", main = "",
                       xlab = "QTR_Year", ylab = "Log (GDP per Capita)",
                       id = c("Kansas", "Nebraska", "Missouri", "Oklahoma", "Colorado", "Arkansas"),
                       theme = "red")



estab_region <- panelview(estabscapita ~ treated,
                   data = kansas, index = c("state", "unit_id"),
                   type = "outcome", main = "",
                   xlab = "QTR_Year", ylab = "Establishments per Capita",
                   id = c("Kansas", "Nebraska", "Missouri", "Oklahoma", "Colorado", "Arkansas"),
                   theme = "red")


employ_region <- panelview(emplvlcapita ~ treated,
                    data = kansas, index = c("state", "unit_id"),
                    type = "outcome", main = "",
                    xlab = "QTR_Year", ylab = "Employment per Capita",
                    id = c("Kansas", "Nebraska", "Missouri", "Oklahoma", "Colorado", "Arkansas"),
                    theme = "red")


week_wages_region <- panelview(avg_wkly_wage ~ treated,
                        data = kansas, index = c("state", "unit_id"),
                        type = "outcome", main = "",
                        xlab = "QTR_Year", ylab = "Employment per Capita",
                        id = c("Kansas", "Nebraska", "Missouri", "Oklahoma", "Colorado", "Arkansas"),
                        theme = "red")




# ---- libs -----

suppressPackageStartupMessages(suppressWarnings({
  
  library("readr")     # read data funcs
  library("dplyr")     # data manipulation
  library("tibble")    # improved data frames
  library("ggplot2")   # plotting
  library("purrr")     # apply functions over lists/vectors
  library("forcats")   # handling categorical variables
  library("tidyr")     # manipulate data to long/short representations
  
  library("mice")      # missing value utilities
  library("car")       # regression utilities
  
  library("lme4")      # linear mixed effects modelling
  library("merTools")  # alows prediction intervals using merMod objects
  library("lmerTest")  # step-wise lmer model selection
  
  library("knitr")     # pretty printing of tables: kable()
  library("gtsummary") # print summary tables of regression mods
  library("lattice")   # diagnostic plots

}))

# ---- consts ----


### plotting characters
# steep x long: "|"
# steep x short:  "/"
# less steep x long: "-"
# less steep x short: "\"
plchs <- 
  c(
    "yes x long" = "|", 
    "yes x short" = "/", 
    "no x long" = "-", 
    "no x short" = "\\"
  )

# plotting colour scheme
col_lohi <- 
  c(
    "Lower IAcc" = "darkorange", 
    "Higher IAcc" = "purple"
  )




# ---- funcs1 ----

# effect size (f^2), R^2 and residual variance functions

get_res_var_lmer <- function(lmer_obj) {
  return(sigma(lmer_obj) ^ 2)
}

get_lmer_r2 <- function(lmer_obj) {
  
  # residual variance of input model
  v_mod <- get_res_var_lmer(lmer_obj)
  
  # get formula for null model (intercept and REs only)
  null_form <- formula(lmer_obj, random.only = TRUE)
  # create null model
  lmer_obj_null <- update(lmer_obj, null_form)
  
  # res var of null mod
  v_null <- get_res_var_lmer(lmer_obj_null)
  
  r2 <- (v_null - v_mod) / v_null
  
  attr(r2, "v_null") <- v_null
  attr(r2, "v_mod") <- v_mod
  
  return(r2)
  
}

rm_terms_lmer <- function(lmer_obj, terms) {
  
  update_form <- as.formula(paste0("~ . -", paste(terms, collapse = " - ")))
  print(update_form)

  lmer_obj_less_term <- update(lmer_obj, update_form)
  
  return(lmer_obj_less_term)
}


eff_size_f2 <- function(lmer_obj, terms) {
  
  r2_full <- get_lmer_r2(lmer_obj)[1]
  r2_less_term <- get_lmer_r2(rm_terms_lmer(lmer_obj, terms))[1]
  
  f2 <- (r2_full - r2_less_term) / (1 - r2_full)
  
  return(f2)
  
}


# ---- funcs2 ----

# function that replaces repeated values in a vector with empty strings
# as the verbose redundancy is too busy in some cases 
rm_rpts <- function(x) {
  x <- as.character(x)
  nx <- length(x)
  rm_ii <- rep(FALSE, nx)
  for (i in 2:nx) {
    if (x[i - 1] == x[i]) 
      rm_ii[i] <- TRUE
  }
  x[rm_ii] <- ""
  return(x)
}


# ---- read -----

dat_col_spec <- 
  cols(
    partic = col_integer(),
    affval = col_integer(),
    perexe = col_integer(),
    int_sens = col_double(),
    block = col_integer(),
    cond = col_character(),
    steep = col_character(),
    dist = col_character()
  )

# read in dataset
hill_dat <- read_csv("dat/seefeel-hill-dat.csv", col_types = dat_col_spec)

# have a peak
hill_dat

# ---- wrangle ----

# make factor variables and default levels
hill_dat <-
  hill_dat %>%
  mutate(
    cond = factor(cond),
    cond = relevel(cond, ref = "flat"),
    steep = factor(steep),
    steep = relevel(steep, ref = "no"),
    dist = factor(dist),
    dist = relevel(dist, ref = "short"),
    block = factor(block)
  ) 

# NAs only present in outcome vars
# md.pattern(hill_dat, rotate.names = TRUE)

### centring the int_sens variable for easier interpretation of model intercept terms
# NOTE: want average of average participant values
isc <-
  hill_dat %>%
  group_by(partic) %>%
  summarise(avg_int_sens = mean(int_sens))

# This is the mean ISC over participants
mean_isc <- isc %>% pull(avg_int_sens) %>% mean(.)
sd_isc <- isc %>% pull(avg_int_sens) %>% sd(.)

# now modify the ISC values in the data
hill_dat <-
  hill_dat %>%
  mutate(int_sens = int_sens - mean_isc)



# ---- predictors_dataset ----


# create a minimal prediction dataset
pred_dat <-
  hill_dat %>% 
  distinct(cond, block) %>%
  # NB: this is the mean ISC as we centred this data previously
  mutate(int_sens = 0, partic = 1L) 

pred_dat_isc_plus_sd <-
  pred_dat %>%
  mutate(int_sens = int_sens + sd_isc) 

pred_dat_isc_less_sd <-
  pred_dat %>%
  mutate(int_sens = int_sens - sd_isc) 

pred_dat <- 
  bind_rows(pred_dat_isc_less_sd, pred_dat, pred_dat_isc_plus_sd) %>%
  as.data.frame(.)



# ---- plot_dat ----


plot_dat <-
  hill_dat %>%
  mutate(
    cond = as.character(cond),
    hill = paste(steep, dist, sep = " x "),
    cond_block_hill = 
      paste(
        as.integer(as.factor(cond)), 
        block, 
        as.integer(steep), 
        as.integer(dist), 
        sep = "_"
      )
  ) %>%
  group_by(partic) %>%
  mutate(
    int_sens = ifelse(row_number() == 1, sprintf("int_sens=%1.3f", int_sens), NA),
    x1 = ifelse(row_number() == 1,  3*4*4/2, NA),
    y1 = ifelse(row_number() == 1,  4, NA)
  ) %>%
  ungroup() 

plot_dat %>%
  distinct(steep, dist, hill)



# ---- affval_exploratory_plot -----

plot_dat %>%
  ggplot(
    ., 
    aes(
      x = cond_block_hill, 
      y = affval, 
      col = as.factor(cond), 
      shape = as.factor(hill)
    )
  ) +
  geom_point(size = 5) + #aes(size = block)
  scale_shape_manual(values = plchs) +
  scale_colour_manual(
    values = c("dodgerblue", "orchid", "firebrick1")
  ) +
  scale_x_discrete(labels = rm_rpts(plot_dat$block)) + 
  geom_text(aes(x = x1, y = y1,label = int_sens), col = "grey60") +
  geom_line(aes(group = partic), col = "black", alpha = 0.25) +
  facet_wrap(~ partic, labeller = "label_both") +
  labs(col = "cond", shape = "hill", x = "block within condition") +
  theme_bw() 


# ---- perexe_exploratory_plot -----



plot_dat %>%
  ggplot(
    ., 
    aes(
      x = cond_block_hill, 
      y = perexe, 
      col = as.factor(cond), 
      shape = as.factor(hill)
    )
  ) +
  geom_point(size = 5) + #aes(size = block)
  scale_shape_manual(values = plchs) +
  scale_colour_manual(
    values = c("dodgerblue", "orchid", "firebrick1")
  ) +
  scale_x_discrete(labels = rm_rpts(plot_dat$block)) + 
  geom_text(aes(x = x1, y = y1,label = int_sens), col = "grey60") +
  geom_line(aes(group = partic), col = "black", alpha = 0.25) +
  facet_wrap(~ partic, labeller = "label_both") +
  labs(col = "cond", shape = "hill", x = "block within condition ") +
  theme_bw() 




# ---- affval_mod_stepwise_sel ----

hill_dat_affmod <- subset(hill_dat, !is.na(affval))

# largest potential model
m1 <- 
  lmer(
    affval ~ 
      int_sens * (cond + steep + dist + block) + 
      (1 | partic), 
    data = hill_dat_affmod, 
    REML = FALSE
  )
# summary(m1)

# elimination of non-significant effects
# partly thanks to code found at:: 
# https://www.rdocumentation.org/packages/lmerTest/versions/2.0-36/topics/step
s1 <- step(m1) # consider optional arguments: test = c("none", "Rao", "LRT", "Chisq", "F")

# look at the model reduction
print(s1)

# plot of post-hoc analysis of the final model
# plot(s1)


# use REML for final fit... see the following links
# https://stats.stackexchange.com/questions/41123/reml-vs-ml-stepaic
# https://stats.stackexchange.com/questions/116770/reml-or-ml-to-compare-two-mixed-effects-models-with-differing-fixed-effects-but
# https://stats.stackexchange.com/questions/414551/forward-selection-with-mixed-model-using-lmer
m1_final <- 
  lmer(
    affval ~ int_sens + cond + block + int_sens:block + 
      (1 | partic), 
    data = hill_dat_affmod, 
    REML = TRUE
  )
summary(m1_final)


# term significance -- Type III Wald chi-square tests
car::Anova(m1_final, type = "III")

# prettier printing of regression model
m1_final %>%
  tbl_regression(
    estimate_fun = function(x) sprintf("%2.2f", x),
    pvalue_fun = function(x) sprintf("%1.8f", x)
  ) %>%
  add_global_p(keep = TRUE) %>%
  as_gt(.)





# ---- affval_mod_diagnostics ----


# look at BLUPs of random effects:
random_ints <- ranef(m1_final)$partic[["(Intercept)"]]
hist(random_ints, main = "BLUPs of random intercepts")

# https://rdrr.io/cran/lme4/man/plot.merMod.html

## standardized residuals versus fitted values by variable
# resid(m1_final, scaled=TRUE)
# fitted(m1_final)
plot(m1_final, resid(., scaled=TRUE) ~ fitted(.) | block, abline = 0) 
plot(m1_final, resid(., scaled=TRUE) ~ fitted(.) | partic, abline = 0) 
plot(m1_final, resid(., scaled=TRUE) ~ fitted(.) | cond, abline = 0) 
plot(m1_final, resid(., scaled=TRUE) ~ fitted(.) | steep, abline = 0) 
plot(m1_final, resid(., scaled=TRUE) ~ fitted(.) | dist, abline = 0) 
## residuals by Subject
plot(m1_final, partic ~ resid(., scaled=TRUE))
qqmath(m1_final, id=0.05)



# ---- affval_effect_sizes ----


### testing and extracting model elements for f^2 calc
# terms(formula(m1_final, fixed.only = TRUE))
# summary(m1_final)
# summary(rm_terms_lmer(m1_final, "cond"))
# summary(rm_terms_lmer(m1_final, c("int_sens", "int_sens:block")))

### test functions
# get_lmer_r2(m1_final)
# get_lmer_r2(rm_terms_lmer(m1_final, "cond"))
### need to consider higher level terms with main effects
# get_lmer_r2(rm_terms_lmer(m1_final, c("int_sens", "int_sens:block")))

# cond eff size
eff_size_f2(m1_final, "cond")

(get_lmer_r2(m1_final) - 
    get_lmer_r2(rm_terms_lmer(m1_final, c("cond")))) /
  (1 - get_lmer_r2(m1_final)) # manual check

# int_sens eff size
eff_size_f2(m1_final, c("int_sens", "int_sens:block"))

(get_lmer_r2(m1_final) - 
    get_lmer_r2(rm_terms_lmer(m1_final, c("int_sens", "int_sens:block")))) /
  (1 - get_lmer_r2(m1_final))  # manual check

# interaction only eff size
eff_size_f2(m1_final, "int_sens:block")

(get_lmer_r2(m1_final) - 
    get_lmer_r2(rm_terms_lmer(m1_final, "int_sens:block"))) /
  (1 - get_lmer_r2(m1_final))  # manual check


# ---- affval_mod_preds ----

# ?predict.merMod
# ?predict

### usage
# predict(
#   object, newdata = NULL, newparams = NULL,
#   re.form = ~0, # or NULL for no REs
#   random.only = FALSE, terms = NULL,
#   type = c("link", "response"), allow.new.levels = FALSE,
#   na.action = na.pass, ...
# )


# test the above centring claim
# hist(hill_dat$int_sens_cont)

pred_est <-
  predict(
    m1_final,
    newdata = pred_dat,
    re.form = ~0
  )


# see:
# https://cran.r-project.org/web/packages/merTools/vignettes/Using_predictInterval.html
# ?merTools::predictInterval

pred_ci_est <-
  merTools::predictInterval(
    merMod = m1_final, 
    newdata = pred_dat,
    which = c("full", "fixed", "random", "all")[4],
    level = 0.95, 
    n.sims = 1000,
    stat = "median", 
    type = "linear.prediction",
    include.resid.var = TRUE, # TRUE for including
    # fix.intercept.variance = TRUE
    seed = 1234567890
  ) %>% 
  dplyr::filter(effect == "fixed") %>%
  arrange(obs)

pred_dat_aff <-
  bind_cols(pred_dat, tibble(fit_analytic = pred_est), pred_ci_est) %>%
  as_tibble() 

# cat("#### Min and max difference between analytic fit and bootstrp median is:\n")
# with(pred_dat_aff, min(fit_analytic - fit))
# with(pred_dat_aff, max(fit_analytic - fit))

# pred_dat_aff %>% 
#   dplyr::select(cond, block, int_sens, fit, lwr, upr) %>%
#   kable(., digits = 2)

pred_dat_aff <-
  pred_dat_aff %>%
  mutate(
    ISC_value = 
      ifelse(
        int_sens < (0 - sd_isc/2), 
        "Lower IAcc", # mean(IAcc) - sd(IAcc)
        ifelse(
          int_sens > (0 + sd_isc/2), 
          "Higher IAcc", # mean(IAcc) + sd(IAcc), 
          "Mean IAcc" # mean(IAcc)
        )
      ),
    ISC_value = factor(ISC_value),
    ISC_value = relevel(ISC_value, ref = "Lower IAcc"),
    cond =  relevel(cond, ref = "downhill")
  )


# change likert scale data recorded as [1, 12] to [-5, 5]
likert_adj <- -6

pred_dat_aff %>%
  dplyr::filter(ISC_value != "Mean IAcc") %>%
  mutate(
    trunc_up = if_else(upr > 11, 11L, as.integer(NA)),
    upr = if_else(!is.na(trunc_up), 11, upr),
    fit = fit + likert_adj,
    lwr = lwr + likert_adj,
    upr = upr + likert_adj
  ) %>%
  ggplot(data = ., aes(x= factor(cond), y = fit, col = ISC_value, group = ISC_value)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr, fill = ISC_value), alpha = 0.2, colour = NA) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ block, ncol = 4, labeller = label_both) +
  theme_bw() +
  theme(text = element_text(family = "serif"), panel.grid = element_blank()) +
  scale_color_manual(values = col_lohi) +
  scale_fill_manual(values = col_lohi) +
  labs(
    y = "Predicted Affective Valence",
    x = "Condition",
    col = "Interoceptive\nAccuracy Value",
    fill = "Interoceptive\nAccuracy Value"
  )


# ---- affval_pred_plot_save ----

ggsave(
  filename = "fig/figure-3.png", 
  width = 2 * 3200 , 
  height = 2 * 1700, 
  units = "px",
  dpi = 2 * 300
)



# ---- perexe_mod_stepwise_sel ----

hill_dat_permod <- subset(hill_dat, !is.na(perexe))

# largest potential model
m2 <- 
  lmer(
    perexe ~ 
      int_sens * (cond + steep + dist + block) + 
      (1 | partic), 
    data = hill_dat_permod, 
    REML = FALSE
  )
# summary(m2)

# elimination of non-significant effects
s2 <- step(m2) 

# look at the model reduction
print(s2)



# use REML for final fit
m2_final <- 
  lmer(
    perexe ~ int_sens + cond + block +
      int_sens:cond + int_sens:block + 
      (1 | partic), 
    data = hill_dat_permod, 
    REML = TRUE
  )
summary(m2_final)

# term significance
car::Anova(m2_final, type = "III")

# prettier printing of regression model
m2_final %>%
  tbl_regression(
    estimate_fun = function(x) sprintf("%2.2f", x),
    pvalue_fun = function(x) sprintf("%1.8f", x)
  ) %>%
  add_global_p(keep = TRUE) %>%
  as_gt(.)





# ---- perexe_mod_diagnostics ----


# look at BLUPs of random effects:
random_ints <- ranef(m2_final)$partic[["(Intercept)"]]
hist(random_ints, main = "BLUPs of random intercepts")


## standardized residuals versus fitted values by variable
plot(m2_final, resid(., scaled=TRUE) ~ fitted(.) | block, abline = 0) 
plot(m2_final, resid(., scaled=TRUE) ~ fitted(.) | partic, abline = 0) 
plot(m2_final, resid(., scaled=TRUE) ~ fitted(.) | cond, abline = 0) 
plot(m2_final, resid(., scaled=TRUE) ~ fitted(.) | steep, abline = 0) 
plot(m2_final, resid(., scaled=TRUE) ~ fitted(.) | dist, abline = 0) 
## residuals by Subject
plot(m2_final, partic ~ resid(., scaled=TRUE))
qqmath(m2_final, id=0.05)




# ---- perexe_effect_sizes ----


### testing and extracting model elements for f^2 calc
# terms(formula(m2_final, fixed.only = TRUE))
# summary(m2_final)
# summary(rm_terms_lmer(m2_final, c("int_sens", "int_sens:block")))


# interaction int_sens:block only eff size
eff_size_f2(m2_final, "int_sens:cond")

(get_lmer_r2(m2_final) - 
    get_lmer_r2(rm_terms_lmer(m2_final, "int_sens:cond"))) /
  (1 - get_lmer_r2(m2_final))  # manual check

# interaction int_sens:block only eff size
eff_size_f2(m2_final, "int_sens:block")

(get_lmer_r2(m2_final) - 
    get_lmer_r2(rm_terms_lmer(m2_final, "int_sens:block"))) /
  (1 - get_lmer_r2(m2_final))  # manual check



# ---- perexe_mod_preds ----



# test the above centring claim
# hist(hill_dat$int_sens_cont)

pred_est <-
  predict(
    m2_final,
    newdata = pred_dat,
    re.form = ~0
  )



pred_ci_est <-
  merTools::predictInterval(
    merMod = m2_final, 
    newdata = pred_dat,
    which = c("full", "fixed", "random", "all")[4],
    level = 0.95, 
    n.sims = 1000,
    stat = "median", 
    type = "linear.prediction",
    include.resid.var = TRUE, # TRUE for including
    # fix.intercept.variance = TRUE
    seed = 1234567890
  ) %>% 
  dplyr::filter(effect == "fixed") %>%
  arrange(obs)

pred_dat_per <-
  bind_cols(pred_dat, tibble(fit_analytic = pred_est), pred_ci_est) %>%
  as_tibble() 

# cat("#### Min and max difference between analytic fit and bootstrp median is:\n")
# with(pred_dat_per, min(fit_analytic - fit))
# with(pred_dat_per, max(fit_analytic - fit))

# pred_dat_per %>% 
#   dplyr::select(cond, block, int_sens, fit, lwr, upr) %>%
#   kable(., digits = 2)

pred_dat_per <-
  pred_dat_per %>%
  mutate(
    ISC_value = 
      ifelse(
        int_sens < (0 - sd_isc/2), 
        "Lower IAcc", # mean(IAcc) - sd(IAcc)
        ifelse(
          int_sens > (0 + sd_isc/2), 
          "Higher IAcc", # mean(IAcc) + sd(IAcc), 
          "Mean IAcc" # mean(IAcc)
        )
      ),
    ISC_value = factor(ISC_value),
    ISC_value = relevel(ISC_value, ref = "Lower IAcc"),
    cond =  relevel(cond, ref = "downhill")
  )


pred_dat_per %>%
  dplyr::filter(ISC_value != "Mean IAcc") %>%
  ggplot(data = ., aes(x= factor(cond), y = fit, col = ISC_value, group = ISC_value)) +
  geom_ribbon(aes(ymax = upr, ymin = lwr, fill = ISC_value), alpha = 0.2, colour = NA) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ block, ncol = 4, labeller = label_both) +
  theme_bw() +
  theme(text = element_text(family = "serif"), panel.grid = element_blank()) +
  scale_color_manual(values = col_lohi) +
  scale_fill_manual(values = col_lohi) +
  labs(
    y = "Predicted Ratings of Perceived Exertion",
    x = "Condition",
    col = "Interoceptive\nAccuracy Value",
    fill = "Interoceptive\nAccuracy Value"
  )


# ---- perexe_pred_plot_save ----

ggsave(
  filename = "fig/figure-4.png", 
  width = 2 * 3200 , 
  height = 2 * 1700, 
  units = "px",
  dpi = 2 * 300
)





# ---- sess_info ----

# for reproducibility
sessionInfo()





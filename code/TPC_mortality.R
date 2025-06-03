library(tidyverse)
library(boot)
library(broom)

# Load fecundity data
dv <- read_csv("/Users/rebecca/Downloads/mortality traitdata.csv") %>%
  mutate(temp = as.numeric(temp), rate = as.numeric(rate)) %>%
  drop_na(temp, rate) %>%
  group_by(species) %>%
  mutate(curve_ID = cur_group_id())

# Fit quadratic model and extract predictions + Tpk
FitModel <- function(ID) {
  local_df <- filter(dv, curve_ID == ID)
  if (nrow(local_df) < 3) return(NULL)
  
  Model <- lm(rate ~ temp + I(temp^2), data = local_df)
  
  temp_seq <- seq(min(local_df$temp), max(local_df$temp), length.out = 100)
  pred_df <- tibble(
    temp = temp_seq,
    fitted = predict(Model, newdata = tibble(temp = temp_seq)),
    species = unique(local_df$species)
  )
  
  n_boot <- 200
  boot_obj <- boot(local_df, statistic = function(d, i) {
    boot_fit <- lm(rate ~ temp + I(temp^2), data = d[i, ])
    predict(boot_fit, newdata = tibble(temp = temp_seq))
  }, R = n_boot)
  
  boot_preds <- boot_obj$t
  pred_bounds <- t(apply(boot_preds, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
  pred_df$lower <- pred_bounds[, 1]
  pred_df$upper <- pred_bounds[, 2]
  
  Tpk <- temp_seq[which.max(pred_df$fitted)]
  list(pred_df = pred_df, Tpk = Tpk)
}

results <- map(unique(dv$curve_ID), ~ FitModel(.x)) %>% compact()
predictions_df <- map_dfr(results, "pred_df")

Tpk_df <- tibble(
  species = map_chr(results, ~ .x$pred_df$species[1]),
  Tpk = map_dbl(results, "Tpk")
)

print(predictions_df)
print(Tpk_df)

ggplot(predictions_df, aes(x = temp, y = fitted)) +
  facet_wrap(~species, scales = 'free_y') +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
  geom_line(color = "red") +
  labs(
    x = "Temperature (°C)",
    y = expression("Mortality rate (" * day^{-1} * ")")
  ) +
  theme_bw()

# Bootstrap Tpk CI
bootstrap_tpk_ci <- function(data, R = 200) {
  temp_seq <- seq(min(data$temp), max(data$temp), length.out = 100)
  boot_result <- boot(data = data, statistic = function(d, i) {
    d_boot <- d[i, ]
    fit <- tryCatch(
      lm(rate ~ temp + I(temp^2), data = d_boot),
      error = function(e) return(NA)
    )
    if (inherits(fit, "lm")) {
      preds <- predict(fit, newdata = tibble(temp = temp_seq))
      temp_seq[which.max(preds)]
    } else {
      NA
    }
  }, R = R)
  
  tpk_vals <- na.omit(boot_result$t)
  if (length(tpk_vals) < 10) return(tibble(lower = NA, median = NA, upper = NA))
  ci <- quantile(tpk_vals, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  tibble(lower = ci[1], median = ci[2], upper = ci[3])
}

Tpk_CI_df <- dv %>%
  group_by(species) %>%
  group_split() %>%
  map_df(~ {
    ci <- bootstrap_tpk_ci(.x, R = 200)
    ci$species <- unique(.x$species)
    ci
  }) %>%
  select(species, lower, median, upper)

print(Tpk_CI_df)

ggplot(Tpk_CI_df, aes(x = reorder(species, median), y = median)) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "grey40") +
  coord_flip() +
  labs(
    x = "Species",
    y = expression(T[pk]*" (°C)")
  ) +
  theme_minimal(base_size = 14)

# Extract full bootstrap Tpk distributions
extract_bootstrap_tpk <- function(data, R = 200) {
  temp_seq <- seq(min(data$temp), max(data$temp), length.out = 100)
  
  boot_result <- boot(data = data, statistic = function(d, i) {
    d_boot <- d[i, ]
    fit <- tryCatch(
      lm(rate ~ temp + I(temp^2), data = d_boot),
      error = function(e) return(NA)
    )
    if (inherits(fit, "lm")) {
      preds <- predict(fit, newdata = tibble(temp = temp_seq))
      temp_seq[which.max(preds)]
    } else {
      NA
    }
  }, R = R)
  
  tibble(Tpk = as.vector(boot_result$t))
}

Tpk_boot_df <- dv %>%
  group_by(species) %>%
  group_split() %>%
  map_df(~ {
    df <- extract_bootstrap_tpk(.x, R = 200)
    df$species <- unique(.x$species)
    df
  })

print(Tpk_boot_df)




Pmax_df <- predictions_df %>%
  group_by(species) %>%
  summarise(Pmax = max(fitted, na.rm = TRUE))

Tpk_df_mor <- results_df %>%
  group_by(species) %>%
  summarise(Tpk = temp[which.max(fitted)])


#Merge with Tpk
Tpk_Pmax_df_fec_fec <- left_join(Tpk_df_fec, Pmax_df_fec, by = "species")
Tpk_Pmax_df_mor <- left_join(Tpk_df, Pmax_df, by = "species")

print(Tpk_Pmax_df_mor)

#Plot Pmax vs Tpk
library(ggplot2)
ggplot(Tpk_Pmax_df_mor, aes(x = Tpk, y = Pmax)) +
  geom_point(size = 3, color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  labs(
    x = expression(T[pk]*" (°C)"),
    y = "Maximum trait value (Pmax)",
    title = "'Hotter-is-better' test: Tpk vs Pmax"
  ) +
  theme_minimal(base_size = 14)

#linear regression
lm_model <- lm(Pmax ~ Tpk, data = Tpk_Pmax_df_mor)
summary(lm_model)

#Non-parametric correlation test
cor_test <- cor.test(Tpk_Pmax_df_mor$Tpk, Tpk_Pmax_df_mor$Pmax, method = "spearman")
print(cor_test)


# bootstrap slope CI function
bootstrap_slope_ci <- function(df, R = 1000) {
  boot_model <- boot(data = df, statistic = function(data, idx) {
    model <- lm(Pmax ~ Tpk, data = data[idx, ])
    coef(model)["Tpk"]
  }, R = R)
  
  ci <- quantile(boot_model$t, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  slope_df <- tibble(
    lower = ci[1],
    median = ci[2],
    upper = ci[3]
  )
  print(slope_df)
  
  #slope distribution
  ggplot(data.frame(slope = boot_model$t), aes(x = slope)) +
    geom_histogram(bins = 30, fill = "#E69F00", color = "white", alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = ci, linetype = "dotted", color = "black") +
    labs(
      x = "Bootstrap slope (Tpk → Pmax)",
      y = "Frequency",
      title = "Bootstrap distribution of regression slope (Mortality)"
    ) +
    theme_minimal(base_size = 14)
}

# Step 6.2: Run it!
bootstrap_slope_ci(Tpk_Pmax_df_mor)

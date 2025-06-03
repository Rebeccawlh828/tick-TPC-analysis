
library(tidyverse)
library(boot)
library(broom)

# Load data
dv <- read_csv("/Users/rebecca/Downloads/development traitdata.csv") %>%
  mutate(temp = as.numeric(temp), rate = as.numeric(rate)) %>%
  drop_na(temp, rate) %>%
  group_by(species) %>%
  mutate(curve_ID = cur_group_id())

# Define simplified model fitting function (quadratic)
FitModel <- function(ID) {
  local_df <- filter(dv, curve_ID == ID)
  message("Processing curve_ID: ", ID, " for species: ", unique(local_df$species))
  if (nrow(local_df) < 3) return(NULL)

  # Fit quadratic model
  Model <- lm(rate ~ temp + I(temp^2), data = local_df)

  # Prediction
  temp_seq <- seq(min(local_df$temp), max(local_df$temp), length.out = 100)
  pred_df <- tibble(
    temp = temp_seq,
    fitted = predict(Model, newdata = tibble(temp = temp_seq)),
    species = unique(local_df$species)
  )

  # Bootstrap for CI
  n_boot <- 200
  boot_obj <- boot(local_df, statistic = function(d, i) {
    boot_fit <- lm(rate ~ temp + I(temp^2), data = d[i, ])
    predict(boot_fit, newdata = tibble(temp = temp_seq))
  }, R = n_boot)

  boot_preds <- boot_obj$t
  pred_bounds <- t(apply(boot_preds, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
  pred_df$lower <- pred_bounds[, 1]
  pred_df$upper <- pred_bounds[, 2]

  # Find Tpk (temp with max predicted rate)
  Tpk <- temp_seq[which.max(pred_df$fitted)]
  list(pred_df = pred_df, Tpk = Tpk)
}

# Run model for all species
results <- map(unique(dv$curve_ID), ~ FitModel(.x)) %>% compact()

# Combine prediction results
predictions_df <- map_dfr(results, "pred_df")

# Extract Tpk per species
Tpk_df <- tibble(
  species = map_chr(results, ~ .x$pred_df$species[1]),
  Tpk = map_dbl(results, "Tpk")
)

# print outputs
print(predictions_df)
print(Tpk_df) 

# Plot
# Plot
ggplot(predictions_df, aes(x = temp, y = fitted)) +
  facet_wrap(~species, scales = 'free_y') +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
  geom_line(color = "blue") +
  labs(
    x = "Temperature (°C)",
    y = expression("Development rate (" * day^{-1} * ")")
  ) +
  theme_bw()

Tpk_df %>% arrange(desc(Tpk))
ggplot(Tpk_df, aes(x = reorder(species, Tpk), y = Tpk)) +
  geom_col(fill = "lightblue") +
  coord_flip() +
  labs(x = "Species", y = expression(T[pk] * " (°C)")) +
  theme_minimal()



# Function to bootstrap Tpk CI using quadratic model
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

# Apply bootstrap Tpk CI extraction to each species
Tpk_CI_df <- dv %>%
  group_by(species) %>%
  group_split() %>%
  map_df(~ {
    ci <- bootstrap_tpk_ci(.x, R = 200)
    ci$species <- unique(.x$species)
    ci
  }) %>%
  select(species, lower, median, upper)

# View results
print(Tpk_CI_df)

ggplot(Tpk_CI_df, aes(x = reorder(species, median), y = median)) +
  geom_point(color = "#1f78b4", size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "grey40") +
  coord_flip() +
  labs(
    x = "Species",
    y = expression(T[pk]*" (°C)")
  ) +
  theme_minimal(base_size = 14)





# Extract full bootstrap Tpk distributions from the quadratic model
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
  
  # Turn into dataframe
  tibble(Tpk = as.vector(boot_result$t))
}

# Combine bootstrap Tpk distributions for all species
Tpk_boot_df <- dv %>%
  group_by(species) %>%
  group_split() %>%
  map_df(~ {
    df <- extract_bootstrap_tpk(.x, R = 200)
    df$species <- unique(.x$species)
    df
  })

# View long-format bootstrap Tpk distribution
print(Tpk_boot_df)

Tpk_boot_df_dev <- Tpk_boot_df


#HOTTER IS BETTER
Pmax_df <- predictions_df %>%
  group_by(species) %>%
  summarise(Pmax = max(fitted, na.rm = TRUE))

# Step 2: Merge with Tpk
Tpk_Pmax_df <- left_join(Tpk_df, Pmax_df, by = "species")

# Step 3: Plot Pmax vs Tpk
library(ggplot2)
ggplot(Tpk_Pmax_df, aes(x = Tpk, y = Pmax)) +
  geom_point(size = 3, color = "blue") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  labs(
    x = expression(T[pk]*" (°C)"),
    y = "Maximum trait value (Pmax)",
    title = "'Hotter-is-better' test: Tpk vs Pmax"
  ) +
  theme_minimal(base_size = 14)

# Step 4: Run linear regression
lm_model <- lm(Pmax ~ Tpk, data = Tpk_Pmax_df)
summary(lm_model)

# Step 5 (Optional): Non-parametric correlation test
cor_test <- cor.test(Tpk_Pmax_df$Tpk, Tpk_Pmax_df$Pmax, method = "spearman")
print(cor_test)




#Significant test fev&dev
sp <- "Haemaphysalis leporispalustris"
fec_boot <- Tpk_boot_df_fec %>% filter(species == sp) %>% pull(Tpk)
dev_boot <- Tpk_boot_df_dev %>% filter(species == sp) %>% pull(Tpk)

diff_boot <- fec_boot - dev_boot
quantile(diff_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)


#Significant test lab&field
sp1 <- "Haemaphysalis longicornis (field)"
sp2 <- "Haemaphysalis longicornis"

boot1 <- Tpk_boot_df_dev %>% filter(species == sp1) %>% pull(Tpk)
boot2 <- Tpk_boot_df_dev %>% filter(species == sp2) %>% pull(Tpk)

boot_diff <- boot1 - boot2
quantile(boot_diff, probs = c(0.025, 0.5, 0.975))


library(ggplot2)
tibble(diff = diff_boot) %>%
  ggplot(aes(x = diff)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Tpk difference: field vs lab (development)",
       x = "Tpk(field) - Tpk(lab)",
       y = "Frequency") +
  theme_minimal()




# 查看 bootstrap Tpk 分布
hist(boot1, main = "Field", xlab = "Tpk", col = "blue", breaks = 20)
hist(boot2, main = "Lab", xlab = "Tpk", col = "red", breaks = 20)

# 或者可视化对比
library(ggplot2)
bind_rows(
  tibble(Tpk = boot1, source = "Field"),
  tibble(Tpk = boot2, source = "Lab")
) %>%
  ggplot(aes(x = Tpk, fill = source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Bootstrap Tpk Distributions", x = "Tpk (°C)", y = "Density") +
  theme_minimal()

############################################################
# Assignment 4 – Wealth and Infant Mortality
############################################################

############################
# Setup
############################

pkgs <- c("haven", "dplyr", "ggplot2",
          "modelsummary", "marginaleffects", "sandwich")

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

PLOT_DIR <- "plots_part2"
if(!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR)

save_plot <- function(name, plot, w=7, h=5){
  ggsave(file.path(PLOT_DIR, name), plot = plot, width = w, height = h, dpi = 300)
}

############################
# 1.1 Data exploration
############################

#DATA_PATH <- file.choose()
df <- read_dta(DATA_PATH)

df <- df |>
  mutate(
    country = as.character(country),
    region  = as.factor(region),
    oil     = as.factor(oil),
    income  = as.numeric(income),
    infant  = as.numeric(infant)
  ) |>
  filter(!is.na(income), !is.na(infant),
         income > 0, infant > 0)

modelsummary::datasummary_skim(df)

n_countries <- nrow(df)
print(n_countries)

# COMMENT (2.1a):
# The dataset contains n_countries countries.

# Histograms
p1 <- ggplot(df, aes(infant)) +
  geom_histogram(bins=30) +
  theme_minimal() +
  labs(title="Distribution of Infant Mortality",
       x="Infant mortality", y="Count")

p2 <- ggplot(df, aes(income)) +
  geom_histogram(bins=30) +
  theme_minimal() +
  labs(title="Distribution of Income",
       x="Income per capita", y="Count")

print(p1); print(p2)
save_plot("2_1b_hist_infant.png", p1)
save_plot("2_1b_hist_income.png", p2)

# COMMENT (2.1b):
# Both variables are right-skewed. Most countries have relatively low infant mortality
# and low-to-moderate income levels, while a small number of observations appear at very high values.

# Scatter level
p3 <- ggplot(df, aes(income, infant, color=region)) +
  geom_point(alpha=.8) +
  theme_minimal() +
  labs(title="Infant Mortality vs Income",
       x="Income", y="Infant mortality")

print(p3)
save_plot("2_1c_scatter_levels.png", p3)

# COMMENT (2.1c):
# The scatter plot shows us strong negative relationship between income and infant mortality.
# Countries with higher income levels tend to have substantially lower infant mortality rates.
# Africans are clustered at lower income levels with higher mortality,
# while european countries appears at higher income levels with much lower mortality.

# Scatter log-log
p4 <- ggplot(df, aes(log(income), log(infant), color=region)) +
  geom_point(alpha=.8) +
  theme_minimal() +
  labs(title="Log-Log Relationship",
       x="log(Income)", y="log(Infant mortality)")

print(p4)
save_plot("2_1d_scatter_loglog.png", p4)

# COMMENT (2.1d):
# The log-log transformation makes the relationship more linear.
# The negative relationship between income and infant mortality becomes clearer,
# and the spread of observations is more evenly distributed across the plot,
# suggesting the log and log specification may be more appropriate for modeling.

############################
# 2.2 Comparing specifications
############################

m1 <- lm(infant ~ income, data=df)
m2 <- lm(log(infant) ~ log(income), data=df)

summary(m1)
summary(m2)

b1 <- coef(m1)["income"]
effect_1000 <- 1000*b1
print(effect_1000)

# COMMENT (2.2c – m1):
# In the level-level model, a $1,000 increase in income is associated with
# a change of effect_1000 in infant mortality, based on the coefficient.

b2 <- coef(m2)["log(income)"]
print(b2)

# COMMENT (2.2c – m2):
# The coefficient on log(income) represents an elasticity.
# A ten percent increase in income is associated with approximately (10*b2)% changes.
# in infant mortality. The negative coefficient indicates that higher income is associated with lower infant mortality.

# Residual plots
p_r1 <- ggplot(data.frame(f=fitted(m1), r=resid(m1)),
               aes(f,r)) + geom_point() + geom_hline(yintercept=0) +
  theme_minimal() +
  labs(title="Residuals vs Fitted (m1)")

p_r2 <- ggplot(data.frame(f=fitted(m2), r=resid(m2)),
               aes(f,r)) + geom_point() + geom_hline(yintercept=0) +
  theme_minimal() +
  labs(title="Residuals vs Fitted (m2)")

print(p_r1); print(p_r2)
save_plot("2_2d_resid_m1.png", p_r1)
save_plot("2_2d_resid_m2.png", p_r2)

# COMMENT (2.2d):
# The residual plot for the level-level model shows greater curvature and
# heteroskedastisity while the log-log model displays a more random scatter
# of residuals around zero. This suggests the log-log specification provides
# a better functional form for the relationship.

############################
# 2.3 Multiple regression with controls
############################

if("Europe" %in% levels(df$region))
  df$region <- relevel(df$region, ref="Europe")

if("no" %in% levels(df$oil))
  df$oil <- relevel(df$oil, ref="no")

m3 <- lm(log(infant) ~ log(income) + region + oil, data=df)
summary(m3)

b3 <- coef(m3)["log(income)"]
print(b3)

# COMMENT (2.3b):
# After controlling for region and oil exporting status, the coefficient on
# log(income) remains negative, indicating that higher income continues to
# associated with lower infant mortality. The effect of the coefficient
# represents the income elasticity of infant mortality after accounting  for regional differences.

if("regionAfrica" %in% names(coef(m3))){
  africa_pct <- (exp(coef(m3)["regionAfrica"])-1)*100
  print(africa_pct)
  # COMMENT (2.3c):
  # The coefficient for africa indicates that, holding income and oil status constant,
  # African countries have substantially higher infant mortality compared to the
  # reference region such as Europe.

ame_m3 <- avg_slopes(m3, variables="income")
print(ame_m3)

############################
# 2.4 Interaction model
############################

m4 <- lm(log(infant) ~ log(income)*oil + region, data=df)
summary(m4)

me_by_oil <- avg_slopes(m4, variables="income", by="oil")
print(me_by_oil)

# COMMENT (2.4c):
# The marginal effects indicate that the relationship between income and infant
# mortality differs between oil-exporting and non-oil countries. In oil-exporting
# countries, increases in income may have a weaker or different effect on reducing
# infant mortality compared to non-oil countries. This reflects differences
# in how resource wealth translates into public health improvements.

# REQUIRED: plot_slopes()

p_slopes <- plot_slopes(
  m4,
  variables = "income",
  condition = "oil"
) +
  theme_minimal() +
  labs(title="Marginal Effect of Income by Oil Status")

print(p_slopes)
save_plot("2_4d_plot_slopes.png", p_slopes, 7,4)

############################
# 2.5 Predicted values
############################

pred <- predictions(
  m3,
  newdata = datagrid(
    income = c(1000,20000,10000),
    region = c("Africa","Europe","Americas"),
    oil    = c("no","no","yes")
  )
)

pred$infant_hat <- exp(pred$estimate)
print(pred)

gap <- pred$infant_hat[1] - pred$infant_hat[2]
print(gap)

# COMMENT (2.5b):
# Predicted infant mortality for a low-income African country is much higher
# than for a high-income European country.

############################
# 2.6 Publication-quality visualization
############################

p_final <- plot_predictions(
  m3,
  condition = c("income","region")
) +
  theme_minimal() +
  labs(
    title="Predicted Infant Mortality by Income and Region",
    x="Income",
    y="Predicted infant mortality"
  )

print(p_final)
save_plot("2_6_prediction_plot.png", p_final,8,5)

# COMMENT (2.6b):
# The prediction plot shows a clear negative relationship between income and infant mortality
# as countries become wealthier, predicted infant mortality declines. The decline is steepest
# at low income levels and flattens as income rises, which means diminishing returns to income
# in improving infant survival. Africa has always consistently higher
# predicted mortality than other regions across the income range, while Europe remains the lowest.
# The Americas and Asia lie in between, indicating persistent regional differences beyond income.
# This pattern likely reflects differences in infrastructure, environments, public
# health etc.

############################
# 2.7 Robust inference
############################

p_diag <- ggplot(data.frame(f=fitted(m3), r=resid(m3)),
                 aes(f,r)) +
  geom_point() +
  geom_hline(yintercept=0) +
  theme_minimal() +
  labs(
    title = "Residuals vs Fitted (m3)",
    x = "Fitted values",
    y = "Residuals"
  )

print(p_diag)

if(exists("p_diag")){
  save_plot("2_7_resid_m3.png", p_diag)
}

modelsummary(
  list("Level"=m1,"Log-Log"=m2,"Controls"=m3,"Interaction"=m4),
  vcov="robust",
  stars=TRUE,
  gof_map=c("r.squared","nobs")
)

modelsummary(list("Default SE"=m3))
modelsummary(list("Robust SE"=m3), vcov="robust")

# COMMENT (2.7):
# Robust standard errors adjust for potential heteroskedasticity in the data.
# While the coefficient estimates remain the same, the standard errors may
# change Using robust standard errors ensures that inference is valid even if the variance of residuals
# is not constant among observations.
############################################################


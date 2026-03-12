## Importing the seminr package
library(seminr)

## Reading the data file
df <- read.csv("all_climate.csv", header = T)

## Creating the measurement model
mm <- constructs(
  composite("RH", c("RH_JFM_1981_2010","RH_DJF_1991_2020","RH_JFM_1991_2020","RH_NDJ_1991_2020"), weights = mode_A),
  composite("TMAX", c("Tx_DJF_1981_2010", "Tx_JFM_1981_2010","Tx_NDJ_1981_2010","Tx_DJF_1991_2020","Tx_JFM_1991_2020","Tx_NDJ_1991_2020"), weights = mode_A),
  composite("TMIN", c("Tn_DJF_1981_2010","Tn_JFM_1981_2010","Tn_NDJ_1981_2010","Tn_DJF_1991_2020","Tn_JFM_1991_2020","Tn_NDJ_1991_2020")),
  composite("WIND", c("WS_DJF_1981_2010","WS_JFM_1981_2010","WS_NDJ_1981_2010","WS_DJF_1991_2020","WS_JFM_1991_2020","WS_NDJ_1991_2020")),
  composite("RAIN", c("RF_DJF_1981_2010","RF_NDJ_1981_2010","RF_DJF_1991_2020","RF_NDJ_1991_2020"))
)
## Creating the relationships
sm <- relationships(
  paths(from = c("TMIN","TMAX","WIND","RAIN"), to = c("RH")),
  paths(from = c("TMIN","TMAX","WIND"), to = c("RAIN")),
  paths(from = c("TMAX"), to = c("TMIN")),
  paths(from = "WIND", to = c("TMIN","TMAX"))
)
## Estimating the model
pls_model <- estimate_pls(data = df, measurement_model = mm, structural_model = sm)
summary_model <- summary(pls_model)
plot(pls_model)
boot_model <- bootstrap_model(seminr_model = pls_model, nboot = 10000, cores = NULL, seed = 123)
summary_boot <- summary(boot_model, alpha = 0.05)
# T = T-stat.
#p_value <- 2 * (1 - pnorm(abs(T)))

## Mediation analysis
specific_effect_significance(boot_model, from = "WIND", through = "TMAX", to = "TMIN")
specific_effect_significance(boot_model, from = "WIND", through = c("TMAX","TMIN"), to = "RAIN")

## Moderation analysis
mm2 <- constructs(
  composite("RH", c("RH_JFM_1981_2010","RH_DJF_1991_2020","RH_JFM_1991_2020","RH_NDJ_1991_2020"), weights = mode_A),
  composite("TMAX", c("Tx_DJF_1981_2010", "Tx_JFM_1981_2010","Tx_NDJ_1981_2010","Tx_DJF_1991_2020","Tx_JFM_1991_2020","Tx_NDJ_1991_2020"), weights = mode_A),
  composite("TMIN", c("Tn_DJF_1981_2010","Tn_JFM_1981_2010","Tn_NDJ_1981_2010","Tn_DJF_1991_2020","Tn_JFM_1991_2020","Tn_NDJ_1991_2020")),
  composite("WIND", c("WS_DJF_1981_2010","WS_JFM_1981_2010","WS_NDJ_1981_2010","WS_DJF_1991_2020","WS_JFM_1991_2020","WS_NDJ_1991_2020")),
  composite("RAIN", c("RF_DJF_1981_2010","RF_NDJ_1981_2010","RF_DJF_1991_2020","RF_NDJ_1991_2020")),
  interaction_term(iv = "WIND", moderator = "RH", method = two_stage),
  interaction_term(iv = "TMAX", moderator = "RH", method = two_stage),
  interaction_term(iv = "TMIN", moderator = "RH", method = two_stage)
)
sm2 <- relationships(
  # Main effects
  paths(from = c("WIND"), to = c("TMAX", "TMIN", "RAIN")),
  paths(from = c("TMAX"), to = c("TMIN", "RAIN")),
  paths(from = c("TMIN"), to = c("RAIN")),
  paths(from = c("RH"), to = c("TMAX","TMIN","RAIN")),
  # Moderation paths (interaction → DV)
  paths(from = "WIND*RH", to = c("TMAX","RAIN","TMIN")),
  paths(from = "TMAX*RH", to = c("TMIN","RAIN")),
  paths(from = "TMIN*RH", to = "RAIN")
)
model_mod <- estimate_pls(data = df, measurement_model = mm2, structural_model = sm2)
summary_model_mod <- summary(model_mod)
plot(model_mod)
boot_model_mod <- bootstrap_model(seminr_model = model_mod, nboot = 10000, cores = parallel::detectCores(), seed = 123)
summary_boot_mod <- summary(boot_model_mod, alpha = 0.05)
slope_analysis(moderated_model = model_mod, dv = "TMIN", moderator = "RH", iv = "WIND", leg_place = "bottomright")

## Predictive analysis
predict_model <- predict_pls(model = pls_model,technique = predict_DA, noFolds = 10, reps = 10)
sum_predict <- summary(predict_model)
##############
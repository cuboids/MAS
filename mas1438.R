install.packages("RCurl")

library(lattice)
library(magrittr)
library(nlme)
library(RCurl)
library(tidyverse)
#library(lme4)

URL <- "https://raw.githubusercontent.com/cuboids/MAS/main/DataMAS.txt"
dataMAS <- getURL(URL)

col_names <- c("dep", "pen_in_dep", "group", "animal_in_pen", "time", "weight")
df <- read_table(dataMAS, col_names = col_names, col_types = "ffffdd")
df %<>% mutate(penID = as.factor(paste0(dep, group, pen_in_dep)),
               animalID = as.factor(paste0(penID, animal_in_pen)),
               logweight = log(weight),
               time_f = as.factor(df$time),
               invsqrt = weight ** -0.5)
df %>% head()
summary(df)

# Advanced model with sqrt transformation with random intercept
model.1 <- lme(sqrt(weight) ~ dep + group + time + group*time, 
               random = ~ 1|penID/animalID, data = df)
summary(model.1)

# Advanced model with log transformation with random intercept
model.2 <- lme(logweight ~ dep + group + poly(time,2) + group*poly(time,2) +
               dep*poly(time, 2), 
               random = ~ time|penID/animalID, data = df)
summary(model.2)

plot(model.2, resid(., type = "n") ~ fitted(.),
  type = c("p", "smooth"), lwd = 3)
plot(model.2, resid(., type = "n") ~ fitted(.) | group,
  type = c("p", "smooth"), lwd = 3)
plot(model.2, resid(., type = "n") ~ time,
  type = c("p", "smooth"), lwd = 3)
plot(model.2, resid(., type = "n") ~ time | group,
  type = c("p", "smooth"), lwd = 3)

# Advanced model with log transformation with random intercept and random slopes
model.3 <- lme(logweight ~ dep + group + poly(time,2) + group*poly(time,2) +
                dep*poly(time,2), 
               random = ~ time|penID/animalID, data = df)
summary(model.3)

# Advanced model with log transformation with random intercept and random slopes - weight = varIdent
model.4 <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) +
                dep*poly(time, 2), 
              weights = varIdent(form = ~time),
               random = ~ time|penID/animalID, data = df)
summary(model.4)

# Advanced model with log transformation with random intercept and random slopes - correlation = AR-1
model.5 <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) +
                dep*poly(time, 2), 
                correlation = corAR1(form = ~ time|penID/animalID),
               random = ~ time|penID/animalID, data = df)
summary(model.5)

plot(model.3, resid(., type = "n") ~ fitted(.),
  type = c("p", "smooth"), lwd = 3)
plot(model.3, resid(., type = "n") ~ fitted(.) | group,
  type = c("p", "smooth"), lwd = 3)
plot(model.3, resid(., type = "n") ~ time,
  type = c("p", "smooth"), lwd = 3)
plot(model.3, resid(., type = "n") ~ time | group,
  type = c("p", "smooth"), lwd = 3)

plot(model.4, resid(., type = "n") ~ fitted(.),
  type = c("p", "smooth"), lwd = 3)
plot(model.4, resid(., type = "n") ~ fitted(.) | group,
  type = c("p", "smooth"), lwd = 3)
plot(model.4, resid(., type = "n") ~ time,
  type = c("p", "smooth"), lwd = 3)
plot(model.4, resid(., type = "n") ~ time | group,
  type = c("p", "smooth"), lwd = 3)

model.2ML <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) + dep*poly(time, 2), 
               random = ~time|penID/animalID, data = df, 
               method = "ML")
summary(model.3ML)

model.3ML <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) +
                dep*poly(time, 2), 
               random = ~1|penID/animalID, data = df, method = "ML")

summary(model.3ML)

model.4ML <- lme(log(weight) ~ dep + group + poly(time,2) + group*poly(time,2) +
                dep*poly(time,2), 
                weights = varIdent(form = ~time),
               random = ~ time|penID/animalID, data = df, method = "ML", control = lmeControl(opt = "optim"))

summary(model.4ML)

"""# Interaction & Correlation"""

# Test significance of interaction terms
ctrl <- lmeControl(opt='optim')

model.no.interaction <- lme(logweight ~ dep + group + poly(time,2), 
               random = list(penID = pdDiag(~time), animalID = pdDiag(~time)), data = df, method = 'ML')
model.group.interaction <- lme(logweight ~ dep + group + poly(time,2) + group*poly(time,2), 
               random = list(penID = pdDiag(~time), animalID = pdDiag(~time)), data = df, method = 'ML')
model.both.interaction <- lme(logweight ~ dep + group + poly(time,2) + group*poly(time,2) + dep*poly(time,2), 
              random = list(penID = pdDiag(~time), animalID = pdDiag(~time)), data = df, method = 'ML')

anova(model.no.interaction, model.group.interaction, model.both.interaction)

# Test correlation matrix
model.unstructured <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) + dep*poly(time, 2), 
                correlation = corSymm(form = ~ 1|penID/animalID),
                random = list(penID = ~ 1, animalID = pdDiag(~time)), data = df,
                control = lmeControl(opt = "optim"))

model.advanced <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) + dep*poly(time, 2), 
                random = list(penID = ~ 1, animalID = pdDiag(~time)), data = df)

model.compound <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) + dep*poly(time, 2), 
                correlation = corCompSymm(form = ~ 1|penID/animalID),
                random = list(penID = ~ 1, animalID = pdDiag(~time)), data = df)

model.AR1 <- lme(logweight ~ dep + group + poly(time, 2) + group*poly(time, 2) + dep*poly(time, 2), 
                correlation = corAR1(form = ~ 1|penID/animalID),
                random = list(penID = ~ 1, animalID = pdDiag(~time)), data = df)

anova(model.unstructured, model.compound)
anova(model.unstructured, model.AR1)
anova(model.unstructured, model.advanced)

"""# Testing for Random Effects

Test random intercepts. RLRT = REML-based Likelihood Ratio Test.
"""

formula = logweight ~ (dep + group) * poly(time, 2)

m_linear = lm(formula, data = df)

m_ri3 = lme(formula, random = ~1 | penID, data = df)
m_ri3_rs3 = lme(formula, random = list(penID = pdDiag(~time)), data = df)

m_ri2 = lme(formula, random = ~1 | animalID, data = df)
m_ri2_rs2 = lme(formula, random = list(animalID = pdDiag(~time)), data = df)

m_ri23 = lme(formula, random = ~1 | penID/animalID, data = df)
m_ri23_rs23 = lme(formula, random = list(penID = pdDiag(~time), animalID = pdDiag(~time)), data = df)

m_ri23_rs2 = lme(formula, random = list(penID = ~1, animalID = pdDiag(~time)), data = df)

m_opt = m_ri23_rs2

compare <- function(m0, mf, chi2_df) {
    RLRTstat <- -2 * as.numeric(logLik(m0, REML = TRUE) - logLik(mf, REML = TRUE))
    p_value <- 0.5 * pchisq(RLRTstat, chi2_df[[1]], lower.tail = FALSE) +
               0.5 * pchisq(RLRTstat, chi2_df[[2]], lower.tail = FALSE)
    print(c("RLRT:" = RLRTstat, "p-value:" = p_value))
} 

compare(m_ri23_rs3, m_ri23_rs23, c(0, 1))

m_test <- update(m_opt, correlation = corAR1(form = ~ time|penID/animalID))
compare(m_opt, m_test, c(1,0))

m_test <- update(m_opt, correlation = corAR1(form = ~ 1|penID/animalID))
compare(m_opt, m_test, c(1,0))


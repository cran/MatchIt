## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE,
                      fig.width=7, fig.height=5)
options(width = 200)

notice <- "Note: if the `optmatch` package is not available, the subsequent lines will not run."
use <- {
  if (requireNamespace("optmatch", quietly = TRUE)) "full"
  else if (requireNamespace("quickmatch", quietly = TRUE)) "quick"
  else "none"
}

me_ok <- requireNamespace("marginaleffects", quietly = TRUE) &&
  requireNamespace("sandwich", quietly = TRUE)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("MatchIt")
data("lalonde")

head(lalonde)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# No matching; constructing a pre-match matchit object
m.out0 <- matchit(treat ~ age + educ + race + married + 
                   nodegree + re74 + re75, data = lalonde,
                 method = NULL, distance = "glm")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Checking balance prior to matching
summary(m.out0)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1:1 NN PS matching w/o replacement
m.out1 <- matchit(treat ~ age + educ + race + married + 
                   nodegree + re74 + re75, data = lalonde,
                 method = "nearest", distance = "glm")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
m.out1

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Checking balance after NN matching
summary(m.out1, un = FALSE)

## ---- fig.alt="Jitter plot of the propensity scores, which shows that no treated unit were dropped, and a large number of control units with low propensity scores were dropped."---------------------
plot(m.out1, type = "jitter", interactive = FALSE)

## ---- fig.alt="Density plots of age, married and re75 in the unmatched and matched samples."----------------------------------------------------------------------------------------------------------
plot(m.out1, type = "density", interactive = FALSE,
     which.xs = ~age + married + re75)

## ---- eval = (use == "full"), include= (use %in% c("full", "none"))-----------------------------------------------------------------------------------------------------------------------------------
# Full matching on a probit PS
m.out2 <- matchit(treat ~ age + educ + race + married + 
                   nodegree + re74 + re75, data = lalonde,
                 method = "full", distance = "glm", link = "probit")
m.out2

## ---- eval = (use == "quick"), include = (use == "quick")---------------------------------------------------------------------------------------------------------------------------------------------
#  # Full matching on a probit PS
#  m.out2 <- matchit(treat ~ age + educ + race + married +
#                     nodegree + re74 + re75, data = lalonde,
#                   method = "quick", distance = "glm", link = "probit")
#  m.out2

## ---- eval = (use != "none")--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Checking balance after full matching
summary(m.out2, un = FALSE)

## ---- eval = (use != "none"), fig.alt = "A love plot with matched dots below the threshold lines, indicaitng good balance after matching, in contrast to the unmatched dots far from the treshold lines, indicating poor balance before matching."----
plot(summary(m.out2))

## ---- eval = (use != "none")--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
m.data <- match.data(m.out2)

head(m.data)

## ---- eval = (use != "none" && me_ok)-----------------------------------------------------------------------------------------------------------------------------------------------------------------
library("marginaleffects")

fit <- lm(re78 ~ treat * (age + educ + race + married + nodegree + 
             re74 + re75), data = m.data, weights = weights)

avg_comparisons(fit,
                variables = "treat",
                vcov = ~subclass,
                newdata = subset(m.data, treat == 1),
                wts = "weights")

## ---- include = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
est <- {
  if (use != "none" && me_ok) {
    avg_comparisons(fit,
                variables = "treat",
                vcov = ~subclass,
                newdata = subset(m.data, treat == 1),
                wts = "weights")
  }
  else data.frame(type = "response", term = "1 - 0", estimate = 2114, 
                  std.error = 646, statistic = 3.27, 
                  p.value = 0.0011, conf.low = 848, 
                  conf.high = 3380)
}


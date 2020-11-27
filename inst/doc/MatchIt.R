## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE,
                      fig.width=7, fig.height=5)
options(width = 200)

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

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(m.out1, type = "jitter", interactive = FALSE)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(m.out1, type = "qq", interactive = FALSE,
     which.xs = c("age", "married", "re75"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Full matching on a probit PS
m.out2 <- matchit(treat ~ age + educ + race + married + 
                   nodegree + re74 + re75, data = lalonde,
                 method = "full", distance = "glm", link = "probit")
m.out2

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Checking balance after full matching
summary(m.out2, un = FALSE)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(summary(m.out2))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
m.data1 <- match.data(m.out1)

head(m.data1)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("lmtest") #coeftest
library("sandwich") #vcovCL

fit1 <- lm(re78 ~ treat + age + educ + race + married + nodegree + 
             re74 + re75, data = m.data1, weights = weights)

coeftest(fit1, vcov. = vcovCL, cluster = ~subclass)

## ---- message = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
m.data2 <- match.data(m.out2)

fit2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + 
             re74 + re75, data = m.data2, weights = weights)

coeftest(fit2, vcov. = vcovCL, cluster = ~subclass)

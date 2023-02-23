## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE,
                      fig.width=7, fig.height=5)
options(width = 200, digits = 4)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("MatchIt")
data("lalonde", package = "MatchIt")

#1:1 NN matching w/ replacement on a logistic regression PS
m.out <- matchit(treat ~ age + educ + race + married + 
                   nodegree + re74 + re75, data = lalonde,
                 replace = TRUE)
m.out

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(m.out, addlvariables = ~ I(age^2) + I(re74==0) + 
          I(re75==0) + educ:race)

## ---- fig.alt="A love plot with most matched dots below the threshold lines, indicaitng good balance after matching, in contrast to the unmatched dots far from the treshold lines, indicating poor balance before matching."----
m.sum <- summary(m.out, addlvariables = ~ I(age^2) + I(re74==0) + 
                   I(re75==0) + educ:race)
plot(m.sum, var.order = "unmatched")

## ---- fig.alt ="eQQ plots of age, nodegree, and re74 in the unmatched and matched samples."-----------------------------------------------------------------------------------------------------------
#eQQ plot
plot(m.out, type = "qq", which.xs = ~age + nodegree + re74)

## ---- fig.alt ="eCDF plots of educ, married, and re75 in the unmatched and matched samples."----------------------------------------------------------------------------------------------------------
#eCDF plot
plot(m.out, type = "ecdf", which.xs = ~educ + married + re75)

## ---- fig.alt ="Density plots of age, educ, and race in the unmatched and matched samples."-----------------------------------------------------------------------------------------------------------
#density plot
plot(m.out, type = "density", which.xs = ~age + educ + race)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Subclassification on a logistic regression PS
s.out <- matchit(treat ~ age + educ + race + married + 
                   nodegree + re74 + re75, data = lalonde,
                 method = "subclass", subclass = 4)
s.out

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(s.out)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(s.out, subclass = TRUE, un = FALSE)

## ---- fig.alt ="Love plot of balance before and after subclassification, with subclass IDs representing balance within each subclass in addition to dots representing balance overall."---------------
s <- summary(s.out, subclass = TRUE)
plot(s, var.order = "unmatched", abs = FALSE)

## ---- fig.alt ="Density plots of educ, married, and re75 in the unmatched sample and in subclass 1."--------------------------------------------------------------------------------------------------
plot(s.out, type = "density", which.xs = ~educ + married + re75,
     subclass = 1)

## ---- include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ok <- requireNamespace("cobalt", quietly = TRUE)

## ---- message = F, eval = ok--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("cobalt")

## ---- eval = ok---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bal.tab(m.out, un = TRUE, stats = c("m", "v", "ks"))

## ---- eval = ok---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Nearest neighbor (NN) matching on the PS
m.out2 <- matchit(treat ~ age + educ + race + married + 
                   nodegree + re74 + re75, data = lalonde)

#Balance on covariates after full and NN matching
bal.tab(treat ~ age + educ + race + married + 
          nodegree + re74 + re75, data = lalonde, 
        un = TRUE, weights = list(full = m.out, nn = m.out2))

## ---- eval = ok, fig.alt ="Minimal love plot of balance before and after matching."-------------------------------------------------------------------------------------------------------------------
love.plot(m.out, binary = "std")

## ---- fig.width=7, eval = ok, fig.alt ="A more elaborate love plot displaying some of the cobalt's capabilities for making publication-ready plots."--------------------------------------------------
love.plot(m.out, stats = c("m", "ks"), poly = 2, abs = TRUE,
          weights = list(nn = m.out2),
          drop.distance = TRUE, thresholds = c(m = .1),
          var.order = "unadjusted", binary = "std",
          shapes = c("circle filled", "triangle", "square"), 
          colors = c("red", "blue", "darkgreen"),
          sample.names = c("Original", "Full Matching", "NN Matching"),
          position = "bottom")

## ---- eval = ok, fig.alt = c("Density plot for educ before and after matching.", "Bar graph for race before and after matching.", "Mirrored histograms of propensity scores before and after matching.")----
#Density plot for continuous variables
bal.plot(m.out, var.name = "educ", which = "both")

#Bar graph for categorical variables
bal.plot(m.out, var.name = "race", which = "both")

#Mirrored histogram
bal.plot(m.out, var.name = "distance", which = "both",
         type = "histogram", mirror = TRUE)


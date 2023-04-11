# -------------------------------------------------------------------------------------------------
# R Codes to reproduce the results on paper:
# glme: An R package for Mixed Effects Model Inference by Generalized Approach

# Created by: Mustafa Cavus in 19/12/2020
# Updated by: Mustafa Cavus in 11/04/2023
# Contact: mustafacavus@eskisehir.edu.tr
# -------------------------------------------------------------------------------------------------
install.packages("glme")
install.packages("nlme")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("BHH2")
install.packages("car")
install.packages("lme4")

library(glme)
library(nlme)
library(dplyr)
library(ggplot2)
library(BHH2) # for retrieving penicillin data
library(car)  # for retrieving blackmore data
library(lme4) # for retrieving sleepystudy data

## 1. USAGE OF glme PACKAGE #######################################################################

# The generic dataset Orthodont is used to show the using of glme(). It is a data frame has 108
# rows and 4 columns of the change in an orthdontic measurement over time for several young
# subjects. Firsty, the following R codes are used to visualize the dataset:
data("Orthodont")

Orthodont |>
  ggplot(aes(x = age, y = distance, color = Sex)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# The glme() is used with intercept, only intercept as in the following.
# run as:
glme(distance ~ age + Sex, data = Orthodont, random = ~ age|Subject, method = "GM")

# without intercept run as:
glme(distance ~ age + Sex, data = Orthodont, random = ~ -1+age|Subject, method = "GM")

# only intercept run as:
glme(distance ~ age + Sex, data = Orthodont, random = ~ 1|Subject, method = "GM")

# It returns the same output of lme() with ML or REML when the method = "ML" or "REML". Also,
# the glme() provides three components are fixed,SD and coefficients.

output <- glme(distance ~ age + Sex, data = Orthodont, random = ~ age | Subject)
output$fixed
output$sd
output$coefficients$fixed
cbind(output$coefficients$random,
      predicted_intercept = output$coefficients$random[,1] + output$coefficients$fixed[1],
      predicted_age = output$coefficients$random[,2] + output$coefficients$fixed[2])

## 2. ILLUSTRATIVE EXAMPLES ########################################################################

# 2.1. Penicillin Data

# Penicillin yield data is the process of manufacturing penicillin under different treatments and
# product blends (Box et al., 2005). The plot of the distribution run vs. yield colored by blend
# is given in Fig 2.
data("penicillin.data")

penicillin.data |>
  ggplot(aes(x = run, y = yield, color = blend)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()

# In this example, yield is explained by the fixed effects of run and treatment and blend is
# random effect with only intercept as in the followings:
glme(yield ~ run + treat, data = penicillin.data, random = ~1 | blend, method = "GM")
glme(yield ~ run + treat, data = penicillin.data, random = ~1 | blend, method = "ML")
glme(yield ~ run + treat, data = penicillin.data, random = ~1 | blend, method = "REML")

# 2.2. Blackmore Data

# The Blackmore data frame has 945 rows and 4 columns. Blackmore and Davis’s data on exercise
# histories of 138 teenaged girls hospitalized for eating disorders and 98 control subjects
# (Fox and Weisberg, 2019). The plot of the distribution age vs. exercise colored by group is
# given in Fig 3.
data("Blackmore")

Blackmore |>
  ggplot(aes(x = age, y = exercise, color = group)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw()

# In this example, exercise is explained by the fixed effects of age and age and subject are
# random effect with only intercept as in the followings:
glme(exercise ~ age, data = Blackmore, ~ 1 + age | subject, method = "GM")
glme(exercise ~ age, data = Blackmore, ~ 1 + age | subject, method = "ML")
glme(exercise ~ age, data = Blackmore, ~ 1 + age | subject, method = "REML")

# 2.3. Sleep Data

# The average reaction time per day for subjects in a sleep deprivation study. On day 0 the
# subjects had their normal amount of sleep. Starting that night they were restricted to 3
# hours of sleep per night. The observations represent the average reaction time on a series
# of tests given each day to each subject (Gregory Belenky et al., 2003). The plot of the
# distribution Days vs. Reaction colored by Subject is given in Fig 4.
data("sleepstudy")

sleepstudy |>
  ggplot(aes(x = Days, y = Reaction, color = Subject)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw()

# In this example, Reaction is explained by the fixed effects of Days, and Days and Subject are
# random effect with only intercept as in the followings:
glme(Reaction ~ Days, data = sleepstudy, ~ 1 + Days | Subject, method = "GM")
glme(Reaction ~ Days, data = sleepstudy, ~ 1 + Days | Subject, method = "ML")
glme(Reaction ~ Days, data = sleepstudy, ~ 1 + Days | Subject, method = "REML")

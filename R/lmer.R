library(papaja)
library(scales)
library(broom)
library(broom.mixed)
library(here)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(merTools)
library(performance)

options(scipen=999)
theme_set(theme_apa(base_size = 14, base_family = "Times New Roman") + 
            theme(axis.title.y = element_text(size = 18), 
                  axis.text.y = element_text(size = 10, 
                                             family = "Arial", colour = "black"),
                  axis.text.x = element_text(size = 12, 
                                             colour = "black", face = "plain"),
                  legend.position = "none", 
                  plot.margin = ggplot2::margin(b = 0.5, l = 0.5, r = 0.4,
                                                t = 0.7, unit = "cm"),
                  axis.line = element_line(size = 0.8),
                  legend.key.size = unit(0.1, "cm"),
                  legend.background=element_blank()))

# load data ---------------------------------------------------------------
path <- here("raw-data", "cing2221.csv")
df <- read_csv(path, na = c("nan", "#DIV/0!"), 
               col_types = "fffiiifffii--ii-iiiiinn-")

# recording ---------------------------------------------------------------
# filter out control and unilateral cingulotomy cases

df_filt <- df %>% dplyr::filter(treatment %in% c("sham", "cingulotomy"))

table(df_filt$opdaybin)

# collapse all pre in 1 bin & all post => 7 in 1 bin = "post7+" 
# because the coverage is too sparse otherwise
df_filt_m <- df_filt %>% dplyr::mutate(opdaybin_r = fct_collapse(opdaybin,
                    pre = c("pre0-1-2", "pre-3-4-5", "pre-6-7-8", "pre-9-10-11"),
                    "post1-3" = "post1+2+3",
                    "post4-6" = "post4+5+6",
                    "post7+" = c("post7+8+9", "post10+11+12", "post13+14+15", 
                               "post16+17+18", "post19+20+21", "post22+23+24",
                               "post28+29+30", "post31+32+33")),
                    treatment = fct_drop(treatment))

# make first to third order orthogonal polynomials of animal's age
# create orthogonal polynomial age variables in data frame
df_filt_m[, paste0("pol_age", 1:3)] <- poly(df_filt_m$pnd, degree = 3)

# analysis ----------------------------------------------------------------
df_filt_m$opdaybin_r <- relevel(df_filt_m$opdaybin_r, "pre")
df_filt_m$treatment <- relevel(df_filt_m$treatment, "sham")

# control for recording bin
model_op <- lmer(percentphee ~ treatment*opdaybin_r + (1|id), data = df_filt_m,
                   REML = FALSE, 
                   control=lmerControl(optimizer="bobyqa")) 

summary(model_op)
nobs(model_op)

# same, but no treatment * op_bin interaction
# fits better than previous model
model_op_noint <- lmer(percentphee ~ treatment + opdaybin_r + (1|id), 
                 data = df_filt_m,
                 REML = FALSE, 
                 control=lmerControl(optimizer="bobyqa")) 

summary(model_op_noint)

# add animal's age to the data
# fits better than previous model
model_op_noint_age <- lmer(percentphee ~ treatment + opdaybin_r + pnd + (1|id), 
                       data = df_filt_m,
                       REML = FALSE, 
                       control=lmerControl(optimizer="bobyqa")) 

summary(model_op_noint_age) 
nobs(model_op_noint_age)

# add animal's `opcondition` to the data
# fits better than previous model
model_age_opcond <- lmer(percentphee ~ treatment*opcondition +
                           pnd + (1|id), 
                           data = df_filt_m,
                           REML = FALSE, 
                           control=lmerControl(optimizer="bobyqa")) 

summary(model_age_opcond) 
nobs(model_age_opcond)

# make animal's age also a random effect
# fits better than previous model
model_age_opcond_rand <- lmer(percentphee ~ treatment*opcondition +
                            pnd + (pnd|id), 
                            data = df_filt_m,
                            REML = FALSE, 
                            control=lmerControl(optimizer="bobyqa")) 

summary(model_age_opcond_rand) 
nobs(model_age_opcond_rand)


# add a second degree polynomials of animal's age
# fits better than previous model
model_age_opcond_pol <- lmer(percentphee ~ treatment*opcondition +
                              pol_age1 + pol_age2 + ((pol_age1 + pol_age2)|id), 
                              data = df_filt_m,
                              REML = FALSE, 
                              control=lmerControl(optimizer="bobyqa")) 

summary(model_age_opcond_pol) 
nobs(model_age_opcond_pol)

# model with no interaction between `opcondition` and treatment
model_age_opcond_pol_noint <- lmer(percentphee ~ opcondition + treatment + 
                                   pol_age1 + pol_age2 +
                              ((pol_age1 + pol_age2)|id), 
                             data = df_filt_m,
                             REML = FALSE, 
                             control=lmerControl(optimizer="bobyqa")) 

summary(model_age_opcond_pol_noint)
nobs(model_age_opcond_pol_noint)

# model performance -------------------------------------------------------
anova(model_op_noint, model_op_noint_age, model_age_opcond, 
      model_age_opcond_rand, model_age_opcond_pol, model_age_opcond_pol_noint)
check_model(model_age_opcond_pol_noint)







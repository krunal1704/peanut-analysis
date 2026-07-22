library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(AICcmodavg)
library(DescTools)
set.seed(1234)

#identifying outliers if any
all_data_for_anova_analysis %>% sample_n_by(group, size = 15)
#all_data_for_anova_analysis [,c(1,7)] %>% sample_n_by(group, size = 15)

all_data_for_anova_analysis %>% group_by(group) %>% get_summary_stats(plantparameter, type = "mean_se") #change me (plantparameter)

ggboxplot(all_data_for_anova_analysis, x = "group", y = "plantparameter") #change me

all_data_for_anova_analysis%>% group_by(group) %>% identify_outliers(plantparameter) #change me

model  <- lm(plantparameter ~ group, data = all_data_for_anova_analysis)

ggqqplot(residuals(model))


#Shairo-Wilk test and levene's test
shapiro_test(residuals(model))

all_data_for_anova_analysis%>% group_by(group) %>% shapiro_test(plantparameter)# change me

ggqqplot(all_data_for_anova_analysis, "plantparameter", facet.by = "group")

plot(model, 1)

all_data_for_anova_analysis %>% levene_test(plantparameter ~ group)

#ANOVA test
res.aov <-all_data_for_anova_analysis%>% anova_test(plantparameter ~ group)
res.aov


#Dunnet-test
model <- aov(plantparameter ~ group, data = all_data_for_anova_analysis)
summary (model)
DunnettTest(x=all_data_for_anova_analysis$plantparameter, g=all_data_for_anova_analysis$group)



#Tukey hsd
res.aov <-all_data_for_anova_analysis%>% anova_test(plantparameter ~ group)
res.aov
model <- aov(plantparameter ~ group, data = all_data_for_anova_analysis)
summary (model)
pwc <-all_data_for_anova_analysis  %>% tukey_hsd(plantparameter ~ group)
pwc%>% print(n=78)


#Graphical presentation
pwc <- pwc %>% add_xy_position(x = "group", ref.group = "control")
ggboxplot(all_data_for_anova_analysis, x = "group", y = "plantparameter") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )



#box plot comparison 
ggbetweenstats(all_data_for_anova_analysis, x = "group", y = "plantparameter", 
               type = "parametric", plot.type = "box", 
               pairwise.comparisons = TRUE, 
               pairwise.display = "significant", centrality.plotting = FALSE, bf.message = FALSE)




# Install packages if necessary
library(readxl)
library(dplyr)
library(FSA)
library(rcompanion)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(AICcmodavg)
library(car)
library(multcompView)
library(DescTools)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(patchwork)
library(writexl)
library(svglite)

set.seed(1234)

#statistical analysis of PGP activity
# Step 1: Import Excel data
data <- read_excel(file.choose(), sheet = 1)

# Rename first two columns if needed
names(data)[1:2] <- c("Isolate", "Value")

data$Isolate <- factor(data$Isolate)
data$Value <- as.numeric(data$Value)

# Step 2: Mean and standard deviation
summary_table <- data %>%
  group_by(Isolate) %>%
  summarise(
    n = n(),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

# Step 3: Kruskal-Wallis test
kw_result <- kruskal.test(Value ~ Isolate, data = data)
print(kw_result)

# Step 4: Dunn's test with Holm correction
dunn_result <- FSA::dunnTest(
  Value ~ Isolate,
  data = data,
  method = "holm"
)

print(dunn_result$res)

# Step 5: Significance letters
letters_result <- rcompanion::cldList(
  P.adj ~ Comparison,
  data = dunn_result$res,
  threshold = 0.05
)

print(letters_result)

# Step 6: Combine mean, SD and letters
final_table <- summary_table %>%
  left_join(
    letters_result,
    by = c("Isolate" = "Group")
  ) %>%
  mutate(
    Mean_SD = sprintf("%.2f ± %.2f", Mean, SD)
  ) %>%
  select(Isolate, n, Mean_SD, Letter)

print(final_table)

# Step 7: Export results
write.csv(
  final_table,
  "Kruskal_Dunn_Significance_Letters.csv",
  row.names = FALSE
)

write.csv(
  dunn_result$res,
  "Dunn_Holm_Pairwise_Results.csv",
  row.names = FALSE
)

### ==================================================###

#optimization of IAA production 
data <- read_excel(file.choose())

data$Isolate <- factor(data$Isolate)

results <- list()

for (temp in unique(data$Temperature)) {
  for (conc in unique(data$Concentration)) {
    
    subset_data <- data %>%
      filter(
        Temperature == temp,
        Concentration == conc
      )
    
    # Overall Kruskal-Wallis test
    kw <- kruskal.test(
      Value ~ Isolate,
      data = subset_data
    )
    
    # Dunn-Holm pairwise comparisons
    dunn <- FSA::dunnTest(
      Value ~ Isolate,
      data = subset_data,
      method = "holm"
    )
    
    # Significance letters
    letters <- rcompanion::cldList(
      P.adj ~ Comparison,
      data = dunn$res,
      threshold = 0.05
    )
    
    # Mean and SD
    summary <- subset_data %>%
      group_by(Isolate) %>%
      summarise(
        Mean = mean(Value),
        SD = sd(Value),
        .groups = "drop"
      ) %>%
      left_join(
        letters,
        by = c("Isolate" = "Group")
      )
    
    key <- paste(temp, conc, sep = "_")
    
    results[[key]] <- list(
      Kruskal_Wallis = kw,
      Dunn_Holm = dunn$res,
      Summary = summary
    )
    
    cat("\nTemperature:", temp,
        "Concentration:", conc, "\n")
    print(kw)
    print(summary)
  }
}

### ==========================================================###

#ACC deaminase activity analysis
# Step 1: Import Excel file
data <- read_excel(file.choose(), sheet = 1)

# Rename first two columns if necessary
names(data)[1:2] <- c("Isolate", "Value")

data$Isolate <- factor(data$Isolate)
data$Value <- as.numeric(data$Value)

# Step 2: Calculate mean and SD
summary_table <- data %>%
  group_by(Isolate) %>%
  summarise(
    n = n(),
    Mean = mean(Value),
    SD = sd(Value),
    .groups = "drop"
  )

print(summary_table)

# Step 3: Perform one-way ANOVA
model <- aov(Value ~ Isolate, data = data)

summary(model)

# Step 4: Test normality of ANOVA residuals
shapiro.test(residuals(model))

# Step 5: Test homogeneity of variances
leveneTest(Value ~ Isolate, data = data)

# Step 6: Inspect residual Q-Q plot
qqnorm(residuals(model),
       main = "ACC deaminase: Residual Q-Q Plot")
qqline(residuals(model), col = "red")

# Step 7: Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(model)

print(tukey_result)

# Step 8: Generate significance letters
letters_result <- multcompView::multcompLetters4(
  model,
  tukey_result
)

letters_df <- data.frame(
  Isolate = names(letters_result$Isolate$Letters),
  Letter = letters_result$Isolate$Letters
)

# Step 9: Combine mean, SD and letters
final_table <- summary_table %>%
  left_join(letters_df, by = "Isolate") %>%
  mutate(
    Mean_SD = sprintf("%.2f ± %.2f", Mean, SD)
  ) %>%
  select(Isolate, n, Mean_SD, Letter)

print(final_table)

# Step 10: Export results
write.csv(
  final_table,
  "ACC_ANOVA_Tukey_Letters.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(tukey_result$Isolate),
  "ACC_Tukey_Pairwise_Results.csv"
)

###===========================================================###


#Pot data analysis
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

#============================================================================================

# 120 DAS pot data 
# Six-panel bar graph with Dunnett's significance

# Enter data
# Each parameter has Mean, SE, and Dunnett significance

Treatment <- c("Control", paste0("T", 1:12))

#replace the data and significance values as per data.
data <- data.frame(
  Treatment = Treatment,
  
  Root_mean = c(
    22.58,41.16,33.79,54.57,43.99,33.12,31.09,
    43.83,39.00,44.96,38.38,44.39,53.36),
  Root_SE = c(
    .754,1.789,.872,3.453,1.00,1.022,1.328,
    1.054,1.101,1.35,3.756,1.823,1.13),
  Root_sig = c(
    "","***","*","***","***","ns","ns",
    "***","***","***","***","***","***"),
  
  Fresh_mean = c(
    5.09,8.71,7.85,8.70,8.20,6.46,5.55,
    7.73,9.57,10.05,7.22,7.65,8.79),
  Fresh_SE = c(
    .316,.621,.559,.239,.200,.846,.247,
    .389,.204,.636,1.721,.405,.162),
  Fresh_sig = c(
    "","***","ns","ns","***","ns","ns",
    "*","*","**","ns","ns","ns"),
  
  Total_mean = c(
    74.60,79.95,88.33,104.54,94.08,80.40,83.81,
    89.44,88.06,92.93,97.93,109.62,99.51),
  Total_SE = c(
    1.86,2.52,2.69,5.22,3.04,1.60,2.72,
    3.28,3.38,2.52,2.91,4.31,1.77),
  Total_sig = c(
    "","ns","ns","ns","**","ns","ns",
    "ns","ns","ns","ns","**","**"),
  
  Dry_mean = c(
    1.32,2.76,2.42,3.02,2.34,1.90,1.54,
    3.07,2.95,2.43,2.20,2.27,3.10),
  Dry_SE = c(
    .139,.175,.102,.116,.284,.059,.035,
    .095,.118,.164,.150,.147,.066),
  Dry_sig = c(
    "","***","**","ns","***","**","ns",
    "***","**","***","*","ns","***"),
  
  Shoot_mean = c(
    53.41,37.44,56.11,48.45,48.70,49.32,54.80,
    44.78,49.64,48.73,58.13,64.89,46.45),
  Shoot_SE = c(
    1.491,2.894,3.613,1.267,1.689,2.556,3.050,
    3.635,3.162,2.098,1.921,2.852,.344),
  Shoot_sig = c(
    "","*","ns","**","**","ns","ns",
    "***","*","*","ns","ns","**"),
  
  Pod_mean = c(
    .976,2.94,1.87,2.56,1.61,.96,1.01,
    2.87,2.84,1.68,1.24,.83,2.91),
  Pod_SE = c(
    .163,.063,.152,.309,.166,.291,.152,
    .281,.180,.085,.126,.219,.203),
  Pod_sig = c(
    "","***","*","ns","***","ns","ns",
    "***","***","ns","ns","ns","***")
)

# Convert data into long format

long_data <- bind_rows(
  lapply(
    c("Root","Fresh","Total","Dry","Shoot","Pod"),
    function(p) {
      data.frame(
        Treatment = data$Treatment,
        Parameter = p,
        Mean = data[[paste0(p, "_mean")]],
        SE = data[[paste0(p, "_SE")]],
        Significance = data[[paste0(p, "_sig")]]
      )
    }
  )
)

# Set panel order and labels

parameter_labels <- c(
  Root = "Root length (cm)",
  Fresh = "Fresh weight (g)",
  Total = "Total length (cm)",
  Dry = "Dry weight (g)",
  Shoot = "Shoot length (cm)",
  Pod = "Pod weight (g)"
)

long_data$Parameter <- factor(
  long_data$Parameter,
  levels = names(parameter_labels),
  labels = parameter_labels
)

long_data$Treatment <- factor(
  long_data$Treatment,
  levels = Treatment
)

# Position significance annotations

long_data <- long_data %>%
  group_by(Parameter) %>%
  mutate(
    panel_max = max(Mean + SE),
    label_y = Mean + SE + 0.025 * panel_max
  ) %>%
  ungroup()

# Create the six-panel figure

p <- ggplot(
  long_data,
  aes(x = Treatment, y = Mean)
) +
  
  geom_col(
    fill = "#2878B5",
    colour = "black",
    linewidth = 0.35,
    width = 0.75
  ) +
  
  geom_errorbar(
    aes(
      ymin = pmax(0, Mean - SE),
      ymax = Mean + SE
    ),
    width = 0.20,
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(
      y = label_y,
      label = Significance
    ),
    size = 2.7,
    fontface = "bold",
    vjust = 0
  ) +
  
  facet_wrap(
    ~Parameter,
    nrow = 3,
    ncol = 2,
    scales = "free_y"
  ) +
  
  # Independent Y-axis scales with room for annotations
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  
  labs(
    title = paste(
      "Effect of PGP rhizobacterial treatments on",
      "peanut growth in non-sterile soil at 120 DAS",
      "(second year)"
    ),
    x = NULL,
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 13,
      margin = margin(b = 15)
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    strip.background = element_blank(),
    
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      size = 8,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.y = element_line(
      colour = "grey85",
      linetype = "dashed",
      linewidth = 0.3
    ),
    
    panel.spacing = grid::unit(1.5, "lines"),
    
    plot.margin = margin(15, 15, 15, 15)
  )

# Display figure

print(p)

# 7. Save high-resolution figure

ggsave(
  filename = "Second_year_nonsterile_120DAS.png",
  plot = p,
  width = 14,
  height = 15,
  units = "in",
  dpi = 300,
  bg = "white"
)

###===============================================###

#PCA plot

# Arrange experimental conditions
panel_order <- c(
  "Sterile, 45 DAS",
  "Non-sterile, 45 DAS",
  "Sterile, 120 DAS",
  "Non-sterile, 120 DAS"
)

df$Panel <- factor(
  df$Panel,
  levels = panel_order
)

# Arrange treatments
df$Treatment <- factor(
  df$Treatment,
  levels = c("C", paste0("T", 1:12))
)

# Check that each panel contains 13 treatments
stopifnot(all(table(df$Panel) == 13))


# Variables for each experimental condition
pca_variables <- list(
"Sterile, 45 DAS" = c(
  "Germination", "VI",
  "RL", "SL", "FW", "DW",
  "Root_bacteria"
),

"Non-sterile, 45 DAS" = c(
  "Germination", "VI",
  "RL", "SL", "FW", "DW"
),

"Sterile, 120 DAS" = c(
  "Germination", "VI",
  "RL", "SL", "FW", "DW", "PW"
),

"Non-sterile, 120 DAS" = c(
  "Germination", "VI",
  "RL", "SL", "FW", "DW", "PW"
)
)
# Function to perform PCA and create biplot

make_pca <- function(data, panel_name, variables, panel_letter) {
  
  # Select data for one condition
  sub_data <- data %>%
    filter(Panel == panel_name) %>%
    arrange(Treatment)
  
  # Select numeric variables
  X <- sub_data %>%
    select(all_of(variables)) %>%
    as.data.frame()
  
  rownames(X) <- sub_data$Treatment
  
  # Check for missing values
  stopifnot(!anyNA(X))
  
  # Perform standardized PCA
  pca <- prcomp(
    X,
    center = TRUE,
    scale. = TRUE
  )
  
  # Calculate explained variance
  variance <- 100 * pca$sdev^2 /
    sum(pca$sdev^2)
  
  # Treatment scores
  scores <- as.data.frame(pca$x[, 1:2])
  
  scores$Treatment <- rownames(scores)
  
  scores$Group <- case_when(
    scores$Treatment == "C" ~ "Control",
    scores$Treatment == "T12" ~ "Consortium",
    TRUE ~ "Individual isolates"
  )
  
  # Variable loadings
  loadings <- as.data.frame(
    pca$rotation[, 1:2]
  )
  
  loadings$Variable <- rownames(loadings)
  
  # Scale arrows for biplot display
  arrow_scale <- 0.75 * min(
    max(abs(scores$PC1)) /
      max(abs(loadings$PC1)),
    max(abs(scores$PC2)) /
      max(abs(loadings$PC2))
  )
  
  loadings$Arrow_PC1 <-
    loadings$PC1 * arrow_scale
  
  loadings$Arrow_PC2 <-
    loadings$PC2 * arrow_scale
  
  # Variable labels
  loadings$Label <- ifelse(
    loadings$Variable == "Root_bacteria",
    "Root bacteria",
    loadings$Variable
  )
  
  # Generate biplot
  p <- ggplot() +
    
    # Reference lines
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey65",
      linewidth = 0.4
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey65",
      linewidth = 0.4
    ) +
    
    # Treatment points
    geom_point(
      data = scores,
      aes(
        x = PC1,
        y = PC2,
        colour = Group,
        shape = Group
      ),
      size = 3
    ) +
    
    # Treatment labels
    geom_text_repel(
      data = scores,
      aes(
        x = PC1,
        y = PC2,
        label = Treatment
      ),
      size = 3.8,
      colour = "black",
      seed = 123,
      max.overlaps = Inf
    ) +
    
    # Loading arrows
    geom_segment(
      data = loadings,
      aes(
        x = 0,
        y = 0,
        xend = Arrow_PC1,
        yend = Arrow_PC2
      ),
      arrow = arrow(
        length = unit(0.22, "cm")
      ),
      colour = "#B22222",
      linewidth = 0.8
    ) +
    
    # Variable labels
    geom_text_repel(
      data = loadings,
      aes(
        x = Arrow_PC1,
        y = Arrow_PC2,
        label = Label
      ),
      colour = "#B22222",
      fontface = "bold",
      size = 3.8,
      seed = 123,
      max.overlaps = Inf
    ) +
    
    # Treatment colours
    scale_colour_manual(
      values = c(
        "Control" = "black",
        "Individual isolates" = "#377EB8",
        "Consortium" = "#2E8B57"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "Control" = 17,
        "Individual isolates" = 16,
        "Consortium" = 18
      )
    ) +
    
    # Axis labels and title
    labs(
      title = paste0(
        panel_letter, ". ", panel_name
      ),
      x = sprintf(
        "PC1 (%.2f%%)", variance[1]
      ),
      y = sprintf(
        "PC2 (%.2f%%)", variance[2]
      ),
      colour = NULL,
      shape = NULL
    ) +
    
    coord_equal() +
    
    theme_classic(base_size = 12) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        colour = "black"
      ),
      legend.position = "bottom",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      plot.margin = margin(
        10, 15, 10, 15
      )
    )
  
  # Explained variance table
  variance_table <- data.frame(
    Panel = panel_name,
    PC1 = variance[1],
    PC2 = variance[2],
    Total = sum(variance[1:2])
  )
  
  # Add panel information
  loadings$Panel <- panel_name
  scores$Panel <- panel_name
  
  # Return all outputs
  return(list(
    plot = p,
    model = pca,
    scores = scores,
    loadings = loadings,
    variance = variance_table
  ))
}

# Run PCA for all four conditions

pca_results <- lapply(
  seq_along(panel_order),
  
  function(i) {
    
    make_pca(
      data = df,
      panel_name = panel_order[i],
      variables = pca_variables[[panel_order[i]]],
      panel_letter = LETTERS[i]
    )
    
  }
)

names(pca_results) <- panel_order

# Extract individual plots
p1 <- pca_results[[1]]$plot
p2 <- pca_results[[2]]$plot
p3 <- pca_results[[3]]$plot
p4 <- pca_results[[4]]$plot

# Arrange as two rows and two columns
final_pca <- (p1 | p2) / (p3 | p4) +
  
  plot_layout(guides = "collect") &
  
  theme(
    legend.position = "bottom"
  )

# Display final figure
print(final_pca)

# Save high-resolution PNG
ggsave(
  filename = "PCA_complete_VI.png",
  plot = final_pca,
  width = 13,
  height = 11,
  dpi = 300,
  bg = "white"
)

# Combine explained variance results
variance_results <- bind_rows(
  lapply(pca_results, function(x) x$variance)
)

# Combine PCA treatment scores
pca_scores <- bind_rows(
  lapply(pca_results, function(x) x$scores)
)

# Combine variable loadings
pca_loadings <- bind_rows(
  lapply(pca_results, function(x) x$loadings)
)

# Display results
print(variance_results)
print(pca_loadings)

# Export CSV files
write.csv(
  variance_results,
  "PCA_explained_variance.csv",
  row.names = FALSE
)

write.csv(
  pca_scores,
  "PCA_treatment_scores.csv",
  row.names = FALSE
)

write.csv(
  pca_loadings,
  "PCA_variable_loadings.csv",
  row.names = FALSE
)

# Export all results into one Excel workbook
write_xlsx(
  list(
    Explained_variance = variance_results,
    Treatment_scores = pca_scores,
    Variable_loadings = pca_loadings
  ),
  "PCA_complete_results.xlsx"
)

#=========================================================================
#=========================================================================

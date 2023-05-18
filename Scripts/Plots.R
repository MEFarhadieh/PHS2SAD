### Preprocessing and Cleaning ADNI data

# Load library and functions 

setwd("Desktop/Neurotract")
library(ADNIMERGE)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(gtsummary)
library(pROC)
library(areaplot)
library(ggridges)
library(ggsignif)

source("https://adni.bitbucket.io/myfunctions.R")
theme_set(theme_bw())
dxpal <- c( "#9CCDBA", "#F1A84A", "#825C9F")
scale_colour_discrete <- function(...) scale_colour_manual(..., values = dxpal)
scale_fill_discrete <- function(...) scale_fill_manual(..., values = dxpal)
load("Data/final_2.RData")


df.bl <- df[df$VISCODE == "bl",]
df.bl$DX <- gsub("Dementia", "AD", df.bl$DX)

# Spearman's correlation PHS vs GAP-43

cor_result_phs_gap <- cor.test(df.bl$PHS, df.bl$GAP.43, method = "spearman")

ggscatter(df.bl, x = "PHS", y = "GAP.43", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "GAP-43")

ggscatter(df.bl, x = "PHS", y = "GAP.43", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "GAP-43", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = FALSE, color = "black")

# Spearman's correlation PHS vs NFL

cor_result_phs_nfl <- cor.test(df.bl$PHS, df.bl$PLASMA_NFL, method = "spearman")

pdf("Results/PHS_NFL_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "PLASMA_NFL", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "NFL")
dev.off()

pdf("Results/PHS_NFL_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "PLASMA_NFL", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "NFL", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = FALSE, color = "black")
dev.off()

# Spearman's correlation GAP.43 vs NFL

cor_result_gap_nfl <- cor.test(df.bl$GAP.43, df.bl$PLASMA_NFL, method = "spearman")

pdf("Results/GAP.43_NFL_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "PLASMA_NFL", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "NFL")
dev.off()

pdf("Results/GAP.43_NFL_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "PLASMA_NFL", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "NFL", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = FALSE, color = "black")
dev.off()

# Spearman's correlation PHS vs ABETA

cor_result_ABETA <- cor.test(df.bl$PHS, df.bl$PLASMA_ABETA, method = "spearman")

pdf("Results/PHS_ABETA_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "PLASMA_ABETA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "ABETA")
dev.off()

pdf("Results/PHS_ABETA_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "PLASMA_ABETA", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "ABETA", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = FALSE, color = "black")
dev.off()

# Spearman's correlation PHS vs TAU

cor_result_PHS_TAU <- cor.test(df.bl$PHS, df.bl$PLASMA_TAU, method = "spearman")

pdf("Results/PHS_TAU_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "TAU.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "TAU")
dev.off()

pdf("Results/PHS_TAU_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "TAU.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "TAU", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = FALSE, color = "black")
dev.off()

# Spearman's correlation PHS vs PTAU

cor_result_PHS_PTAU <- cor.test(df.bl$PHS, df.bl$PTAU, method = "spearman")

pdf("Results/PHS_PTAU_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "PTAU.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "PTAU")
dev.off()

pdf("Results/PHS_PTAU_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PHS", y = "PTAU.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PHS", ylab = "PTAU", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "grey") + 
  labs(title = "PHS vs PTAU", color = "Baseline Diagnosis")
dev.off()

# Spearman's correlation GAP.43 vs PTAU

cor_result_GAP.43_PTAU <- cor.test(df.bl$GAP.43, df.bl$PTAU, method = "spearman")

pdf("Results/GAP.43_PTAU_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "PTAU.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "PTAU")
dev.off()

pdf("Results/GAP.43_PTAU_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "PTAU.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "PTAU", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "grey") + 
  labs(title = "GAP.43 vs PTAU", color = "Baseline Diagnosis")
dev.off()

# Spearman's correlation GAP.43 vs TAU

cor_result_GAP.43_TAU <- cor.test(df.bl$GAP.43, df.bl$TAU, method = "spearman")

pdf("Results/GAP.43_TAU_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "TAU.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "TAU")
dev.off()

pdf("Results/GAP.43_TAU_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "TAU.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "TAU", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "grey") + 
  labs(title = "GAP.43 vs TAU", color = "Baseline Diagnosis")
dev.off()

# Spearman's correlation GAP.43 vs ABETA

cor_result_GAP.43_ABETA <- cor.test(df.bl$GAP.43, df.bl$ABETA, method = "spearman")

pdf("Results/GAP.43_ABETA_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "ABETA.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "ABETA")
dev.off()

pdf("Results/GAP.43_ABETA_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "GAP.43", y = "ABETA.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "GAP.43", ylab = "ABETA", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "grey") + 
  labs(title = "GAP.43 vs ABETA", color = "Baseline Diagnosis")
dev.off()

# Spearman's correlation PLASMA_NFL vs ABETA

cor_result_PLASMA_NFL_ABETA <- cor.test(df.bl$PLASMA_NFL, df.bl$ABETA, method = "spearman")

pdf("Results/PLASMA_NFL_ABETA_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PLASMA_NFL", y = "ABETA.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PLASMA_NFL", ylab = "ABETA")
dev.off()

pdf("Results/PLASMA_NFL_ABETA_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PLASMA_NFL", y = "ABETA.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PLASMA_NFL", ylab = "ABETA", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "grey") + 
  labs(title = "PLASMA_NFL vs ABETA", color = "Baseline Diagnosis")
dev.off()

# Spearman's correlation PLASMA_NFL vs PTAU

cor_result_PLASMA_NFL_PTAU <- cor.test(df.bl$PLASMA_NFL, df.bl$PTAU, method = "spearman")

pdf("Results/PLASMA_NFL_PTAU_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PLASMA_NFL", y = "PTAU.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PLASMA_NFL", ylab = "PTAU")
dev.off()

pdf("Results/PLASMA_NFL_PTAU_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PLASMA_NFL", y = "PTAU.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PLASMA_NFL", ylab = "PTAU", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "grey") + 
  labs(title = "PLASMA_NFL vs PTAU", color = "Baseline Diagnosis")
dev.off()

# Spearman's correlation PLASMA_NFL vs TAU

cor_result_PLASMA_NFL_TAU <- cor.test(df.bl$PLASMA_NFL, df.bl$TAU, method = "spearman")

pdf("Results/PLASMA_NFL_TAU_Cor_plot_BW.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PLASMA_NFL", y = "TAU.bl", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PLASMA_NFL", ylab = "TAU")
dev.off()

pdf("Results/PLASMA_NFL_TAU_Cor_plot.pdf",
    width = 6, height = 5)
ggscatter(df.bl, x = "PLASMA_NFL", y = "TAU.bl", color = "DX", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "PLASMA_NFL", ylab = "TAU", 
          palette = c("#825C9F", "#9CCDBA", "#F1A84A")) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "grey") + 
  labs(title = "PLASMA_NFL vs TAU", color = "Baseline Diagnosis")
dev.off()


# Subset data for CN and AD groups
df_roc <- subset(df.bl, DX %in% c("CN", "AD"))

# Create an empty data frame to store ROC results
roc_data <- data.frame(variable = character(),
                       sensitivity = numeric(),
                       specificity = numeric())

# Loop over the variables
variables <- c("PHS", "ABETA.bl", "TAU.bl", "PTAU.bl", "GAP.43", "PLASMA_NFL")
colors <- c("#9370DB", "#87CEFA", "#66CD00", "#FFD700", "#FFA500", "#CD0000")

for (i in seq_along(variables)) {
  # Create ROC curve
  roc <- roc(df_roc$DX, df_roc[, variables[i]], levels = c("CN", "AD"))
  
  # Store ROC results in the data frame
  roc_data <- rbind(roc_data, data.frame(variable = variables[i],
                                         sensitivity = roc$sensitivities,
                                         specificity = roc$specificities))
}

# Create ROC plot
roc_plot <- ggplot(roc_data, aes(1 - specificity, sensitivity, color = variable)) +
  geom_line(size = 1.2, aes(group = variable), alpha = 0.8, lineend = "round") +
  geom_abline(color = "gray40", linetype = "dashed", size = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "ROC Curves for CN vs AD",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Variable") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",  # Legend on the right side
        legend.title = element_blank(),  # Remove legend title
        legend.key.height = unit(0.5, "cm"),  # Set the height of legend keys
        legend.key.width = unit(1.5, "cm"),  # Set the width of legend keys
        legend.margin = margin(10, 10, 10, 10))  # Set margin around the legend

# Save plot
pdf("Results/ROC_CN_AD.pdf",
    width = 6, height = 4)
roc_plot + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(color = guide_legend(title = "Variable"))
dev.off()

# Subset data for CN and MCI groups
df_roc <- subset(df.bl, DX %in% c("CN", "MCI"))

# Create an empty data frame to store ROC results
roc_data <- data.frame(variable = character(),
                       sensitivity = numeric(),
                       specificity = numeric())

# Loop over the variables
variables <- c("PHS", "ABETA.bl", "TAU.bl", "PTAU.bl", "GAP.43", "PLASMA_NFL")
colors <- c("#9370DB", "#87CEFA", "#66CD00", "#FFD700", "#FFA500", "#CD0000")

for (i in seq_along(variables)) {
  # Create ROC curve
  roc <- roc(df_roc$DX, df_roc[, variables[i]], levels = c("CN", "MCI"))
  
  # Store ROC results in the data frame
  roc_data <- rbind(roc_data, data.frame(variable = variables[i],
                                         sensitivity = roc$sensitivities,
                                         specificity = roc$specificities))
}

# Create ROC plot
roc_plot <- ggplot(roc_data, aes(1 - specificity, sensitivity, color = variable)) +
  geom_line(size = 1.2, aes(group = variable), alpha = 0.8, lineend = "round") +
  geom_abline(color = "gray40", linetype = "dashed", size = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "ROC Curves for CN vs MCI",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Variable") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",  # Legend on the right side
        legend.title = element_blank(),  # Remove legend title
        legend.key.height = unit(0.5, "cm"),  # Set the height of legend keys
        legend.key.width = unit(1.5, "cm"),  # Set the width of legend keys
        legend.margin = margin(10, 10, 10, 10))  # Set margin around the legend

# Save plot
pdf("Results/ROC_CN_MCI.pdf",
    width = 6, height = 4)
roc_plot + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(color = guide_legend(title = "Variable"))
dev.off()

# Subset data for MCI and AD groups
df_roc <- subset(df.bl, DX %in% c("MCI", "AD"))

# Create an empty data frame to store ROC results
roc_data <- data.frame(variable = character(),
                       sensitivity = numeric(),
                       specificity = numeric())

# Loop over the variables
variables <- c("PHS", "ABETA.bl", "TAU.bl", "PTAU.bl", "GAP.43", "PLASMA_NFL")
colors <- c("#9370DB", "#87CEFA", "#66CD00", "#FFD700", "#FFA500", "#CD0000")

for (i in seq_along(variables)) {
  # Create ROC curve
  roc <- roc(df_roc$DX, df_roc[, variables[i]], levels = c("MCI", "AD"))
  
  # Store ROC results in the data frame
  roc_data <- rbind(roc_data, data.frame(variable = variables[i],
                                         sensitivity = roc$sensitivities,
                                         specificity = roc$specificities))
}

# Create ROC plot
roc_plot <- ggplot(roc_data, aes(1 - specificity, sensitivity, color = variable)) +
  geom_line(size = 1.2, aes(group = variable), alpha = 0.8, lineend = "round") +
  geom_abline(color = "gray40", linetype = "dashed", size = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "ROC Curves for MCI vs AD",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Variable") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",  # Legend on the right side
        legend.title = element_blank(),  # Remove legend title
        legend.key.height = unit(0.5, "cm"),  # Set the height of legend keys
        legend.key.width = unit(1.5, "cm"),  # Set the width of legend keys
        legend.margin = margin(10, 10, 10, 10))  # Set margin around the legend

# Save plot
pdf("Results/ROC_MCI_AD.pdf",
    width = 6, height = 4)
roc_plot + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(color = guide_legend(title = "Variable"))
dev.off()

# Subset data for CN and AD&MCI groups
df.bl_C.MA <- df.bl
df.bl_C.MA$DX <- ifelse(
  df.bl_C.MA$DX%in% c("AD", "MCI"), "AD&MCI", df.bl_C.MA$DX)
df.bl_C.MA$DX <- gsub("2", "CN", df.bl_C.MA$DX)
unique(df.bl_C.MA$DX)
df_roc <- subset(df.bl_C.MA, DX %in% c("CN", "AD&MCI"))

# Create an empty data frame to store ROC results
roc_data <- data.frame(variable = character(),
                       sensitivity = numeric(),
                       specificity = numeric())

# Loop over the variables
variables <- c("PHS", "ABETA.bl", "TAU.bl", "PTAU.bl", "GAP.43", "PLASMA_NFL")
colors <- c("#9370DB", "#87CEFA", "#66CD00", "#FFD700", "#FFA500", "#CD0000")

for (i in seq_along(variables)) {
  # Create ROC curve
  roc <- roc(df_roc$DX, df_roc[, variables[i]], levels = c("CN", "AD&MCI"))
  
  # Store ROC results in the data frame
  roc_data <- rbind(roc_data, data.frame(variable = variables[i],
                                         sensitivity = roc$sensitivities,
                                         specificity = roc$specificities))
}

# Create ROC plot
roc_plot <- ggplot(roc_data, aes(1 - specificity, sensitivity, color = variable)) +
  geom_line(size = 1.2, aes(group = variable), alpha = 0.8, lineend = "round") +
  geom_abline(color = "gray40", linetype = "dashed", size = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "ROC Curves for CN vs AD&MCI",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Variable") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",  # Legend on the right side
        legend.title = element_blank(),  # Remove legend title
        legend.key.height = unit(0.5, "cm"),  # Set the height of legend keys
        legend.key.width = unit(1.5, "cm"),  # Set the width of legend keys
        legend.margin = margin(10, 10, 10, 10))  # Set margin around the legend

# Save plot
pdf("Results/ROC_CN_AD&MCI.pdf",
    width = 6, height = 4)
roc_plot + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(color = guide_legend(title = "Variable"))
dev.off()

# Subset data for AD and CN&MCI groups
df.bl_A.CM <- df.bl
df.bl_A.CM$DX <- ifelse(
  df.bl_A.CM$DX%in% c("CN", "MCI"), "CN&MCI", df.bl_A.CM$DX)
df.bl_A.CM$DX <- gsub("2", "AD", df.bl_A.CM$DX)
unique(df.bl_A.CM$DX)
df_roc <- subset(df.bl_A.CM, DX %in% c("AD", "CN&MCI"))

# Create an empty data frame to store ROC results
roc_data <- data.frame(variable = character(),
                       sensitivity = numeric(),
                       specificity = numeric())

# Loop over the variables
variables <- c("PHS", "ABETA.bl", "TAU.bl", "PTAU.bl", "GAP.43", "PLASMA_NFL")
colors <- c("#9370DB", "#87CEFA", "#66CD00", "#FFD700", "#FFA500", "#CD0000")

for (i in seq_along(variables)) {
  # Create ROC curve
  roc <- roc(df_roc$DX, df_roc[, variables[i]], levels = c("AD", "CN&MCI"))
  
  # Store ROC results in the data frame
  roc_data <- rbind(roc_data, data.frame(variable = variables[i],
                                         sensitivity = roc$sensitivities,
                                         specificity = roc$specificities))
}

# Create ROC plot
roc_plot <- ggplot(roc_data, aes(1 - specificity, sensitivity, color = variable)) +
  geom_line(size = 1.2, aes(group = variable), alpha = 0.8, lineend = "round") +
  geom_abline(color = "gray40", linetype = "dashed", size = 0.5) +
  scale_color_manual(values = colors) +
  labs(title = "ROC Curves for AD vs CN&MCI",
       x = "1 - Specificity",
       y = "Sensitivity",
       color = "Variable") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",  # Legend on the right side
        legend.title = element_blank(),  # Remove legend title
        legend.key.height = unit(0.5, "cm"),  # Set the height of legend keys
        legend.key.width = unit(1.5, "cm"),  # Set the width of legend keys
        legend.margin = margin(10, 10, 10, 10))  # Set margin around the legend

# Save plot
pdf("Results/ROC_AD_CN&MCI.pdf",
    width = 6, height = 4)
roc_plot + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(color = guide_legend(title = "Variable"))
dev.off()

## Find cutoff threshold of PHS by Youden index
# CN vs AD
df_roc <- subset(df.bl, DX %in% c("CN", "AD"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "AD"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# CN vs MCI
df_roc <- subset(df.bl, DX %in% c("CN", "MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# AD vs MCI
df_roc <- subset(df.bl, DX %in% c("AD", "MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("AD", "MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# CN vs AD&MCI
df_roc <- subset(df.bl_C.MA, DX %in% c("CN", "AD&MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "AD&MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# AD vs CN&MCI
df_roc <- subset(df.bl_A.CM, DX %in% c("AD", "CN&MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("AD", "CN&MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

### Progression comparison Analysis 

# Check VISCODE distrbution 
viscount <- as.data.frame(table(df$VISCODE))
viscount$Var1 <- factor(viscount$Var1, 
                   levels = viscount$Var1[order(
                     viscount$Freq, decreasing = T)])
ggplot(viscount, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Individuals") +
  ylab("VISCODE") +
  ggtitle("Distrbution of VISCODEs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the stacked histogram of PHS
df_subset <- subset(df.bl, DX %in% c("AD", "CN", "MCI"))
ggplot(df_subset, aes(x = PHS, fill = DX)) +
  geom_histogram(breaks = seq(min(df_subset$PHS), max(df_subset$PHS), length.out = 21),
                 position = "stack", color = "black") +
  scale_fill_manual(values = c("#825C9F", "#9CCDBA", "#F1A84A"),
                    labels = c("AD", "CN", "MCI")) +
  labs(x = "PHS",
       y = "Number of Individuals",
       title = "Stacked Histogram of PHS by Group",
       fill = "Baseline Diagnosis") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Create the density plot of PHS
ggplot(df.bl, aes(x = PHS, y = DX, fill = DX)) +
  geom_density_ridges(scale = 0.9) +
  scale_fill_manual(values = c(AD = "#825C9F", MCI = "#F1A84A", CN = "#9CCDBA")) +
  xlab("PHS") +
  ylab("Density") +
  ggtitle("Distribution of basline diagnosis based on PHS") + 
  theme_set(theme_ridges())

## Subset df based on 24m and 48m after basline
df_m24 <- df[df$VISCODE == "m24",]
df_m48 <- df[df$VISCODE == "m48",]

write.table(df_m24, "Data/df_m24.csv", sep = ",",
            row.names = F, col.names = T, quote = F)
write.table(df_m48, "Data/df_m48.csv", sep = ",",
            row.names = F, col.names = T, quote = F)

# Remove rows with NA in DX column
df_m24 <- df_m24[!is.na(df_m24$DX), ]
df_m48 <- df_m48[!is.na(df_m48$DX), ]

# Subset data for AD and CN&MCI groups
df_m24_A.CM <- df_m24
df_m24_A.CM$DX <- ifelse(
  df_m24_A.CM$DX%in% c("CN", "MCI"), "CN&MCI", df_m24_A.CM$DX)
df_m24_A.CM$DX <- gsub("2", "AD", df_m24_A.CM$DX)
unique(df_m24_A.CM$DX)

df_m24_C.MA <- df_m24
df_m24_C.MA$DX <- ifelse(
  df_m24_C.MA$DX%in% c("AD", "MCI"), "AD&MCI", df_m24_C.MA$DX)
df_m24_C.MA$DX <- gsub("2", "CN", df_m24_C.MA$DX)
unique(df_m24_C.MA$DX)

df_m48_A.CM <- df_m48
df_m48_A.CM$DX <- ifelse(
  df_m48_A.CM$DX%in% c("CN", "MCI"), "CN&MCI", df_m48_A.CM$DX)
df_m48_A.CM$DX <- gsub("2", "AD", df_m48_A.CM$DX)
unique(df_m48_A.CM$DX)

df_m48_C.MA <- df_m48
df_m48_C.MA$DX <- ifelse(
  df_m48_C.MA$DX%in% c("AD", "MCI"), "AD&MCI", df_m48_C.MA$DX)
df_m48_C.MA$DX <- gsub("2", "CN", df_m48_C.MA$DX)
unique(df_m48_C.MA$DX)


## Find cutoff threshold of PHS by Youden index for m24 data
# CN vs AD
df_roc <- subset(df_m24, DX %in% c("CN", "AD"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "AD"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# CN vs MCI
df_roc <- subset(df_m24, DX %in% c("CN", "MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# AD vs MCI
df_roc <- subset(df_m24, DX %in% c("AD", "MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("AD", "MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# CN vs AD&MCI
df_roc <- subset(df_m24_C.MA, DX %in% c("CN", "AD&MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "AD&MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# AD vs CN&MCI
df_roc <- subset(df_m24_A.CM, DX %in% c("AD", "CN&MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("AD", "CN&MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# Create the stacked histogram of PHS
df_subset <- subset(df_m24, DX %in% c("AD", "CN", "MCI"))
ggplot(df_subset, aes(x = PHS, fill = DX)) +
  geom_histogram(breaks = seq(min(df_subset$PHS), max(df_subset$PHS), length.out = 21),
                 position = "stack", color = "black") +
  scale_fill_manual(values = c("#825C9F", "#9CCDBA", "#F1A84A"),
                    labels = c("AD", "CN", "MCI")) +
  labs(x = "PHS",
       y = "Number of Individuals",
       title = "Stacked Histogram of PHS by Group",
       fill = "24m Diagnosis") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Create the density plot of PHS
ggplot(df_m24, aes(x = PHS, y = DX, fill = DX)) +
  geom_density_ridges(scale = 0.9) +
  scale_fill_manual(values = c(AD = "#825C9F", MCI = "#F1A84A", CN = "#9CCDBA")) +
  xlab("PHS") +
  ylab("Density") +
  ggtitle("Distribution of 24m diagnosis based on PHS") + 
  theme_set(theme_ridges())

## Find cutoff threshold of PHS by Youden index for m48 data
# CN vs AD
df_roc <- subset(df_m48, DX %in% c("CN", "AD"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "AD"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# CN vs MCI
df_roc <- subset(df_m48, DX %in% c("CN", "MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# AD vs MCI
df_roc <- subset(df_m48, DX %in% c("AD", "MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("AD", "MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# CN vs AD&MCI
df_roc <- subset(df_m48_C.MA, DX %in% c("CN", "AD&MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("CN", "AD&MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# AD vs CN&MCI
df_roc <- subset(df_m48_A.CM, DX %in% c("AD", "CN&MCI"))
roc <- roc(df_roc$DX, df_roc$PHS, levels = c("AD", "CN&MCI"))
coords(roc, "best", ret = "threshold", best.method = "youden")

# Create the stacked histogram of PHS
df_subset <- subset(df_m48, DX %in% c("AD", "CN", "MCI"))
ggplot(df_subset, aes(x = PHS, fill = DX)) +
  geom_histogram(breaks = seq(min(df_subset$PHS), max(df_subset$PHS), length.out = 21),
                 position = "stack", color = "black") +
  scale_fill_manual(values = c("#825C9F", "#9CCDBA", "#F1A84A"),
                    labels = c("AD", "CN", "MCI")) +
  labs(x = "PHS",
       y = "Number of Individuals",
       title = "Stacked Histogram of PHS by Group",
       fill = "48m Diagnosis") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Create the density plot of PHS
ggplot(df_m48, aes(x = PHS, y = DX, fill = DX)) +
  geom_density_ridges(scale = 0.9) +
  scale_fill_manual(values = c(AD = "#825C9F", MCI = "#F1A84A", CN = "#9CCDBA")) +
  xlab("PHS") +
  ylab("Density") +
  ggtitle("Distribution of 48m diagnosis based on PHS") + 
  theme_set(theme_ridges())

### Longitudinal analysis 
# Select bl to 48m data by each 12 month
unique(df$VISCODE)
df <- df[df$VISCODE %in% c("bl", "m12", "m24","m36","m48"),]
df$VISCODE <- gsub("bl", "0", df$VISCODE)
df$VISCODE <- gsub("m12", "12", df$VISCODE)
df$VISCODE <- gsub("m24", "24", df$VISCODE)
df$VISCODE <- gsub("m36", "36", df$VISCODE)
df$VISCODE <- gsub("m48", "48", df$VISCODE)

## Separate PHS based on the cutoff threshold 
# Subset the data and calculate PHS_Group
df_subset <- df %>%
  mutate(PHS_Group = ifelse(PHS <= 0.652, "PHS Negative", "PHS Positive")) %>%
  dplyr::select(VISCODE, CDRSB, PHS_Group)

# Calculate the mean and 95% confidence interval (CI) for CDRSB score
mean_CI <- df_subset %>%
  group_by(VISCODE, PHS_Group) %>%
  summarise(Mean = mean(CDRSB),
            CI_lower = t.test(CDRSB)$conf.int[1],
            CI_upper = t.test(CDRSB)$conf.int[2]) %>%
  ungroup()

# Perform Mann-Whitney test for comparing PHS groups at each VISCODE
mw_test <- df_subset %>%
  group_by(VISCODE) %>%
  do(test = wilcox.test(.$CDRSB ~ .$PHS_Group)) %>%
  mutate(p_value = test$p.value) %>%
  dplyr::select(VISCODE, p_value)

# Rearrange data for connecting points with lines
line_data <- mean_CI %>%
  group_by(PHS_Group) %>%
  arrange(as.numeric(VISCODE)) %>%
  mutate(line_group = row_number()) %>%
  ungroup()

# Sort line_data by line_group
line_data <- line_data[order(line_data$line_group), ]

# Create the plot with connected lines
plot <- ggplot() +
  geom_line(data = line_data, aes(x = as.numeric(VISCODE), y = Mean, color = PHS_Group, group = line_group)) +
  geom_point(data = mean_CI, aes(x = as.numeric(VISCODE), y = Mean, color = PHS_Group)) +
  geom_errorbar(data = mean_CI, aes(x = as.numeric(VISCODE), ymin = CI_lower, ymax = CI_upper, color = PHS_Group), width = 0.1) +
  scale_color_manual(values = c("#9CCDBA", "#825C9F")) +
  labs(x = "Month", y = "CDR-SB Score") +
  theme_minimal()

# Add asterisks to indicate significance level
p_values <- mw_test$p_value
significance <- ifelse(p_values < 0.001, "***", ifelse(p_values < 0.01, "**", ifelse(p_values < 0.05, "*", "")))
plot <- plot + geom_text(data = mw_test, aes(x = as.numeric(VISCODE), y = max(mean_CI$Mean), label = significance))

# Display the plot
plot

# Create a new grouping variable based on PHS threshold
subset_df$Group <- ifelse(subset_df$PHS <= 0.652, "PHS Negative", "PHS Positive")

# Split the data into PHS Negative and PHS Positive groups
negative_group <- subset_df %>%
  filter(Group == "PHS Negative")
positive_group <- subset_df %>%
  filter(Group == "PHS Positive")

# Calculate the Interquartile Range for each VISCODE
iq_range <- subset_df %>%
  group_by(VISCODE) %>%
  summarize(Q1 = quantile(CDRSB, 0.25, na.rm = TRUE),
            Q3 = quantile(CDRSB, 0.75, na.rm = TRUE), .groups = "drop")


# Perform Mann-Whitney test between PHS Negative and PHS Positive groups
p_value <- wilcox.test(CDRSB ~ Group, data = subset_df)$p.value


# Merge subset_df and iq_range data frames
merged_df <- merge(subset_df, iq_range, by = "VISCODE")

# Create the plot
plot <- ggplot(merged_df, aes(x = VISCODE, y = CDRSB, color = PHS > 0.652)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.1, size = 0.7) +
  geom_signif(y_position = 30, xmin = 1, xmax = 2, annotation = "p < 0.05", tip_length = 0.03) +
  geom_signif(y_position = 32, xmin = 1, xmax = 2, annotation = "p < 0.01", tip_length = 0.03) +
  geom_signif(y_position = 34, xmin = 1, xmax = 2, annotation = "p < 0.001", tip_length = 0.03) +
  scale_color_manual(values = c("#9CCDBA", "#825C9F"), labels = c("PHS Negative", "PHS Positive")) +
  labs(x = "Month", y = "CDR-SB Score") +
  theme_minimal()

# Display the plot
plot

write.table(df,"Data/df_long.csv", sep = ",", quote = F,
            row.names = F, col.names = T)

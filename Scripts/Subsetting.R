### The goal is generation subsets cohort for longitudinal analysis 

library(dplyr)
library(tidyverse)
library(gtsummary)

setwd("Desktop/Neurotract/")
df.bl <- read.csv("df_bl.csv")
df_m24 <- read.csv("Data/df_m24.csv")
df_m48 <- read.csv("Data/df_m48.csv")
df_long <- read.csv("Data/df_long.csv")


# Select bl to 24m data by each 6 month
df <- read.delim("Data/Final_2.tsv")
unique(df$VISCODE)
df <- df[df$VISCODE %in% c("bl", "m12", "m24","m06","m18"),]
df$VISCODE <- gsub("bl", "0", df$VISCODE)
df$VISCODE <- gsub("m12", "12", df$VISCODE)
df$VISCODE <- gsub("m24", "24", df$VISCODE)
df$VISCODE <- gsub("m06", "6", df$VISCODE)
df$VISCODE <- gsub("m18", "18", df$VISCODE)
unique(df$VISCODE)

m18 <- df[df$VISCODE == "18",]

sub_24 <- df[df$RID %in% m18$RID, ]

write.table(sub_24,"Data/sub_24.csv", sep = ",", quote = F,
            row.names = F, col.names = T)

sub_24[sub_24$VISCODE == "0",] %>% 
  select(c(AGE,PTGENDER,PTEDUCAT,APOE4,PHS,ABETA.bl,TAU.bl,PTAU.bl,
           GAP.43,PLASMA_NFL,MMSE.bl,CDRSB.bl,ADAS11.bl,ADAS13.bl,
           ADASQ4.bl,MOCA.bl,DX)) %>% 
  tbl_summary(by = DX)

# Select bl to 48m data by each 12 month
df <- read.delim("Data/Final_2.tsv")
df <- df[df$VISCODE %in% c("bl", "m12", "m24","m36","m48"),]
df$VISCODE <- gsub("bl", "0", df$VISCODE)
df$VISCODE <- gsub("m12", "12", df$VISCODE)
df$VISCODE <- gsub("m24", "24", df$VISCODE)
df$VISCODE <- gsub("m36", "36", df$VISCODE)
df$VISCODE <- gsub("m48", "48", df$VISCODE)
unique(df$VISCODE)

m36 <- df[df$VISCODE == "36",]
sub_48 <- df[df$RID %in% m36$RID, ]

write.table(sub_48,"Data/sub_48.csv", sep = ",", quote = F,
            row.names = F, col.names = T)

sub_48[sub_48$VISCODE == "0",] %>% 
  select(c(AGE,PTGENDER,PTEDUCAT,APOE4,PHS,ABETA.bl,TAU.bl,PTAU.bl,
           GAP.43,PLASMA_NFL,MMSE.bl,CDRSB.bl,ADAS11.bl,ADAS13.bl,
           ADASQ4.bl,MOCA.bl,DX)) %>% 
  tbl_summary(by = DX)


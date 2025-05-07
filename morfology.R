library(tidyverse)
library(ggpubr)
library(rstatix)

data <-
  read.csv("C:/Users/mient/Desktop/biology/phD/2025_04_14_1744653553.csv")

data <- data %>%
  select(-any_of(
    c(
      "Operator",
      "Age",
      "Tel.",
      "Case.ID",
      "Send.time",
      "Auditor",
      "Sender",
      "Sampling.time",
      "Pat..Name",
      "Master",
      "Time",
      "Mode",
      "WBC.time",
      "RBC.time"
    )
  )) %>%
  select(where( ~ !all(is.na(.)))) %>%
  mutate(
    Group = case_when(
      ID %in% paste0("'", 101:106, "'") ~ "VEH_OR",
      ID %in% paste0("'", 201:206, "'") ~ "GEN_OR",
      ID %in% paste0("'", 301:306, "'") ~ "VEH_IP",
      ID %in% paste0("'", 401:406, "'") ~ "17A_IP",
      ID %in% paste0("'", 501:506, "'") ~ "TUB_IP",
      ID %in% paste0("'", 601:606, "'") ~ "COMB_IP",
      ID %in% paste0("'", 701:706, "'") ~ "MP_IP",
      TRUE ~ "NAIVE"
    ),
    Gender = "female"
  ) %>%
  mutate(across(WBC:last_col(), as.character)) %>%
  filter(ID != "'0000000'")

data_long <- data %>%
  pivot_longer(cols = WBC:last_col(),
               names_to = "Analyte",
               values_to = "Result") %>%
  mutate(Result = str_remove_all(Result, "[^0-9]"),
         Result = as.integer(Result)) %>%
  filter(!is.na(Result)) %>%
  mutate(
    Analyte = case_when(
      Analyte == "LYM." ~ "LYM_Abs",
      Analyte == "LYM..1" ~ "LYM_Perc",
      Analyte == "MON." ~ "MONO_Abs",
      Analyte == "MON..1" ~ "MONO_Perc",
      Analyte == "NEU." ~ "NEU_Abs",
      Analyte == "NEU..1" ~ "NEU_Perc",
      Analyte == "EOS." ~ "EOS_Abs",
      Analyte == "EOS..1" ~ "EOS_Perc",
      Analyte == "BASO." ~ "BASO_Abs",
      Analyte == "BASO..1" ~ "BASO_Perc",
      Analyte == "NRBC." ~ "NRBC_Abs",
      Analyte == "NRBC..1" ~ "NRBC_Perc",
      Analyte == "ALY." ~ "ALY_Abs",
      Analyte == "ALY..1" ~ "ALY_Perc",
      Analyte == "LIC." ~ "LIC_Abs",
      Analyte == "LIC..1" ~ "LIC_Perc",
      TRUE ~ Analyte
    )
  )

### ANALIZA IP ###
ip_groups <- c("VEH_IP", "17A_IP", "TUB_IP", "COMB_IP", "MP_IP")
comparison_pairs <-
  lapply(ip_groups[-1], function(gr)
    c("VEH_IP", gr))

data_long_ip <- data_long %>%
  filter(Group %in% ip_groups) %>%
  mutate(Group = factor(Group, levels = ip_groups))

pval_df <- data.frame()

for (analyte_value in unique(data_long_ip$Analyte)) {
  for (groups in comparison_pairs) {
    df_pair <-
      data_long_ip %>% filter(Analyte == analyte_value, Group %in% groups)
    
    if (nrow(df_pair) > 1 && n_distinct(df_pair$Group) == 2) {
      normality <- df_pair %>%
        group_by(Group) %>%
        summarise(shapiro_p = shapiro.test(Result)$p.value,
                  .groups = "drop")
      
      levene <- df_pair %>% levene_test(Result ~ Group)
      equal_var <- levene$p > 0.05
      
      if (all(normality$shapiro_p > 0.05)) {
        test_result <-
          t.test(Result ~ Group, data = df_pair, var.equal = equal_var)
      } else {
        test_result <- wilcox.test(Result ~ Group, data = df_pair)
      }
      
      pval <- test_result$p.value
      label <- ifelse(pval < 0.05, "*", NA)
      y_max <- max(df_pair$Result, na.rm = TRUE)
      
      pval_df <- rbind(
        pval_df,
        data.frame(
          Analyte = analyte_value,
          group1 = groups[1],
          group2 = groups[2],
          p.adj = label,
          y.position = y_max * 1.05
        )
      )
    }
  }
}

plot_stat_ip <- ggplot(data_long_ip, aes(x = Group, y = Result)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 1,
    color = "black"
  ) +
  facet_wrap( ~ Analyte, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(
    pval_df,
    label = "p.adj",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 2,
    steps.increase = 0.5
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    strip.text = element_text(size = 8),
    legend.position = "top"
  ) +
  labs(y = "Value", x = "Group", fill = "Group")

print(plot_stat_ip)


### ANALIZA OR ###

or_groups <- c("NAIVE", "VEH_OR", "GEN_OR")
comparison_pairs_or <- list(c("VEH_OR", "NAIVE"),
                            c("VEH_OR", "GEN_OR"))

data_long_or <- data_long %>%
  filter(Group %in% or_groups) %>%
  mutate(Group = factor(Group, levels = or_groups))

pval_df_or <- data.frame()

for (analyte_value in unique(data_long_or$Analyte)) {
  for (groups in comparison_pairs_or) {
    df_pair <- data_long_or %>%
      filter(Analyte == analyte_value, Group %in% groups)
    
    if (nrow(df_pair) > 1 && n_distinct(df_pair$Group) == 2) {
      normality <- df_pair %>%
        group_by(Group) %>%
        summarise(
          shapiro_p = if (n() >= 3)
            shapiro.test(Result)$p.value
          else
            NA_real_,
          .groups = "drop"
        )
      
      levene <- df_pair %>% levene_test(Result ~ Group)
      equal_var <- levene$p > 0.05
      
      if (all(normality$shapiro_p > 0.05, na.rm = TRUE)) {
        test_result <-
          t.test(Result ~ Group, data = df_pair, var.equal = equal_var)
      } else {
        test_result <- wilcox.test(Result ~ Group, data = df_pair)
      }
      
      pval <- test_result$p.value
      label <- ifelse(pval < 0.05, "*", NA)
      y_max <- max(df_pair$Result, na.rm = TRUE)
      
      pval_df_or <- rbind(
        pval_df_or,
        data.frame(
          Analyte = analyte_value,
          group1 = groups[1],
          group2 = groups[2],
          p.adj = label,
          y.position = y_max * 1.05
        )
      )
    }
  }
}

plot_stat_or <- ggplot(data_long_or, aes(x = Group, y = Result)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size = 1,
    color = "black"
  ) +
  facet_wrap( ~ Analyte, scales = "free_y", ncol = 6) +
  stat_pvalue_manual(
    pval_df_or,
    label = "p.adj",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 2,
    steps.increase = 0.5
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    strip.text = element_text(size = 8),
    legend.position = "top"
  ) +
  labs(y = "Value", x = "Group", fill = "Group")

print(plot_stat_or)

ggsave("plot_stat_ip_morf.png", plot = plot_stat_ip, width = 12, height = 8)
ggsave("plot_stat_or_morf.png", plot = plot_stat_or, width = 12, height = 8)

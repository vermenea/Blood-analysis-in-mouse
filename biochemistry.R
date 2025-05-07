library(dplyr)
library(readr)
library(stringr)
library(rstatix)
library(ggplot2)
library(ggpubr)

extract_analytes <- function(file_path) {
  rows <- list()
  header_found <- FALSE
  medical_id <- NA
  animal_name <- NA
  flag_hem <- FALSE
  
  lines <- readLines(file_path, encoding = "UTF-8", warn = FALSE)
  
  if (any(grepl("Abnormal sample: HEM+", lines))) {
    flag_hem = TRUE
  }
  
  for (line in lines) {
    if (is.na(line))
      next
    line <- str_squish(as.character(line)) #spaces & tabs cleanup
    if (line == "")
      next
    
    if (str_detect(line, "^Medical ID:")) {
      medical_id <- str_replace(line, "^Medical ID:\\s*", "")
      medical_id <- str_remove_all(medical_id, "^,|,$")
      medical_id <- str_trim(medical_id)
      next
    }
    if (str_detect(line, "^Animal:")) {
      animal_name <- str_replace(line, "^Animal:\\s*", "")
      animal_name <- str_remove_all(animal_name, "^,|,$")
      animal_name <- str_trim(animal_name)
      next
    }
    
    if (!header_found && str_detect(line, "^(Assay|Analyte),")) {
      header_found <- TRUE
      next
    }
    
    if (header_found) {
      if (str_starts(line, "-") || str_starts(line, "WAT")) {
        break
      }
      parts <- str_split(line, ",")[[1]] %>% str_trim()
      
      if (length(parts) < 2 || parts[1] == "") {
        next
      }
      
      analyte <- parts[1]
      result <- parts[2]
      ref <- NA
      unit <- NA
      
      if (length(parts) == 5) {
        ref <- parts[4]
        unit <- parts[5]
      } else if (length(parts) == 4) {
        ref <- parts[3]
        unit <- parts[4]
      } else if (length(parts) == 3) {
        ref <- parts[3]
      }
      
      rows <- append(rows, list(
        data.frame(
          ID = medical_id,
          Animal = animal_name,
          Analyte = analyte,
          Result = result,
          Ref = ref,
          Unit = unit,
          Flags = ifelse(flag_hem, "HEM+", NA),
          stringsAsFactors = FALSE
        )
      ))
    }
  }
  
  if (length(rows) == 0) {
    return(
      data.frame(
        ID = character(),
        Animal = character(),
        Analyte = character(),
        Result = character(),
        Ref = character(),
        Unit = character()
      )
    )
  }
  
  bind_rows(rows)
}


FOLDER <- "C:/Users/mient/Desktop/biology/phD/biochemia"
files <- list.files(FOLDER, pattern = "^Us.*\\.csv$", full.names = TRUE)


frames <- lapply(files, extract_analytes)
df_all <- bind_rows(frames)


df_all <- df_all %>%
  mutate(
    ID = as.numeric(ID),
    Group = case_when(
      ID >= 100 & ID < 200 ~ "VEH_OR",
      ID >= 200 & ID < 300 ~ "GEN_OR",
      ID >= 300 & ID < 400 ~ "VEH_IP",
      ID >= 400 & ID < 500 ~ "17A_IP",
      ID >= 500 & ID < 600 ~ "TUB_IP",
      ID >= 600 & ID < 700 ~ "COMB_IP",
      ID >= 700 & ID < 800 ~ "MP_IP",
      TRUE ~ NA_character_
    ),
    Result = na_if(Result, "--"),
    Result = str_remove(Result, "[<>]\\s*"),
    Result = as.numeric(Result),
    hundreds = ID %% 1000 %/% 100,
    tens = ID %% 100
  ) %>%
  arrange(hundreds, tens) %>%
  select(ID, Group, Animal, Analyte, Result, Ref, Unit, Flags)


View(df_all)

## ANALIZA IP ##

ip_groups <- c("VEH_IP", "17A_IP", "TUB_IP", "COMB_IP", "MP_IP")
comparison_pairs <- lapply(ip_groups[-1], function(gr)
  c("VEH_IP", gr))

data_long_ip <- df_all %>%
  filter(Group %in% ip_groups) %>%
  mutate(Group = factor(Group, levels = ip_groups))

pval_df <- data.frame()


for (analyte_value in unique(data_long_ip$Analyte)) {
  for (groups in comparison_pairs) {
    df_pair <- data_long_ip %>%
      filter(Analyte == analyte_value, Group %in% groups)
    
    group_counts <- df_pair %>% count(Group)
    
    
    if (nrow(df_pair) > 1 &&
        n_distinct(df_pair$Group) == 2 &&
        all(group_counts$n >= 3)) {
      df_pair <- df_pair %>% drop_na(Result)
      
      
      results_normality <- df_pair %>%
        group_by(Group) %>%
        filter(n() >= 3) %>%
        summarise(shapiro_p = tryCatch(
          shapiro.test(Result)$p.value,
          error = function(e)
            NA_real_
        ), .groups = "drop") %>%
        filter(!is.na(shapiro_p))
      
      
      levene_test_result <- tryCatch(
        levene_test(Result ~ Group, data = df_pair),
        error = function(e)
          NA
      )
      
      
      if (!inherits(levene_test_result, "try-error") &&
          !is.na(levene_test_result$p) &&
          !is.nan(levene_test_result$p)) {
        if (levene_test_result$p < 0.05) {
          t_test_result <- t.test(Result ~ Group,
                                  data = df_pair,
                                  var.equal = FALSE)
        } else {
          t_test_result <- t.test(Result ~ Group,
                                  data = df_pair,
                                  var.equal = TRUE)
        }
        
        pval <- t_test_result$p.value
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


or_groups <- c("VEH_OR", "GEN_OR")
comparison_pairs_or <- lapply(or_groups[-1], function(gr)
  c("VEH_OR", gr))

## ANALIZA OR ##
data_long_or <- df_all %>%
  filter(Group %in% or_groups) %>%
  mutate(Group = factor(Group, levels = or_groups))

pval_df_or <- data.frame()


for (analyte_value in unique(data_long_or$Analyte)) {
  for (groups in comparison_pairs_or) {
    df_pair <- data_long_or %>%
      filter(Analyte == analyte_value, Group %in% groups)
    
    group_counts <- df_pair %>% count(Group)
    
    
    if (nrow(df_pair) > 1 &&
        n_distinct(df_pair$Group) == 2 &&
        all(group_counts$n >= 3)) {
      df_pair <- df_pair %>% drop_na(Result)
      
      
      results_normality <- df_pair %>%
        group_by(Group) %>%
        filter(n() >= 3) %>%
        summarise(shapiro_p = tryCatch(
          shapiro.test(Result)$p.value,
          error = function(e)
            NA_real_
        ), .groups = "drop") %>%
        filter(!is.na(shapiro_p))
      
      
      levene_test_result <- tryCatch(
        levene_test(Result ~ Group, data = df_pair),
        error = function(e)
          NA
      )
      
      
      if (!inherits(levene_test_result, "try-error") &&
          !is.na(levene_test_result$p) &&
          !is.nan(levene_test_result$p)) {
        if (levene_test_result$p < 0.05) {
          t_test_result <- t.test(Result ~ Group,
                                  data = df_pair,
                                  var.equal = FALSE)
        } else {
          t_test_result <- t.test(Result ~ Group,
                                  data = df_pair,
                                  var.equal = TRUE)
        }
        
        
        pval <- t_test_result$p.value
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

data_long_ip_no_hem <- data_long_ip %>% filter(is.na(Flags))
data_long_ip_hem <- data_long_ip %>% filter(Flags == "HEM+")

data_long_or_no_hem <- data_long_or %>% filter(is.na(Flags))
data_long_or_hem <- data_long_or %>% filter(Flags == "HEM+")

make_plot <- function(data, title_suffix) {
  ggplot(data, aes(x = Group, y = Result)) +
    geom_boxplot(
      outlier.shape = 16,
      outlier.size = 1,
      color = "black"
    ) +
    facet_wrap( ~ Analyte, scales = "free_y", ncol = 6) +
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
    labs(
      y = "Value",
      x = "Group",
      fill = "Group",
      title = paste("Analytes", title_suffix)
    )
}

### saving plots

# without flags
ggsave(
  "plot_stat_ip.png",
  plot = plot_stat_ip,
  width = 12,
  height = 8
)
ggsave(
  "plot_stat_or.png",
  plot = plot_stat_or,
  width = 12,
  height = 8
)

# IP
plot_ip_no_hem <- make_plot(data_long_ip_no_hem, "(IP, without HEM+)")
plot_ip_hem <- make_plot(data_long_ip_hem, "(IP, HEM+)")

ggsave(
  "plot_stat_ip_no_hem.png",
  plot = plot_ip_no_hem,
  width = 12,
  height = 8,
  dpi = 300
)
ggsave(
  "plot_stat_ip_hem.png",
  plot = plot_ip_hem,
  width = 12,
  height = 8,
  dpi = 300
)

# OR
plot_or_no_hem <- make_plot(data_long_or_no_hem, "(OR, without HEM+)")
plot_or_hem <- make_plot(data_long_or_hem, "(OR, HEM+)")

ggsave(
  "plot_stat_or_no_hem.png",
  plot = plot_or_no_hem,
  width = 12,
  height = 8,
  dpi = 300
)
ggsave(
  "plot_stat_or_hem.png",
  plot = plot_or_hem,
  width = 12,
  height = 8,
  dpi = 300
)

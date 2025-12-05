#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(scales)
  library(emmeans) })

input_dir <- "KRAKEN_VALID_FILES"
fastq_files <- list.files(input_dir, pattern = ".gz$", full.names = TRUE)
count_fastq_reads <- function(filepath) {
  # count lines beginning with '@' or total lines divided by 4
  cmd <- sprintf("zgrep -c '^@' '%s' || zcat '%s' | wc -l", filepath, filepath)
  line_count <- suppressWarnings(as.numeric(system(cmd, intern = TRUE)[1]))
  reads <- round(line_count / 4)
  return(reads) }

# Parallel counting (optional, if parallel package available)
read_counts <- sapply(fastq_files, count_fastq_reads)

# 🧩 Build dataframe
data <- data.frame(
  filepath = fastq_files,
  filename = basename(fastq_files),
  reads = read_counts,
  stringsAsFactors = FALSE )

# 🧬 Extract metadata
data <- data %>%
  mutate(
    species = ifelse(str_detect(filename, "mouse"), "mouse", "VEEV"),
    # clean sample name (remove pair number + extension)
    sample_base = filename %>% # Control_RNA_1_S3_1.mouse.fq.gz
      str_replace(".mouse", "") %>%
      str_replace(".VEEV", "") %>% 
      str_replace(".fq.gz", "") %>%
      str_replace("_RNA_", "_")  )

# 🧩 Merge read pairs (1 and 2)
data_merged <- data %>%
  group_by(sample_base, species) %>%
  summarise(total_reads = sum(reads), .groups = "drop")

veev_mouse_df <- data_merged %>%
  pivot_wider(names_from = species, values_from = total_reads) %>%
  mutate(
    VEEV = ifelse(is.na(VEEV), 0, VEEV),
    mouse = ifelse(is.na(mouse), 0, mouse),
    veev_mouse_ratio = ifelse(mouse > 0, VEEV / mouse, NA),
    percent_veev = ifelse(VEEV + mouse > 0, 100 * VEEV / (VEEV + mouse), NA)  )

# 🧩 Assign groups (Control / WT / 4X)
veev_mouse_df <- veev_mouse_df %>%
  mutate(
    group = case_when(
      str_detect(sample_base, "Control") ~ "Control",
      str_detect(sample_base, "WT") ~ "WT",
      str_detect(sample_base, "4X") ~ "4X",
      TRUE ~ "Other"    ),
    replicate = str_extract(sample_base, "_(\\d+)") %>% str_remove("_"),
    replicate = as.factor(replicate)  )

# 💾 Save count_host.txt
write.table(  veev_mouse_df %>%
    select(sample_id = sample_base, viral_reads = VEEV, mouse_reads = mouse, percent_viral = percent_veev),
  file = "count_host.txt",  sep = "\t", quote = FALSE, row.names = FALSE )

# 🧮 Summary stats per group
summary_stats <- veev_mouse_df %>%  group_by(group) %>%
  summarise(    n = n(),
    mean_ratio = round(mean(veev_mouse_ratio, na.rm = TRUE), 4),
    median_ratio = round(median(veev_mouse_ratio, na.rm = TRUE), 4),
    sd_ratio = round(sd(veev_mouse_ratio, na.rm = TRUE), 4),
    mean_percent = round(mean(percent_veev, na.rm = TRUE), 2),
    .groups = "drop"  )

print(summary_stats)
write.table(summary_stats, "VEEV_mouse_summary.txt", sep = "\t", quote = FALSE, row.names = F)

# 🎨 Plot setup
cb_palette <- c("Control" = "#999999", "WT" = "#0072B2", "4X" = "#E69F00")

p1 <- ggplot(veev_mouse_df, aes(x = group, y = percent_veev, color = group)) +
  geom_boxplot(alpha = 0.6, width = 0.6 ) +
  geom_jitter(width = 0.2, size = 4, alpha = 0.8) +
  scale_color_manual(values = cb_palette) +
  labs(    x = "", y = "VEEV reads (%)", color = "Group" ) + 
  theme_minimal(base_size = 15) +
  theme(    axis.text.x = element_text(angle =90, hjust = 1),
    legend.position = "right"  )

ggsave("VEEV_mouse_boxplot.pdf", p1, width =4, height = 5)

p1 <- ggplot(veev_mouse_df, aes(x = group, y = percent_veev, color = group)) +
  geom_boxplot(alpha = 0.6, width = 0.6 ) +
  geom_jitter(width = 0.2, size = 4, alpha = 0.8) +
  scale_color_manual(values = cb_palette) + ylim(0,15.6) + 
  labs(    x = "", y = "VEEV reads (%)", color = "Group" ) + 
  theme_minimal(base_size = 15) +
  theme(    axis.text.x = element_text(angle =90, hjust = 1),
    legend.position = "right"  )

ggsave("VEEV_mouse_boxplot.ylim.pdf", p1, width =4, height = 5)

# 🧪 ANOVA and post-hoc
anova_model <- aov(percent_veev ~ group + replicate, data = veev_mouse_df)
summary(anova_model)
print("emmeans")
em <- emmeans(anova_model, ~ group)
pairwise_comparisons <- pairs(em)
print(pairwise_comparisons)
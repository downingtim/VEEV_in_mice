library(tidyverse)
library(ggplot2)

# Read all RDAF files (excluding _7141_)
files <- list.files(pattern = '_rdaf.csv') %>%
  str_subset('_7141_', negate = TRUE)

rdaf_data <- map_df(files, ~read_csv(.x, show_col_types = FALSE) %>%
                      mutate(filename = .x))

# Extract condition and replicate from filename
rdaf_data <- rdaf_data %>%
  mutate(
    condition = case_when(
      grepl('Control', filename) ~ 'Control',
      grepl('VEEV_WT', filename) ~ 'WT',
      grepl('VEEV_4X', filename) ~ '4X'
    ),
    replicate = case_when(
      grepl('_1_', filename) ~ '1',
      grepl('_2_', filename) ~ '2',
      grepl('_3_', filename) ~ '3'
    ),
    sample_label = paste0(condition, '_', replicate)
  ) %>%
  filter(!is.na(condition))

# Create color palette (9 distinct colors)
colors <- c(
  'Control_1' = '#A9A9A9',
  'Control_2' = '#808080',
  'Control_3' = '#505050',
  'WT_1' = '#B3D9FF',
  'WT_2' = '#4292C6',
  'WT_3' = '#08519C',
  '4X_1' = '#FFE699',
  '4X_2' = '#FE9929',
  '4X_3' = '#D95F0E'
)

rdaf_data <- rdaf_data %>%  filter(rdaf > 0.05)

# Create histogram with overlapping distributions
p_hist <- ggplot(rdaf_data, aes(x = rdaf, fill = sample_label)) +
  geom_histogram(binwidth = 0.05, alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = colors) +
  labs(x = 'Read depth allele frequency',
       y = 'Count',       fill = 'Sample') +
  theme_minimal() +
    theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
      panel.grid.major = element_line(color = 'grey70', linewidth = 0.2),
      panel.grid.minor = element_line(color = 'grey90', linewidth = 0.04),
        legend.position = 'right')
ggsave('rdaf_histogram.pdf', p_hist, width =8, height =5)

# Optional: density plot instead
p_density <- ggplot(rdaf_data, aes(x = rdaf, fill = sample_label)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = colors) +
  labs(x = 'Read depth allele frequency',
       y = 'Frequency',
       fill = 'Sample') +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'right')
ggsave('rdaf_density.pdf', p_density, width = 10, height = 6)
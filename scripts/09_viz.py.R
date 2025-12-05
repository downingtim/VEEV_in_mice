library(tidyverse)
library(ggplot2)
library(ggExtra)

# Read all RDAF files (excluding _7141_)
files <- list.files(pattern = 'daf.csv') %>% str_subset('_7141_', negate=T)

rdaf_data <- map_df(files, ~read_csv(.x, show_col_types = FALSE) %>% 
                      mutate(filename = .x))

# Extract condition and replicate from filename
rdaf_data <- rdaf_data %>%
  mutate(    condition = case_when(
      grepl('Control', filename) ~ 'Control',
      grepl('VEEV_WT', filename) ~ 'WT',
      grepl('VEEV_4X', filename) ~ '4X'    ),
    replicate = case_when(
      grepl('_1_', filename) ~ '1',
      grepl('_2_', filename) ~ '2',
      grepl('_3_', filename) ~ '3'    ),
    sample_label = paste0(condition, '_', replicate)
  ) %>%  filter(!is.na(condition))

# Determine genome length
genome_length <- 11442  # 68U201_v vs 14525

# Create complete pos x sample grid
all_samples <- unique(rdaf_data$sample_label)
complete_grid <- expand_grid(
  pos = 1:genome_length,
  sample_label = all_samples ) %>%
  left_join(rdaf_data %>% select(pos, sample_label, depth, rdaf, condition),
            by = c('pos', 'sample_label'))

# Fill in missing depth values as 0
complete_grid <- complete_grid %>% mutate(depth = replace_na(depth, 0))

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
  '4X_3' = '#D95F0E')
  # pos,ref,depth,rdaf,A,C,G,T

p2 <- ggplot(complete_grid, aes(x = pos, y = (depth), color = sample_label)) +
      geom_smooth(se = FALSE, span = 0.15, method = 'loess', linewidth=1.8, alpha=0.6) +
    scale_x_continuous(limits = c(7500, 8400)) +
  scale_y_continuous(limits = c(0, 80))+ scale_color_manual(values = colors) +
  labs(x = 'Genome Position (bp)',
       y = 'Number of reads', color = 'Sample') +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = 'grey80', linewidth = 0.4),
        panel.grid.minor = element_line(color = 'grey95', linewidth = 0.2),
        legend.position = 'right')
ggsave('rdaf_depth_combined.py2.pdf', p2, width = 9, height =4)

p22 <- ggplot(complete_grid, aes(x = pos, y = (depth), color = sample_label)) +
      geom_smooth(se = FALSE, span = 0.2, method = 'loess', linewidth=1.8, alpha=0.6) +
    scale_x_continuous(limits = c(7540, 7900)) +
  scale_y_continuous(limits = c(0, 80))+ scale_color_manual(values = colors) +
  labs(x = 'Genome Position (bp)',
       y = 'Number of reads', color = 'Sample') +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = 'grey80', linewidth = 0.4),
        panel.grid.minor = element_line(color = 'grey95', linewidth = 0.2),
        legend.position = 'right')
ggsave('rdaf_depth_combined.py3.pdf', p22, width = 9, height =4)


# Define genes with adjusted coordinates (-16)
genes <- tibble(
  xmin = c(7562-16, 8387-16, 8564-16, 9833-16, 10001-16),
  xmax = c(8386-16, 8563-16, 9832-16, 10000-16, 11326-16),
  gene = c("capsid", "E3", "E2", "6K", "E1") )

p3 <- ggplot(complete_grid, aes(x = pos, y = depth, color = sample_label)) +
geom_rect(data = genes, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill = 'grey90', alpha = 0.8, color = 'black', linewidth = 0.5, inherit.aes = FALSE) +
	    geom_smooth(se = FALSE, span = 0.1, method = 'loess', linewidth = 1.8, alpha = 0.6) +
  geom_text(data = genes, aes(x = (xmin + xmax) / 2, y = 70, label = gene),
            color = 'black', size = 5, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(7500, 11400)) +
  scale_y_continuous(limits = c(0, 80)) +
  scale_color_manual(values = colors) +
  labs(x = 'Genome Position (bp)',
       y = 'Number of reads', color = 'Sample') +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = 'grey80', linewidth = 0.4),
        panel.grid.minor = element_line(color = 'grey95', linewidth = 0.2),
        legend.position = 'right')

ggsave('rdaf_depth_combined.py.pdf', p3, width = 9, height = 4)

p2 <- ggplot(complete_grid, aes(x = pos, y = (rdaf), color = sample_label)) +
   geom_point(alpha=0.6) + 
#      geom_smooth(se = FALSE, span = 0.15, method = 'loess', linewidth=1.8, alpha=0.6) +
    scale_x_continuous(limits = c(7500, 8400)) +
    scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = colors) +
  labs(x = 'Genome Position (bp)',
       y = 'Read depth allele frequency', color = 'Sample') +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = 'grey80', linewidth = 0.4),
        panel.grid.minor = element_line(color = 'grey95', linewidth = 0.2),
        legend.position = 'right')
ggsave('rdaf_combined.py2.pdf', p2, width = 9, height =4)

p22 <- ggplot(complete_grid, aes(x = pos, y = (rdaf), color = sample_label)) +
#      geom_smooth(se = FALSE, span = 0.1, method = 'loess', linewidth=1.8, alpha=0.6) +
    scale_x_continuous(limits = c(1, 11440)) +
    scale_y_continuous(limits = c(0,1)) +
    geom_point(alpha=0.6) +
  scale_color_manual(values = colors) +
  labs(x = 'Genome Position (bp)',
       y = 'Read depth allele frequency', color = 'Sample') +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = 'grey80', linewidth = 0.4),
        panel.grid.minor = element_line(color = 'grey95', linewidth = 0.2),
        legend.position = 'right')
ggsave('rdaf_combined.py_genome.pdf', p22, width = 9, height =4)


# Define genes with adjusted coordinates (-16)
genes <- tibble(
  xmin = c(7562-16, 8387-16, 8564-16, 9833-16, 10001-16),
  xmax = c(8386-16, 8563-16, 9832-16, 10000-16, 11326-16),
  gene = c("capsid", "E3", "E2", "6K", "E1") )

p3 <- ggplot(complete_grid, aes(x = pos, y = rdaf, color = sample_label)) +
   geom_rect(data = genes, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill = 'grey90', alpha = 0.8, color = 'black', linewidth = 0.5, inherit.aes = FALSE) +
geom_point(alpha=0.6 ) +
#	    geom_smooth(se = FALSE, span = 0.1, method = 'loess', linewidth = 1.8, alpha = 0.6) +
  geom_text(data = genes, aes(x = (xmin + xmax) / 2, y = 0.9, label = gene),
            color = 'black', size = 5, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(7500, 11400)) +
    scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = colors) +
  labs(x = 'Genome Position (bp)',
       y = 'Read depth allele frequency', color = 'Sample') +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = 'grey80', linewidth = 0.4),
        panel.grid.minor = element_line(color = 'grey95', linewidth = 0.2),
        legend.position = 'right')
ggsave('rdaf_combined.py.pdf', p3, width = 9, height = 4)

## histogram


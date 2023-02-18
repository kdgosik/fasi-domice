#' @title Figure 2A Cell Proportions
#' @author Kirk Gosik
#' @description
#'
#'

library(data.table)
library(tidyverse)

figure_path <- "results/figures/"
data_path <- "data/"
source(paste0(figure_path, "helpers.R"))
plot_df <- fread(paste0(data_path, "manuscript-plot-data.csv"), data.table = FALSE)


df_to_plot <- plot_df %>%
  mutate(BestCall = factor(BestCall, levels = sample(unique(BestCall)))) %>%
  filter(grepl("ILC", called_cell_types_new)) %>% 
  count(BestCall, called_cell_types_new) %>%
  group_by(BestCall) %>%
  mutate(percent = n / sum(n),
         order = dplyr::first(percent))


p1 <- ggplot(df_to_plot, aes(reorder(BestCall, order), n, fill = called_cell_types_new)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          legend.position = "top") +
    labs(title = "ILC Proportion by Mouse",
         x = "Mouse ID",
         y = "Proportion",
         fill = "Cell Type")

ggsave(filename = paste0(figure_path, "Figure-2A-ilc-proportion-bars.png"),
       plot = p1,
       dpi = 330)

ggsave(filename = paste0(figure_path, "Figure-2A-ilc-proportion-bars.pdf"),
       plot = p1,
       dpi = 330)


p2_density <- ggplot(df_to_plot, aes(percent, fill = called_cell_types_new, group = called_cell_types_new)) + 
  geom_density(alpha = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "top") +
  labs(title = "ILC Proportion by Mouse",
       x = "Percent",
       y = "Density",
       fill = "Cell Type")

ggsave(filename = paste0(figure_path, "Figure-2A-ilc-proportion-density.png"),
       plot = p2_density,
       dpi = 330)

ggsave(filename = paste0(figure_path, "Figure-2A-ilc-proportion-density.pdf"),
       plot = p2_density,
       dpi = 330)

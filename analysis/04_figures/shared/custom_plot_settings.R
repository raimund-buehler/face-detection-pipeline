# Manuscript section: Shared figure helper
# Analysis family: shared helper
# Original source path: scripts/publication/format_markdown/custom_plot_settings.R
# Primary input dataset(s): figure scripts across manuscript analyses
# Primary output(s): reusable ggplot formatting helpers
# Known TODOs: copied helper is not yet wired into copied analysis scripts; original helper path remains the active source path
# Scientific logic note: copied from source without changing scientific logic

library(ggplot2)
library(grid)

# Simplified Custom Settings Function without Aesthetic Detection
apply_custom_settings <- function(
  values = c("#FDC010", "#0B7A6B"),
  base_size = 10,
  axis_text_size = 10,
  axis_title_size = 12,
  title_size = axis_title_size,
  legend_text_size = 10,
  strip_text_size = axis_text_size,
  panel_spacing = unit(1, "lines")
) {
  list(
    theme_classic(base_size = base_size) +
      theme(
        legend.position = "top",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = title_size),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black", size = strip_text_size),
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.title = element_blank(),

        # Ensure axis lines are shown across all facets
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),

        # Remove panel border to keep only classic axis lines
        panel.border = element_blank(),

        # Adjust facet spacing to create a balanced layout
        panel.spacing = panel_spacing
      ),

    # Customize color and fill scales
    scale_color_manual(values = values),
    scale_fill_manual(values = values),

    # Customize legend guides
    guides(
      color = guide_legend(override.aes = list(
        fill = values,
        color = values,   # Remove border color
        alpha = 1,
        shape = 15,
        size = 5
      )),
      fill = guide_legend(override.aes = list(
        fill = values,
        color = NA,   # Remove border color
        alpha = 1,
        shape = 15,
        size = 5
      ))
    )
  )
}



# Manuscript section: Shared figure helper
# Analysis family: shared helper
# Original source path: scripts/publication/format_markdown/custom_theme_settings_publication.R
# Primary input dataset(s): publication-style figure scripts
# Primary output(s): reusable theme settings
# Known TODOs: determine whether this helper is still actively used by manuscript scripts
# Scientific logic note: copied from source without changing scientific logic

#' Publication plotting helpers
#'
#' Utility functions to standardise theming and palettes for publication figures.
#'
#' @param base_size Base font size passed to ggplot2::theme_classic().
#' @param font_face Font face applied to textual elements.
#' @param values Character vector of colour hex codes.
#' @param position Legend position (passed to ggplot2::theme()).
#' @param text_size Legend text size.
#' @param aesthetics A vector of aesthetics that should receive manual palettes.
#' @param name Optional legend title.
#' @param drop_guides Logical indicating whether guides should be stripped.
#' @param keep_guides A character vector of aesthetics whose guides should be retained.
#' @param plot A ggplot object to modify.
#'
#' @return A ggplot theme, scale, or modified plot.
#' @keywords internal
NULL

publication_panel_theme <- function(base_size = 7.2, font_face = "plain") {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text                = ggplot2::element_text(face = font_face),
      plot.title.position = "plot",
      plot.title          = ggplot2::element_text(
        size   = 8,
        face   = font_face,
        hjust  = 0.6,
        margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
      ),
      axis.title.x        = ggplot2::element_text(
        size   = 8,
        face   = font_face,
        margin = ggplot2::margin(t = 1)
      ),
      axis.text.x         = ggplot2::element_text(
        size   = 6,
        face   = font_face,
        margin = ggplot2::margin(t = 0.5)
      ),
      axis.title.y        = ggplot2::element_text(
        size   = 8,
        face   = font_face,
        margin = ggplot2::margin(t = 1)
      ),
      axis.text.y         = ggplot2::element_text(
        size   = 6,
        face   = font_face,
        margin = ggplot2::margin(t = 1, r = 1)
      ),
      strip.text          = ggplot2::element_text(size = base_size, face = font_face),
      strip.background    = ggplot2::element_blank(),
      panel.border        = ggplot2::element_blank(),
      panel.spacing       = grid::unit(1, "pt"),
      plot.margin         = ggplot2::margin(t = 1, r = 1, b = 1, l = 1)
    )
}

publication_palette <- function(values = c("#FDC010", "#0B7A6B"),
                                 aesthetics = c("fill", "colour"),
                                 name = NULL) {
  aesthetics <- match.arg(aesthetics, several.ok = TRUE)

  scales <- list()

  if ("fill" %in% aesthetics) {
    scales <- c(scales, list(ggplot2::scale_fill_manual(values = values, name = name)))
  }

  if ("colour" %in% aesthetics || "color" %in% aesthetics) {
    # treat US/UK spelling interchangeably
    scales <- c(scales, list(ggplot2::scale_colour_manual(values = values, name = name)))
  }

  guide_overrides <- list(
    fill = ggplot2::guide_legend(
      title          = name,
      nrow           = 1,
      byrow          = TRUE,
      label.position = "right",
      override.aes   = list(shape = 22, size = 4, colour = NA, alpha = 1)
    )
  )

  if ("colour" %in% aesthetics || "color" %in% aesthetics) {
    guide_overrides$colour <- ggplot2::guide_legend(
      title        = name,
      override.aes = list(fill = values, colour = values, alpha = 1, shape = 22, size = 4)
    )
  }

  scales <- c(scales, list(do.call(ggplot2::guides, guide_overrides)))
  scales
}

strip_all_guides <- function(plot, keep_guides = NULL) {
  drop_guides <- setdiff(c("colour", "color", "fill", "linetype", "shape", "alpha", "size", "stroke"), keep_guides)
  if (length(drop_guides) == 0) {
    return(plot)
  }

  guide_args <- rep(list("none"), length(drop_guides))
  names(guide_args) <- drop_guides
  plot + do.call(ggplot2::guides, guide_args)
}

style_publication_panel <- function(plot,
                                    base_size = 7.2,
                                    font_face = "plain",
                                    drop_guides = FALSE,
                                    keep_guides = NULL) {
  styled_plot <- plot + publication_panel_theme(base_size = base_size, font_face = font_face)

  if (isTRUE(drop_guides)) {
    styled_plot <- strip_all_guides(styled_plot, keep_guides = keep_guides)
  }

  styled_plot
}

publication_legend_theme <- function(position = "bottom", text_size = 8, direction = "horizontal") {
  ggplot2::theme(
    legend.position       = position,
    legend.direction      = direction,
    legend.box            = "horizontal",
    legend.background     = ggplot2::element_blank(),
    legend.box.background = ggplot2::element_blank(),
    legend.key            = ggplot2::element_blank(),
    legend.box.spacing    = grid::unit(2, "pt"),
    legend.margin         = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
    legend.title          = ggplot2::element_blank(),
    legend.text           = ggplot2::element_text(size = text_size, face = "plain")
  )
}

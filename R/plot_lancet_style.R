# ============================================================
# LANCET STYLE HELPERS
# Reusable utilities for all figures in the thesis
# ============================================================

library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 1. Middle-dot decimal formatter (Lancet house style)
#    Converts "." to "·" in axis labels
#    Usage: scale_y_continuous(labels = lancet_num())
# ------------------------------------------------------------
lancet_num <- function(accuracy = NULL, ...) {
  function(x) {
    lbl <- scales::number(x, accuracy = accuracy, ...)
    gsub("\\.", "\u00b7", lbl)
  }
}

# Percentage version (e.g. 25% stays 25%, but 25.5% -> 25·5%)
lancet_pct <- function(accuracy = 1, ...) {
  function(x) {
    lbl <- scales::percent(x, accuracy = accuracy, ...)
    gsub("\\.", "\u00b7", lbl)
  }
}

# ------------------------------------------------------------
# 2. Lancet-style theme
#    Minimalist, clean, classic axes, no grid
# ------------------------------------------------------------
theme_lancet <- function(base_size = 10, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Text
      plot.title       = element_text(size = rel(1.05), face = "plain",
                                      hjust = 0, margin = ggplot2::margin(b = 6)),
      plot.subtitle    = element_text(size = rel(0.95), colour = "grey30",
                                      hjust = 0, margin = ggplot2::margin(b = 8)),
      axis.title       = element_text(size = rel(1.00), colour = "grey10"),
      axis.title.x     = element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y     = element_text(margin = ggplot2::margin(r = 6)),
      axis.text        = element_text(size = rel(0.90), colour = "grey20"),

      # Axes
      axis.line        = element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks       = element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),

      # No grid (Lancet prefers clean background)
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),

      # Legend
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = rel(0.90)),
      legend.key       = element_rect(fill = "white", colour = NA),
      legend.key.width = unit(0.9, "cm"),
      legend.key.height= unit(0.4, "cm"),
      legend.margin    = ggplot2::margin(0, 0, 2, 0),
      legend.box.margin= ggplot2::margin(0, 0, 0, 0),

      # Facets (for panel figures)
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(size = rel(0.95), colour = "grey10",
                                      margin = ggplot2::margin(4, 4, 4, 4)),

      # Plot area
      plot.margin      = ggplot2::margin(10, 12, 10, 10)
    )
}

# ------------------------------------------------------------
# 3. Lancet-inspired colour palette
#    Based on ggsci::pal_lancet() with slight tuning
#    10 colours — enough for all 10 models
# ------------------------------------------------------------
lancet_palette <- c(
  "#00468B",  # deep blue      — DRF
  "#ED0000",  # red            — QRF
  "#42B540",  # green          — RF
  "#0099B4",  # teal           — Ensemble DRF
  "#925E9F",  # purple         — Hybrid NN-DRF
  "#FDAF91",  # salmon         — XGBoost
  "#AD002A",  # dark red       — Deep Ensemble
  "#ADB6B6",  # grey           — MC Dropout
  "#1B1919",  # near-black     — BNN
  "#E18727"   # orange         — DDNN
)

# Named version — if you want each model to always get the same colour
# across all figures, set the name = model key used in your plotting data
model_colours <- c(
  "DRF"          = "#00468B",
  "QRF"          = "#ED0000",
  "RF"           = "#42B540",
  "Ensemble-DRF" = "#0099B4",
  "Hybrid"       = "#925E9F",
  "XGBoost"      = "#FDAF91",
  "ENS"          = "#AD002A",
  "MCD"          = "#ADB6B6",
  "BNN"          = "#1B1919",
  "DDNN"         = "#E18727"
)

# ------------------------------------------------------------
# 4. Scale wrappers
# ------------------------------------------------------------
scale_colour_lancet <- function(...) {
  ggplot2::discrete_scale("colour", "lancet",
    palette = function(n) lancet_palette[seq_len(n)], ...)
}
scale_fill_lancet <- function(...) {
  ggplot2::discrete_scale("fill", "lancet",
    palette = function(n) lancet_palette[seq_len(n)], ...)
}

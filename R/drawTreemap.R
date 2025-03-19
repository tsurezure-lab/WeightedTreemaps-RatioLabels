#' drawTreemap
#'
#' Draws the treemap object that was obtained by running \code{\link{voronoiTreemap}} or
#' \code{\link{sunburstTreemap}}. Many graphical parameters can be customized, but some
#' settings that determine the appearance of treemaps (e.g., cell size and initial shape)
#' are set during treemap generation.
#'
#' @param treemap (treemapResult) Either a \code{voronoiResult} or \code{sunburstResult}
#'   object from \code{\link{voronoiTreemap}} or \code{\link{sunburstTreemap}}.
#' @param levels (numeric) Vector of hierarchical levels to draw (default: all levels).
#' @param color_type (character) One of "categorical", "cell_size", "both", or "custom_color".
#'   Determines how cells are colored (see Details).
#' @param color_level (numeric) Level for cell coloring (default: lowest level for Voronoi, all for Sunburst).
#' @param color_palette (character) Vector of colors for filling cells (default: \code{\link[colorspace]{rainbow_hcl}}).
#' @param color_steps (numeric) Number of steps for color gradient when \code{color_type = "cell_size"} (default: 10).
#' @param border_level (numeric) Levels for drawing cell borders (default: all levels).
#' @param border_size (numeric) Line width for borders (default: 6, scaled by level for Voronoi).
#' @param border_color (character) Color for borders (default: light grey).
#' @param label_level (numeric) Levels for drawing labels (default: deepest level).
#' @param label_size (numeric) Relative size of labels (default: 1, or vector per level).
#' @param label_color (character) Color for labels (default: light grey).
#' @param label_fontfamily (character) Font family for labels and title (default: "sans").
#' @param label_fontweight (character) Font weight for labels ("normal" or "bold", default: "normal").
#' @param label_autoscale (logical) Automatically scale labels by cell area (default: TRUE).
#' @param label_line_spacing (numeric or numeric vector) Line spacing multiplier for multi-line labels (default: 0.5, range 0.1–1.0). If a numeric vector, specify values for each level.
#' @param label_ratio (logical or numeric) Either a logical value or a numeric vector. If TRUE, draw ratio labels for all levels; if FALSE, draw no ratio labels; if a numeric vector, draw ratio labels only for the specified levels (default: FALSE).
#' @param label_ratio_color (character or character vector) Color for ratio labels (default: light grey). If a character vector, specify colors for each level.
#' @param label_ratio_size (numeric) Size for ratio labels (default: 1).
#' @param label_ratio_fontfamily (character) Font family for ratio labels (default: "sans").
#' @param label_ratio_fontweight (character) Font weight for ratio labels ("normal" or "bold", default: "normal").
#' @param title (character) Optional title for the treemap (default: NULL).
#' @param title_size (numeric) Size of the title (default: 1).
#' @param title_color (character) Color for the title (default: grey(0.5)).
#' @param title_fontweight (character) Font weight for the title ("normal" or "bold", default: "normal").
#' @param legend (logical) Draw a color legend (default: FALSE).
#' @param legend_position (character) Position of the legend ("left", "right", "top", "bottom", default: "left").
#' @param legend_size (numeric) Relative size of the legend (0 to 1, default: 0.1).
#' @param custom_range (numeric) Vector of length 2 to rescale custom color values (default: NULL).
#' @param width (numeric) Width of the viewport (0 to 0.9, default: 0.9).
#' @param height (numeric) Height of the viewport (0 to 0.9, default: 0.9).
#' @param layout (numeric) Vector of length 2 for rows and columns in the plot area (default: c(1, 1)).
#' @param position (numeric) Vector of length 2 for the treemap position (default: c(1, 1)).
#' @param add (logical) Add to existing plot instead of creating a new page (default: FALSE).
#' @param y_offset (numeric) Offset for vertical positioning of the treemap (default: 0, negative for downward movement, e.g., -0.1 for 10% downward, positive for upward movement. Values should typically be between -0.5 and 0.5 to stay within the viewport).
#'
#' @return NULL (creates a grid viewport and draws the treemap).
#'
#' @details
#' - `color_type` options:
#'   - "categorical": Colors cells by category (may repeat colors).
#'   - "cell_size": Colors by relative cell area.
#'   - "both": Combines categorical coloring with lightness based on area.
#'   - "custom_color": Uses a custom color index from \code{custom_color} in \code{voronoiTreemap}.
#' - Use \code{\link{voronoiTreemap}} to generate the treemap input.
#'
#' @examples
#' library(WeightedTreemaps)
#' df <- data.frame(
#'   A = rep(c("abcd", "efgh"), each = 4),
#'   B = letters[1:8],
#'   size = c(37, 52, 58, 27, 49, 44, 34, 45),
#'   primary_cluster_ratio = c(0.1, 0.2, 0.15, 0.05, 0.25, 0.1, 0.05, 0.1),
#'   secondary_cluster_ratio = c(0.3, 0.4, 0.2, 0.1, 0.5, 0.2, 0.1, 0.2)
#' )
#' tm <- voronoiTreemap(
#'   data = df, levels = c("A", "B"), cell_size = "size",
#'   shape = "circle", positioning = "regular", seed = 123,
#'   label_ratios = c("primary_cluster_ratio", "secondary_cluster_ratio")
#' )
#' drawTreemap(tm, label_size = 1, color_type = "categorical",
#'   label_ratio = c(1, 2), label_ratio_color = c("red", "blue"),
#'   label_fontfamily = "zenmaru", label_fontweight = "bold",
#'   label_line_spacing = c(0.3, 0.7), label_ratio_fontweight = "bold",
#'   title_fontweight = "normal", y_offset = -0.1)
#'
#' @importFrom dplyr %>%
#' @importFrom grid grid.newpage grid.text grid.polygon grid.draw grid.layout gpar unit viewport pushViewport popViewport
#' @importFrom colorspace rainbow_hcl
#' @importFrom scales rescale
#' @importFrom lattice draw.colorkey
#' @importFrom grDevices colorRampPalette grey
#' @importFrom stats setNames
#' @importFrom utils tail
#'
#' @export drawTreemap
#'
drawTreemap <- function(
  treemap,
  levels = 1:length(treemap@call$levels),
  color_type = "categorical",
  color_level = NULL,
  color_palette = NULL,
  color_steps = 10,
  border_level = levels,
  border_size = 6,
  border_color = grey(0.9),
  border_alpha = 0.6,
  label_level = max(levels),
  label_size = 1,
  label_color = grey(0.9),
  label_fontfamily = "sans",
  label_fontweight = "normal",
  label_autoscale = TRUE,
  label_line_spacing = 0.5,
  label_ratio = FALSE,
  label_ratio_color = grey(0.9),
  label_ratio_size = 1,
  label_ratio_fontfamily = "sans",
  label_ratio_fontweight = "normal",
  title = NULL,
  title_size = 1,
  title_color = grey(0.5),
  title_fontweight = "normal",
  legend = FALSE,
  legend_position = "left",
  legend_size = 0.1,
  custom_range = NULL,
  width = 0.9,
  height = 0.9,
  layout = c(1, 1),
  position = c(1, 1),
  add = FALSE,
  y_offset = 0  # 新しいパラメータ: Y 軸方向のオフセット（負で下方移動）
) {
  # 入力データの検証
  validate_treemap(treemap,
    width, height, layout, position, add,
    levels, color_level, border_level,
    label_level, color_palette,
    border_color, label_color,
    custom_range, title)

  # 色付けのレベルを決定（Sunburst/Voronoiに応じて）
  if (is.null(color_level)) {
    if (inherits(treemap, "sunburstResult")) {
      color_level = levels
    } else {
      color_level = min(levels)
    }
  }

  # label_ratio の処理
  if (is.logical(label_ratio)) {
    if (label_ratio) {
      label_ratio_levels <- levels
    } else {
      label_ratio_levels <- numeric(0)
    }
  } else if (is.numeric(label_ratio)) {
    label_ratio_levels <- label_ratio
    if (!all(label_ratio_levels %in% levels)) {
      warning("Some values in label_ratio are not valid levels. Ignoring invalid levels.")
      label_ratio_levels <- label_ratio_levels[label_ratio_levels %in% levels]
    }
  } else {
    stop("label_ratio must be a logical value or a numeric vector.")
  }

  # label_line_spacing の処理
  if (length(label_line_spacing) == 1) {
    label_line_spacing_per_level <- rep(label_line_spacing, max(levels))
  } else if (is.numeric(label_line_spacing)) {
    if (length(label_line_spacing) != max(levels)) {
      stop("length of label_line_spacing must match the number of levels.")
    }
    label_line_spacing_per_level <- label_line_spacing
  } else {
    stop("label_line_spacing must be a numeric value or vector.")
  }

  # label_ratio_color の処理
  if (length(label_ratio_color) == 1) {
    label_ratio_color_per_level <- rep(label_ratio_color, max(levels))
  } else if (is.character(label_ratio_color)) {
    if (length(label_ratio_color) != max(levels)) {
      stop("length of label_ratio_color must match the number of levels.")
    }
    label_ratio_color_per_level <- label_ratio_color
  } else {
    stop("label_ratio_color must be a character or character vector.")
  }

  # 新しいページを開始（add = FALSE の場合）
  if (!add) {
    grid::grid.newpage()
  }

  # メインのグリッドビューを作成（layoutで分割可能）
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(nrow = layout[1], ncol = layout[2])
    )
  )

  # ツリーマップを描画するサブビューを作成
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = position[1],
      layout.pos.col = position[2]
    )
  )

  # ツリーマップの描画エリアを定義（マージンや凡例、タイトルを考慮）
  if (!legend) { legend_position <- "none" }
  key_offsets <- switch(
    legend_position,
    left = c(0.5 + legend_size/2, 0.5, width - legend_size, height, legend_size/2, 0.5, legend_size, height),
    right = c(0.5 - legend_size/2, 0.5, width - legend_size, height, 1 - legend_size/2, 0.5, legend_size, height),
    top = c(0.5, 0.5 - legend_size/2, width, height - legend_size, 0.5, 1 - legend_size/2, width, legend_size),
    bottom = c(0.5, 0.5 + legend_size/2, width, height - legend_size, 0.5, legend_size/2, width, legend_size),
    none = c(0.5, 0.5, width, height, 0.5, 0.5, 0, 0)
  )

  if (!is.null(title)) {
    title_offset <- 0.05
  } else {
    title_offset <- 0
  }

  # y_offset を適用してビューポートを移動
  grid::pushViewport(
    grid::viewport(
      x = key_offsets[1],
      y = key_offsets[2] - title_offset/2 + y_offset,  # y_offset を適用
      width = key_offsets[3],
      height = key_offsets[4] - title_offset,
      xscale = c(0, 2000),
      yscale = c(0, 2000)
    )
  )

  # ポリゴンの描画（色付け）
  treemap <- add_color(treemap, color_palette, color_type,
    color_level, color_steps, custom_range)
  lapply(treemap@cells, function(tm_slot) {
    if (tm_slot$level %in% color_level) {
      drawPoly(tm_slot$poly, tm_slot$name, fill = tm_slot$color, lwd = NA, col = NA)
    }
    if (color_type == "both" && tm_slot$level == max(levels) && !(tm_slot$level %in% color_level)) {
      drawPoly(tm_slot$poly, tm_slot$name, fill = tm_slot$color, lwd = NA, col = NA)
    }
  }) %>% invisible()

  # 境界線の描画
  if (!is.null(border_color) && !is.null(border_size)) {
    lapply(treemap@cells, function(tm_slot) {
      if (tm_slot$level %in% border_level) {
        if (length(border_size) == 1) {
          if (inherits(treemap, "sunburstResult")) {
            border_lwd <- border_size
          } else {
            border_lwd <- border_size / tm_slot$level
          }
        } else {
          border_lwd <- border_size[tm_slot$level]
        }

        if (length(border_color) > 1) {
          border_col <- border_color[tm_slot$level]
        } else {
          border_col <- border_color
        }

        drawPoly(tm_slot$poly, tm_slot$name,
          fill = NA, lwd = border_lwd, col = border_col)
      }
    }) %>% invisible()
  }

  # ラベルの描画
  if (!is.null(label_level) && !is.null(label_size) && !is.null(label_color)) {
    if (inherits(treemap, "sunburstResult")) {
      stop("Sunburst treemaps are not fully implemented in this example.")
    } else {
      # 構成比ラベルの準備
      ratio_labels <- list()
      cat("Debug: label_ratio_levels =", label_ratio_levels, "\n")
      cat("Debug: treemap@label_ratios exists =", !is.null(treemap@label_ratios), "\n")
      cat("Debug: label_ratios names =", names(treemap@label_ratios), "\n")
      if (length(label_ratio_levels) > 0 && !is.null(treemap@label_ratios) && 
          all(c("secondary_cluster_ratio", "primary_cluster_ratio") %in% names(treemap@label_ratios))) {
        for (cell_name in names(treemap@cells)) {
          level <- treemap@cells[[cell_name]]$level
          cell_name_raw <- treemap@cells[[cell_name]]$name
          cell_name_clean <- sub("LEVEL[12]_", "", cell_name)
          cat(sprintf("Processing cell: %s (Raw: %s, Clean: %s), Level: %d\n", 
                      cell_name, cell_name_raw, cell_name_clean, level))
          
          if (level == 1) {
            matched_indices <- which(treemap@data$secondary_cluster_name == cell_name_clean)
            if (length(matched_indices) > 0) {
              matched_name <- treemap@data$secondary_cluster_name[matched_indices[1]]
              ratio <- treemap@label_ratios$secondary_cluster_ratio[matched_indices[1]]
              cat("Level 1 - Name:", cell_name_clean, 
                  "Matched:", matched_name, 
                  "Ratio:", ratio, "\n")
              if (!is.na(ratio)) {
                ratio_labels[[cell_name]] <- sprintf("%.1f%%", ratio)
              } else {
                cat("Debug: No match or NA ratio for", cell_name_clean, "at Level 1\n")
              }
            } else {
              cat("Debug: No match found for", cell_name_clean, "at Level 1\n")
            }
          } else if (level == 2) {
            matched_indices <- which(treemap@data$primary_cluster_name == cell_name_clean)
            if (length(matched_indices) > 0) {
              matched_name <- treemap@data$primary_cluster_name[matched_indices[1]]
              ratio <- treemap@label_ratios$primary_cluster_ratio[matched_indices[1]]
              cat("Level 2 - Name:", cell_name_clean, 
                  "Matched:", matched_name, 
                  "Ratio:", ratio, "\n")
              if (!is.na(ratio)) {
                ratio_labels[[cell_name]] <- sprintf("%.1f%%", ratio)
              } else {
                cat("Debug: No match or NA ratio for", cell_name_clean, "at Level 2\n")
              }
            } else {
              cat("Debug: No match found for", cell_name_clean, "at Level 2\n")
            }
          }
        }
        print("Ratio labels prepared:")
        print(ratio_labels)
      } else {
        cat("Debug: Condition failed - label_ratio_levels:", length(label_ratio_levels), 
            "treemap@label_ratios:", !is.null(treemap@label_ratios), 
            "names check:", all(c("secondary_cluster_ratio", "primary_cluster_ratio") %in% names(treemap@label_ratios)), "\n")
      }

      draw_label_voronoi <- function(cells, label_level, label_size, label_color, label_autoscale,
                                     ratio_labels, label_fontfamily, label_fontweight,
                                     label_line_spacing, label_ratio_levels, label_ratio_color,
                                     label_ratio_size, label_ratio_fontfamily, label_ratio_fontweight) {
        lapply(cells, function(tm_slot) {
          if (tm_slot$level %in% label_level) {
            lab_size <- if (length(label_size) == 1) label_size else label_size[which(label_level == tm_slot$level)]
            lab_color <- if (length(label_color) == 1) label_color else label_color[which(label_level == tm_slot$level)]

            if (label_autoscale) {
              cell_area <- tm_slot$area
              mean_area <- mean(sapply(cells, function(c) c$area))
              lab_size <- lab_size * sqrt(cell_area / mean_area)
            }

            # クラスタ名ラベルの描画（y_offset はビューポートで処理済み）
            cat("Drawing label for:", tm_slot$name, "at x:", tm_slot$site[1], "y:", tm_slot$site[2] + (lab_size * 10), "\n")
            grid::grid.text(
              label = tm_slot$name,
              x = unit(tm_slot$site[1], "native"),
              y = unit(tm_slot$site[2] + (lab_size * 10), "native"),  # y_offset をここでは適用しない
              gp = grid::gpar(
                cex = lab_size,
                col = lab_color,
                fontfamily = label_fontfamily,
                fontface = label_fontweight
              ),
              check.overlap = TRUE
            )

            # 構成比ラベルの描画（y_offset はビューポートで処理済み）
            cell_id <- names(cells)[match(list(tm_slot), cells)]
            if (is.null(cell_id)) {
              cat("Debug: Could not find cell_id for tm_slot with name:", tm_slot$name, "\n")
              return(invisible(NULL))
            }
            cat("Debug: cell_id =", cell_id, "for tm_slot$name =", tm_slot$name, "\n")

            if (tm_slot$level %in% label_ratio_levels && !is.null(ratio_labels[[cell_id]]) && !is.na(ratio_labels[[cell_id]])) {
              current_spacing <- label_line_spacing[tm_slot$level]
              current_color <- label_ratio_color[tm_slot$level]
              base_offset <- 40
              scaling_factor <- 20
              offset <- base_offset + (lab_size * scaling_factor * current_spacing)
              offset <- min(offset, 80)

              cat("Drawing ratio label for:", tm_slot$name, "Value:", ratio_labels[[cell_id]], 
                  "at x:", tm_slot$site[1], "y:", tm_slot$site[2] - offset, "\n")
              grid::grid.text(
                label = ratio_labels[[cell_id]],
                x = unit(tm_slot$site[1], "native"),
                y = unit(tm_slot$site[2] - offset, "native"),  # y_offset をここでは適用しない
                gp = grid::gpar(
                  cex = lab_size,
                  col = current_color,
                  fontfamily = label_ratio_fontfamily,
                  fontface = label_ratio_fontweight
                ),
                check.overlap = TRUE
              )
            } else {
              cat("Debug: No ratio label drawn for", tm_slot$name, 
                  "label_ratio_levels includes level:", tm_slot$level %in% label_ratio_levels, 
                  "ratio_labels exists:", !is.null(ratio_labels[[cell_id]]), 
                  "NA check:", !is.na(ratio_labels[[cell_id]]), "\n")
            }
          }
        }) %>% invisible()
      }

      draw_label_voronoi(
        treemap@cells, label_level, label_size, label_color, label_autoscale,
        ratio_labels = ratio_labels,
        label_fontfamily = label_fontfamily,
        label_fontweight = label_fontweight,
        label_line_spacing = label_line_spacing_per_level,
        label_ratio_levels = label_ratio_levels,
        label_ratio_color = label_ratio_color_per_level,
        label_ratio_size = label_ratio_size,
        label_ratio_fontfamily = label_ratio_fontfamily,
        label_ratio_fontweight = label_ratio_fontweight
      )
    }
  }

  # タイトルの描画（オプション）
  if (!is.null(title)) {
    grid::popViewport()
    grid::pushViewport(
      grid::viewport(
        x = key_offsets[1],
        y = 1 - title_offset/2 + y_offset,  # y_offset を適用
        width = key_offsets[3],
        height = title_offset
      )
    )
    grid::grid.text(
      title,
      y = 0.5,
      gp = grid::gpar(
        cex = title_size,
        col = title_color,
        fontfamily = label_fontfamily,
        fontface = title_fontweight
      )
    )
  }

  # 凡例の描画（オプション）
  if (legend) {
    grid::popViewport()
    grid::pushViewport(
      grid::viewport(
        x = key_offsets[5],
        y = key_offsets[6] - title_offset/2 + y_offset,  # y_offset を適用
        width = key_offsets[7],
        height = key_offsets[8] - title_offset
      )
    )
    pal <- treemap@call$palette
    colorkey <- list(
      space = legend_position,
      col = pal,
      at = 0:length(pal),
      labels = c(names(pal), ""),
      width = 0.8,
      axis.line = list(alpha = 1, col = border_color, lwd = 1, lty = 1),
      axis.text = list(alpha = 1, cex = 0.8, col = title_color, font = 1, lineheight = 1)
    )
    grid::grid.draw(lattice::draw.colorkey(key = colorkey))
  }

  # ビューを元に戻す
  grid::popViewport(ifelse(!is.null(title) || legend, 3, 2))
}
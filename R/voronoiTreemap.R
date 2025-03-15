#' voronoiTreemap
#'
#' Create nested additively weighted Voronoi treemaps.
#'
#' @param data (data.frame) A data.frame with one column for each hierarchical level
#' @param levels (character) Character vector indicating the column names to be used.
#' @param fun (function) Function to be used to aggregate cell sizes of parental cells. Default is `sum`.
#' @param sort (logical) Should the columns of the data.frame be sorted before treemap generation? Default is `TRUE`.
#' @param filter (numeric) Filter threshold for small cells. Default is `0`.
#' @param cell_size (character) The name of the column used to control cell size.
#' @param custom_color (character) An optional column to control cell color. Default is `NULL`.
#' @param shape (list or character) Set the initial shape of the treemap. Default is `"rectangle"`.
#' @param maxIteration (numeric) Maximum number of iterations. Default is `100`.
#' @param error_tol (numeric) The allowed maximum error tolerance of a cell. Default is `0.01`.
#' @param convergence (character) Convergence speed: "slow", "intermediate", or "fast". Default is `"intermediate"`.
#' @param seed (integer) Seed for reproducible cell arrangement. Default is `NULL`.
#' @param positioning (character) Algorithm for positioning of starting coordinates. Default is `"regular"`.
#' @param verbose (logical) If `TRUE`, print progress messages. Default is `FALSE`.
#' @param debug (logical) If `TRUE`, draw iterations for visual inspection. Default is `FALSE`.
#'
#' @return `voronoiTreemap` returns an object of the formal class `voronoiResult`.
#'   It contains the following slots:
#'     \item{cells}{`list` of polygons for drawing a treemap, with each cell potentially containing a `ratio` slot for primary or secondary cluster ratios}
#'     \item{data}{`data.frame`, the original data that was supplied to calling `voronoiTreemap`}
#'     \item{call}{`list` of arguments used to call `voronoiTreemap`}
#'
#' @export voronoiTreemap
voronoiTreemap <- function(
  data,
  levels,
  fun = sum,
  sort = TRUE,
  filter = 0,
  cell_size = NULL,
  custom_color = NULL,
  shape = "rectangle",
  maxIteration = 100,
  error_tol = 0.01,
  convergence = "intermediate",
  seed = NULL,
  positioning = "regular",
  verbose = FALSE,
  debug = FALSE
) {

  # validate input data and parameters
  data <- validate_input(
    data, levels, fun,
    sort, filter, cell_size,
    custom_color, verbose)

  # levels の長さが2であることを確認（構成比表示のため）
  if (length(levels) != 2) {
    stop("This function expects exactly 2 levels (secondary_cluster_name and primary_cluster_name).")
  }

  # in debug mode, open a viewport to draw iterations
  if (debug) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      width = 0.9,
      height = 0.9,
      xscale = c(0, 2000),
      yscale = c(0, 2000)
    ))
  }

  # CORE FUNCTION (RECURSIVE)
  voronoi_core <- function(level, df, parent = NULL, output = list()) {
    counter = 1
    repeat {
      # 1. Define the boundary polygon
      if (level == 1) {
        if (is.list(shape)) {
          ParentPoly <- poly_transform_shape(shape)
        } else {
          if (shape == "rectangle") {
            ParentPoly <- list(x = c(0, 0, 2000, 2000, 0), y = c(0, 2000, 2000, 0, 0))
          } else if (shape == "circle") {
            ParentPoly <- list(x = sin(seq(0, 2, 2/50)*pi) * 1000 + 1000, y = cos(seq(0, 2, 2/50)*pi) * 1000 + 1000)
          } else if (shape == "hexagon") {
            ParentPoly <- list(x = sin(seq(0, 2, 2/6)*pi) * 1000 + 1000, y = cos(seq(0, 2, 2/6)*pi) * 1000 + 1000)
          } else if (shape == "rounded_rect") {
            if (!exists("rounded_rect")) {
              stop("rounded_rect is not defined. Ensure it is loaded with data(rounded_rect).")
            }
            ParentPoly <- list(x = rounded_rect[[1]], y = rounded_rect[[2]])
          } else {
            stop("shape is not a coordinate list, nor one of 'rectangle', 'rounded_rect', 'circle', or 'hexagon'.")
          }
        }
        sfpoly <- to_sfpoly(ParentPoly)
      } else {
        stopifnot(!is.null(parent))
        sfpoly <- parent
        ParentPoly <- list(x = parent[[1]][, 1], y = parent[[1]][, 2])
      }

      # 2. Generate starting coordinates
      ncells <- tibble::deframe(dplyr::count(df, get(levels[level])))
      positioning <- ifelse(length(positioning) == 1, positioning, positioning[level])
      if (length(ncells) != 1) {
        sampledPoints <- samplePoints(ParentPoly = ParentPoly, n = length(ncells), seed = seed, positioning = positioning)
      }

      # 3. Generate weights
      if (is.null(cell_size)) {
        weights <- ncells / sum(ncells)
      } else {
        stopifnot(is.numeric(df[[cell_size]]))
        weights <- df %>%
          dplyr::group_by(get(levels[level])) %>%
          dplyr::summarise(fun(get(cell_size)))
        weights <- convertInput(weights[[2]])
      }
      if (length(ncells) != 1 && positioning %in% c("regular_by_area", "clustered_by_area")) {
        sampledPoints <- sampledPoints[order(order(weights)), ]
      }

      # 4. Generate custom color values
      if (!is.null(custom_color)) {
        color_value <- df %>%
          dplyr::group_by(get(levels[level])) %>%
          dplyr::summarise(fun(get(custom_color)))
        color_value <- color_value[[2]]
        color_value <- setNames(color_value, names(ncells))
      }

      # 5. Generate treemap
      if (length(ncells) == 1) {
        treemap <- list(list(
          name = names(ncells),
          poly = sfpoly,
          site = poly_centroid(ParentPoly[[1]], ParentPoly[[2]]),
          weight = weights,
          area = sf::st_area(sfpoly),
          target = weights,
          count = 0
        ))
        names(treemap) <- names(ncells)[1]
      } else {
        treemap <- allocate(
          names = names(ncells),
          s = list(x = sampledPoints[, 1], y = sampledPoints[, 2]),
          w = weights,
          target = weights,
          maxIteration = maxIteration,
          error_tol = error_tol,
          convergence = convergence,
          outer = sfpoly,
          debug = debug
        )

        if (is.null(treemap) & counter < 10) {
          if (!is.null(seed)) { seed = seed + 1 }
          counter = counter + 1
          message("Iteration failed, randomising positions...")
          next
        } else if (is.null(treemap) & counter >= 10) {
          stop("Iteration failed after 10 randomisation trials, try to rerun treemap with new seed")
        }

        if (debug || verbose) {
          tessErr <- sapply(treemap, function(tm) tm$area)
          tessErr <- abs(tessErr / sum(tessErr) - weights)
          message("Level ", level, " tesselation: ",
                  round(mean(tessErr) * 100, 2), " % mean error, ",
                  round(max(tessErr) * 100, 2), " % max error, ",
                  treemap[[1]]$count, " iterations.")
        }
      }

      # Add level and custom color info
      for (i in names(ncells)) {
        treemap[[i]]$level <- level
        treemap[[i]]$custom_color <- if (!is.null(custom_color)) color_value[[i]] else NA
      }

      # Recursive call
      if (level != length(levels)) {
        res <- lapply(1:length(ncells), function(i) {
          voronoi_core(
            level = level + 1,
            df = subset(df, get(levels[level]) %in% names(ncells)[i]),
            parent = treemap[[i]]$poly,
            output = {
              output[[paste0("LEVEL", level, "_", names(ncells)[i])]] <- treemap[[i]]
              output
            }
          )
        }) %>% unlist(recursive = FALSE)
        return(res)
      } else {
        names(treemap) <- paste0("LEVEL", level, "_", names(ncells))
        return(c(output, treemap))
      }
    }
  }

  # MAIN FUNCTION CALL
  tm <- voronoi_core(level = 1, df = data)
  tm <- tm[!duplicated(tm)]
  tm <- tm[names(tm) %>% order]
  if (debug || verbose) {
    message("Treemap successfully created.")
  }

  # Set S4 class BEFORE adding ratios
  tm <- voronoiResult(
    cells = tm,
    data = data,
    call = list(
      levels = levels,
      fun = fun,
      sort = sort,
      filter = filter,
      cell_size = cell_size,
      custom_color = custom_color,
      shape = shape,
      maxIteration = maxIteration,
      error_tol = error_tol,
      seed = seed,
      positioning = positioning
    )
  )

  # 構成比の追加（S4クラス変換後）
  if (all(c("primary_cluster_ratio", "secondary_cluster_ratio") %in% names(data))) {
    level1 <- levels[1]  # "secondary_cluster_name"
    level2 <- levels[2]  # "primary_cluster_name"
    
    # primary_cluster_nameごとにprimary_cluster_ratioを一意にマッピング（全体に対する割合）
    primary_ratios <- data %>%
      group_by(!!sym(level2)) %>%
      summarise(ratio = first(primary_cluster_ratio),  # 最初の値を採用
                total_count = sum(!!sym(cell_size))) %>%
      distinct()
    ratio_map_primary <- setNames(primary_ratios$ratio, primary_ratios[[level2]])
    
    # secondary_cluster_name内でprimary_cluster_nameの割合を再計算（secondary_cluster_nameの合計に対する割合）
    secondary_ratios <- data %>%
      group_by(!!sym(level1), !!sym(level2)) %>%
      summarise(sub_count = sum(!!sym(cell_size))) %>%
      group_by(!!sym(level1)) %>%
      mutate(total_count = sum(sub_count),
             ratio = (sub_count / total_count) * 100) %>%
      distinct() %>%
      select(!!sym(level1), ratio)
    ratio_map_secondary <- setNames(secondary_ratios$ratio, secondary_ratios[[level1]])
    
    # 各セルに構成比を追加
    tm@cells <- lapply(tm@cells, function(tm_slot) {
      cell_name <- tm_slot$name
      if (tm_slot$level == 1) {
        tm_slot$ratio <- ratio_map_secondary[cell_name]
      } else if (tm_slot$level == 2) {
        tm_slot$ratio <- ratio_map_primary[cell_name]
      }
      if (is.na(tm_slot$ratio)) {
        warning(paste("No ratio found for", cell_name, "- setting to NA"))
      }
      tm_slot
    })
  }

  return(tm)
}

# #' drawTreemap
# #'
# #' Draw a nested additively weighted Voronoi treemap.
# #'
# #' @param x (voronoiResult) An object of class `voronoiResult` as returned by `voronoiTreemap`.
# #' @param label_level (numeric) Vector of levels to be labelled. Default is `1`.
# #' @param label_size (numeric) Size of labels. Default is `3`.
# #' @param label_color (character) Color of labels. Default is `"black"`.
# #' @param label_autoscale (logical) Should label size be scaled by cell area? Default is `TRUE`.
# #' @param label_ratio_size (numeric) Size of ratio labels (relative to `label_size`). Default is `NULL` (0.75 times `label_size`).
# #' @param label_ratio_color (character) Color of ratio labels. Default is `NULL` (same as `label_color` but with transparency).
# #' @param legend (logical) Should a legend be drawn? Default is `TRUE`.
# #' @param title (character) Title of the treemap. Default is `""`.
# #' @param color_type (character) Type of coloring: "categorical" or "continuous". Default is `"categorical"`.
# #' @param color_level (numeric) Level to color by. Default is `1`.
# #' @param border_size (numeric) Size of cell borders. Default is `1`.
# #' @param border_color (character) Color of cell borders. Default is `"black"`.
# #' @param palette (character) Color palette for categorical coloring. Default is `NULL` (uses `RColorBrewer`).
# #'
# #' @export drawTreemap
# drawTreemap <- function(
#   x,
#   label_level = 1,
#   label_size = 3,
#   label_color = "black",
#   label_autoscale = TRUE,
#   label_ratio_size = NULL,
#   label_ratio_color = NULL,
#   legend = TRUE,
#   title = "",
#   color_type = "categorical",
#   color_level = 1,
#   border_size = 1,
#   border_color = "black",
#   palette = NULL
# ) {
#   if (!inherits(x, "voronoiResult")) {
#     stop("x must be an object of class 'voronoiResult'")
#   }

#   # Set default label_ratio_size and label_ratio_color if not provided
#   if (is.null(label_ratio_size)) {
#     label_ratio_size <- label_size * 0.75  # クラスタ名の0.75倍
#   }
#   if (is.null(label_ratio_color)) {
#     label_ratio_color <- adjustcolor(label_color, alpha.f = 0.7)  # 透明度を追加
#   }

#   # Prepare colors
#   if (color_type == "categorical") {
#     if (is.null(palette)) {
#       palette <- RColorBrewer::brewer.pal(12, "Set3")[1:length(unique(sapply(x@cells, function(cell) cell$custom_color)))]
#     }
#     colors <- sapply(x@cells, function(cell) {
#       if (cell$level == color_level) {
#         if (!is.na(cell$custom_color)) {
#           palette[which(unique(sapply(x@cells, function(c) c$custom_color)) %in% cell$custom_color)[1]]
#         } else {
#           "grey"
#         }
#       } else {
#         adjustcolor("grey", alpha.f = 0.2)
#       }
#     })
#   } else if (color_type == "continuous") {
#     colors <- sapply(x@cells, function(cell) {
#       if (cell$level == color_level) {
#         if (!is.na(cell$custom_color)) {
#           scales::rescale(cell$custom_color, to = c(0, 1))
#         } else {
#           0.5
#         }
#       } else {
#         adjustcolor("grey", alpha.f = 0.2)
#       }
#     })
#     colors <- colorRampPalette(c("blue", "red"))(100)[round(colors * 100)]
#   }

#   # Draw polygons
#   grid::grid.newpage()
#   for (i in seq_along(x@cells)) {
#     cell <- x@cells[[i]]
#     poly <- cell$poly[[1]]
#     grid::grid.polygon(
#       x = poly[, 1],
#       y = poly[, 2],
#       default.units = "native",
#       gp = grid::gpar(
#         fill = colors[i],
#         col = border_color,
#         lwd = border_size
#       )
#     )
#   }

#   # Draw labels
#   draw_label_voronoi(
#     x@cells,
#     label_level = label_level,
#     label_size = label_size,
#     label_color = label_color,
#     label_autoscale = label_autoscale,
#     label_ratio_size = label_ratio_size,
#     label_ratio_color = label_ratio_color
#   )

#   # Draw title
#   if (nchar(title) > 0) {
#     grid::grid.text(
#       title,
#       x = 0.5,
#       y = 0.95,
#       gp = grid::gpar(fontsize = 18)
#     )
#   }

#   # Draw legend (if requested)
#   if (legend && color_type == "categorical") {
#     unique_colors <- unique(colors[!is.na(colors) & colors != "grey"])
#     unique_names <- unique(sapply(x@cells[colors %in% unique_colors], function(cell) cell$name))
#     n <- length(unique_colors)
#     x_pos <- seq(0.05, 0.95, length.out = n + 1)[-c(1, n + 1)]
#     for (i in seq_along(unique_colors)) {
#       grid::grid.rect(
#         x = x_pos[i],
#         y = 0.05,
#         width = 0.05,
#         height = 0.05,
#         gp = grid::gpar(fill = unique_colors[i], col = border_color)
#       )
#       grid::grid.text(
#         unique_names[i],
#         x = x_pos[i] + 0.03,
#         y = 0.025,
#         just = "left",
#         gp = grid::gpar(fontsize = 8)
#       )
#     }
#   }
# }

# #' draw_label_voronoi
# #'
# #' Draw labels for Voronoi treemap cells.
# #'
# #' @param cells (list) List of cell objects from a `voronoiResult`.
# #' @param label_level (numeric) Vector of levels to be labelled.
# #' @param label_size (numeric) Size of labels.
# #' @param label_color (character) Color of labels.
# #' @param label_autoscale (logical) Should label size be scaled by cell area?
# #' @param label_ratio_size (numeric) Size of ratio labels (relative to `label_size`).
# #' @param label_ratio_color (character) Color of ratio labels.
# #'
# #' @importFrom grid grid.text
# draw_label_voronoi <- function(
#   cells,
#   label_level,
#   label_size,
#   label_color,
#   label_autoscale,
#   label_ratio_size = NULL,
#   label_ratio_color = NULL
# ) {
#   if (is.null(label_ratio_size)) {
#     label_ratio_size <- label_size * 0.75  # クラスタ名の0.75倍
#   }
#   if (is.null(label_ratio_color)) {
#     label_ratio_color <- label_color
#   }

#   for (tm_slot in rev(cells)) {
#     if (tm_slot$level %in% label_level) {
#       if (label_autoscale) {
#         label_cex <- sqrt(tm_slot$area) / (100 * nchar(tm_slot$name)) %>% round(1)
#       } else {
#         label_cex <- 0.5
#       }

#       if (length(label_size) == 1) {
#         label_cex <- label_cex * label_size
#       } else {
#         label_cex <- label_cex * label_size[which(label_level %in% tm_slot$level)]
#       }
#       if (length(label_color) == 1) {
#         label_col <- label_color
#       } else {
#         label_col <- label_color[which(label_level %in% tm_slot$level)]
#       }

#       if (length(label_ratio_size) == 1) {
#         ratio_cex <- label_cex * label_ratio_size
#       } else {
#         ratio_cex <- label_cex * label_ratio_size[which(label_level %in% tm_slot$level)]
#       }
#       if (length(label_ratio_color) == 1) {
#         ratio_col <- label_ratio_color
#       } else {
#         ratio_col <- label_ratio_color[which(label_level %in% tm_slot$level)]
#       }

#       # クラスタ名を描画（1行目）
#       grid::grid.text(
#         tm_slot$name,
#         tm_slot$site[1],
#         tm_slot$site[2],
#         default = "native",
#         gp = gpar(cex = label_cex, col = label_col)
#       )

#       # 構成比ラベルを描画（2行目）
#       if (!is.na(tm_slot$ratio)) {
#         ratio_text <- sprintf("%.1f%%", tm_slot$ratio)
#         y_offset <- -label_cex * 0.5  # フォントサイズに応じたオフセット
#         grid::grid.text(
#           ratio_text,
#           tm_slot$site[1],
#           tm_slot$site[2] + y_offset,
#           default = "native",
#           gp = gpar(cex = ratio_cex, col = ratio_col)
#         )
#       }
#     }
#   }
# }

# #' @importFrom Rcpp evalCpp
# #' @importFrom grid grid.newpage
# #' @importFrom grid pushViewport
# #' @importFrom grid viewport
# #' @importFrom dplyr %>%
# #' @importFrom dplyr mutate_if
# #' @importFrom dplyr group_by
# #' @importFrom dplyr summarise
# #' @importFrom dplyr count
# #' @importFrom tibble deframe
# #' @importFrom scales rescale
# #' @importFrom sf st_polygon
# #' @importFrom sf st_area
# #' @importFrom sp Polygon
# #' @importFrom sp spsample
# #'
# #' @useDynLib WeightedTreemaps, .registration = TRUE

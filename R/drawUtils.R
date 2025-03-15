#' @importFrom grid gpar
#' @importFrom grid grid.polygon
#' @importFrom grid grid.text
#' @importFrom grid grid.lines
#' @importFrom sf st_area
#' @importFrom scales rescale
#' @importFrom methods as
#' @importFrom methods new
#' @importFrom stats median
#' @importFrom colorspace lighten
#' @importFrom colorspace rainbow_hcl

# function to coerce and rescale different types of input to
# numeric range between 1 and 100 (for color coding)
convertInput <- function(x, from = NULL, to = c(1, 100)) {
  if (is.character(x)) {
    if (all(!is.na(suppressWarnings(as.numeric(x))))) {
      x = as.numeric(x)
    } else {
      x = as.factor(x) %>%
        as.numeric
    }
  }
  if (is.numeric(x)) {
    res <- scales::rescale(x,
      from = {if (!is.null(from)) from else range(x)},
      to = to) %>% round
    res <- replace(res, res > to[2], to[2])
    res <- replace(res, res < to[1], to[1])
    res
  } else {
    stop("Input data is not of type numeric, factor, or character. Color-coding impossible.")
  }
}

drawPoly <- function(sfpoly, name, fill, lwd, col) {
  if (length(sfpoly)) {
    pts <- to_coords(sfpoly)
    grid::grid.polygon(
      pts$x,
      pts$y,
      default = "native",
      gp = gpar(col = col, lwd = lwd, fill = fill),
      name = name)
  }
}

polyRangeX <- function(sfpoly) {
  if (length(sfpoly)) {
    pts <- to_coords(sfpoly)
    range(pts$x)
  } else {
    NA
  }
}

polyRangeY <- function(sfpoly) {
  if (length(sfpoly)) {
    pts <- to_coords(sfpoly)
    range(pts$y)
  } else {
    NA
  }
}

drawRegions <- function(
  result,
  debug = FALSE,
  label = TRUE,
  label.col = grey(0.5),
  lwd = 2, col = grey(0.8),
  fill = NA)
{
  names <- result$names
  k <- result$k
  sites <- result$s

  # draw polygon, pass graphical parameters to drawPoly function
  mapply(drawPoly, k, names, fill = fill,
    SIMPLIFY = FALSE,
    MoreArgs = list(lwd = lwd, col = col)
  )

  if (label) {
    # function to determine label sizes for each individual cell
    # based on cell dimension and label character length
    cex = sqrt(unlist(result$a)) * 0.01 / nchar(names) %>%
      round(1)
    grid::grid.text(names,
      sites$x,
      sites$y,
      default = "native",
      gp = gpar(cex = cex, col = label.col)
    )
  }
}

# calculate sector polygon from boundary input
draw_sector <- function(
  level,
  lower_bound,
  upper_bound,
  diameter_inner,
  diameter_sector,
  name,
  custom_color) {

  # compute_sector from lower and upper bounds and diameter arguments
  segment <- c(lower_bound, upper_bound) * 2 * pi
  a <- diameter_inner + (diameter_sector * (level - 1))
  z <- seq(segment[1], segment[2], by = pi/400)
  xx <- c(a * cos(z), rev((a + diameter_sector) * cos(z)))
  yy <- c(a * sin(z), rev((a + diameter_sector) * sin(z)))
  # rescale for canvas dimensions [0, 2000] and convert into sfpoly polygon
  poly = to_sfpoly(list(x = (xx+1)*1000, y = (yy+1)*1000))

  # return list of polygon properties
  list(
    name = name,
    poly = poly,
    area = sf::st_area(poly),
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    level = level,
    custom_color = custom_color
  )
}

# draw_label_voronoi 関数（Voronoiツリーマップ用のラベル描画）
draw_label_voronoi <- function(cells, levels, size, color, autoscale, label_ratio_format = "%.1f%%", fontfamily = "sans") {
  grid::grid.newpage()
  grid::pushViewport(viewport(xscale = c(0, 2000), yscale = c(0, 2000)))

  lapply(cells, function(tm_slot) {
    cat("Processing cell:", tm_slot$name, "Level:", tm_slot$level, "\n")
    if (tm_slot$level %in% levels) {
      cat("Drawing label for:", tm_slot$name, "Level:", tm_slot$level, "\n")
      level_idx <- which(levels == tm_slot$level)
      label_size <- if (length(size) > 1) size[level_idx] else size
      label_color <- if (length(color) > 1) color[level_idx] else color

      # ポリゴンデータの存在確認
      if (is.null(tm_slot$poly) || !inherits(tm_slot$poly, "sfg")) {
        warning("Invalid polygon data for cell: ", tm_slot$name)
        return(NULL)
      }

      # サイトデータの存在確認
      if (is.null(tm_slot$site) || length(tm_slot$site) < 2) {
        warning("Invalid site data for cell: ", tm_slot$name)
        return(NULL)
      }

      # ラベル名
      label <- tm_slot$name
      if (!is.null(label_ratio_format) && !is.null(tm_slot$ratio) && !is.na(tm_slot$ratio)) {
        ratio_text <- sprintf(label_ratio_format, tm_slot$ratio * 100)
        label <- paste(label, ratio_text, sep = "\n")
      }

      # 重心計算
      if (inherits(tm_slot$poly, "sfg")) {
        coords <- sf::st_coordinates(tm_slot$poly)[, c("X", "Y")]
        centroid <- c(mean(coords[, "X"]), mean(coords[, "Y"]))
      } else {
        centroid <- tm_slot$site
      }
      x <- centroid[1] / 2000
      y <- centroid[2] / 2000

      # オートスケーリング
      if (autoscale) {
        # to_coords を使用して座標を抽出
        poly_coords <- to_coords(tm_slot$poly)
        if (length(poly_coords$x) > 1) {
          width <- (max(poly_coords$x) - min(poly_coords$x)) / 2000
          height <- (max(poly_coords$y) - min(poly_coords$y)) / 2000
          char_width <- 0.05
          char_height <- 0.08
          label_aspect_ratio <- nchar(label) * char_width / (2 * char_height)
          cell_aspect_ratio <- width / height

          label_size <- min(
            width / (nchar(label) * char_width),
            height / (2 * char_height),
            max_font_size / 300
          ) * label_size

          cat("Cell:", tm_slot$name, "Width:", width, "Height:", height, "Label Size:", label_size, "\n")
        } else {
          warning("Insufficient polygon data for autoscale in cell: ", tm_slot$name)
          label_size <- min(max(label_size, 1), 10)
        }
      } else {
        label_size <- min(max(label_size, 1), 10)
      }

      cat("Label drawn for:", tm_slot$name, "at (", x, ",", y, ") with size", label_size, "\n")

      # ラベルの描画
      grid::grid.text(
        label,
        x = unit(x, "npc"),
        y = unit(y, "npc"),
        just = "center",
        hjust = 0.5,
        vjust = 0.5,
        gp = grid::gpar(
          fontsize = label_size * 300,
          col = label_color,
          fontfamily = fontfamily,
          lineheight = 0.8
        )
      )
    } else {
      cat("Skipping cell:", tm_slot$name, "Level:", tm_slot$level, "not in", levels, "\n")
    }
  }) %>% invisible

  grid::popViewport()
}

# draw_label_sunburst 関数（Sunburstツリーマップ用のラベル描画）
draw_label_sunburst <- function(
  cells,
  label_level,
  label_size,
  label_color,
  diameter,
  label_ratio_format = "%.1f%%",
  fontfamily = "zenmaru"
) {
  lapply(cells, function(tm_slot) {
    if (tm_slot$level %in% label_level) {
      # ラベルサイズと色の決定
      if (length(label_size) > 1) {
        label_cex <- label_size[1]
        warning("'label_size' should only have length 1. Using first argument.")
      } else {
        label_cex <- label_size
      }

      if (length(label_color) > 1) {
        label_col <- label_color[1]
        warning("'label_color' should only have length 1. Using first argument.")
      } else {
        label_col <- label_color
      }

      # セクターの角度計算
      segment <- c(tm_slot$lower_bound, tm_slot$upper_bound) * 2 * pi
      z <- seq(segment[1], segment[2], by = pi/400)
      if (diameter * cos(stats::median(z)) >= 0) side = 1 else side = -1
      sinz <- sin(stats::median(z))
      cosz <- cos(stats::median(z))
      d1 <- diameter + 0.02
      d2 <- diameter + 0.05
      d3 <- diameter + 0.10

      # ラベルテキスト（クラスター名と構成比を2行表示）
      label_text <- tm_slot$name
      if (!is.null(label_ratio_format) && !is.null(tm_slot$ratio) && !is.na(tm_slot$ratio)) {
        ratio_text <- sprintf(label_ratio_format, tm_slot$ratio * 100)  # パーセントに変換
        label_text <- paste(substr(label_text, 1, 18), ratio_text, sep = "\n")  # 2行表示、名前を18文字に制限
      }

      # ラベル用の円弧と線の描画
      z <- z[-c(1, length(z))]
      grid::grid.lines(
        (c(d1 * cos(z[1]), d2 * cos(z), d1 * cos(tail(z, 1))) + 1) * 1000,
        (c(d1 * sin(z[1]), d2 * sin(z), d1 * sin(tail(z, 1))) + 1) * 1000,
        default.units = "native",
        gp = gpar(lwd = label_cex, col = label_col)
      )

      grid::grid.lines(
        x = (c(d2 * cosz, d2 * cosz + 0.15 * cosz * abs(sinz), d3 * side) + 1) * 1000,
        y = (c(d2 * sinz, d2 * sinz + 0.15 * sinz * abs(sinz),
               d2 * sinz + 0.15 * sinz * abs(sinz)) + 1) * 1000,
        default.units = "native",
        gp = gpar(lwd = label_cex, col = label_col)
      )

      # ラベルの描画（2行表示を考慮）
      grid::grid.text(
        label = label_text,
        x = ((d3 + 0.02) * side + 1) * 1000,
        y = ((d2 * sinz + 0.15 * sinz * abs(sinz)) + 1) * 1000,
        just = ifelse(side == 1, "left", "right"),
        default.units = "native",
        gp = gpar(
          cex = label_cex,
          col = label_col,
          fontfamily = fontfamily,
          lineheight = 0.8  # 2行表示時の行間調整
        )
      )
    }
  }) %>% invisible
}

# function to add colors to a treemap object
add_color <- function(treemap, color_palette = NULL,
  color_type = "categorical", color_level = 1,
  color_steps = 10, custom_range = NULL) {

  # CASE 1: CATEGORICAL
  if (color_type %in% c("categorical", "both")) {
    # determine number of required colors
    if (length(color_level) == 1) {
      color_list <- unique(treemap@data[[treemap@call$levels[color_level]]])
    } else {
      color_list <- apply(treemap@data[treemap@call$levels[color_level]], 2, unique) %>%
        unlist
    }
  }

  # CASE 2: CELL AREA
  # determine total area per level
  level_areas <- list()
  for (lvl in 1:length(treemap@call$levels)) {
    total_area <- lapply(treemap@cells, function(tm_slot) {
      if (tm_slot$level == lvl) tm_slot$area
    }) %>% unlist %>% sum
    level_areas[[lvl]] <- total_area
  }
  treemap@call$level_areas <- level_areas  # レベルごとの総面積を保存

  # determine number of required colors
  if (color_type == "cell_size") {
    cell_sizes <- lapply(treemap@cells, function(tm_slot) {
      if (tm_slot$level %in% color_level) tm_slot$area / level_areas[[tm_slot$level]]
    }) %>% unlist
    color_list <- cell_sizes %>% pretty(n = color_steps)
  }

  # CASE 3: CUSTOM COLOR
  if (color_type == "custom_color") {
    color_list <- lapply(treemap@cells, function(tm_slot) {
        if (tm_slot$level %in% color_level) tm_slot$custom_color
      }) %>% unlist %>% pretty(n = 10)
  }

  # DEFINE PALETTE
  if (!is.null(custom_range) & !(color_type %in% c("categorical", "both"))) {
    color_list <- custom_range %>% pretty(n = 10)
  }
  if (is.null(color_palette)) {
    pal <- colorspace::rainbow_hcl(length(color_list), start = 60)
  } else {
    pal <- colorRampPalette(color_palette)(length(color_list))
  }
  pal <- setNames(pal, color_list)

  # ADD COLORS TO TREEMAP OBJECT
  treemap@cells <- lapply(treemap@cells, function(tm_slot) {
    if (tm_slot$level %in% color_level) {
      if (color_type %in% c("categorical", "both")) {
        tm_slot$color <- pal[[tm_slot$name]]
      } else if (color_type == "cell_size") {
        area <- tm_slot$area / level_areas[[tm_slot$level]]
        tm_slot$color <- pal[[findInterval(area, as.numeric(names(pal)))]]
      } else if (color_type == "custom_color") {
        if (tm_slot$custom_color < as.numeric(names(pal))[[1]]) {
          tm_slot$color <- pal[[1]]
        } else {
          tm_slot$color <- pal[[findInterval(tm_slot$custom_color, as.numeric(names(pal)))]]
        }
      }
    }
    tm_slot
  })

  # SPECIAL CASE "BOTH": DARKEN OR LIGHTEN LOWEST CELL LEVEL
  if (color_type == "both") {
    cell_area <- lapply(treemap@cells, function(tm_slot) {
      if (tm_slot$level == length(treemap@call$levels)) tm_slot$area
    }) %>% unlist
    treemap@cells <- lapply(treemap@cells, function(tm_slot) {
      if (tm_slot$level == length(treemap@call$levels)) {
        area <- tm_slot$area / level_areas[[tm_slot$level]]
        corr_factor <- scales::rescale(area, from = range(cell_area / level_areas[[tm_slot$level]]), to = c(-0.2, 0.2))
        if (tm_slot$level %in% color_level) {
          tm_slot$color <- colorspace::lighten(tm_slot$color, corr_factor)
        } else {
          tm_slot$color <- grey(0.5 + (2 * corr_factor), alpha = (1.5 * abs(corr_factor)))
        }
      }
      tm_slot
    })
  }

  # return treemap with colors and palette
  treemap@call$palette <- pal
  treemap
}

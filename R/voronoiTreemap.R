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
    level1 <- levels[1]
    level2 <- levels[2]
    ratio_map_primary <- setNames(data$primary_cluster_ratio, data[[level2]])
    ratio_map_secondary <- setNames(data$secondary_cluster_ratio, data[[level1]])
    tm@cells <- lapply(tm@cells, function(tm_slot) {
      cell_name <- tm_slot$name
      if (tm_slot$level == 1) {
        tm_slot$ratio <- ratio_map_primary[cell_name]
      } else if (tm_slot$level == 2) {
        tm_slot$ratio <- ratio_map_secondary[cell_name]
      }
      if (is.na(tm_slot$ratio)) {
        warning(paste("No ratio found for", cell_name, "- setting to NA"))
      }
      tm_slot
    })
  }

  return(tm)
}

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
#' @param shape (list or character) Set the initial shape of the treemap. Default is `"rounded_rect"`.
#' @param maxIteration (numeric) Maximum number of iterations. Default is `400`.
#' @param error_tol (numeric) The allowed maximum error tolerance of a cell. Default is `0.2`.
#' @param positioning (character) Algorithm for positioning of starting coordinates. Default is `"regular"`.
#' @param verbose (logical) If `TRUE`, print progress messages. Default is `FALSE`.
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
  cell_size,
  fun = sum,  # 追加
  sort = TRUE,  # 追加
  filter = 0,  # 追加
  custom_color = NULL,  # 追加
  shape = "rounded_rect",
  positioning = "regular",
  error_tol = 0.2,
  maxIteration = 400,
  verbose = FALSE
) {
  # 入力データの検証
  data <- validate_input(data, levels, fun, sort, filter, cell_size, custom_color, verbose)

  # levels の長さが2であることを確認
  if (length(levels) != 2) {
    stop("This function expects exactly 2 levels (secondary_cluster_name and primary_cluster_name).")
  }

  # データを変換
  input <- convertInput(data, levels, cell_size)

  # ツリーマップの生成
  result <- cropped_voronoi(input, shape, positioning, error_tol, maxIteration, verbose)

  # 結果をクラスに格納
  result <- voronoiResult(
    cells = result$cells,
    levels = levels,
    cell_size = cell_size,
    shape = shape,
    positioning = positioning,
    error = result$error,
    iterations = result$iterations
  )

  # ratio を cells に追加
  if (all(c("primary_cluster_ratio", "secondary_cluster_ratio") %in% names(data))) {
    level1 <- levels[1]  # "secondary_cluster_name"（二次クラスタ名）
    level2 <- levels[2]  # "primary_cluster_name"（一次クラスタ名）
    
    # データフレームから ratio をキー付きで取得
    ratio_map_primary <- setNames(data$primary_cluster_ratio, data[[level2]])
    ratio_map_secondary <- setNames(data$secondary_cluster_ratio, data[[level1]])
    
    # cells に ratio を追加（ベクトル化処理）
    result@cells <- lapply(result@cells, function(tm_slot) {
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

  # CORE FUNCTION (RECURSIVE)
  voronoi_core <- function(level, df, parent = NULL, output = list()) {

    # set counter for number of maximum tries to not get stuck
    # in repeat loop
    counter = 1

    repeat {

      # CREATE VORONOI TREEMAP OBJECT
      #
      # 1. define the boundary polygon
      # either predefined rectangular bounding box for 1st level
      if (level == 1) {

        if (is.list(shape)) {
          ParentPoly <- poly_transform_shape(shape)
        } else {
          if (shape == "rectangle") {
            ParentPoly <- list(
              x = c(0, 0, 2000, 2000, 0),
              y = c(0, 2000, 2000, 0, 0)
            )
          } else if (shape == "circle") {
            ParentPoly <- list(
              x = sin(seq(0, 2, 2/50)*pi) * 1000 + 1000,
              y = cos(seq(0, 2, 2/50)*pi) * 1000 + 1000
            )
          } else if (shape == "hexagon") {
            ParentPoly <- list(
              x = sin(seq(0, 2, 2/6)*pi) * 1000 + 1000,
              y = cos(seq(0, 2, 2/6)*pi) * 1000 + 1000
            )
          } else if (shape == "rounded_rect") {
            if (!exists("rounded_rect")) {
              stop("rounded_rect is not defined. Ensure it is loaded with data(rounded_rect).")
            }
            ParentPoly <- list(
              x = rounded_rect[[1]],
              y = rounded_rect[[2]]
            )
          } else {
            stop("shape is not a coordinate list, nor one of 'rectangle', 'rounded_rect', 'circle', or 'hexagon'.")
          }
        }

        # turn boundary polygon into sf polygon object for treemap generation
        sfpoly <- to_sfpoly(ParentPoly)

      } else {

        # or the parental polygon in case of all lower levels > 1
        stopifnot(!is.null(parent))
        sfpoly <- parent
        ParentPoly <- list(x = parent[[1]][, 1], y = parent[[1]][, 2])

      }

      # 2. generate starting coordinates within the boundary polygon
      # using sp package's spsample function.
      ncells <- tibble::deframe(dplyr::count(df, get(levels[level])))

      # positioning can be defined globally or for each level independently
      positioning <- ifelse(
        length(positioning) == 1,
        positioning,
        positioning[level]
      )

      if (length(ncells) != 1) {
        sampledPoints <- samplePoints(
          ParentPoly = ParentPoly,
          n = length(ncells),
          seed = seed,
          positioning = positioning
        )
      }

      # 3. generate the weights, these are the (aggregated) scaling factors
      # supplied by the user or simply the n members per cell
      if (is.null(cell_size)) {
        # average cell size by number of members, if no function is given
        weights <- ncells / sum(ncells)
      } else {
        # average cell size by user defined function, e.g. sum of expression values
        # the cell size is calculated as aggregated relative fraction of total
        stopifnot(is.numeric(df[[cell_size]]))
        weights <- df %>%
          dplyr::group_by(get(levels[level])) %>%
          dplyr::summarise(fun(get(cell_size)))
        weights <- weights[[2]] / sum(weights[[2]])
      }
      # reorder starting coordinate positions by weights (= target cell areas)
      # if sorting by area is toggled
      if (length(ncells) != 1 &&
          positioning %in% c("regular_by_area", "clustered_by_area")) {
        sampledPoints <- sampledPoints[order(order(weights)), ]
      }

      # 4. generate custom color values for each cell that can be used
      # with different palettes when drawing;
      if (!is.null(custom_color)) {
        color_value <- df %>%
          dplyr::group_by(get(levels[level])) %>%
          dplyr::summarise(fun(get(custom_color)))
        color_value <- color_value[[2]]
        color_value <- setNames(color_value, names(ncells))
      }

      # 5. generate additively weighted voronoi treemap object;
      # the allocate function returns a list of polygons to draw,
      # among others.
      # if the parent has only 1 child, skip map generation
      # and make pseudo treemap object instead
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
          s = list(
            x = sampledPoints[, 1],
            y = sampledPoints[, 2]),
          w = weights,
          target = weights,
          maxIteration = maxIteration,
          error_tol = error_tol,
          convergence = "slow",
          outer = sfpoly,
          debug = debug
        )

        # error handling in case of failed tesselation:
        # try up to ten new random starting positions before finally giving up
        if (is.null(treemap) & counter < 10) {
          if (!is.null(seed)) { seed = seed + 1 }
          counter = counter + 1
          message("Iteration failed, randomising positions...")
          next
        } else if (is.null(treemap) & counter >= 10) {
          stop("Iteration failed after 10 randomisation trials, try to rerun treemap with new seed")
        }

        # print summary of cell tesselation
        if (debug || verbose) {
          tessErr <- sapply(treemap, function(tm) tm$area)
          if (any(tessErr < 0)) {
            warning("Negative area detected. Adjusting weights and retrying may help.")
          }
          tessErr <- abs(tessErr / sum(tessErr) - weights)
          message("Level ", level, " tesselation: ",
                  round(mean(tessErr) * 100, 2), " % mean error, ",
                  round(max(tessErr) * 100, 2), " % max error, ",
                  treemap[[1]]$count, " iterations.")
        }
      }

      # add level and custom color info to treemap
      for (i in names(ncells)) {
        treemap[[i]]$level <- level
        treemap[[i]]$custom_color <- if (!is.null(custom_color)) color_value[[i]] else NA
      }

      # CALL CORE FUNCTION RECURSIVELY
      if (level != length(levels)) {
        # iterate through all possible sub-categories,
        # these are the children of the parental polygon
        # and pass the children's polygon as new parental
        # also add current tesselation results to output list
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
        }) %>%
          unlist(recursive = FALSE)
        return(res)
      } else {
        names(treemap) <- paste0("LEVEL", level, "_", names(ncells))
        return(c(output, treemap))
      }
    }
  }

  # MAIN FUNCTION CALL
  # ------------------
  # iterate through all levels,
  # collect results in list, remove duplicated polygons
  # and order by hierarchical level
  tm <- voronoi_core(level = 1, df = data)
  tm <- tm[!duplicated(tm)]
  tm <- tm[names(tm) %>% order]
  if (debug || verbose) {
    message("Treemap successfully created.")
  }

  # set S4 class and return result
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

  return(tm)
}

#' voronoiTreemap
#'
#' Create nested additively weighted Voronoi treemaps.
#'
#' This is a recursive wrapper function, making use of the original implementation
#' of the voronoi tesselation from Paul Murrell, University of Auckland.
#' The original functions were obtained and slightly modified from
#' \url{https://www.stat.auckland.ac.nz/~paul/Reports/VoronoiTreemap/voronoiTreeMap.html}
#' This function returns a treemap object instead of a plot. In order
#' to actually draw the treemap, use \code{\link{drawTreemap}}.
#'
#'
#' @param data (data.frame) A data.frame with one column for each hierarchical level
#' @param levels (character) Character vector indicating the column names to
#'   be used. The order of names must correspond to the hierarchical levels,
#'   going from broad to fine
#' @param fun (function) Function to be used to aggregate cell sizes of parental cells
#' @param sort (logical) Should the columns of the data.frame be sorted before treemap generation?
#' @param filter (numeric) Filter the supplied data frame to remove very small
#'   cells that may not be visible. The default is to remove cells with a
#'   relative target area below a threshold of zero (no negative values allowed).
#'   Computation time can increase when many small cells are present. For example,
#'   a threshold of 0.01 filters out all observations/cells below 1 \% of the total area.
#' @param cell_size (character) The name of the column used to control cell size.
#'   Can be one of \code{levels} or any other column with numerical data. NA or
#'   values equal or less than zero are not allowed as the cell area needs to be positive.
#'   The values in this column are aggregated by the function specified by \code{fun}.
#'   If \code{cell_size = NULL}, cell area is simply computed by the number of members
#'   for the respective cell (corresponding to rows in the data.frame).
#' @param custom_color (character) An optional column that can be specified to
#'   control cell color. Cell colors are determined when drawing the treemap
#'   using \code{\link{drawTreemap}}, but the default is to use one of
#'   \code{levels} or \code{cell size}. Any other data source that shall be used
#'   instead has to be included in the treemap generation and explicitly
#'   specified here. The default value is \code{NULL}.
#' @param shape (list or character) Set the initial shape of the treemap. Currently
#'   supported are the keywords "rectangle", "rounded_rect", "circle" or "hexagon".
#'   Alternatively the user can supply a named list with coordinates for a custom polygon.
#'   The slots of the list must be labeled 'x' and 'y'. The coordinates are not tested
#'   for validity, use on your own risk.
#' @param maxIteration (numeric) Force algorithm to stop at this number of iterations
#'   for each parent cell. The algorithm usually converges to an acceptable
#'   solution fairly quickly, so it seems reasonable to restrict this number
#'   in order to save computation time. However, more iterations give higher
#'   accuracy.
#' @param error_tol (numeric) The allowed maximum error tolerance of a cell.
#'   The algorithm will stop when all cells have lower error than this value.
#'   It is calculated as the absolute difference of a cell's area to its target
#'   area. The default is 0.01 (or 1 \%) of the total parental area. Note: this
#'   is is different from a relative per-cell error, where 1 \% would be more
#'   strict.
#' @param convergence (character) One of "slow", "intermediate", or "fast".
#'   Intermediate (default) and fast try to adjust cell weights stronger such
#'   that the algorithm converges faster towards the final size of the cell.
#'   However this comes at the price of stability, with a larger number of
#'   polygons possibly being misformed, e.g. by having self-intersections.
#'   Set convergence to "slow" if you experience problems to calculate treemaps
#'   with very unequal cell sizes or very large treemaps.
#' @param seed (integer) The default seed is NULL, which will lead to a new
#'   random sampling of cell coordinates for each tesselation. If you want
#'   a reproducible arrangement of cells, set seed to an arbitrary number.
#' @param positioning (character) Algorithm for positioning of starting
#'   coordinates of child cells in the parental cell using \code{spsample()};
#'   "random" for completely random positions, "regular" for cells aligned
#'   to a grid sorted from bottom to top by name, "clustered" with regular
#'   positions of cells but sorted by name from inside out. Two variants
#'   "regular_by_area" and "clustered_by_area" will work as their counterparts
#'   but will sort by cell target area instead of cell name. \code{positioning}
#'   can be a single character or a vector of \code{length(levels)} to allow
#'   different positioning algorithms for each level.
#' @param verbose (logical) If verbose is TRUE (default is FALSE), messages
#'   with statistics for each iteration of a treemap as well as a success message
#'   are printed to the console.
#' @param debug (logical) If debug is TRUE (default is FALSE), the solution
#'   for each iteration is drawn to the viewport to allow some visual
#'   inspection. The weights, target area, and difference are printed to the
#'   console. It is not recommended to set this option to TRUE unless you know
#'   what you are doing, as it makes treemap generation much slower.
#'
#' @return `voronoiTreemap` returns an object of the formal class `voronoiResult`.
#'   It is essentially a list of objects related to the graphical
#'   representation of the treemap (polygons, labels, cell data) as well as data from the call
#'   of the function. It contains the following slots:
#'     \item{cells}{`list` of polygons for drawing a treemap}
#'     \item{data}{`data.frame`, the original data that was supplied to calling `voronoiTreemap`}
#'     \item{call}{`list` of arguments used to call `voronoiTreemap`}
#'
#' @seealso \code{\link{drawTreemap}} for drawing the treemap.
#'
#' @examples
#' # load package
#' library(WeightedTreemaps)
#'
#' # generate dummy data
#' df <- data.frame(
#'   A = rep(c("abcd", "efgh"), each = 4),
#'   B = letters[1:8],
#'   size = c(37, 52, 58, 27, 49, 44, 34, 45)
#' )
#'
#' # compute treemap
#' tm <- voronoiTreemap(
#'   data = df,
#'   levels = c("B"),
#'   cell_size = "size",
#'   shape = "circle",
#'   positioning = "regular",
#'   seed = 123
#' )
#'
#' # plot treemap with each cell colored by name (default)
#' drawTreemap(tm, label_size = 1, color_type = "categorical")
#'
#' # plot treemap with each cell colored by name, but larger cells
#' # lighter and smaller cells darker
#' drawTreemap(tm, label_size = 1, color_type = "both")
#'
#' # plot treemap with different color palette and style
#' drawTreemap(tm, label_size = 1, label_color = grey(0.3),
#'             border_color = grey(0.3), color_palette = heat.colors(6)
#' )
#'
#' @importFrom Rcpp evalCpp
#' @importFrom grid grid.newpage
#' @importFrom grid pushViewport
#' @importFrom grid viewport
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate_if
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr count
#' @importFrom tibble deframe
#' @importFrom scales rescale
#' @importFrom sf st_polygon
#' @importFrom sf st_area
#' @importFrom sp Polygon
#' @importFrom sp spsample
#'
#' @useDynLib WeightedTreemaps, .registration = TRUE
#'
#' @export voronoiTreemap
#'
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

  # 入力データのバリデーション（既存のまま）
  data <- validate_input(
    data, levels, fun,
    sort, filter, cell_size,
    custom_color, verbose)

  # デバッグモードの設定（既存のまま）
  if (debug) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      width = 0.9,
      height = 0.9,
      xscale = c(0, 2000),
      yscale = c(0, 2000)
    ))
  }

  # コア関数（再帰処理）
  voronoi_core <- function(level, df, parent = NULL, output = list()) {

    # 最大試行回数のカウンタ（既存のまま）
    counter = 1

    repeat {

      # 境界ポリゴンの定義（既存のまま）
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
            ParentPoly <- list(
              x = rounded_rect[[1]],
              y = rounded_rect[[2]]
            )
          } else {
            stop("shape is not a coordinate list, nor one of 'rectangle', 'rounded_rect', circle', or 'hexagon'.")
          }
        }
        sfpoly <- to_sfpoly(ParentPoly)
      } else {
        stopifnot(!is.null(parent))
        sfpoly <- parent
        ParentPoly <- list(x = parent[[1]][, 1], y = parent[[1]][, 2])
      }

      # 開始座標の生成（既存のまま）
      ncells <- tibble::deframe(dplyr::count(df, get(levels[level])))

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

      # ウェイトの生成（既存のまま）
      if (is.null(cell_size)) {
        weights <- ncells / sum(ncells)
      } else {
        stopifnot(is.numeric(df[[cell_size]]))
        weights <- df %>%
          dplyr::group_by(get(levels[level])) %>%
          dplyr::summarise(fun(get(cell_size)))
        weights <- weights[[2]]/sum(weights[[2]])
      }

      if (length(ncells) != 1 &
          positioning %in% c("regular_by_area", "clustered_by_area")) {
        sampledPoints <- sampledPoints[order(order(weights)), ]
      }

      # カスタムカラーの生成（既存のまま）
      if (!is.null(custom_color)) {
        color_value <- df %>%
          dplyr::group_by(get(levels[level])) %>%
          dplyr::summarise(fun(get(custom_color)))
        color_value <- color_value[[2]]
        color_value <- setNames(color_value, names(ncells))
      }

      # ボロノイテッセレーションの生成（既存のまま）
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
          convergence = convergence,
          outer = sfpoly,
          debug = debug
        )

        if (is.null(treemap) & counter < 10) {
          if (!is.null(seed)) {seed = seed + 1}
          counter = counter + 1
          message("Iteration failed, randomising positions...")
          next
        } else if (is.null(treemap) & counter >= 10) {
          stop("Iteration failed after 10 randomisation trials, try to rerun treemap with new seed")
        }

        if (debug || verbose) {
          tessErr <- sapply(treemap, function(tm) tm$area)
          tessErr <- abs(tessErr/sum(tessErr) - weights)
          message("Level ", level, " tesselation: ",
            round(mean(tessErr) * 100, 2), " % mean error, ",
            round(max(tessErr) * 100, 2), " % max error, ",
            treemap[[1]]$count, " iterations."
          )
        }
      }

      # 構成比の追加（ここが新たに追加される部分）
      for (i in names(ncells)) {
        treemap[[i]]$level <- level
        treemap[[i]]$custom_color <- {if (!is.null(custom_color))
          color_value[[i]] else NA}
        
        # 構成比の保存
        cluster_name <- names(ncells)[i]
        if (level == 1) {  # レベル1はsecondary_cluster_nameに対応
          treemap[[i]]$ratio <- df$secondary_cluster_ratio[df$secondary_cluster_name == cluster_name][1]
        } else if (level == 2) {  # レベル2はprimary_cluster_nameに対応
          treemap[[i]]$ratio <- df$primary_cluster_ratio[df$primary_cluster_name == cluster_name][1]
        } else {
          treemap[[i]]$ratio <- NA  # それ以降のレベルはNA（必要に応じて調整）
        }
      }

      # 再帰呼び出し（既存のまま）
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
        }) %>%
        unlist(recursive = FALSE)
        return(res)
      } else {
        names(treemap) <- paste0("LEVEL", level, "_", names(ncells))
        return(c(output, treemap))
      }
    }
  }

  # メイン関数呼び出し（既存のまま）
  tm <- voronoi_core(level = 1, df = data)
  tm <- tm[!duplicated(tm)]
  tm <- tm[names(tm) %>% order]
  if (debug || verbose) {
    message("Treemap successfully created.")
  }

  # S4クラスの設定と結果の返却（既存のまま）
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

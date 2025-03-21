# /Users/hdymacuser/R_date_analysis/WeightedTreemaps/R/test2.r

# 必要なライブラリ
library(sysfonts)
library(showtext)
library(WeightedTreemaps)
library(tidyverse)
library(ggplot2)
library(grid)
library(IRdisplay)

# フォントの追加と確認
font_add_google("Zen Maru Gothic", "zenmaru")
showtext_auto()

# フォントが利用可能か確認
font_families()  # "zenmaru" がリストに含まれているか確認

# ファイルパス
file_path <- "/Users/hdymacuser/R_date_analysis/WeightedTreemaps/data/clustered_texts_2stage (3).csv" # ご自身のファイルパスに変更

# ファイルの存在確認
cat("### ファイルの存在確認 ###\n")
if (!file.exists(file_path)) {
  stop("エラー: ファイル '", file_path, "' が見つかりません。")
} else {
  cat("ファイル '", file_path, "' が見つかりました。\n")
}

# ファイルサイズの確認
file_size <- file.info(file_path)$size
cat("ファイルサイズ:", file_size, "バイト\n")
if (file_size == 0) {
  stop("エラー: ファイルは空です。")
}

# CSVファイルの読み込み
cat("### CSVファイルの読み込み ###\n")
tryCatch({
  # エンコーディングを明示的に指定しないバージョンで読み込む
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  cat("エンコーディング自動判定で正常に読み込まれました。\n")
  cat("読み込まれた行数:", nrow(df), "\n")
  cat("読み込まれた列数:", ncol(df), "\n")
  cat("列名:", paste(names(df), collapse = ", "), "\n")
  str(df) # データフレームの構造を確認
  print(head(df))  # 最初の数行を表示
}, error = function(e) {
  cat("エンコーディング自動判定でエラーが発生しました。\n")
  
  tryCatch({
      df <- read.csv(file_path, fileEncoding = "CP932", stringsAsFactors = FALSE)
      cat("エンコーディング 'CP932' で正常に読み込まれました。\n")
      cat("読み込まれた行数:", nrow(df), "\n")
      cat("読み込まれた列数:", ncol(df), "\n")
      cat("列名:", paste(names(df), collapse = ", "), "\n")
       str(df)
       print(head(df))
  },error = function(e){
      tryCatch({
        df <- read.csv(file_path, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
        cat("エンコーディング 'UTF-8' で正常に読み込まれました。\n")
        cat("読み込まれた行数:", nrow(df), "\n")
        cat("読み込まれた列数:", ncol(df), "\n")
        cat("列名:", paste(names(df), collapse = ", "), "\n")
        str(df)
        print(head(df))
    },error = function(e){
       stop("エンコーディング 'UTF-8' と 'CP932'と自動判断でエラーが発生しました。他のエンコーディングをお試しくださいまたはファイル形式を再度ご確認ください。\n", e$message)
    })
  })

})

# 生データの確認
cat("### 生データの確認 ###\n")
cat("df の行数:", nrow(df), "\n")
cat("df の列名:", names(df), "\n")
print(head(df))

# データが空の場合の処理
if (nrow(df) == 0) {
  stop("エラー: CSVファイルにデータが含まれていません。ファイルの内容を確認してください。")
}

# 必要な列の存在確認
required_cols <- c("primary_cluster_name", "secondary_cluster_name", "primary_cluster_ratio", "secondary_cluster_ratio")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("エラー: 必要な列が欠けています: ", paste(missing_cols, collapse = ", "))
}


# データの前処理
cat("### 前処理の確認 ###\n")

# 初期データ
cat("初期データ df の行数:", nrow(df), "\n")


# "未定義" を含む行を除外
cat("primary_cluster_name の分布:\n")
print(table(df$primary_cluster_name))
cat("secondary_cluster_name の分布:\n")
print(table(df$secondary_cluster_name))

if (any(df$primary_cluster_name == "未定義") || any(df$secondary_cluster_name == "未定義")) {
  df <- df %>%
    filter(primary_cluster_name != "未定義",
           secondary_cluster_name != "未定義")
  cat("未定義を除外後の行数:", nrow(df), "\n")
} else {
  cat("未定義は含まれていません。フィルタリングをスキップします。\n")
}

# 各クラスタのコメント数をカウント、ratioも保持
df_counts <- df %>%
  group_by(secondary_cluster_name, primary_cluster_name, primary_cluster_ratio, secondary_cluster_ratio) %>% # .drop = FALSE は不要
  summarise(count = n(), .groups = 'drop')
cat("df_counts の行数（未定義フィルタリング直後）:", nrow(df_counts), "\n")
print(head(df_counts))

# count が0または負のものを除外(一時的にコメントアウト)
# df_counts <- df_counts %>% filter(count > 0)
cat("count > 0 フィルタ後の行数:", nrow(df_counts), "\n")
print(head(df_counts))

# データが空の場合の処理
if (nrow(df_counts) == 0) {
  stop("エラー: データが空です。フィルタリング条件を見直してください。")
}

# サブセット化をスキップし、すべての行を使用
cat("### すべての行を使用します ###\n")
cat("処理対象の行数:", nrow(df_counts), "\n")

# Voronoiツリーマップの生成
# 計算時間を計測
cat("### Treemap の生成開始 ###\n")
flush.console()  # 出力を強制的にフラッシュ
system.time({
  tm <- voronoiTreemap(
    data = df_counts,  # サブセットをスキップし、すべての行を使用
    levels = c("secondary_cluster_name", "primary_cluster_name"),
    cell_size = "count",
    shape = "rounded_rect",
    positioning = "clustered",
    error_tol = 0.01,  # 許容誤差をさらに大きくして計算負荷を軽減
    maxIteration = 100,  # 最大反復回数をさらに減らして計算負荷を軽減
    verbose = TRUE  # 進捗状況を表示
  )
}) -> treemap_time

cat("### Treemap の生成時間 ###\n")
print(treemap_time)


# 画像サイズとファイル名の指定
width_px <- 1200
height_px <- 675
dpi <- 300
output_filename <- "treemap.png"

# ツリーマップの描画時に使用するパラメータを事前に設定
label_level_val <- c(1, 2)
label_color_val <- c(adjustcolor(grey(1.0), alpha.f = 0.85), adjustcolor(grey(0.9), alpha.f = 0.80))
label_size_val <- c(8, 8) # デフォルトからサイズ調整
border_size_val <- 1.0

# ツリーマップのgrobを取得、空文字列でデフォルトタイトルを無効化し、カスタムタイトルを追加
cat("### Treemap の描画開始 ###\n")
system.time({
    treemap_grob <- grid.grabExpr({
        grid::grid.newpage()
        
        # drawTreemap関数の前にテキストを描画するのは間違いです。
        # drawTreemap関数内で適切な位置に描画します。

        drawTreemap(tm,
                    label_level = label_level_val,
                    label_color = label_color_val,
                    label_size = label_size_val,
                    label_autoscale = TRUE,
                    legend = FALSE,
                    title = "",  # 空文字列を使用してデフォルトタイトルを無効化
                    color_type = "categorical",
                    color_level = 1,
                    border_size = border_size_val)

          # カスタムタイトルを「Zen Maru Gothic」で追加
        grid::grid.text(
        "「米津玄師 Kenshi Yonezu - BOW AND ARROW」コメントクラスター", # タイトル
        x = unit(0.5, "npc"),
        y = unit(0.95, "npc"),
        gp = gpar(fontsize = 18, fontfamily = "zenmaru", col = grey(0))
    )
    }, wrap.grobs = TRUE)
}) -> draw_time

cat("### Treemap の描画時間 ###\n")
print(draw_time)


# コメント総数の計算
total_reviews <- nrow(df)

# キャプションの作成、「Zen Maru Gothic」を適用
caption_grob <- textGrob(
  paste0("YouTube動画の", format(total_reviews, big.mark = ","), "件のコメントより徒然研究室が2025年3月10日作成"),
  x = unit(0.01, "npc"),
  y = unit(0.01, "npc"),
  hjust = 0,
  vjust = 0,
  gp = gpar(fontsize = 8, fontfamily = "zenmaru", col = grey(0))
)

# 最終プロットの作成
final_plot <- ggplot() +
  annotation_custom(treemap_grob) +
  annotation_custom(caption_grob) +
  theme_void() +
  theme(plot.background = element_rect(fill = "#F9F4EE", color = NA))

# PNGファイルとして保存し、表示
ggsave(output_filename, final_plot, width = width_px, height = height_px, units = "px", dpi = dpi)
display_png(file = output_filename, width = width_px, height = height_px)


# 全体の処理が終了したことを確認
cat("### 処理が正常に終了しました ###\n")

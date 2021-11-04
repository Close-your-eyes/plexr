#' Standard report on magpix measurement
#'
#'
#' @param df the data frame returned by read_magpix; filter this df to focus on certain analytes or samples only
#' @param count_threshold minimal bead count to include in the
#' @param meta_data a data frame of meta_data; should contain a column named well_id for joining and a column named group_col (next argument) for grouping observations
#' @param group_col column name in meta_data with values to group observations
#' @param out_folder path to the folder to save plots and tables to
#'
#' @return a data frame of plots and tables may also be saved to out_folder if provided
#' @export
#'
#' @examples
#' \dontrun{
#' # get the path to the parent folder of a script
#' wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
#' report_magpix(df = read_magpix('path_to_your_file.csv'),
#' count_threshold = 30,
#' meta_data = openxlsx::read.xlsx('path_to_meta_data_xlsx'),
#' out_folder = file.path(wd, "report"))
#' }
report_magpix <- function(df,
                          count_threshold = 30,
                          meta_data,
                          group_col = "group",
                          out_folder) {

  ggplot2::theme_set(ggplot2::theme_bw())
  ggplot2::theme_update(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), text = ggplot2::element_text(family = "Courier"))

  var1 <- paste0("n_analytes_below_", count_threshold)
  var2 <- paste0("bead_count_below_", count_threshold)

  df_count_summary <-
    df %>%
    dplyr::group_by(Sample, well_id, well_row, well_col) %>%
    dplyr::summarise(mean_count_per_analyte = mean(Count), "{var1}" := sum(Count < count_threshold), .groups = "drop") %>%
    dplyr::group_by(well_col, well_row) %>%
    dplyr::mutate("{var2}" := any(!!dplyr::sym(var1) > 0)) %>%
    dplyr::ungroup()
  p1 <- ggplot2::ggplot(df_count_summary, ggplot2::aes(x = well_col, y = well_row, fill = !!dplyr::sym(var2), label = !!dplyr::sym(var1))) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::geom_text() +
    ggplot2::scale_fill_manual(values = c("#999999", "#E69F00"))

  p2 <- ggplot2::ggplot(df %>% dplyr::mutate("{var2}" := Count < count_threshold), ggplot2::aes(x = well_col, y = well_row, fill = !!dplyr::sym(var2))) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::scale_fill_manual(values = c("#999999", "#E69F00")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::facet_wrap(dplyr::vars(analyte))

  df_cor <-
    df %>%
    dplyr::filter(Count >= count_threshold) %>%
    dplyr::filter(sample_type == "Unknown") %>%
    dplyr::select(-Count, -Result) %>%
    tidyr::pivot_wider(names_from = analyte, values_from = Median) %>%
    dplyr::select(unique(df$analyte))
  cor_mat <- stats::cor(df_cor, df_cor, use = "complete.obs")
  cor_mat_upper <- cor_mat
  cor_mat_upper[upper.tri(cor_mat_upper)] <- NA
  cor_df_save <-
    as.data.frame(cor_mat_upper) %>%
    tibble::rownames_to_column("analyte1") %>%
    tidyr::pivot_longer(cols = -analyte1, names_to = "analyte2", values_to = "pearson_corr") %>%
    tidyr::drop_na(pearson_corr) %>%
    dplyr::filter(analyte1 != analyte2)

  cor_df_plot <-
    as.data.frame(cor_mat) %>%
    tibble::rownames_to_column("analyte1") %>%
    tidyr::pivot_longer(cols = -analyte1, names_to = "analyte2", values_to = "pearson_corr") %>%
    dplyr::filter(analyte1 != analyte2)
  p5 <- ggplot2::ggplot(cor_df_plot, ggplot2::aes(x = analyte1, y = analyte2, fill = pearson_corr)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")), interpolate = "linear")(100), na.value = "white") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))


  if (!missing(meta_data)) {
    df <- dplyr::left_join(df, meta_data)
    df_medians <-
      df %>%
      dplyr::filter(Count >= 30) %>%
      dplyr::filter(sample_type == "Unknown") %>%
      dplyr::group_by(analyte) %>%
      dplyr::mutate(Median_zscore = scale(Median)) %>%
      dplyr::mutate(is.outlier = !dplyr::between(Median_zscore, stats::median(Median_zscore) - 3 * stats::mad(Median_zscore), stats::median(Median_zscore) + 3 * stats::mad(Median_zscore))) %>%
      dplyr::mutate(is.extreme.outlier = !dplyr::between(Median_zscore, stats::median(Median_zscore) - 5 * stats::mad(Median_zscore), stats::median(Median_zscore) + 5 * stats::mad(Median_zscore))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(analyte, group) %>%
      dplyr::mutate(avg_group_analyte_Median_zscore = mean(Median_zscore)) %>%
      dplyr::group_by(analyte) %>%
      dplyr::mutate(rank = dplyr::dense_rank(avg_group_analyte_Median_zscore)) %>%
      dplyr::ungroup()

    p3 <- ggplot2::ggplot(df_medians, ggplot2::aes(x = group, y = Median_zscore)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_boxplot(data = df_medians %>% dplyr::group_by(analyte) %>% dplyr::filter(rank == max(rank)), outlier.shape = NA, color = "red") +
      ggplot2::geom_point(ggplot2::aes(color = is.outlier)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_text(data = df_medians %>% dplyr::filter(is.outlier), ggplot2::aes(label = well_id), hjust = "top") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = c("#999999", "#E69F00")) +
      ggpubr::stat_compare_means(method = "kruskal.test", vjust = 1, label = "p.format") +
      ggplot2::facet_wrap(dplyr::vars(analyte), scales = "free_y")

    df_medians2 <-
      df %>%
      dplyr::filter(Count >= 30) %>%
      dplyr::filter(sample_type == "Unknown") %>%
      dplyr::group_by(analyte) %>%
      dplyr::mutate(Median_zscore = scale(Median)) %>%
      dplyr::mutate(is.outlier = !dplyr::between(Median_zscore, stats::median(Median_zscore) - 3 * stats::mad(Median_zscore), stats::median(Median_zscore) + 3 * stats::mad(Median_zscore))) %>%
      dplyr::mutate(is.extreme.outlier = !dplyr::between(Median_zscore, stats::median(Median_zscore) - 5 * stats::mad(Median_zscore), stats::median(Median_zscore) + 5 * stats::mad(Median_zscore))) %>%
      dplyr::filter(!is.extreme.outlier) %>%
      dplyr::mutate(Median_zscore = scale(Median)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(analyte, group) %>%
      dplyr::mutate(avg_group_analyte_Median_zscore = mean(Median_zscore)) %>%
      dplyr::group_by(analyte) %>%
      dplyr::mutate(rank = dplyr::dense_rank(avg_group_analyte_Median_zscore)) %>%
      dplyr::ungroup()

    p4 <- ggplot2::ggplot(df_medians2 , ggplot2::aes(x = group, y = Median_zscore)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_boxplot(data = df_medians2 %>% dplyr::group_by(analyte) %>% dplyr::filter(rank == max(rank)), outlier.shape = NA, color = "red") +
      ggplot2::geom_point(ggplot2::aes(color = is.outlier)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::geom_text(data = df_medians2 %>% dplyr::filter(is.outlier), ggplot2::aes(label = well_id), hjust = "top") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = c("#999999", "#E69F00")) +
      ggplot2::ggtitle("extreme outliers excluded") +
      ggpubr::stat_compare_means(method = "kruskal.test", vjust = 1, label = "p.format") +
      ggplot2::facet_wrap(dplyr::vars(analyte), scales = "free_y")
  } else {
    p3 <- NULL
    p4 <- NULL
  }

  if (!missing(out_folder)) {
    dir.create(out_folder, recursive = T)
    ggplot2::ggsave(plot = p1, filename = "p1.png", device = "png", path = out_folder, width = 6, height = 6)
    ggplot2::ggsave(plot = p2, filename = "p2.png", device = "png", path = out_folder, width = 16, height = 16)
    if (!missing(meta_data)) {
      ggplot2::ggsave(plot = p3, filename = "p3.png", device = "png", path = out_folder, width = 16, height = 16)
      ggplot2::ggsave(plot = p4, filename = "p4.png", device = "png", path = out_folder, width = 16, height = 16)
    }
    ggplot2::ggsave(plot = p5, filename = "p5.png", device = "png", path = out_folder, width = 8, height = 8)
    openxlsx::write.xlsx(cor_df_save, file.path(out_folder, "correlations.xlsx"))
  }

   return(list(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,cor_df=cor_df_save))
}



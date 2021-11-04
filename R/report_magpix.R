report_magpix <- function(df, count_threshold = 30, meta_data, group_col = "group") {


  var1 <- paste0("n_analytes_below_", count_threshold)
  var2 <- paste0("bead_count_below_", count_threshold)

  df_count_summary <-
    df %>%
    dplyr::group_by(Sample, well_id, well_row, well_col) %>%
    dplyr::summarise(mean_count_per_analyte = mean(Count), "{var1}" := sum(Count < count_threshold), .groups = "drop") %>%
    dplyr::group_by(well_col, well_row) %>%
    dplyr::mutate("{var2}" := any(!!sym(var1) > 0)) %>%
    dplyr::ungroup()
  p1 <- ggplot(df_count_summary, aes(x = well_col, y = well_row, fill = !!sym(var2), label = !!sym(var1))) +
    geom_tile(color = "black") +
    geom_text() +
    scale_fill_manual(values = c("#999999", "#E69F00"))

  p2 <- ggplot(df %>% dplyr::mutate("{var2}" := Count < count_threshold), aes(x = well_col, y = well_row, fill = !!sym(var2))) +
    geom_tile(color = "black") +
    scale_fill_manual(values = c("#999999", "#E69F00")) +
    facet_wrap(vars(analyte))

  df <-
    df %>%
    dplyr::left_join(meta_data)


  df_medians <-
    df %>%
    dplyr::filter(Count >= 30) %>%
    dplyr::filter(sample_type == "Unknown") %>%
    dplyr::group_by(analyte) %>%
    dplyr::mutate(Median_zscore = scale(Median)) %>%
    dplyr::mutate(is.outlier = !dplyr::between(Median_zscore, median(Median_zscore) - 3 * mad(Median_zscore), median(Median_zscore) + 3 * mad(Median_zscore))) %>%
    dplyr::mutate(is.extreme.outlier = !dplyr::between(Median_zscore, median(Median_zscore) - 5 * mad(Median_zscore), median(Median_zscore) + 5 * mad(Median_zscore))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(analyte, group) %>%
    dplyr::mutate(avg_group_analyte_Median_zscore = mean(Median_zscore)) %>%
    dplyr::group_by(analyte) %>%
    dplyr::mutate(rank = dense_rank(avg_group_analyte_Median_zscore)) %>%
    dplyr::ungroup()

  p3 <- ggplot(df_medians, aes(x = group, y = Median_zscore)) +
    geom_boxplot(outlier.shape = NA) +
    geom_boxplot(data = df_medians %>% dplyr::group_by(analyte) %>% dplyr::filter(rank == max(rank)), outlier.shape = NA, color = "red") +
    geom_point(aes(color = is.outlier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(data = df_medians %>% dplyr::filter(is.outlier), aes(label = well_id), hjust = "top") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank()) +
    scale_color_manual(values = c("#999999", "#E69F00")) +
    ggpubr::stat_compare_means(method = "kruskal.test", vjust = 1, label = "p.format") +
    facet_wrap(vars(analyte), scales = "free_y")
  p3

  df_medians2 <-
    df %>%
    dplyr::filter(Count >= 30) %>%
    dplyr::filter(sample_type == "Unknown") %>%
    dplyr::group_by(analyte) %>%
    dplyr::mutate(Median_zscore = scale(Median)) %>%
    dplyr::mutate(is.outlier = !dplyr::between(Median_zscore, median(Median_zscore) - 3 * mad(Median_zscore), median(Median_zscore) + 3 * mad(Median_zscore))) %>%
    dplyr::mutate(is.extreme.outlier = !dplyr::between(Median_zscore, median(Median_zscore) - 5 * mad(Median_zscore), median(Median_zscore) + 5 * mad(Median_zscore))) %>%
    dplyr::filter(!is.extreme.outlier) %>%
    dplyr::mutate(Median_zscore = scale(Median)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(analyte, group) %>%
    dplyr::mutate(avg_group_analyte_Median_zscore = mean(Median_zscore)) %>%
    dplyr::group_by(analyte) %>%
    dplyr::mutate(rank = dense_rank(avg_group_analyte_Median_zscore)) %>%
    dplyr::ungroup()

  p4 <- ggplot(df_medians2 , aes(x = group, y = Median_zscore)) +
    geom_boxplot(outlier.shape = NA) +
    geom_boxplot(data = df_medians2 %>% dplyr::group_by(analyte) %>% dplyr::filter(rank == max(rank)), outlier.shape = NA, color = "red") +
    geom_point(aes(color = is.outlier)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(data = df_medians2 %>% dplyr::filter(is.outlier), aes(label = well_id), hjust = "top") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank()) +
    scale_color_manual(values = c("#999999", "#E69F00")) +
    ggtitle("extreme outliers excluded") +
    ggpubr::stat_compare_means(method = "kruskal.test", vjust = 1, label = "p.format") +
    facet_wrap(vars(analyte), scales = "free_y")
  p4


}

pheno_data <- read_magpix("/Volumes/AG_Hiepe-1/Christopher.Skopnik/202010-202109_Experiments/20211025_Stefan_Frischbutter_ZytokinDaten/Zytokindaten/210625_45plex_Neutroderm_Jana.csv")
df <- pheno_data[["df"]]

meta_data <-
  openxlsx::read.xlsx("/Volumes/AG_Hiepe-1/Christopher.Skopnik/202010-202109_Experiments/20211025_Stefan_Frischbutter_ZytokinDaten/Zytokindaten/Neutroderm_Serum Aliquots.xlsx", sheet = 2) %>%
  dplyr::mutate(group = tolower(group))

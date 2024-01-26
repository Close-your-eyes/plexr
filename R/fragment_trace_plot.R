#' Standard plot of nucleic acid trace data from Fragment- or Bio-Analyzer
#'
#' Create a standardized plot for nucleid acid trace data. The electropherogram
#' table is joined with quality table. Hence, RQN values or concentration or R28S_18S (ratio 28S/18S) can be
#' used as a facetting variable. Combinations are possible: RQN_R28S_18S,
#' RQN_conc, R28S_18S_conc, RQN_R28S_18S_conc.
#'
#' @param FA_data list of nucleic acid trace data read with plexr::read_fragmentanalyzer
#' @param electropherogram_entry name of the electropherogram list entry
#' @param quality_entry name of the quality list entry
#' @param x_breaks breaks on x-axis
#' @param theme ggplot2 theme
#' @param theme_args arguments to modify the theme
#' @param x_lims x-axis limits
#' @param facetting how to facet the plot; see example
#'
#' @return a ggplot2 plot object
#' @export
#'
#' @examples
#' \dontrun{
#' FA_files <- list.files(FA_folder, pattern = "\\.csv", full.names = T)
#' FA_data <- plexr::read_fragmentanalyzer(FA_files)
#' fragment_trace_plot(FA_data, facetting = ggplot2::facet_wrap(vars(chemical, Buffer, DTT_mM, RQN), scales = "free_y"))
#' }
fragment_trace_plot <- function(FA_data,
                                electropherogram_entry = "electropherogram",
                                quality_entry = "quality",
                                x_breaks = c(15,500,1500,3000,6000),
                                theme = ggplot2::theme_bw(),
                                theme_args = list(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                                  panel.grid.minor = element_blank(),
                                                  strip.background = element_rect(fill = "cornsilk")),
                                facetting = NULL,
                                x_lims = c(5,NA)) {

  if (!is.list(FA_data)) {
    stop("FA_data has to be a list.")
  }

  if (!electropherogram_entry %in% names(FA_data)) {
    stop(electropherogram_entry, " not found in FA_data.")
  }

  if (!quality_entry %in% names(FA_data)) {
    stop(quality_entry, " not found in FA_data.")
  }

  FA_data_plot <-
    FA_data[[electropherogram_entry]] %>%
    dplyr::left_join(FA_data[[quality_entry]] %>%
                       dplyr::mutate(RQN = paste0("RQN = ", RQN), conc = paste0(conc, " ", unit), R28S_18S = paste0("28S/18S = ", R28S_18S)) %>%
                       dplyr::select(RQN, sample_ID, conc, R28S_18S) %>%
                       dplyr::mutate(RQN_R28S_18S = paste0(RQN, ", ", R28S_18S)) %>%
                       dplyr::mutate(RQN_conc = paste0(RQN, ", ", conc)) %>%
                       dplyr::mutate(R28S_18S_conc = paste0(R28S_18S, ", ", conc)) %>%
                       dplyr::mutate(RQN_R28S_18S_conc = paste0(RQN, ", ", R28S_18S, ", ", conc)))

  FA_plot <- ggplot(FA_data_plot, aes(x = size_nt, y = signal)) +
    geom_line(color = "tomato2") +
    scale_x_log10(breaks = x_breaks, limits = x_lims)

  FA_plot <- FA_plot + theme
  FA_plot <- FA_plot + do.call(ggplot2::theme, args = theme_args)

  if (!is.null(facetting)) {
    FA_plot <- FA_plot + facetting
  }

  return(FA_plot)

}

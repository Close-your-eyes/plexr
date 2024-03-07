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
#' @param line_args arguments to ggplot2::geom_line
#' @param text_args arguments to ggrepel::geom_text_repel; pass empty list to omit
#' plotting text annotation; pass a list of lists to plot multiple annotations;
#' within ggplot2::aes it is selected which variable is being plotted
#' @param text_fun which function to use for plotting text annotation;
#' pass a list of functions for every list element in text_args or one function
#' for every element in text_args
#' @param sample_ID_col_name column name of sample_IDs, used for joining data tables from FA
#'
#' @return a ggplot2 plot object
#' @export
#'
#' @examples
#' \dontrun{
#' FA_files <- list.files(FA_folder, pattern = "\\.csv", full.names = T)
#' FA_data <- plexr::read_fragmentanalyzer(FA_files)
#' # color = "tomato2" is passed to geom_line
#' fragment_trace_plot(FA_data, facetting = ggplot2::facet_wrap(vars(chemical, Buffer, DTT_mM), scales = "free_y"))
#' # use the deprecated ggplot2::aes_q to define color by multiple columns
#' fragment_trace_plot(FA_data, line_args = list(mapping = ggplot2::aes_q(color = ~paste(chemical, Buffer, DTT_mM))))
#' # add RQN value as text annotation
#' fragment_trace_plot(FA_data, facetting = ggplot2::facet_grid(cols = vars(chemical, Buffer), rows = vars(DTT_mM), scales = "free_y"), text_args = list(data = FA_data[["quality"]], mapping = ggplot2::aes(label = RQN), x = -Inf, y = Inf))
#' }
fragment_trace_plot <- function(FA_data,
                                electropherogram_entry = "electropherogram",
                                quality_entry = "quality",
                                sample_ID_col_name = "sample_ID_FA",
                                x_breaks = c(15,500,1500,3000,6000),
                                theme = ggplot2::theme_bw(),
                                theme_args = list(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                                                  panel.grid.minor = ggplot2::element_blank(),
                                                  strip.background = ggplot2::element_rect(fill = "cornsilk")),
                                facetting = NULL,
                                x_lims = c(5,NA),
                                text_fun = list(ggrepel::geom_label_repel, ggrepel::geom_text_repel),
                                text_args = list(list(data = FA_data[[quality_entry]],
                                                      mapping = ggplot2::aes(label = format(RQN, nsmall = 1)),
                                                      x = -Inf, y = Inf,
                                                      size = 4,
                                                      min.segment.length = Inf),
                                                 list(data = FA_data[[quality_entry]],
                                                      mapping = ggplot2::aes(label = chemical_code),
                                                      x = -Inf, y = -Inf,
                                                      size = 3,
                                                      min.segment.length = Inf)),
                                line_args = list(color = "tomato2")) {

  # setting x or y to Inf or -Inf fails to work in ggrepel 0.9.5; but works in ggrepel 0.9.4
  # https://github.com/slowkow/ggrepel/issues/252
  # but works again with the latest, fixed version, which enabled x or y to be INF
  # devtools::install_github("slowkow/ggrepel")

  if (!is.list(FA_data)) {
    stop("FA_data has to be a list.")
  }
  if (is.list(text_args) && any(!sapply(text_args, is.list))) {
    # this tests if the list is a list of lists
    # wraps a list around text_args if text_args is an immediate list of text arguments
    text_args <- list(text_args)
  }

  if (length(text_fun) < length(text_args) && length(text_fun) != 1) {
    stop("length of text_fun unequal to length of text_args. provide equal length of both or text_fun of length 1.")
  }

  if (!is.list(text_fun)) {
    text_fun <- list(text_fun)
  }

  if (length(text_fun) < length(text_args)) {
    text_fun <- rep(text_fun, length(text_args))
  }

  if (!electropherogram_entry %in% names(FA_data)) {
    stop(electropherogram_entry, " not found in FA_data.")
  }

  if (!quality_entry %in% names(FA_data)) {
    stop(quality_entry, " not found in FA_data.")
  }

  if (!sample_ID_col_name %in% names(FA_data[[electropherogram_entry]]) && !sample_ID_col_name %in% names(FA_data[[quality_entry]])) {
    stop(sample_ID_col_name, " not found in ", electropherogram_entry, " and ", quality_entry, " slots of FA_data.")
  }


  if (any(grepl("RQN", names(facetting[["params"]][["facets"]])))) {
    message("Did you put RQN or RQN2 into facetting vars? You may want to avoid that.")
  }
  if (any(grepl("chemical_code", names(facetting[["params"]][["facets"]])))) {
    message("Did you put chemical_code into facetting vars? You may want to avoid that.")
  }


  FA_data_plot <-
    FA_data[[electropherogram_entry]] %>%
    dplyr::left_join(FA_data[[quality_entry]] %>%
                       dplyr::mutate(RQN = paste0("RQN = ", RQN), conc = paste0(conc, " ", unit), R28S_18S = paste0("28S/18S = ", R28S_18S)) %>%
                       dplyr::select(RQN, !!rlang::sym(sample_ID_col_name), conc, R28S_18S) %>%
                       dplyr::mutate(RQN_R28S_18S = paste0(RQN, ", ", R28S_18S)) %>%
                       dplyr::mutate(RQN_conc = paste0(RQN, ", ", conc)) %>%
                       dplyr::mutate(R28S_18S_conc = paste0(R28S_18S, ", ", conc)) %>%
                       dplyr::mutate(RQN_R28S_18S_conc = paste0(RQN, ", ", R28S_18S, ", ", conc)),
                     by = dplyr::join_by(!!rlang::sym(sample_ID_col_name)))

  FA_plot <- ggplot2::ggplot(FA_data_plot, ggplot2::aes(x = size_nt, y = signal)) +
    do.call(ggplot2::geom_line, args = line_args) +
    ggplot2::scale_x_log10(breaks = x_breaks, limits = x_lims)

  if (length(text_args) > 0) {
    for (i in seq_along(text_args)) {
      FA_plot <-
        FA_plot +
        do.call(text_fun[[i]], args = text_args[[i]])
    }
  }

  FA_plot <- FA_plot + theme
  FA_plot <- FA_plot + do.call(ggplot2::theme, args = theme_args)

  if (!is.null(facetting)) {
    FA_plot <- FA_plot + facetting
  }

  return(FA_plot)

}

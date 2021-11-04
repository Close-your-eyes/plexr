#' Read csv file from MAGPIX assay
#'
#' Read the csv file which is the standard output of MAGPIX.
#'
#' @param file_path character, path to the csv file
#'
#' @return list of all table parts and a processed data frame of Count, Median and Result
#' @export
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' read_magpix('path_to_your_file.csv')
#' }
read_magpix <- function(file_path) {

  if (rev(strsplit(basename(file_path), "\\.")[[1]])[1] != "csv") {
    stop("file_path has to be a csv file.")
  }
  file <- suppressMessages(as.data.frame(readr::read_csv(file_path, skip_empty_rows = F, show_col_types = F)))

  lines <- suppressWarnings(which(file[,1] == "DataType:") + 1)
  names(lines) <- file[lines-1,2]
  read_part <- function(l1, l2, url)  {
    as.data.frame(readr::read_csv(url, skip = l1, show_col_types = F, n_max = l2))
  }

  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*10)
  files <- mapply(read_part, lines, c(lines[-1], Inf) - lines - 3, SIMPLIFY = F, MoreArgs = list(url = file_path))

  ## ---- process_file --------

  df <- purrr::reduce(lapply(c("Count", "Median", "Result"), function(x) {
    files[[x]] %>%
      tidyr::pivot_longer(cols = -c(Location, Sample, `Total Events`), names_to = "analyte", values_to = x)
  }), dplyr::left_join, by = c("Location", "Sample", "Total Events", "analyte"))

  df <- tryCatch({
    df %>%
      dplyr::mutate(well_number = as.numeric(stringr::str_extract(Location, "^[:digit:]{1,}"))) %>%
      dplyr::mutate(well_id = stringr::str_extract(Location, "[:alpha:]{1}[:digit:]{1,}")) %>%
      dplyr::mutate(well_row = stringr::str_extract(well_id, "[:alpha:]")) %>%
      dplyr::mutate(well_col = stringr::str_extract(well_id, "[:digit:]{1,}")) %>%
      dplyr::mutate(well_col = factor(well_col, levels = as.character(sort(as.numeric(unique(well_col)))))) %>%
      dplyr::mutate(sample_type = stringr::str_extract(Sample, "[:alpha:]{1,}"))
  }, error = function(e) {
    print("df could not be further processed")
  })

  tryCatch({
    df$well_id <- factor(df$well_id, levels = df %>% dplyr::arrange(well_row, well_col) %>% dplyr::distinct(well_id) %>% dplyr::pull(well_id))
    df$Sample <- factor(df$Sample, levels = unique(df$Sample)[order(as.numeric(stringr::str_extract(unique(df$Sample), "[:digit:]{1,}")))])
  }, error = function(e) {
    print("df could not be further processed")
  })

  return(list(table_parts = files, df = df))

}




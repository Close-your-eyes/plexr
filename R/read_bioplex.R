#' Read data from Bioplex Managers standard output xlsx
#'
#' Get standard values, sample values and parameters from standard curve as calculated by Bioplex Manager.
#' This function was written based a file from Bioplex Manager version 4.
#' Every sheet in the xlsx represents one analyte from the assay.
#' Reading the xlsx file is done by openxlsx::read.xlsx.
#'
#' Conc column is derived from Exp Conc in case of standard and blank; and from Obs Conc in case
#' of samples. OOR values are made NA; extrapolated values (with one asterisks (*)) are simply used
#' as they are but asterisks are removed. Warning that NAs are created by conversion
#' to numeric is normal. Columns to document previous OOR or *values are added. FI - Bkgd column
#' is renamed to FI_bkgd_corr.
#' NLS model is simply created based on values from file.
#'
#' @param file path to the xlsx file (export_per_analyte)
#' @param sheets names of sheets to consider; if NULL all sheets are considered
#' @param rows rows which contain non-aggregated measured values including the header;
#' e.g. 59:160 or c(59, 160); if NULL the range is guessed by occurrence of "Analyte"
#' in the first column of worksheets
#' @param list_names one of c("analyte", "region") or both; both will leave names as they occur in file
#' @param calc_alt_5pl_model calculate an alternative 5pl regression with plexr::alt_5pl_model
#' @param ... arguments passed to plexr::alt_5pl_model
#'
#' @importFrom magrittr "%>%"
#'
#' @return
#' @export
#'
#' @examples
read_bioplex <- function(file,
                         sheets = NULL,
                         rows = NULL,
                         list_names = "analyte",
                         calc_alt_5pl_model = F,
                         ...) {

  options(warn = 1)
  if (!file.exists(file)) {
    stop("file not found, is the path correct?")
  }

  list_names <- match.arg(list_names, c("analyte", "region"), several.ok = T)

  if (is.null(sheets)) {
    sheets <- openxlsx::getSheetNames(file)
  } else {
    if (!any(sheets %in% openxlsx::getSheetNames(file))) {
      stop(paste0("Sheets ", paste(sheets[!which(sheets %in% openxlsx::getSheetNames(file))], collapse = ", "), " not found in file."))
    }
  }

  list <- sapply(sheets, read_bioplex_file, file = file, rows = rows, simplify = F, USE.NAMES = T)

  if (calc_alt_5pl_model) {
    list <- lapply(list, alt_5pl_model, model_name = "alt_nls_model", ...)
  }

  if (length(list_names) == 1) {
    if (list_names == "analyte") {
      names(list) <- stringr::str_replace(names(list), " \\([:digit:]{1,}\\)", "")
    }
    if (list_names == "region") {
      names(list) <- stringr::str_extract(names(list), "([:digit:]{1,})")
    }
  }

  return(list)
}


read_bioplex_file <- function(sheet, file, rows) {

  print(sheet)
  temp <- openxlsx::read.xlsx(xlsxFile = file, sheet = sheet, colNames = F, skipEmptyRows = F)
  if (is.null(rows)) {
    start <- max(which(grepl("Analyte", temp[,1])))
    end <- which(is.na(temp[,1]))[which(is.na(temp[,1])) > start] - 1
    rows <- c(start, end)
  }

  std <- openxlsx::read.xlsx(xlsxFile = file, sheet = sheet, rows = min(rows):max(rows), colNames = T, skipEmptyRows = F) %>%
    dplyr::filter(grepl("S|B", Type)) %>%
    dplyr::mutate(Type = ifelse(Type == "B", "B1", Type)) %>%
    dplyr::mutate(sample_type = ifelse(grepl("S", Type), "standard", "blank")) %>%
    dplyr::mutate(sample_num = stringr::str_extract(Type, "[:digit:]{1,}")) %>%
    dplyr::select(sample_type, sample_num, Well, Exp.Conc, FI, `FI.-.Bkgd`, Type, Dilution, Sampling.Errors) %>%
    dplyr::rename("FI_bkgd_corr" = `FI.-.Bkgd`) %>%
    dplyr::mutate(Exp.Conc = ifelse(is.na(Exp.Conc), 0, Exp.Conc)) %>%
    dplyr::mutate(FI = as.numeric(FI), Exp.Conc = as.numeric(Exp.Conc), FI_bkgd_corr = as.numeric(FI_bkgd_corr)) %>%
    dplyr::rename("Conc" = Exp.Conc) %>%
    dplyr::mutate(sheet = sheet) %>%
    dplyr::mutate(analyte = stringr::str_replace(sheet, " \\([:digit:]{1,}\\)", "")) %>%
    dplyr::mutate(region = stringr::str_extract(sheet, "([:digit:]{1,})")) %>%
    dplyr::arrange(FI) %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(replicate = seq(1:n())) %>%
    dplyr::ungroup()

  smp <- openxlsx::read.xlsx(xlsxFile = file, sheet = sheet, rows = min(rows):max(rows), colNames = T, skipEmptyRows = F) %>%
    dplyr::filter(grepl("X", Type)) %>%
    dplyr::mutate(sample_type = "sample") %>%
    dplyr::mutate(sample_num = stringr::str_extract(Type, "[:digit:]{1,}")) %>%
    dplyr::select(sample_type, sample_num, Well, Obs.Conc, FI, `FI.-.Bkgd`, Type, Dilution, Sampling.Errors) %>%
    dplyr::rename("FI_bkgd_corr" = `FI.-.Bkgd`) %>%
    dplyr::mutate(OOR = ifelse(grepl("OOR", Obs.Conc), T, F)) %>%
    dplyr::mutate(expol = ifelse(grepl("\\*", Obs.Conc), T, F)) %>%
    dplyr::mutate(FI = as.numeric(FI), Obs.Conc = as.numeric(gsub("\\*", "", Obs.Conc)), FI_bkgd_corr = as.numeric(FI_bkgd_corr)) %>%
    dplyr::rename("Conc" = Obs.Conc) %>%
    dplyr::mutate(sheet = sheet) %>%
    dplyr::mutate(analyte = stringr::str_replace(sheet, " \\([:digit:]{1,}\\)", "")) %>%
    dplyr::mutate(region = stringr::str_extract(sheet, "([:digit:]{1,})")) %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(replicate = seq(1:n())) %>%
    dplyr::ungroup()


  # this regex matches all numbers, with or without decimals and optional minus as prefix
  all_nums <- "(-)?[:digit:]{1,}(\\.[:digit:]{1,})?"
  std_curve_pars <- stats::setNames(as.numeric(stringr::str_extract_all(temp[which(grepl("Std. Curve", temp[,1])),1], all_nums)[[1]][-c(3,4)]),
                             c("dd", "aa", "cc", "bb", "gg"))


  nls_model <- nls.multstart::nls_multstart(formula = FI ~ dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg,
                                            data = std,
                                            iter = 1,
                                            start_lower = std_curve_pars,
                                            start_upper = std_curve_pars,
                                            lower = std_curve_pars,
                                            upper = std_curve_pars)

  formula = 'FI = dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg'

  return(list(standard = std, samples = smp, std_curve_pars = std_curve_pars, std_curve_formula = formula, nls_model = nls_model))

}


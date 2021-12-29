#' Read data from Bioplex Managers standard output xlsx
#'
#' Get standard values, sample values and parameters from standard curve as calculated by Bioplex Manager.
#' This function was written based a file from Bioplex Manager version 4.
#' Every sheet in the xlsx represents one analyte from the assay.
#' Reading the xlsx file is done by openxlsx::read.xlsx.
#'
#' @param file path to the xlsx file (export_per_analyte)
#' @param sheets names of sheets to consider; if NULL all sheets are considered
#' @param rows rows which contain non-aggregated measured values including the header;
#' e.g. 59:160 or c(59, 160); if NULL the range is guessed by occurrence of "Analyte"
#' in the first column of worksheets
#' @param list_names one of c("analyte", "region") or both; both will leave names as they occur in file
#'
#' @return
#' @export
#'
#' @examples
read_bioplex <- function(file, sheets = NULL, rows = NULL, list_names = "analyte") {

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

  list <- sapply(sheets, read_file, file = file, rows = rows, simplify = F, USE.NAMES = T)

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


read_file <- function(sheet, file, rows) {

  temp <- openxlsx::read.xlsx(xlsxFile = file, sheet = sheet, colNames = F, skipEmptyRows = F)
  if (is.null(rows)) {
    start <- max(which(grepl("Analyte", temp[,1])))
    end <- which(is.na(temp[,1]))[which(is.na(temp[,1])) > start] - 1
    rows <- c(start, end)
  }

  std <- openxlsx::read.xlsx(xlsxFile = file, sheet = sheet, rows = min(rows):max(rows), colNames = T, skipEmptyRows = F) %>%
    dplyr::filter(grepl("S|B", Type)) %>%
    dplyr::select(Exp.Conc, FI, `FI.-.Bkgd`, Type, Dilution, Sampling.Errors) %>%
    dplyr::mutate(Exp.Conc = ifelse(is.na(Exp.Conc), 0, Exp.Conc)) %>%
    dplyr::mutate(FI = as.numeric(FI), Exp.Conc = as.numeric(Exp.Conc)) %>%
    dplyr::rename("Conc" = Exp.Conc) %>%
    dplyr::mutate(sheet = sheet) %>%
    dplyr::mutate(analyte = stringr::str_replace(sheet, " \\([:digit:]{1,}\\)", "")) %>%
    dplyr::mutate(region = stringr::str_extract(sheet, "([:digit:]{1,})"))

  smp <- openxlsx::read.xlsx(xlsxFile = file, sheet = sheet, rows = min(rows):max(rows), colNames = T, skipEmptyRows = F) %>%
    dplyr::filter(grepl("X", Type)) %>%
    dplyr::select(Obs.Conc, FI, `FI.-.Bkgd`, Type, Dilution, Sampling.Errors) %>%
    dplyr::mutate(OOR = ifelse(grepl("OOR", Obs.Conc), T, F)) %>%
    dplyr::mutate(expol = ifelse(grepl("\\*", Obs.Conc), T, F)) %>%
    suppressWarnings(dplyr::mutate(FI = as.numeric(FI), Obs.Conc = as.numeric(gsub("\\*", "", Obs.Conc)))) %>%
    dplyr::rename("Conc" = Obs.Conc) %>%
    dplyr::mutate(sheet = sheet) %>%
    dplyr::mutate(analyte = stringr::str_replace(sheet, "\\([:digit:]{1,}\\)", "")) %>%
    dplyr::mutate(region = stringr::str_extract(sheet, "([:digit:]{1,})"))

  std_curve_pars <- setNames(as.numeric(stringr::str_extract_all(temp[which(grepl("Std. Curve", temp[,1])),1], "(-)?[:digit:]{1,}(\\.[:digit:]{1,})?")[[1]][-c(3,4)]),
                             c("dd", "aa", "cc", "bb", "gg"))

  formula = 'FI = dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg'

  return(list(standard = std, samples = smp, std_curve_pars = std_curve_pars, std_curve_formula = formula))

}


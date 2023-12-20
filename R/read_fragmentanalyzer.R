#' Read csv output files from Fragment Analyzer into data frames
#'
#' @param files full path to all files
#' @param Electropherogram_file_pattern pattern to identify electropherogram files within files
#' @param QualityTable_file_pattern pattern to identify Quality Table files within files
#' @param PeakTable_file_pattern pattern to identify Peak Table files within files
#' @param filter_names vector of sample name patterns to remove, can be regular expressions; set to NULL or empty vector to not filter any sample
#' @param return which type to return, list of tables for each file (list); or one data.frame for all (data.frame)
#'
#' @return list of lists or list of data.frames
#' @export
#'
#' @examples
#' \dontrun{
#' # get all file paths
#' files <- list.files("parent_dir", recursive = T, full.names = T, pattern = "\\.csv")
#' FA_data <- plexr::read_fragmentanalyzer(files = files)
#' }
read_fragmentanalyzer <- function(files,
                                  Electropherogram_file_pattern = "Electropherogram",
                                  QualityTable_file_pattern = "Quality Table",
                                  PeakTable_file_pattern = "Peak Table",
                                  filter_names = c("Ladder", "Samp.[[:digit:]]{1,2}"),
                                  return = c("data.frame", "list")) {

  return <- match.arg(return, c("data.frame", "list"))

  if (length(Electropherogram_file_pattern) > 1) {
    stop("Electropherogram_file_pattern should be length 1.")
  }
  if (length(QualityTable_file_pattern) > 1) {
    stop("QualityTable_file_pattern should be length 1.")
  }
  if (length(PeakTable_file_pattern) > 1) {
    stop("PeakTable_file_pattern should be length 1.")
  }

  ## read Electropherograms
  pherogram_files <- grep(Electropherogram_file_pattern, files, value = T, ignore.case = T)
  if (length(pherogram_files) > 0) {
    if (any(temp <- !file.exists(pherogram_files))) {
      message("Electropherogram files")
      message(paste(pherogram_files[temp], collapse = "\n"))
      message("not found.")
      pherogram_files <- pherogram_files[which(!temp)]
    }

    pherograms <- tryCatch({
      read_electropherograms(subfiles = pherogram_files,
                             filter_names = filter_names)
    }, error = function(e) {
      message("Error in reading Electropherogram files.")
      return(NULL)
    })

  } else {
    message("No Electropherogram files found.")
  }

  ## read QualityTables
  quality_files <- grep(QualityTable_file_pattern, files, value = T, ignore.case = T)
  if (length(quality_files) > 0) {
    if (any(temp <- !file.exists(quality_files))) {
      message("Quality Table files")
      message(paste(quality_files[temp], collapse = "\n"))
      message("not found.")
      quality_files <- quality_files[which(!temp)]
    }

    qualities <- tryCatch({
      read_qualities(subfiles = quality_files,
                     filter_names = filter_names)
    }, error = function(e) {
      message("Error in reading Quality Table files.")
      return(NULL)
    })

  } else {
    message("No Quality Table files found.")
  }

  ## read PeakTables
  peak_files <- grep("Peak Table", files, value = T)
  if (length(peak_files) > 0) {
    if (any(temp <- !file.exists(peak_files))) {
      message("Peak Table files")
      message(paste(peak_files[temp], collapse = "\n"))
      message("not found.")
      peak_files <- peak_files[which(!temp)]
    }

    peak_list <- tryCatch({
      read_peaks(subfiles = peak_files,
                 filter_names = filter_names)
    }, error = function(e) {
      message("Error in reading Peak Table files.")
      return(NULL)
    })

  } else {
    message("No Peak Table files found.")
  }


  if (return == "data.frame") {
    return(list(electropherogram = tibble::as_tibble(pherograms),
                quality = tibble::as_tibble(qualities),
                peaks = tibble::as_tibble(peak_list[["peaks"]]),
                peaks_summary = tibble::as_tibble(peak_list[["summary"]])))
  } else if (return == "list") {
    return(list(electropherogram = split(pherograms,pherograms$file),
                quality = split(qualities,qualities$file),
                peaks = split(peak_list[["peaks"]],peak_list[["peaks"]]$file),
                peaks_summary = split(peak_list[["summary"]], peak_list[["summary"]]$file)))
  }

}


read_electropherograms <- function(subfiles,
                                   filter_names = NULL) {
  pherograms <- purrr::map(stats::setNames(subfiles, basename(subfiles)), read.csv, check.names = F)
  pherograms <- purrr::map_dfr(pherograms, function(x) {
    names(x) <- sapply(strsplit(names(x), ": "), "[", 2)
    names(x)[1] <- "size_nt"
    if (!is.null(filter_names) && length(filter_names) > 1) {
      x <- x[,!apply(sapply(filter_names, grepl, x = names(x), simplify = T), 1, any)]
    }
    x <- tidyr::pivot_longer(x, cols = -size_nt, names_to = "sample_ID", values_to = "signal")
    return(x)
  }, .id = "file")

  pherograms <- pherograms[,c("sample_ID", "size_nt", "signal", "file")]
  return(pherograms)
}

read_qualities <- function(subfiles,
                           filter_names = NULL) {
  qualities <- purrr::map_dfr(stats::setNames(subfiles, basename(subfiles)), read.csv, check.names = F, .id = "file")
  unit <- strsplit(grep("Conc", names(qualities), value = T), " ")[[1]][2]
  unit <- gsub("[\\(\\)]", "", unit)
  if (length(unit) == 0) {
    unit <- "unit not recognized"
  }
  names(qualities) <- c("file", "well", "sample_ID", "conc", "RQN", "R28S_18S")
  if (!is.null(filter_names) && length(filter_names) > 1) {
    qualities <- qualities[!apply(sapply(filter_names, grepl, x = qualities[,"sample_ID",drop=T], simplify = T), 1, any),]
  }
  qualities$unit <- unit
  qualities <- qualities[,c("well", "sample_ID", "conc", "unit", "RQN", "R28S_18S", "file")]
  return(qualities)
}

read_peaks <- function(subfiles,
                       filter_names = NULL) {

  peaks_list <- purrr::map(stats::setNames(subfiles, basename(subfiles)), read.csv, check.names = F, header = F)

  ## loop over files
  peaks_list2 <- purrr::map(peaks_list, function(x) {
    #splits <- which(diff(trimws(x[,1,drop=T]) == "") == -1)
    splits <- which(trimws(x[,2,drop=T]) == "Total Conc.:") + 1
    splits <- splits[-length(splits)]
    split_rows <- brathering::seq2(c(1,splits+1), c(splits, nrow(x)))
    sampleIDs <- x[c(1,splits+1),2]
    x_split <- split(x, rep(sampleIDs, lengths(split_rows)))

    ## loop over samples
    sample_tables <- purrr::map(x_split, function(y) {
      names(y) <- unname(unlist(y[2,,drop=T]))
      y <- y[-c(1:2),]
      split_rows2 <- which(diff(trimws(y[,1]) == "") == 1)
      if (length(split_rows2) == 0) {
        # empty tables sometimes for empty lanes
        return(list(peak_table = NULL, summary_table = NULL))
      } else {
        names(y) <- gsub("/", "_", names(y))
        names(y)[1:3] <- c("peak_ID", "size_nt", "conc_pct")
        z1 <- y[1:split_rows2,]
        rownames(z1) <- seq(1, nrow(z1))
        z1$size_nt <- factor(z1$size_nt, levels = z1$size_nt)
        for (i in 3:6) {
          z1[,i] <- as.numeric(z1[,i,drop=T])
        }

        z2 <- y[(split_rows2+3):(nrow(y)-1),]
        z2 <- z2[,2:4]
        names(z2) <- c("parameter", "amount", "unit")
        z2$parameter <- gsub(":", "", z2$parameter)
        z2$parameter <- gsub("Total Conc.", "Total_Conc", z2$parameter)
        z2$amount <- as.numeric(z2$amount)
        z2$unit <- gsub("/", "_", z2$unit)
        return(list(peak_table = z1, summary_table = z2))
      }
    })
    if (!is.null(filter_names) && length(filter_names) > 1) {
      sample_tables <- sample_tables[!apply(sapply(filter_names, grepl, x = names(sample_tables), simplify = T), 1, any)]
    }
    peak_sample_tables <- dplyr::bind_rows(stats::setNames(sapply(sample_tables, "[", 1), names(sample_tables)), .id = "sample_ID")
    summary_sample_tables <- dplyr::bind_rows(stats::setNames(sapply(sample_tables, "[", 2), names(sample_tables)), .id = "sample_ID")

    return(list(peak_sample_tables = peak_sample_tables, summary_sample_tables = summary_sample_tables))
  })

  # put file column as last
  peak_df <- dplyr::bind_rows(stats::setNames(sapply(peaks_list2, "[", 1), names(peaks_list2)), .id = "file")
  summary_df <- dplyr::bind_rows(stats::setNames(sapply(peaks_list2, "[", 2), names(peaks_list2)), .id = "file")
  peak_df <- peak_df[,c(names(peak_df)[which(names(peak_df) != "file")], "file")]
  summary_df <- summary_df[,c(names(summary_df)[which(names(summary_df) != "file")], "file")]

  return(list (peaks = peak_df, summary = summary_df))
}

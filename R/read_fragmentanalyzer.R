#' Read csv output files from Fragment Analyzer into data frames
#'
#' It may be advisable to provide a set of files from one measurement only at a time.
#' So, if there is more than one folder with csv files from the machine, loop
#' over these folders and provide the respective contained files (their paths) to this function.
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
                                  filter_names = c("Ladder", "Samp.[[:digit:]]{1,2}", "blank"),
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


  sample_name_replace <- NULL # if there is separator in header strings, this will be assigned from read_electropherograms
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
                     filter_names = filter_names,
                     sample_name_replace = sample_name_replace)
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
                 filter_names = filter_names,
                 sample_name_replace = sample_name_replace)
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
                                   filter_names = NULL,
                                   sep = ",") {

  if (any(purrr::map_lgl(subfiles, check_header_for_commas, sep = sep))) {
    message("Greater number of separator '", sep, "' found in header than in row 1 of electropherogram csv file.")
    message("Will try remove every second separator.")
    read_file_fun <- read_csv_while_replacing_sep_in_header
  } else {
    read_file_fun <- read.csv
  }

  pherograms <- purrr::map(stats::setNames(subfiles, basename(subfiles)), read_file_fun, check.names = F)
  pherograms <- purrr::map_dfr(pherograms, function(x) {
    names(x) <- sapply(strsplit(names(x), ": "), "[", 2)
    names(x)[1] <- "size_nt"
    if (!is.null(filter_names) && length(filter_names) > 1) {
      x <- x[,!apply(sapply(filter_names, grepl, x = names(x), simplify = T, ignore.case = T), 1, any)]
    }
    x <- tidyr::pivot_longer(x, cols = -size_nt, names_to = "sample_ID", values_to = "signal")
    return(x)
  }, .id = "file")

  pherograms <- pherograms[,c("sample_ID", "size_nt", "signal", "file")]
  return(pherograms)
}

read_qualities <- function(subfiles,
                           filter_names = NULL,
                           sample_name_replace = NULL) {

  if (!is.null(sample_name_replace)) {
    qualities <- purrr::map_dfr(stats::setNames(subfiles, basename(subfiles)), read_quality_csv_with_sample_name_replacement, check.names = F, sample_name_replace = sample_name_replace, .id = "file")
  } else {
    qualities <- purrr::map_dfr(stats::setNames(subfiles, basename(subfiles)), read.csv, check.names = F, .id = "file")
  }

  unit <- strsplit(grep("Conc", names(qualities), value = T), " ")[[1]][2]
  unit <- gsub("[\\(\\)]", "", unit)
  if (length(unit) == 0) {
    unit <- "unit not recognized"
  }
  names(qualities) <- c("file", "well", "sample_ID", "conc", "RQN", "R28S_18S")
  if (!is.null(filter_names) && length(filter_names) > 1) {
    qualities <- qualities[!apply(sapply(filter_names, grepl, x = qualities[,"sample_ID",drop=T], simplify = T, ignore.case = T), 1, any),]
  }

  qualities$unit <- unit
  qualities <- qualities[,c("well", "sample_ID", "conc", "unit", "RQN", "R28S_18S", "file")]
  qualities$RQN2 <- paste0("RQN = ", qualities$RQN)
  qualities$conc2 <- paste0(qualities$conc, " ", qualities$unit)
  qualities$R28S_18S2 <- paste0("28S/18S = ", qualities$R28S_18S)
  return(qualities)
}

read_peaks <- function(subfiles,
                       filter_names = NULL,
                       sample_name_replace = NULL) {

  if (!is.null(sample_name_replace)) {
    peaks_list <- purrr::map(stats::setNames(subfiles, basename(subfiles)), read_peak_csv_with_sample_name_replacement, check.names = F, header = F, sample_name_replace = sample_name_replace)
  } else {
    peaks_list <- purrr::map(stats::setNames(subfiles, basename(subfiles)), read.csv, check.names = F, header = F)
  }

  ## loop over files
  peaks_list2 <- purrr::map(peaks_list, function(x) {
    #splits <- which(diff(trimws(x[,1,drop=T]) == "") == -1)
    splits <- which(trimws(x[,2,drop=T]) == "Total Conc.:") + 1
    splits <- splits[-length(splits)]
    split_rows <- seq2(c(1,splits+1), c(splits, nrow(x)))
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
      } else if (length(split_rows2) == 1) {
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
      } else {
        message("Peak file reading: More than one split_rows2 found. Skipping.")
        return(list(peak_table = NULL, summary_table = NULL))
      }
    })
    if (!is.null(filter_names) && length(filter_names) > 1) {
      sample_tables <- sample_tables[!apply(sapply(filter_names, grepl, x = names(sample_tables), simplify = T, ignore.case = T), 1, any)]
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

str_count2 <- function(string, pattern = ",") {
  nchar(string) - nchar(gsub(pattern, "", string))
}

check_header_for_commas <- function(filepath, sep = ",") {
  lines <- vroom::vroom_lines(filepath, n_max = 2)
  str_count2(lines[1], sep) != str_count2(lines[2], sep)
}

read_csv_while_replacing_sep_in_header <- function(filepath, sep = ",", check.names = F) {

  lines <- vroom::vroom_lines(filepath)
  # replace every second separator
  sep_positions_in_header <- gregexpr(",", lines[1])[[1]]
  sep_to_rm <- sep_positions_in_header[seq(2,length(sep_positions_in_header), 2)]
  sep_for_split <- sep_positions_in_header[seq(1,length(sep_positions_in_header), 2)]
  temp <- strsplit(lines[1], "")[[1]]
  lines[1] <- paste(temp[-sep_to_rm], collapse = "")

  # this is for read peaks below

  sample_name_replace <- stats::setNames(gsub("^[[:alpha:]][[:digit:]]{1,2}: ", "", strsplit(lines[1], ",")[[1]]),
                                         nm = gsub("^[[:alpha:]][[:digit:]]{1,2}: ", "", split_string_at_indices(input_string = paste(temp, collapse = ""), indices = sep_for_split)))
  assign("sample_name_replace", sample_name_replace, envir = parent.frame(4)) # walk 4 environment upwards to set variable there, which is the read_fragmentanalyzer function environment

  if (str_count2(lines[1], sep) == str_count2(lines[2], sep)) {
    lines <- strsplit(lines, ",")
    df <- as.data.frame(do.call(rbind, lines[-1]))
    names(df) <- lines[[1]]
    if (check.names) {
      names(df) <- make.names(names(df))
    }
    df[] <- lapply(df, as.numeric)
    return(df)
  } else {
    stop("Replacing every second separator did not work. Try to fix manually.")
  }
}

read_quality_csv_with_sample_name_replacement <- function(filepath, sample_name_replace,
                                                          check.names = F) {
  lines <- vroom::vroom_lines(filepath)
  for (i in seq_along(lines)) {
    for (j in seq_along(sample_name_replace)) {
      lines[i] <- gsub(names(sample_name_replace)[j], sample_name_replace[j], lines[i])
    }
  }

  lines <- strsplit(lines, ",")
  df <- as.data.frame(do.call(rbind, lines[-1]))
  names(df) <- lines[[1]]
  df[3:5]  <- lapply(df[3:5], as.numeric)
  return(df)
}

read_peak_csv_with_sample_name_replacement <- function(filepath, sample_name_replace,
                                                       check.names = F, header = F) {


  lines <- vroom::vroom_lines(filepath)
  for (i in seq_along(lines)) {
    for (j in seq_along(sample_name_replace)) {
      lines[i] <- gsub(names(sample_name_replace)[j], sample_name_replace[j], lines[i])
    }
  }

  lines <- strsplit(lines, ",")
  for (i in which(lengths(lines) < max(lengths(lines)))) {
    lines[[i]] <- c(lines[[i]], "")
  }
  #unname(sapply(lines, str_count2))
  #lines2[which(lengths(lines2) < max(lengths(lines2)))]

  df <- as.data.frame(do.call(rbind, lines))
  return(df)
}

split_string_at_indices <- function(input_string, indices) {
  # this removes the character at indices positions
  # this could be changed by an additional fun argument
  return_string <- character(length = length(indices))
  indices <- c(0, indices, nchar(input_string)+1) # if clause to check for 1 at start?
  for (i in seq_along(indices[-1])) {
    return_string[i] <- substr(input_string, indices[i]+1, indices[i+1]-1)
  }

  return(return_string)
}


seq2_default <- Vectorize(seq.default, vectorize.args = c("from", "to"))
seq2 <- function(from = 1, to = 1) {
  x <- seq2_default(from = from, to = to)
  # make sure always a list is returned
  if (is.matrix(x)) {
    x <- unname(as.list(as.data.frame(x)))
  }
  return(x)
}

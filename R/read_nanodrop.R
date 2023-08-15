#' Read reports and spectra files from NanoDrop into tidy tables
#'
#' The tsv files from Nanodrop are not very handy in R. Several conversions
#' are necessary (column names, character to numeric etc.) which is time consuming.
#' Also, in order to plot spectra the respective table from Nandrop requires quite
#' some wrangling. This function tries to solve these tasks elegantly.
#'
#' @param reports vector of paths to report files (tsv)
#' @param spectra vector of paths to spectra files (tsv)
#'
#' @return a list of length 2; index one is report only, index two is report
#' with spectra values added
#' @export
#'
#' @examples
read_nanodrop <- function(reports, spectra) {

  # to do: DNA and protein report from Nanodrop

  repo <- do.call(rbind,lapply(reports, function(x) {
    re <- utils::read.csv(x, sep = "\t", header = T, row.names = 1, check.names = F)
    re$Unit <- stringr::str_replace(re$Unit, "\xb5", "u")
    #re$Unit <- gsub("\\\\", "", re$Unit)
    #re$Unit <- gsub("xb5", "u", re$Unit)
    names(re)[which(names(re) == "260/280")] <- "R260_280"
    names(re)[which(names(re) == "260/230")] <- "R260_230"
    names(re)[which(names(re) == "Nucleic Acid Conc.")] <- "Conc"
    names(re) <- gsub(" ", "_", names(re))
    re[,"Date"] <- sapply(strsplit(re[,"Date_and_Time"], " "), "[", 1)
    re[,"Time"] <- sapply(strsplit(re[,"Date_and_Time"], " "), "[", 2)
    re <- re[,-which(names(re) == "Date_and_Time")]

    to_num <- c("Conc", "A260", "A280", "R260_280", "R260_230", "Factor")
    re[to_num] <- sapply(re[to_num], function(y) {
      as.numeric(gsub(",",".",y))
    })
    re[,"sheet"] <- basename(x)
    re[,"Sample_ID"] <- as.character(re[,"Sample_ID"])
    return(re)
  }))

  spec <- do.call(rbind, lapply(spectra, function (x) {
    sp <- utils::read.csv(x, sep = "\t", header = F, row.names = NULL, blank.lines.skip = F)
    lines <- which(sp[,1] == "Wavelength (nm)")
    lines <- c(lines,nrow(sp))
    wl_df <- do.call(rbind, lapply(rev(rev(seq_along(lines))[-1]), function(z) {
      wl <- sp[(lines[z]+1):(lines[z+1]-5),]
      names(wl) <- sp[lines[z],c(1,2)]
      names(wl) <- gsub(" \\(nm\\)", "_nm", names(wl))
      wl[names(wl)] <- sapply(wl[names(wl)], function(y) {
        as.numeric(gsub(",",".",y))
      })
      wl[,"Date"] <- strsplit(sp[(lines[z]-1),1], " ")[[1]][1]
      wl[,"Time"] <- strsplit(sp[(lines[z]-1),1], " ")[[1]][2]
      wl[,"Sample_ID"] <- as.character(sp[(lines[z]-2),1])
      return(wl)
    }))
    return(wl_df)
  }))

  repo_spec <- merge(repo, spec)
  if (nrow(repo_spec) > nrow(spec)) {
    warning("Merging reports and spectra by Date, Time and Sample_ID was not unambiougsly.")
  }
  return(list(report = repo, spectra = repo_spec))
}


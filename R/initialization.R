# Copyright 2020-2022 Bundesanstalt für Gewässerkunde
# This file is part of shewhart4hrms
# shewhart4hrms is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any 
# later version.
# 
# shewhart4hrms is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along 
# with shewhart4hrms. If not, see <https://www.gnu.org/licenses/>.



  
#' Create directory for shewhart4hrms
#' 
#' Creates a new shewhart4hrms directory tree and example files
#'
#' @param newDir Path to new directory
#'
#' @return TRUE if successful
#' @export
#'
#' @examples shewhart4hrms::create_dir("~/linsensuppe")
create_dir <- function(path) {
  if (!dir.exists(path))
    dir.create(path) else stop("You must create a new directory")
  
  dir.create(file.path(path, "mzXML-data-files"))
  dir.create(file.path(path, "results"))
  dir.create(file.path(path, "settings"))
  
  # in the settings folder, place IS-Table examples
  file.copy(
    system.file("example-settings", "is-table-pos.csv", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  file.copy(
    system.file("example-settings", "is-table-neg.csv", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  # place an example settings.yml file
  
  file.copy(
    system.file("example-settings", "settings.yml", package = "shewhart4hrms"),
    file.path(path, "settings")
  )
  
  # place bat file in directory for starting the app
  
  ptr <- normalizePath(file.path(R.home(), "bin", "R.exe"))
  pth <- normalizePath(path, winslash = "/")
  batpth <- file.path(path, sprintf("shewhart4hrms-%s.bat", basename(path)))
  cat(sprintf("\"%s\" -e \"shewhart4hrms::viewshewhart('%s')\"\npause", ptr, pth), 
      file = batpth)
  
  # place readme file in the directory
  readme <- sprintf("
    ---- Using shewhart4hrms ----
    Step 1: Convert your data into mzXML using Proteowizard
    Step 2: Save data into the %s/mzXML-data-files directory
    Step 3: Open %s
    Step 4: Click 'Refresh' to process the new files
    ", pth, basename(batpth)
  )
  
  cat(readme, file = file.path(path, "README.txt"))
  
  message("Edit the IS tables in ", file.path(path, "settings"))
  message("Copy mzXML files to ", file.path(path, "mzXML-data-files"))
  
  TRUE
}


#' Initiate processing of mzXML files
#'
#' @param path Path to shewhart4hrms directory
#'
#' @return TRUE if it runs to the end
#' @export
#'
#' @examples shewhart4hrms::initiate("~/linsensuppe")
initiate <- function(path) {
  # path <- "~/shewhart-test"
  
  # make some checks ####
  if (!all(c("mzXML-data-files", "results", "settings") %in% list.dirs(path, full.names = F)))
    stop("Incorrect directory in ", path)
  
  # check for IS tables
  
  if (!file.exists(file.path(path, "settings", "is-table-pos.csv")))
    stop("No is-table-pos file found in ", file.path(path, "settings"))
  if (!file.exists(file.path(path, "settings", "is-table-neg.csv")))
    stop("No is-table-neg file found in ", file.path(path, "settings"))
  
  # check for settings
  
  if (!file.exists(file.path(path, "settings", "settings.yml")))
    stop("No settings file found in ", file.path(path, "settings"))
  sings <- yaml::read_yaml(file.path(path, "settings", "settings.yml"))
  
  # TODO check for mzXML data with pos/neg in name
  
  # define function
  
  initiate_pol <- function(pol) {  # pol <- "pos"
    message("Initiating ", pol)
    
    shewhart <- Report$new()
    
    details <- file.info(get(paste0(pol, "Files")))
    first_file <- rownames(details[with(details, order(as.POSIXct(mtime))), ])[1]
    
    shewhart$addRawFiles(F, first_file)
    
    shewhart$addIS(F, file.path(path, "settings", sprintf("is-table-%s.csv", pol)))
    shewhart$changeSettings("area_threshold", sings[[paste0("area_threshold_", pol)]])
    shewhart$changeSettings("use_int_threshold","area")
    shewhart$changeSettings("EIC_extraction", sings[[paste0("eic_extraction_width_", pol)]])
    shewhart$changeSettings("ISmztol", sings[[paste0("mz_tol_", pol)]])
    shewhart$changeSettings("ISrttolm", sings[[paste0("rt_tol_", pol)]])
    shewhart$changeSettings("pol", pol)
    
    shewhart$process_all()
    
    #shewhart$view()
    
    # check that IS were found 
    if (nrow(shewhart$ISresults) == 0)
      stop("No IS found in ", pol, " mode in file ", first_file)
    
    shewhart$clearAndSave(F, file.path(path, "results", sprintf("shewhart-%s.Report", pol)))
    
    # create a table with relevant IS information
    
    tab <- shewhart$ISresults
    tab$mtime <- as.integer(file.mtime(first_file))
    
    write.csv(tab, file = file.path(path, "results", paste0(pol, "_res.csv")), row.names = F)
    
    #prepare pos_res_graph for visualisation
    is_pol <- read.csv(file.path(path, "settings", sprintf("is-table-%s.csv", pol)))
    is_pol$mass <- mapply(
      shewhart4hrms::get_mass, 
      formula = is_pol$formula, 
      adduct = is_pol$adduct, 
      MoreArgs = list(charge = switch(pol, pos = 1, neg = -1))
    )
    is_pol <- is_pol[, c("name", "mass", "rt")]
    colnames(is_pol) <- c("IS", "mz_calc", "rt_known")
    tab <- merge(tab, is_pol)
    tab$delta_mz_mDa <- abs(tab$mz - tab$mz_calc) * 1000 
    tab$delta_rt_min <- abs(tab$rt - tab$rt_known)
    
    # compute peak width
    tab$peak_width_min <- tab$peak_end - tab$peak_start
    
    tab$time <- as.POSIXct(tab$mtime, origin = "1970-01-01 00:00.00 UTC", tz = "Europe/Berlin")
    write.csv(tab, file = file.path(path, "results", paste0(pol,"_res_graph.csv")), row.names = F)
    message("Completed initialization ", pol)
  }
  
  
  # initiate pos and neg ####
  posFiles <- list.files(
    file.path(path, "mzXML-data-files"), 
    full.names = TRUE, 
    pattern=".*pos.*\\.mzXML$"
  )
  
  negFiles <- list.files(
    file.path(path, "mzXML-data-files"), 
    full.names = TRUE, 
    pattern=".*neg.*\\.mzXML$"
  )
  
  if (length(posFiles) == 0 && length(negFiles) == 0)
    stop("No mzXML files with 'pos' or 'neg' in name found in ", file.path(path, "mzXML-data-files"))
  
  
  if (length(posFiles) > 0) {
    initiate_pol("pos")
  }
  
  if (length(negFiles) > 0) {
    initiate_pol("neg")
  }
  
  TRUE
}

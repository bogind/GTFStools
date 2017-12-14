#' @export downloadGTFSil
downloadGTFSil <- function(directory){
  directory <- choose.dir()
  date <- Sys.Date()
  date <- as.character(as.POSIXct(strptime(date, format = "%Y-%m-%d")))
  file_name <- as.character(paste("GTFS",date,".zip", sep = ""))
  download.file("ftp://gtfs.mot.gov.il/israel-public-transportation.zip", destfile = paste0(directory,file_name), method = "libcurl")}

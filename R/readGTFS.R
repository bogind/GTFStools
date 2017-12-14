#' @importFrom utils installed.packages
#' @importFrom data.table fread
#' @importFrom easycsv Identify.OS choose_dir
#' @export
readGTFS = function(directory = NULL,
                        extension = "BOTH",
                        sep="auto",
                        nrows=-1L,
                        header="auto",
                        na.strings="NA",
                        stringsAsFactors=FALSE,
                        verbose=getOption("datatable.verbose"),
                        autostart=1L,
                        skip=0L,
                        drop=NULL,
                        colClasses=NULL,
                        integer64=getOption("datatable.integer64"),
                        dec=if (sep!=".") "." else ",",
                        check.names=FALSE,
                        encoding="unknown",
                        quote="\"",
                        strip.white=TRUE,
                        fill=FALSE,
                        blank.lines.skip=FALSE,
                        key=NULL,
                        prefix="GTFS",
                        minimal=FALSE,
                        showProgress=getOption("datatable.showProgress"),
                        data.table=FALSE
){
  if ("data.table" %in% rownames(installed.packages()) == FALSE) {
    stop("data.table needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(is.null(directory)){
    os = Identify.OS()
    if(tolower(os) == "windows"){
      directory <- utils::choose.dir()
      if(tolower(os) == "linux" | tolower(os) == "macosx"){
        directory <- choose_dir()
      }
    }else{
      stop("Please supply a valid local directory")
    }

  }

  directory = paste(gsub(pattern = "\\", "/", directory,
                         fixed = TRUE))



  endings = list()

  if(tolower(extension) == "txt"){
    endings[1] =  "*\\.txt$"
  }
  if(tolower(extension) == "csv"){
    endings[1] =  "*\\.csv$"

  }
  if(tolower(extension) == "both"){
    endings[1] =  "*\\.txt$"
    endings[2] =  "*\\.csv$"
  }
  if((tolower(extension) %in% c("txt","csv","both")) == FALSE){
    stop("Pleas supply a valid value for 'extension',\n
         allowed values are: 'TXT','CSV','BOTH'.")
  }
  tempfiles = list()
  temppath = list()
  num = 1
  for(i in endings){
    temppath = paste(directory,list.files(path = directory, pattern=i), sep = "/")
    tempfiles = list.files(path = directory, pattern=i)
    num = num +1
    if(length(temppath) < 1 | length(tempfiles) < 1){
      num = num+1
    }
    else{
      temppath = unlist(temppath)
      tempfiles = unlist(tempfiles)
      count = 0

      minvec = c("agency","routes","stop_times","stops","trips","calendar")

      g = 0
      f = vector()
      if(minimal == TRUE){
        for(file in tempfiles){
          filename = paste0(gsub(i, "", file))
          inside = length(grep(filename,minvec, ignore.case = TRUE)) > 0

          g = g+1
          f[g] = inside
        }
      }else{
        f = rep(TRUE,length(temppath))
      }

      temppath = temppath[f]
      for(t in 1:length(temppath)){
        count = count+1
        tbl = temppath[count]
        DTname1 = paste0(gsub(directory, "", tbl))
        DTname2 = paste0(gsub("/", "", DTname1))
        DTname3 = paste0(gsub(i, "", DTname2))

        if(!is.null(prefix) && is.character(prefix)){
          DTname4 = paste(prefix,DTname3, sep = "")
        }else{
          print("Warning: without the `GTFS` prefix you will need to name the tables manually")
          DTname4 = DTname3
        }

        DTable <- data.table::fread(input = tbl,
                                    sep=sep,
                                    nrows=nrows,
                                    header=header,
                                    na.strings=na.strings,
                                    stringsAsFactors=stringsAsFactors,
                                    verbose = verbose,
                                    autostart=autostart,
                                    skip=skip,
                                    drop=drop,
                                    colClasses=colClasses,
                                    dec=if (sep!=".") "." else ",",
                                    check.names=check.names,
                                    encoding=encoding,
                                    quote=quote,
                                    strip.white=strip.white,
                                    fill=fill,
                                    blank.lines.skip=blank.lines.skip,
                                    key=key,
                                    showProgress=getOption("datatable.showProgress"),
                                    data.table=data.table
        )
        assign_to_global <- function(pos=1){
          assign(x = DTname4,value = DTable, envir=as.environment(pos) )
        }
        assign_to_global()

        rm(DTable)
      }
    }
  }
  print("Finished reading all GTFS tables into Global environment")
  }





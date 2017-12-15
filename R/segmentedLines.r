#' @importFrom SIRItoGTFS stopstoSP
#' @importFrom sp CRS spTransform SpatialPointsDataFrame SpatialLines
#' @importFrom plyr join
#' @importFrom rgeos gDistance
#' @importFrom rgdal make_EPSG
#' @export

segmentedLines = function(route,
                          projcode,
                          shapes = GTFSshapes,
                          routes = GTFSroutes,
                          trips = GTFStrips,
                          stop_times = GTFSstop_times,
                          stops = GTFSstops){


  # Also uses rgeos and sp, but SIRItoGTFS is dependent on them

  proj = make_EPSG()
  crswgs84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  crs1 <- CRS(proj$prj4[projcode])


  # Subset and order the shape and the stops within it.

  rt = GTFSroutes[GTFSroutes$route_id %in% route,]
  trip = GTFStrips[GTFStrips$route_id %in% rt$route_id,]
  shp = GTFSshapes[GTFSshapes$shape_id %in% trips$shape_id,]
  shp = shp[order(shp$shape_pt_sequence),]
  st = GTFSstop_times[GTFSstop_times$trip_id %in% trip$trip_id,]
  stops1 = GTFSstops[GTFSstops$stop_id %in% st$stop_id,]

  # Join for the stop sequence, needed for each route/trip
  stops1 = join(stops1,st,by = "stop_id", type = "left", match = "first")
  stops1 = stops1[order(stops1$stop_sequence),]

  # convert stops to ILTM and visualize them
  stops2 = stopstoSP(stops1,4326,useSIRI = FALSE)
  stops2 = spTransform(stops2,crs1)

  # create an spdf from the shapes
  startstop = st$stop_id[st$stop_sequence == 1]
  endstop = st$stop_id[st$stop_sequence == max(st$stop_sequence)]
  shp2 = SpatialPointsDataFrame(coords = shp[,3:2], data = shp, proj4string = crswgs84)
  shp3 = spTransform(shp2,crs1)

  distable = gDistance(stops2, shp3,byid = TRUE)

  shrtdst = apply(distable,2,which.min)
  shrtdst[1] = ifelse(1 %in% shp$shape_pt_sequence, 1, shrtdst[1])

  # Subset points for the new lines
  seqslist = list()
  seqdf = data.frame()
  for(i in 1:(length(shrtdst)-1)){
    if(i < (length(shrtdst)-2) ){
      shp4 = shp3[shp3$shape_pt_sequence >= shrtdst[i] & shp3$shape_pt_sequence < shrtdst[i+1],]
      shp4@data$sequence_number = paste0(shp4@data$shape_id,"_",i)
      seqslist[[i]] = shp4
      seqdf = rbind(seqdf,shp4@data)
    }
    else{
      shp4 = shp3[shp3$shape_pt_sequence >= shrtdst[i] & shp3$shape_pt_sequence <= max(shp3$shape_pt_sequence),]
      shp4@data$sequence_number = paste0(shp4@data$shape_id,"_",i)
      seqslist[[i]] = shp4
      seqdf = rbind(seqdf,shp4@data)
    }

  }

  # Create a list of Lines s4 objects and convert them to SpatialLines
  lineslist = lapply(seqslist,FUN = function(i){Lines(Line(i@coords), ID=i@data$sequence_number[1])})
  sl = SpatialLines(lineslist, proj4string = crs1)

  return(sl)

}

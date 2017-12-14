#' @export GTFSbyshp
GTFSbyshp <- function(geom = NULL,
                      proj4string = NULL,
                      stops = GTFSstops,
                      routes = GTFSroutes,
                      trips = GTFStrips,
                      stop.times = GTFSstop_times,
                      calender = GTFScalendar){
  require(sp)
  require(rgdal)
  require(rgeos)
  if(class(geom) == "SpatialPointsDataFrame" | class(geom) == "SpatialPoints"){
    geombox = gConvexHull(geom)
  }
    if(class(geom) == "SpatialPolygonsDataFrame" | class(geom) == "SpatialPolygons"){
    geombox = geom
    }
    else{stop("'geom' must be a SpatialPolygonsDataFrame,SpatialPointsDataFrame,SpatialPolygons or SpatialPolygons")}
  proj = rgdal::make_EPSG()
  wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  crs1 <- sp::CRS(wgs84)
  if(geombox@proj4string@projargs != crs1@projargs){
    geombox = spTransform(geombox, crs1)
  }
  s <- sp::SpatialPointsDataFrame(coords = data.frame(GTFSstops$stop_lon,GTFSstops$stop_lat),data = GTFSstops,proj4string = crs1)
  s1 <- s@data[gContains(geombox,s),]
  GTFSstops <<- GTFSstops[GTFSstops$stop_id %in% s1$stop_id,]
  GTFSstop_times <<- GTFSstop_times[GTFSstop_times$stop_id %in% GTFSstops$stop_id,]
  GTFStrips <<- GTFStrips[GTFStrips$trip_id %in% GTFSstop_times$trip_id,]
  GTFSroutes <<- GTFSroutes[GTFSroutes$route_id %in% GTFStrips$route_id,]
  GTFScalendar <<- GTFScalendar[GTFScalendar$service_id %in% GTFStrips$service_id,]
  GTFSshapes <<- GTFSshapes[GTFSshapes$shape_id %in% GTFStrips$shape_id,]


}

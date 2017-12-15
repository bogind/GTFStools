# GTFStools
current version 0.0.1

An R package that provieds a set of tools for using and manipulating GTFS tables.

current functions:    
readGTFS - reads multiple tables from csv or txt files and assignes them to the Global environment with a *GTFS* prefix.   
segmentedLines - creates SpatialLines objects from a selected route, the line is segmented between each of the route's stops.

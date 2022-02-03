## -----------------------------------------------------------------
##
## Script name: funcs.R
##
## Purpose of script: Usefull functions for handeling EO Data 
##
## Author: Moritz Roesch
##
## Date: 2020-04-09
##
## -----------------------------------------------------------------


##---- Clip Raster to Shape ----

# Clips a Raster* object to an SpatialPolygons*, SpatialLines*, or SpatialPoints* object

clip_raster2shp<- function(raster, shape, inv = FALSE){
  
  library(raster) #load package dependencies
  library(sf)
  
  crop1<- crop(raster, shape) # crop raster to extent of shape
  ra_shp<- rasterize(shape, crop1) # rasterize shape to raster extent
  if (inv == TRUE){
    mask(crop1, ra_shp, inverse = TRUE)
  } else {
    mask(crop1, ra_shp) # mask/clip raster extent to rasterized shape
  }
}





##---- Mosaic Planetscope ----

mosaic_from_list <- function(filelist, func){
  # Mosaics all raster files in filelist with defined function using the raster::mosaic function
  # filelist:   list, list of raster files
  # func:       string, function for overlapping values 
  
  require(raster)
  
  names(filelist) <- NULL # resets list names
  filelist$fun <- func # adds mosaic function to list
  mos <- do.call(mosaic, filelist) # mosaics all rasters in list with defined func
  return(mos)
}





##---- NDVI Planetscope ----

ndvi <- function(raster, nir = 4, red = 3, add = TRUE, write_raster = FALSE, filename){
  
  require(raster)
  
  ndvi <- (raster[[nir]] - raster[[red]]) / (raster[[nir]] + raster[[red]]) # calculates NDVI based on selected bands
  names(ndvi) <- "NDVI" # renames layer
  if (add == TRUE){
    raster <- stack(raster, ndvi) # stacks new NDVI layer onto raster scene
    return(raster)
  } else {
    return(ndvi)
  }
  if (write_raster == TRUE){
    if (add == TRUE){
      writeRaster(raster, filename = filename, format = "GTiff") # writes stacked PS scene
    } else {
      writeRaster(ndvi, filename = filename, format = "GTiff") # write only NDVI band
    }
  }
}





##---- NDVI ----

ndvi <- function(raster, nir, red, add = TRUE, write_raster = FALSE, filename){
  
  require(raster)
  
  ndvi <- (raster[[nir]] - raster[[red]]) / (raster[[nir]] + raster[[red]]) # calculates NDVI based on Planetscope bands
  names(ndvi) <- "ndvi" # renames layer
  if (add == TRUE){
    raster <- stack(raster, ndvi) # stacks new NDVI layer onto Planetscope scene
    return(raster)
  } else {
    return(ndvi)
  }
  if (write_raster == TRUE){
    if (add == TRUE){
      writeRaster(raster, filename = filename, format = "GTiff") # writes stacked PS scene
    } else {
      writeRaster(ndvi, filename = filename, format = "GTiff") # write only NDVI band
    }
  }
}



##---- Fit extent of similar raster objects ----

fit_extent <- function(ra_list, ext, overwrite = FALSE){
  # Fits rasters with slightly different extents (caused by resampling or reprojecting the data) to a given extent based on their intersection with the extent object (ext)
  # ra_list:  List of rasters with slightly different extents. The rasters in the list are going to be aligned by their intersection.
  # ext:      Extent object of the extent which alle raster should be aligned to.
  
  require(raster)
  
  
  output_list <- list()
  for (i in 1:length(ra_list)){
    if (ext == extent(ra_list[[i]])){
      print(paste("Extent of", names(ra_list)[i], "is equal, no need to crop"), sep = " ")
      if (overwrite == FALSE){
        output_list[i] <- ra_list[[i]] # stores unchanged raster new list
      } else {
        ra_list[i] <- ra_list[[i]] #stores unchanged raster in input list
      }
      
    } else {
      # calculate overlap between the two datasets
      overlap <- intersect(ext, extent(ra_list[[i]]))
      # now let's crop both datasets to the overlap region
      output <- crop(ra_list[[i]], overlap)
      print(paste("Extent of", names(ra_list)[i], "is different, cropping data"), sep = " ")
      if (overwrite == FALSE){
        output_list[i] <- output # stores new raster new list
      } else {
        ra_list[i] <- output #stores new raster in input list
      }
    }
  }
  if (overwrite == FALSE){
    return(output_list)
  } else{
    return(ra_list)
  }
}


#---- Process Sentinel-2 data from raw datafolder ----
preproc_s2 <- function(path_to_dir, dir_name, resolution, add_aux = FALSE,
                       stack_order = NULL, crs = NULL,
                       clip_to_aoi = FALSE, aoi, inverse_clip = FALSE,
                       write_raster = FALSE, filename){
  # Processes Sentinel-2 Level 2 data from raw data folder to a stack. Selection of bands
  # to stack and resolution is possible. Additionally, reprojecting the datastack
  # to another CRS. To reduce dataset to an AOI clipping of datastack is possible.
  
  # path_to_dir:    string,
  #                 Path to Sentinel-2 data folder
  # dir_name:       string,
  #                 Name of Sentinel-2 data folder
  # resolution:     int,
  #                 Select resolution (10, 20 or 60) of Sentinel-2 bands. For each
  #                 resolution all bands in this resolution are included. If 20 is
  #                 selected, all bands with 10m resolution are downsampled to 20.
  #                 If 60m is selected, all bands with 10 or 20m are downsampled.
  # add_aux:        bool, default = FALSE,
  #                 If TRUE, all auxillary datasets of Sentinel (e.g. scene class-
  #                 ification layer (SCL), etc.) are added to stack. See documentation
  #                 for all available aux-datasets.
  # stack_order:    string/vector, default = NULL,
  #                 Select bands and order them in stack. Bands must be provided in
  #                 Sentinel-2 naming convention (e.g. B02). Define order by defining
  #                 a vector with desired bands (e.g. c("B02", "B03", "B04", "B8A"))
  # crs:            string, default = NULL,
  #                 To reproject the created raster stack, provide CRS details in
  #                 PROJ.4 format.
  # clip_to_aoi:    bool, default = FALSE,
  #                 If TRUE, layer stack will be clipped to area of `aoi`
  # aoi:            sf object/geodataframe
  #                 Sf object or geodataframe with polygons which are used for clip.
  # inverse_clip:   bool, default = FALSE,
  #                 If TRUE, inverse clipping will be done (areas outside of polygons)
  # write_raster:   bool, default = FALSE,
  #                 If TRUE, writeRaster function writes raster to directory,
  # filename:       string
  #                 Define path and filename for raster layer which is used in
  #                 writeRaster function
  
  
  require(raster)
  require(sf)
  require(stringr)
  
  gran_dir <- list.files(str_c(path_to_dir, dir_name, "/", dir_name, ".SAFE/GRANULE")) # name of folder in GRANULE
  path_to_bands <- str_c(path_to_dir, dir_name, "/", dir_name, ".SAFE/GRANULE/",
                         gran_dir, "/IMG_DATA/", "R", as.character(resolution), "m") # path to band files in defined resolution
  fl <- list.files(path_to_bands, full.names = TRUE) #file list
  
  # select only satelite bands or auxillary bands
  if (add_aux == FALSE){
    fl_bands <- fl[str_detect(fl, pattern = "B\\d{2}|B8A")] # select only raw data bands
  } else {
    fl_bands <- fl
  }
  
  # load data
  data <- sapply(fl_bands, raster)
  
  # rename bands in list
  for (i in 1:length(data)){
    name <- str_extract(names(data[[i]]), pattern = "B\\d{2}|B8A|WVP|TCI|AOT|SCL") # extract useful band name
    names(data[[i]]) <- name
    names(data)[i] <- name
  }
  
  # create layer stack
  # stack bands
  
  if (is.null(stack_order)){ # if no layer order is provided
    datastack <- stack(data)
  } else {
    data_ordered <- list()
    for (i in stack_order){
      data_ordered[i] <- data[i] # stores layers in order 
    }
    datastack <- stack(data_ordered)
  }
  
  # reproject raster
  if (!is.null(crs)){ # if crs is defined do the reprojection
    datastack <- projectRaster(datastack, crs = crs)
  }
  
  # clip to vector
  if (clip_to_aoi == TRUE){
    crop1<- crop(datastack, aoi) # crop raster to extent of shape
    ra_shp<- rasterize(aoi, crop1) # rasterize shape to raster extent
    if (inverse_clip == TRUE){
      datastack <- mask(crop1, ra_shp, inverse = TRUE)
    } else {
      datastack <- mask(crop1, ra_shp) # mask/clip raster extent to rasterized shape
    }
  }
  
  # write raster
  if (write_raster == TRUE){
    writeRaster(datastack, filename = filename, format = "GTiff", overwrite = TRUE)
  }
  return(datastack)
}





#---- Process Landsat 8 data from raw datafolder ----
preproc_l8 <- function(path_to_dir, dir_name, add_aux = FALSE,
                       stack_order = NULL, crs = NULL,
                       clip_to_aoi = FALSE, aoi, inverse_clip = FALSE,
                       write_raster = FALSE, filename){
  # Processes Landsat 8 Level 2 data from raw data folder to a stack. Define selection of
  # desired bands for stacking is possible. Additionally, reprojecting the datastack
  # to another CRS. To reduce dataset to an AOI clipping of datastack is possible.
  
  # path_to_dir:    string,
  #                 Path to Landsat 8 Level 2 data folder
  # dir_name:       string,
  #                 Name of Landsat 8 Level 2 data folder
  # add_aux:        bool, default = FALSE,
  #                 If TRUE, all auxillary datasets of Landsat (e.g. pixel quality,
  #                 radiometric saturation, areosol) are added to stack. See
  #                 documentation for all available aux-datasets.
  # stack_order:    string/vector, default = NULL,
  #                 Select bands and order them in stack. Bands must be provided in
  #                 Landsat 8 naming convention (e.g. B2). Define order by defining
  #                 a vector with desired bands (e.g. c("B2", "B3", "B4", "B7"))
  # crs:            string, default = NULL,
  #                 To reproject the created raster stack, provide CRS details in
  #                 PROJ.4 format.
  # clip_to_aoi:    bool, default = FALSE,
  #                 If TRUE, layer stack will be clipped to area of `aoi`
  # aoi:            sf object/geodataframe
  #                 Sf object or geodataframe with polygons which are used for clip.
  # inverse_clip:   bool, default = FALSE,
  #                 If TRUE, inverse clipping will be done (areas outside of polygons)
  # write_raster:   bool, default = FALSE,
  #                 If TRUE, writeRaster function writes raster to directory,
  # filename:       string
  #                 Define path and filename for raster layer which is used in
  #                 writeRaster function
  
  
  require(raster)
  require(sf)
  require(stringr)
  
  path_to_bands <- str_c(path_to_dir, "/", dir_name)
  fl <- list.files(path_to_bands, full.names = TRUE) #file list
  
  # select only satelite bands or auxillary bands
  if (add_aux == FALSE){
    fl_bands <- fl[str_detect(fl, pattern = "B\\d{1}.TIF")] # select only raw data bands
  } else {
    fl_bands <- fl[str_detect(fl, pattern = "(B\\d{1}|QA_\\w+).TIF")] # adds pixel quality, radiometric and aerosol layer
  }
  
  # load data
  data <- sapply(fl_bands, raster)
  
  # rename bands in list
  for (i in 1:length(data)){
    name <- str_extract(names(data[[i]]), pattern = "B\\d{1}|QA_PIXEL|QA_RADSAT|QA_AEROSOL") # extract useful band name
    names(data[[i]]) <- name
    names(data)[i] <- name
  }
  
  # create layer stack
  # stack bands
  
  if (is.null(stack_order)){ # if no layer order is provided
    datastack <- stack(data)
  } else {
    data_ordered <- list()
    for (i in stack_order){
      data_ordered[i] <- data[i] # stores layers in order 
    }
    datastack <- stack(data_ordered)
  }
  
  # reproject raster
  if (!is.null(crs)){ # if crs is defined do the reprojection
    datastack <- projectRaster(datastack, crs = crs)
  }
  
  # clip to vector
  if (clip_to_aoi == TRUE){
    crop1<- crop(datastack, aoi) # crop raster to extent of shape
    ra_shp<- rasterize(aoi, crop1) # rasterize shape to raster extent
    if (inverse_clip == TRUE){
      datastack <- mask(crop1, ra_shp, inverse = TRUE)
    } else {
      datastack <- mask(crop1, ra_shp) # mask/clip raster extent to rasterized shape
    }
  }
  
  # write raster
  if (write_raster == TRUE){
    writeRaster(datastack, filename = filename, format = "GTiff", overwrite = TRUE)
  }
  return(datastack)
}
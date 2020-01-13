install.packages("oceanmap")
install.packages("ncdf4")
install.packages("raster")
install.packages("viridis")
install.packages("MODISTools")

library(oceanmap)
library(ncdf4)
library(raster)
library(viridis)
library(fishualize)
library(rvest)     
library(stringr)

#############################2010
single2010 = ('Data/sst_data/2010/AQUA_MODIS.20100102.L3m.DAY.SST.sst.4km.nc')
single_nc2010 = nc_open(single2010)
single_raster2010 = nc2raster(single_nc2010, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2010 = raster::crop(single_raster2010, extent(c(48, 58, 22, 32))) 

files2010 = list.files("Data/sst_data/2010/",pattern='*.nc',full.names=TRUE)


for (i in files2010) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2010 <- stack(stacked_rasters2010, raster_i_cropped)
}

max_sst2010 <- stackApply(mystack, indices = rep(1, nlayers(mystack)), fun = max)
writeRaster(max_sst2010, 
            filename="Data/sst_data/max_sst2010.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
min_sst2010 <- stackApply(stacked_rasters2010, indices = rep(1, nlayers(stacked_rasters2010)), fun = min)
writeRaster(min_sst2010, 
            filename="Data/sst_data/min_sst2010.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)

#############################2010
single2011 = ('Data/sst_data/2011/AQUA_MODIS.20110101.L3m.DAY.SST.sst.4km.nc')
single_nc2011 = nc_open(single2011)
single_raster2011 = nc2raster(single_nc2011, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2011 = raster::crop(single_raster2011, extent(c(48, 58, 22, 32))) 

files2011 = list.files("Data/sst_data/2011/",pattern='*.nc',full.names=TRUE)


for (i in files2011) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2011 <- stack(stacked_rasters2011, raster_i_cropped)
}


max_sst2011 <- stackApply(stacked_rasters2011, indices = rep(1, nlayers(stacked_rasters2011)), fun = max)
writeRaster(max_sst2011, 
            filename="Data/sst_data/max_sst2011.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)

min_sst2011 <- stackApply(stacked_rasters2011, indices = rep(1, nlayers(stacked_rasters2011)), fun = min)
writeRaster(min_sst2011, 
            filename="Data/sst_data/min_sst2011.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
#######################################2012

single2012 = ('Data/sst_data/2012/AQUA_MODIS.20120101.L3m.DAY.SST.sst.4km.nc')
single_nc2012 = nc_open(single2012)
single_raster2012 = nc2raster(single_nc2012, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2012 = raster::crop(single_raster2012, extent(c(48, 58, 22, 32))) 

files2012 = list.files("Data/sst_data/2012",pattern='*.nc',full.names=TRUE)


for (i in files2012) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2012 <- stack(stacked_rasters2012, raster_i_cropped)
}



max_sst2012 <- stackApply(stacked_rasters2012, indices = rep(1, nlayers(stacked_rasters2012)), fun = max)
writeRaster(max_sst2012, 
            filename="Data/sst_data/max_sst2012.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
min_sst2012 <- stackApply(stacked_rasters2012, indices = rep(1, nlayers(stacked_rasters2012)), fun = min)
writeRaster(min_sst2012, 
            filename="Data/sst_data/min_sst2012.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)

#############################2013
single2013 = ('Data/sst_data/2013/AQUA_MODIS.20130101.L3m.DAY.SST.sst.4km.nc')
single_nc2013 = nc_open(single2013)
single_raster2013 = nc2raster(single_nc2013, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2013 = raster::crop(single_raster2013, extent(c(48, 58, 22, 32))) 

files2013 = list.files("Data/sst_data/2013/",pattern='*.nc',full.names=TRUE)


for (i in files2013) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2013 <- stack(stacked_rasters2013, raster_i_cropped)
}

max_sst2013 <- stackApply(stacked_rasters2013, indices = rep(1, nlayers(stacked_rasters2013)), fun = max)
writeRaster(max_sst2013, 
            filename="Data/sst_data/max_sst2013.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)

min_sst2013 <- stackApply(stacked_rasters2013, indices = rep(1, nlayers(stacked_rasters2013)), fun = min)
writeRaster(min_sst2013, 
            filename="Data/sst_data/min_sst2013.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)



#############################2014
single2014 = ('Data/sst_data/2014/AQUA_MODIS.20140101.L3m.DAY.SST.sst.4km.nc')
single_nc2014 = nc_open(single2014)
single_raster2014 = nc2raster(single_nc2014, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2014 = raster::crop(single_raster2014, extent(c(48, 58, 22, 32))) 

files2014 = list.files("Data/sst_data/2014/",pattern='*.nc',full.names=TRUE)


for (i in files2014) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2014 <- stack(stacked_rasters2014, raster_i_cropped)
}


max_sst2014 <- stackApply(stacked_rasters2014, indices = rep(1, nlayers(stacked_rasters2014)), fun = max)
writeRaster(max_sst2014, 
            filename="Data/sst_data/max_sst2014.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
min_sst2014 <- stackApply(stacked_rasters2014, indices = rep(1, nlayers(stacked_rasters2014)), fun = min)
writeRaster(min_sst2014, 
            filename="Data/sst_data/min_sst2014.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)



#############################2015
single2015 = ('Data/sst_data/2015/AQUA_MODIS.20150101.L3m.DAY.SST.sst.4km.nc')
single_nc2015 = nc_open(single2015)
single_raster2015 = nc2raster(single_nc2015, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2015 = raster::crop(single_raster2015, extent(c(48, 58, 22, 32))) 

files2015 = list.files("Data/sst_data/2015/",pattern='*.nc',full.names=TRUE)


for (i in files2015) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2015 <- stack(stacked_rasters2015, raster_i_cropped)
}


max_sst2015 <- stackApply(stacked_rasters2015, indices = rep(1, nlayers(stacked_rasters2015)), fun = max)
writeRaster(max_sst2015, 
            filename="Data/sst_data/max_sst2015.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
min_sst2015 <- stackApply(stacked_rasters2015, indices = rep(1, nlayers(stacked_rasters2015)), fun = min)
writeRaster(min_sst2015, 
            filename="Data/sst_data/min_sst2015.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)


#############################2016
single2016 = ('Data/sst_data/2016/AQUA_MODIS.20160101.L3m.DAY.SST.sst.4km.nc')
single_nc2016 = nc_open(single2016)
single_raster2016 = nc2raster(single_nc2016, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2016 = raster::crop(single_raster2016, extent(c(48, 58, 22, 32))) 

files2016 = list.files("Data/sst_data/2016/",pattern='*.nc',full.names=TRUE)


for (i in files2016) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2016 <- stack(stacked_rasters2016, raster_i_cropped)
}

save(stacked_rasters2016,file="Data/sst_data/stacked_rasters2016.RData")

max_sst2016 <- stackApply(stacked_rasters2016, indices = rep(1, nlayers(stacked_rasters2016)), fun = max)
writeRaster(max_sst2016, 
            filename="Data/sst_data/max_sst2016.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
min_sst2016 <- stackApply(stacked_rasters2016, indices = rep(1, nlayers(stacked_rasters2016)), fun = min)
writeRaster(min_sst2016, 
            filename="Data/sst_data/min_sst2016.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)

#############################2017
single2017 = ('Data/sst_data/2017/AQUA_MODIS.20170101.L3m.DAY.SST.sst.4km.nc')
single_nc2017 = nc_open(single2017)
single_raster2017 = nc2raster(single_nc2017, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2017 = raster::crop(single_raster2017, extent(c(48, 58, 22, 32))) 

files2017 = list.files("Data/sst_data/2017/",pattern='*.nc',full.names=TRUE)


for (i in files2017) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2017 <- stack(stacked_rasters2017, raster_i_cropped)
}

max_sst2017 <- stackApply(stacked_rasters2017, indices = rep(1, nlayers(stacked_rasters2017)), fun = max)
writeRaster(max_sst2017, 
            filename="Data/sst_data/max_sst2017.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
min_sst2017 <- stackApply(stacked_rasters2017, indices = rep(1, nlayers(stacked_rasters2017)), fun = min)
writeRaster(min_sst2017, 
            filename="Data/sst_data/min_sst2017.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)

#############################2017
single2018 = ('Data/sst_data/2018/AQUA_MODIS.20180101.L3m.DAY.SST.sst.4km.nc')
single_nc2018 = nc_open(single2018)
single_raster2018 = nc2raster(single_nc2018, "sst", lonname="lon", latname="lat", date=T)
stacked_rasters2018 = raster::crop(single_raster2018, extent(c(48, 58, 22, 32))) 

files2018 = list.files("Data/sst_data/2018/",pattern='*.nc',full.names=TRUE)


for (i in files2018) {
  nc_i <- nc_open(i)
  raster_i = nc2raster(nc_i, "sst", lonname="lon", latname="lat", date=T)
  raster_i_cropped <- raster::crop(raster_i, extent(c(48, 58, 22, 32)))
  stacked_rasters2018 <- stack(stacked_rasters2018, raster_i_cropped)
}


max_sst2018 <- stackApply(stacked_rasters2018, indices = rep(1, nlayers(stacked_rasters2018)), fun = max)
writeRaster(max_sst2018, 
            filename="Data/sst_data/max_sst2018.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)

min_sst2018 <- stackApply(stacked_rasters2018, indices = rep(1, nlayers(stacked_rasters2018)), fun = min)
writeRaster(min_sst2018, 
            filename="Data/sst_data/min_sst2018.tif", 
            options="INTERLEAVE=BAND", overwrite=TRUE)
#######################################
#######################################
#######################################
max2010 = stack("Data/sst_data/max_sst2010.tif")
max2011 = stack("Data/sst_data/max_sst2011.tif")
max2012 = stack("Data/sst_data/max_sst2012.tif")
max2013 = stack("Data/sst_data/max_sst2013.tif")
max2014 = stack("Data/sst_data/max_sst2014.tif")
max2015 = stack("Data/sst_data/max_sst2015.tif")
max2016 = stack("Data/sst_data/max_sst2016.tif")
max2017 = stack("Data/sst_data/max_sst2017.tif")
max2018 = stack("Data/sst_data/max_sst2018.tif")

#get max estimates of all years
max.sst.all <- stack(max2010,
                     max2011,
                     max2012,
                     max2013,
                     max2014, 
                     max2015, 
                     max2016, 
                     max2017, 
                     max2018)


max.max.sst.all <- max(max.sst.all)

vpal = fish(100, alpha = 1, begin = 0.1, end = 1, option = "Trimma_lantana")

v(max.max.sst.all, cbpos = "l", pal = vpal, zlim = c(28, 35), grid = T,  Save = T, plotname = "max.max.sst", fileformat = "png")
v(max.max.sst.all, cbpos = "r", pal = vpal, zlim = c(28, 35), grid = T,  Save = T, plotname = "max.max.sst.legend", fileformat = "png")
###minimum

min2010 = stack("Data/sst_data/min_sst2010.tif")
min2011 = stack("Data/sst_data/min_sst2011.tif")
min2012 = stack("Data/sst_data/min_sst2012.tif")
min2013 = stack("Data/sst_data/min_sst2013.tif")
min2014 = stack("Data/sst_data/min_sst2014.tif")
min2015 = stack("Data/sst_data/min_sst2015.tif")
min2016 = stack("Data/sst_data/min_sst2016.tif")
min2017 = stack("Data/sst_data/min_sst2017.tif")
min2018 = stack("Data/sst_data/min_sst2018.tif")

min.sst.all <- stack(min2010, 
                     min2011, 
                     min2012, 
                     min2013, 
                     min2014, 
                     min2015, 
                     min2016, 
                     min2017, 
                     min2018)

min.min.sst.all <- min(min.sst.all)

vpalmin = fish(100, alpha = 1, begin = 0, end = 0.8, option = "Coryphaena_hippurus")

v(min.min.sst.all, cbpos = "l", pal = vpalmin, zlim = c(16, 22), grid = T,  Save = T, plotname = "min.min.sst", fileformat = "png")
v(min.min.sst.all, cbpos = "r", pal = vpalmin, zlim = c(16, 22), grid = T,  Save = T, plotname = "min.min.sst.legend", fileformat = "png")


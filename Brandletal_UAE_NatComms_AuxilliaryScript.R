# Auxilliary script that accomplished data processing for sea surface temperatures. Files were scraped from the MODIS Aqua Database (https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/Daily/4km/sst) and placed in folders with respective year labels.

library(oceanmap)
library(ncdf4)
library(raster)
library(viridis)
library(fishualize)
library(rvest)     
library(stringr)
library(RStoolbox)

# function for getting quantiles of temperature ranges
treat.tempdat <- function(data){
  fort.dat <- fortify(data) %>%
    pivot_longer(cols = -c(1:2), names_to = "day", values_to = "temp") %>%
    mutate(lat = round(y, digits = 1), long = round(x, digits = 1)) %>%
    inner_join(loccoord) %>%
    group_by(Site) %>%
    mutate(id = row_number()) %>%
    drop_na(temp)
  
  sum.dat <- fort.dat %>%
    group_by(Site) %>% 
    do(data.frame(quant.vals = quantile(.$temp, probs = c(0.025, 0.50, 0.975)))) %>%
    ungroup() %>%
    mutate(quantile = rep(c("low", "med", "high"), 6)) %>%
    pivot_wider(id_cols = Site, names_from = quantile, values_from = quant.vals) %>%
    inner_join(fort.dat) %>%
    filter(temp < low | temp > high) %>%
    mutate(quant = case_when(temp < low ~ "low",
                             TRUE ~ "high"))
  return(sum.dat)
}


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

temp2010 <- treat.tempdat(stacked_rasters2010)
write.csv(temp2010, file = "Data/temp_summaries/temp2010.csv")

max_sst2010 <- stackApply(stacked_rasters2010, indices = rep(1, nlayers(stacked_rasters2010)), fun = max)
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

temp2011 <- treat.tempdat(stacked_rasters2011)
write.csv(temp2011, file = "Data/temp_summaries/temp2011.csv")

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

temp2012 <- treat.tempdat(stacked_rasters2012)
write.csv(temp2012, file = "Data/temp_summaries/temp2012.csv")


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

temp2013 <- treat.tempdat(stacked_rasters2013)
write.csv(temp2013, file = "Data/temp_summaries/temp2013.csv")

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


temp2014 <- treat.tempdat(stacked_rasters2014)
write.csv(temp2014, file = "Data/temp_summaries/temp2014.csv")

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

temp2015 <- treat.tempdat(stacked_rasters2015)
write.csv(temp2015, file = "Data/temp_summaries/temp2015.csv")

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

temp2016 <- treat.tempdat(stacked_rasters2016)
write.csv(temp2016, file = "Data/temp_summaries/temp2016.csv")

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

temp2017 <- treat.tempdat(stacked_rasters2017)
write.csv(temp2017, file = "Data/temp_summaries/temp2017.csv")

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

temp2018 <- treat.tempdat(stacked_rasters2018)
write.csv(temp2018, file = "Data/temp_summaries/temp2018.csv")

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


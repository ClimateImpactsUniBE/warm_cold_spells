
require(ncdf4); require(lubridate)

data_path <- 'path-to-data'

compDoY <- function(yy,mm,dd){
  if ((yy %% 4)==0){
    prevDays <- c(0,31,60,91,121,152,182,213,244,274,305,335)
  } else{
    prevDays <- c(0,31,59,90,120,151,181,212,243,273,304,334)
  }
  return(prevDays[mm]+dd)
}

# Load daily temperature data
file <- nc_open(paste(data_path,'t2m/ERA5_daily_2m_temperature_',index,'.nc',sep=''))
lon <- ncvar_get(file,'lon')
lat <- ncvar_get(file,'lat')
date_vec <- as.Date(ncvar_get(file,'time'),origin='1900-01-01 00:00:00')
t2m <- ncvar_get(file,'t2m')
nc_close(file)

doyVec <- mapply(compDoY,year(date_vec),month(date_vec),day(date_vec))

# Normalisation; 30-day window
t2m_rescaled <- 0*t2m
for (tt in 1:length(tasVec)){
  inds <- which(((abs(doyVec-doyVec[tt])<=15)|(abs(doyVec-doyVec[tt])>=350)) & (abs(date_vec-date_vec[tt])<=1100))
  t_temp <- t2m[,,inds]
  t2m_rescaled[i,j,tt] <- (t_temp-apply(t_temp,c(1,2),mean))/apply(t_temp,c(1,2),sd)
}

londim <- ncdim_def("lon","degrees_east",lon)
latdim <- ncdim_def("lat","degrees_north",lat)
tdim <- ncdim_def("time","days since 1979-01-01",1:(dim(t2m)[3])-1)
var_def <- ncvar_def('t2m_rescaled','float',list(londim,latdim,tdim),NA,'Rescaled daily 2-meter temperature',prec="float")
ncfname <- paste(data_path,'t2m/ERA5_daily_2m_temperature_rescaled_',index,'.nc',sep='')
ncout <- nc_create(ncfname,list(var_def1),force_v4=T)
ncvar_put(ncout,var_def,t2m_rescaled)
nc_close(ncout)




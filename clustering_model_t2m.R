
library(glmnet); library(quantreg); library(WRTDStidal)
library(clusterCrit); library(dbscan); library(cluster); library(inaparc)

data_path <-'data-directory'
save_path <- 'save-directory'

dERA5 <- seq.Date(as.Date('1979-1-1'),as.Date('2020-12-31'),by='day')
mth_seq <- list(c(12,1,2),3:5,6:8,9:11)

comp.EOF <- function(lon,lat,data,n,anomaly=T,weight=T){
  
  # lon, lat: longitude/latitude vectors (required for weighting)
  # n: number of EOFs to calculate
  # anomaly: compute from anomalies? (default is True)
  # weight: apply latitudinal weight? (default is True)
  
  grid <- expand.grid(x=lon,y=lat)
  steps <- dim(data)[3]
  numCells <- length(lon)*length(lat)
  
  data.df = data.frame(x=rep(grid$x,steps),y=rep(grid$y,steps),
                       timestep=rep(c(1:steps),each=numCells),
                       d=as.vector(data))
  
  # Weigh data
  if (weight==T){
    data.df$d <- cos(pi*data.df$y/180)*data.df$d 
  }
  
  data.EOF <- acast(data.df,x+y~timestep,value.var='d')
  
  if (anomaly==T){
    data.EOF <- sweep(data.EOF,1,rowMeans(data.EOF),"-")
  }
  
  EOF <- prcomp(t(data.EOF),scale.=F,center=T,tol=0.05)
  
  return(list(EOF$x[,1:n],(cumsum(EOF$sdev^2)/sum(EOF$sdev^2))[1:n]))
}

################################################################################
## REGRESSION MODEL ##
################################################################################

# ---------------------------- #
'CHOOSE PARAMETERS'
# ---------------------------- #

# Timescale selection
L <- 21
NN <- floor(length(dERA5)/L)
stepVec <- c(rep(1:NN,each=L),rep(NN+1,length(dERA5)%%L))
df_temp <- data.frame(step=stepVec,month=month(dERA5))
monthVec <- aggregate(month~step,df_temp,median)$month

# Choose hemisphere

hm <- 'NH'
# hm <- 'SH'

indexList <- switch(hm,NH=c('5N','25N','45N','65N','85N'),SH=c('55S','35S','15S','5N'))
lat_range <- switch(hm,NH='0-85N',SH='80-0S')
num_cores <- length(indexList)

# Choose temperature quantile to model

QQ <- '0.95'
# QQ <- '0.05'

# ---------------------------- #
'LOAD ERA5 GRIDDED DATA'
# ---------------------------- #

file <- nc_open(paste0(data_path,"Z500/ERA5_",L,"day_average_Z500_",lat_range,"_1979-2020.nc"))
lonB <- ncvar_get(file,'lon')
latB <- ncvar_get(file,'lat')
meanZ <- ncvar_get(file,'z')
nc_close(file)

meanZ <- meanZ[,order(latB),]
latB <- latB[order(latB)]
gc()

# First 25 PCs (seasonal; ~85-95% of variance explained)
coeffs <- vector(mode="list",length=4)
for (sindex in 1:4){
  mth_sea <- mth_seq[[sindex]]
  eof <- comp.EOF(lonB,latB,meanZ[,,monthVec%in%mth_sea],25)
  coeffs[[sindex]] <- eof[[1]]
}

# ---------------------------- #
'REGRESSION MODEL'
# ---------------------------- #

warn <- getOption("warn")
options(warn=-1)
registerDoParallel(cores=num_cores)
foreach(index=indexList) %dopar% {
  
  print(paste(index,' began'),sep='')
  
  file <- nc_open(paste(data_path,'t2m/ERA5_avg_daily_temperature_rescaled_',L,'days_1deg_',index,'_1979-2020.nc',sep=''))
  lon <- ncvar_get(file,'lon')
  lat <- ncvar_get(file,'lat')
  t2m <- ncvar_get(file,'t2m')
  nc_close(file); gc()
  
  # NH
  if (hm=='NH'){
    t2m <- t2m[,lat>=0,]
    lat <- lat[lat>=0]
  }else{
    t2m <- t2m[,lat<=0,]
    lat <- lat[lat<=0]
  }
  
  matCoeff_pc <- array(0,c(nrow(t2m),ncol(t2m),4,ncol(coeffs[[1]])))
  matDR <- array(0,c(nrow(t2m),ncol(t2m),4))
  
  for (i in 1:nrow(t2m)){
    print(i)
    for (j in 1:ncol(t2m)){
      for (sindex in 1:4){
        
        tryCatch({ # sometimes rq() doesn't converge and throws an error
          
          mth_sea <- mth_seq[[sindex]]
          
          tVec <- t2m[i,j,monthVec%in%mth_sea]
          
          ### QUANTILE REGRESSION ###
          
          df.train <- scale(cbind(tVec,coeffs[[sindex]]))
          df.train <- as.data.frame(df.train)
          colnames(df.train) <- c('t',colnames(coeffs[[sindex]]))
          
          # Stepwise covariate selection (somehow LASSO doesn't performs as well)
          
          listAIC <- vector()
          modelVars <- vector()
          remainingVars <- 2:ncol(df.train)
          currentAIC <- AIC(rq(t~1,tau=as.numeric(QQ),data=df.train))[1]
          listAIC[1] <- AIC(rq(t~1,tau=as.numeric(QQ),data=df.train))[1]
          foundBetterModel <- TRUE
          
          while (foundBetterModel){
            
            # Loop on all remaining variables
            foundBetterModel <- F
            minAIC <- tail(listAIC,1) # current minimum AIC value
            
            for (var in remainingVars){
              
              for (varBis in modelVars){
                
                fit <- rq(t~.,tau=as.numeric(QQ),data=df.train[,c(1,modelVars[!modelVars==varBis],var)])
                if (AIC(fit)[1]<minAIC){
                  foundBetterModel <- T
                  whichModel <- c(modelVars[!modelVars==varBis],var)
                  minAIC <- AIC(fit)[1]
                }
              }
              fit <- rq(t~.,tau=as.numeric(QQ),data=df.train[,c(1,modelVars,var)])
              if (is.infinite(AIC(fit)[1])){
                break
              }
              # Select variable with best improvement
              if (AIC(fit)[1]<minAIC){
                foundBetterModel <- T
                whichModel <- c(modelVars,var)
                selVarAIC <- var
                minAIC <- AIC(fit)[1]
              }
            }
            
            # Update
            if (foundBetterModel){
              listAIC <- c(listAIC,minAIC)
              modelVars <- whichModel
              remainingVars <- remainingVars[!remainingVars%in%whichModel]
            } else{
              break
            }
          }
          preds <- colnames(df.train)[modelVars]
          if (length(preds)==0){
            model.formula <- as.formula("t~1")
            fit.train <- rq(formula=model.formula,data=df.train,tau=as.numeric(QQ))
          }else{
            model.formula <- as.formula(paste("t~",paste(preds,collapse="+"),sep=""))
            fit.train <- rq(formula=model.formula,data=df.train,tau=as.numeric(QQ))
          }
          x <- fit.train$coefficients[-1]
          ii <- 1
          cVec <- 0*1:(ncol(df.train)-1)
          for (pp in preds){
            cVec[which(colnames(df.train)==pp)-1] <- x[ii]
            ii <- ii+1
          }
          
          # Save coefficients
          matCoeff_pc[i,j,sindex,] <- cVec
          matDR[i,j,sindex] <- goodfit(resid(fit.train),resid(rq(t~1,tau=as.numeric(QQ),data=df.train)),as.numeric(QQ))
          
        },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
  }
  londim <- ncdim_def("lon","degrees_east",lon)
  latdim <- ncdim_def("lat","degrees_north",lat)
  sdim <- ncdim_def("season","integer",1:4)
  pcdim <- ncdim_def("pc","integer",1:25)
  var_def1 <- ncvar_def('coeff_pc','float',list(londim,latdim,sdim,pcdim),NA,'NAO coefficient',prec="float")
  var_def2 <- ncvar_def('dev_ratio','float',list(londim,latdim,sdim),NA,'Deviance ratio',prec="float")
  ncfname <- paste0(save_path,'Coefficients_rescaled_z500_t2m_',L,'day_avg_stepwise_',QQ,'_1deg_',hm,'_',index,'.nc')
  ncout <- nc_create(ncfname,list(var_def1,var_def2),force_v4=T)
  ncvar_put(ncout,var_def1,matCoeff_pc)
  ncvar_put(ncout,var_def2,matDR)
  nc_close(ncout)
  
  print(paste(index,' ended'),sep='')
}
options(warn=warn)


# ---------------------------- #
'MERGE RESULTS'
# ---------------------------- #

Llat <- switch(hm,NH=85,SH=66)
matCoeff_pc <- array(0,c(360,Llat,4,25))
matDR <- array(0,c(360,Llat,4))
lat <- 0*1:Llat
index <- 1
for (ll in indexList){
  file <- nc_open(paste0(save_path,'Coefficients_rescaled_z500_t2m_',L,'day_avg_stepwise_',QQ,'_1deg_',hm,'_',ll,'.nc'))
  lon <- ncvar_get(file,'lon')
  LL <- length(ncvar_get(file,'lat'))
  lat[index:(index+LL-1)] <- ncvar_get(file,'lat')[LL:1]
  matCoeff_pc[,index:(index+LL-1),,] <- ncvar_get(file,'coeff_pc')[,LL:1,,]
  matDR[,index:(index+LL-1),] <- ncvar_get(file,'dev_ratio')[,LL:1,]
  nc_close(file)
  index <- index+LL
}
londim <- ncdim_def("lon","degrees_east",lon)
latdim <- ncdim_def("lat","degrees_north",lat)
sdim <- ncdim_def("season","integer",1:4)
pcdim <- ncdim_def("pc","integer",1:25)
var_def1 <- ncvar_def('coeff_pc','float',list(londim,latdim,sdim,pcdim),NA,'PC coefficients',prec="float")
var_def2 <- ncvar_def('dev_ratio','float',list(londim,latdim,sdim),NA,'Deviance ratio',prec="float")
ncfname <- paste0(save_path,'Coefficients_rescaled_z500_t2m_',L,'day_avg_stepwise_',QQ,'_1deg_',hm,'.nc')
ncout <- nc_create(ncfname,list(var_def1,var_def2),force_v4=T)
ncvar_put(ncout,var_def1,matCoeff_pc)
ncvar_put(ncout,var_def2,matDR)
nc_close(ncout)

for (ll in indexList){
  file.remove(paste0(save_path,'Coefficients_rescaled_z500_t2m_',L,'day_avg_stepwise_',QQ,'_1deg_',hm,'_',ll,'.nc'))
}

################################################################################
## SPATIAL CLUSTERING ##
################################################################################

kVec <- 3:50

# Hemisphere
hm <- 'NH'
# hm <- 'SH'

# Quantile
QQ <- '0.95'
# QQ <- '0.05'

file <- nc_open(paste(save_path,'Coefficients_rescaled_z500_t2m_',L,'day_avg_stepwise_',QQ,'_1deg_',hm,'.nc',sep=''))
lon <- ncvar_get(file,'lon')
lat <- ncvar_get(file,'lat')
matCoeff_pc <- ncvar_get(file,'coeff_pc')
matDR <- ncvar_get(file,'dev_ratio')
nc_close(file)

grid <- expand.grid(x=lon,y=lat)

# Loop on seasons
for (sindex in c(1,3)){
  
  print(sindex)
  season <- c('DJF','MAM','JJA','SON')[sindex]
  df <- data.frame(lon=rep(grid$x,dim(matCoeff_pc)[4]),
                   lat=rep(grid$y,dim(matCoeff_pc)[4]),
                   pc=rep(c(1:dim(matCoeff_pc)[4]),each=nrow(grid)),
                   dr=rep(as.vector(matDR[,,sindex]),dim(matCoeff_pc)[4]),
                   c=as.vector(matCoeff_pc[,,sindex,]))
  df <- df[df$dr>=0.4,]

  # dbscan filtering before clustering
  data <- cbind(df$lon[df$pc==1],df$lat[df$pc==1])
  dGeo <- as.dist(geosphere::distm(data))
  db <- dbscan(dGeo,eps=500E3,minPts=10)
  df$cluster <- rep(db$cluster,dim(matCoeff_pc)[4])
  df <- df[df$cluster>0,]
  
  A <- matrix(df$c,nrow=nrow(df)/dim(matCoeff_pc)[4],ncol=dim(matCoeff_pc)[4])
  # B <- matrix(df$lon,nrow=nrow(df)/dim(matCoeff_pc)[4],ncol=dim(matCoeff_pc)[4])
  # C <- matrix(df$lat,nrow=nrow(df)/dim(matCoeff_pc)[4],ncol=dim(matCoeff_pc)[4])
  W <- dist(A,method='euclidean'); gc()
  
  # Hierarchical clustering
  K <- 26
  hc <- hclust(W, method="ward.D" )
  cl <- cutree(hc,k=K)
  pam.clusters <- data.frame(lon=B[,1],lat=C[,1],g=cl)
  
  # No dbscan filtering
  df_res <- matrix(0,nrow=length(pam.clusters$g),ncol=3)
  index <- 1
  clus_ind <- 1
  for (uu in 1:max(pam.clusters$g)){
    
    data <- cbind(pam.clusters$lon[pam.clusters$g==uu],pam.clusters$lat[pam.clusters$g==uu])
    
    # Save
    if (nrow(data)>0){
      df_res[index:(index+nrow(data)-1),1:2] <- data
      df_res[index:(index+nrow(data)-1),3] <- clus_ind
      index <- index+nrow(data)
      clus_ind <- clus_ind+1
    }
  }
  df_res <- df_res[1:(index-1),]
  df_res <- as.data.frame(df_res)
  colnames(df_res) <- c('lon','lat','cluster')
  write.table(df_res,row.names=F,
              file=paste0(save_path,'clusters_t2m_',L,'day_',season,'_z500_stepwise_',QQ,'_1deg_',hm,'_HC_K=',K,'_no_dbscan.txt'))
  
  # Apply dbscan filtering
  df_res <- matrix(0,nrow=length(pam.clusters$g),ncol=3)
  index <- 1
  clus_ind <- 1
  for (uu in 1:max(pam.clusters$g)){
    
    data <- cbind(pam.clusters$lon[pam.clusters$g==uu],pam.clusters$lat[pam.clusters$g==uu])
    dGeo <- as.dist(geosphere::distm(data))
    db <- dbscan(dGeo,eps=500E3,minPts=10)
    data <- data[db$cluster>0,]
    
    # Save
    if (nrow(data)>0){
      df_res[index:(index+nrow(data)-1),1:2] <- data
      df_res[index:(index+nrow(data)-1),3] <- clus_ind
      index <- index+nrow(data)
      clus_ind <- clus_ind+1
    }
  }
  df_res <- df_res[1:(index-1),]
  df_res <- as.data.frame(df_res)
  colnames(df_res) <- c('lon','lat','cluster')
  write.table(df_res,row.names=F,
              file=paste0(save_path,'clusters_t2m_',L,'day_',season,'_z500_stepwise_',QQ,'_1deg_',hm,'_HC_K=',K,'.txt'))
}

################################################################################
## MERGE CLUSTERS ##
################################################################################

L <- 21

hm <- 'NH'
# hm <- 'SH'

QQ <- '0.95'
# QQ <- '0.05'

file <- nc_open(paste(save_path,'Coefficients_z500_rescaled_t2m_',L,'day_avg_stepwise_',QQ,'_1deg_',hm,'.nc',sep=''))
lon <- ncvar_get(file,'lon')
lat <- ncvar_get(file,'lat')
nc_close(file)

matClusters <- array(NA,c(length(lon),length(lat),4))
df <- read.table(paste0(save_path,'clusters_t2m_',L,'day_DJF_z500_stepwise_',QQ,'_1deg_',hm,'_PAM_best.txt'),
                 header=T)
for (i in 1:nrow(df)){
  matClusters[which(lon==df$lon[i]),which(lat==df$lat[i]),1] <- df$cluster[i]
}
df <- read.table(paste0(save_path,'clusters_t2m_',L,'day_MAM_z500_stepwise_',QQ,'_1deg_',hm,'_PAM_best.txt'),
                 header=T)
for (i in 1:nrow(df)){
  matClusters[which(lon==df$lon[i]),which(lat==df$lat[i]),2] <- df$cluster[i]
}
df <- read.table(paste0(save_path,'clusters_t2m_',L,'day_JJA_z500_stepwise_',QQ,'_1deg_',hm,'_PAM_best.txt'),
                 header=T)
for (i in 1:nrow(df)){
  matClusters[which(lon==df$lon[i]),which(lat==df$lat[i]),3] <- df$cluster[i]
}
df <- read.table(paste0(save_path,'clusters_t2m_',L,'day_SON_z500_stepwise_',QQ,'_1deg_',hm,'_PAM_K=20.txt'),
                 header=T)
for (i in 1:nrow(df)){
  matClusters[which(lon==df$lon[i]),which(lat==df$lat[i]),4] <- df$cluster[i]
}
londim <- ncdim_def("lon","degrees_east",lon)
latdim <- ncdim_def("lat","degrees_north",lat)
sdim <- ncdim_def("season","integer",1:4)
var_def <- ncvar_def('cluster','integer',list(londim,latdim,sdim),NA,'Z500 cluster',prec="float")
ncfname <- paste0(save_path,'clusters_t2m_',L,'day_stepwise_',QQ,'_1deg_',hm,'_PAM.nc')
ncout <- nc_create(ncfname,list(var_def),force_v4=T)
ncvar_put(ncout,var_def,matClusters)
nc_close(ncout)

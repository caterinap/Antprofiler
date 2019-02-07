
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                     Build SDMs for Ants
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                     Charge library
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

library(biomod2)
library(raster)
library(rgdal)
library(beepr)
library(rasterVis)
library(sdmpredictors)

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                     Environmental variables
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

predictors<- load_layers(c("WC_bio14","WC_bio3","WC_bio2","WC_bio15","WC_bio8",
                           "WC_bio18","WC_bio13","WC_bio19"))# datadir = tempdir())

names(predictors) <- c("bio14","bio3","bio2","bio15","bio8","bio18","bio13","bio19")
plot(predictors)


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                    Ocuurence data
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

load("GPS_Lepisiota_canescens.RData")
df[, 1] <- as.character(df[,1])
  
z <- strsplit(df[1,1], "_")
sp <- paste(z[[1]][1], z[[1]][2], sep = ".")
  
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                    Current distribution
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

myRespName <- sp
myRespXY <- df[, c(3, 2)]
myResp <-  rep(1, dim(df)[1])
  
myExpl <- predictors

myBiomodData<-BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 1000,
                                     PA.strategy = 'disk',
                                     PA.dist.min= 10000)
  
myBiomodOptions <- BIOMOD_ModelingOptions()
  
myBiomodModelOut <- BIOMOD_Modeling ( myBiomodData, 
                                        models = c('GLM','GBM','ANN', 'FDA','MARS','RF',"MAXENT"), #GMB = BRT
                                        models.options = myBiomodOptions, 
                                        NbRunEval=4, 
                                        DataSplit=70, 
                                        VarImport=1, 
                                        models.eval.meth = c('TSS','ROC'),
                                        Prevalence = 0.5, # Equal weightings were given to presences and PAs
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE, 
                                        do.full.models = TRUE, # run the final model with 100% of the data
                                        modeling.id= paste("mod_", df[1,1], sep = ""))
  
  
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                         chosen.models = 'all',
                                         em.by = 'all',
                                         eval.metric = 'TSS',
                                         eval.metric.quality.threshold = 0.6,
                                         models.eval.meth = c('TSS','ROC'),
                                         prob.mean.weight = TRUE, 
                                         prob.mean.weight.decay = 'proportional',
                                         VarImport = 10) 
  
  
a <- get_evaluations(myBiomodEM)[2]

myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut, 
                                            new.env = myExpl, 
                                            proj.name = "current",
                                            selected.models = 'all',
                                            binary.meth = 'TSS',
                                            compress = 'gzip',
                                            build.clamping.mask = FALSE)


myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection, 
                                             binary.meth = 'TSS',
                                             EM.output = myBiomodEM)
    
    
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                    Plot maps
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

myCurrentProj <- get_predictions(myBiomodProjection)
my.colors = colorRampPalette(c("cornsilk", "darkorchid4")) 
plot(myCurrentProj ,col=my.colors(1000),axes=FALSE, 
     box=FALSE, main="Lepisiota canescens")


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                     Future distribution
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

myExplFuture <- stack(future_predictors)
head(myExplFuture)
plot(myExplFuture)
names(myExpl)
names(myExplFuture)

myBiomodProjectionFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                              new.env = myExplFuture,
                                              proj.name = '2050',
                                              selected.models = 'all',
                                              binary.meth = 'TSS',
                                              compress = FALSE,
                                              build.clamping.mask = TRUE)

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#                    Plot maps
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

my.colors = colorRampPalette(c("cornsilk", "darkorchid4")) 
plot(myBiomodProjectionFuture,col=my.colors(1000),axes=FALSE, 
     box=FALSE, main="Lepisiota canescens Future")




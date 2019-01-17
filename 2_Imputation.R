#### Second part of the analysis: calibrate the imputation and impute the dataset ####
library(missForest)
library(ggplot2)
library(reshape2)
library(Amelia)
library(data.table)

### Calibration: Select optimal number of random forest trees (ntree) 

df <- dfall[rownames(dfall) %in% sp1002, c(1:24, 31:170)] #select target species and traits

dfntree <- data.frame(matrix(NA, nrow = 500, ncol = 6))
colnames(dfntree) <- c("ntree", "Ubiquitous","NestingType","DisturbanceSpecialist","IndependentFoundation","SuperColonial")

for (n in (1:500)){
  dfntree[n, 1] <- n
    o <- missForest(df, maxiter = 10, ntree = n, variablewise = TRUE) 
    dfntree[n,2] <- o$OOBerror[1]; dfntree[n,3] <- o$OOBerror[2];    #save OOBerror for target traits only
    dfntree[n,4] <- o$OOBerror[3]; dfntree[n,5] <- o$OOBerror[7]; 
    dfntree[n,6] <- o$OOBerror[24]
    save(dfntree, file ="./data/data_antprofiler/missForest/error_ntre_selection_dfall.RData")
  }

dfntree

### Calibration: Select optimal number of variables randomly sampled at each random forest split (argument mtry) 

df <- dfall[rownames(dfall) %in% sp1002, c(1:24, 31:170)] #select target species and traits

dfmtry <- data.frame(matrix(NA, nrow = 25, ncol = 6))
colnames(dfmtry) <- c("mtry", "Ubiquitous","NestingType","DisturbanceSpecialist","IndependentFoundation","SuperColonial")

for (n in (1:25)){
  dfmtry[n, 1] <- n
  o <- missForest(df, maxiter = 15, mtry = n, ntree = 100 , variablewise = TRUE) 
  dfmtry[n,2] <- o$OOBerror[1]; dfmtry[n,3] <- o$OOBerror[2]; 
  dfmtry[n,4] <- o$OOBerror[3]; dfmtry[n,5] <- o$OOBerror[7];  #save OOBerror for target traits only
  dfmtry[n,6] <- o$OOBerror[24]; dfmtry[counter,7] <- o$OOBerror[23]
  save(dfmtry, file ="./data/data_antprofiler/missForest/error_mtry_selection_dfall.RData")
}

### Calibration: Select best number of eigenvectors (k) per variable

df <- dfall[rownames(dfall) %in% sp1002, c(1:24, 31:170)] #select target species and traits
names(df)

dfk <- data.frame(matrix(NA, nrow = 1000, ncol = 7))
colnames(dfk) <- c("k","IndependentFoundation","SuperColonial","Ubiquitous", "ColonyFoundation", "DisturbanceSpecialist","NestingType")
counter=1
for (i in 1:10) {
for (n in (1:100)){
  dfk[counter, 1] <- n
  dfimp <- df[, 1: ( 24+n)]
  o <- missForest(dfimp, maxiter = 15, mtry = 10, ntree = 100 , variablewise = TRUE) 
  dfk[counter,2] <- o$OOBerror[1]; dfk[counter,3] <- o$OOBerror[2]; #save OOBerror for target traits only
  dfk[counter,4] <- o$OOBerror[3]; dfk[counter,5] <- o$OOBerror[7]; 
  dfk[counter,6] <- o$OOBerror[24];dfk[counter,7] <- o$OOBerror[23]
  save(dfk, file ="./data/data_antprofiler/missForest/error_k_selection_dfall.RData")
  counter=counter+1
}
 }

summary(dfk)
#average number of eigenvectors of 10 runs
dfk2<-aggregate(.~k,data=dfk,FUN=mean)

#find optimal number of eigenvectors (the one with lower error values)
bestk <- vector()

for(i in 1:6){
  bestk[i] <- which(dfk2[, i+1] == min(dfk2[, i+1])) [1]
}

#save optimal number of eigenvectors per trait
save(bestk, file ="./data/data_antprofiler/missForest/error_bestk_selection.RData") #bestk


### Comparison between OOB and true imputation error

load("./data/data_antprofiler/dfall.RData")

# prepare empty dataset
dferror1002 <- data.frame(matrix(NA, nrow = 5, ncol = 102))
rownames(dferror1002) <- c("IndependentFoundation","SuperColonial","Ubiquitous", "DisturbanceSpecialist","NestingType")
colnames(dferror1002) <- c(paste(rep ("rand", 100), 1:100, sep = ""), "moy", "se")
dferror1002OOB <- dferror1002

# load information on optimal number of eigenvectors
load("./data/data_antprofiler/missForest/error_bestk_selection.RData") #bestk
df <- dfall[rownames(dfall) %in% sp1002,c("IndependentFoundation","SuperColonial","Ubiquitous", "DisturbanceSpecialist","NestingType")]

tauxNA1002 <- c(62.60, 62.50, 5.44, 52.52,5.74) #this is the proportion of NA's per target trait

for (i in 1:5){ # i is the trait column
  l <- which(is.na(df[, i]) == TRUE)
  dfimp <- df[-l, 1: (24+bestk[i])]
  N <- round(tauxNA1002[i]*dim(dfimp)[1]/100)
  
  #randomly include NAs in the dataset (same proportion as original % of NAs
  for (r in 1 : 100){  
    l <- sample(1:dim(dfimp)[1], N)
    dfmiss <- dfimp
    dfmiss[l,i] <- NA
    o <- missForest(dfmiss, maxiter = 15, mtry = 10, ntree = 100 , variablewise = TRUE) # OOB error
    dferror1002[i,r] <- mixError(data.frame(o$ximp[,i]),data.frame(dfmiss[,i]),data.frame(dfimp[,i])) # "true" error
    dferror1002OOB[i,r] <- o$OOBerror[i]
    save(dferror1002, file ="./data/data_antprofiler/missForest/ManualError_dferror1002NRMSE.RData")
    save(dferror1002OOB, file ="./data/data_antprofiler/missForest/ManualError_dferror1002OOB.RData")
  }
}


# calulate mean and sd of both OOB and "true" error
dferror1002$moy<-apply(dferror1002,1,mean,na.rm=T)
dferror1002$se<-apply(dferror1002,1,sd,na.rm=T)

dferror1002OOB$moy<-apply(dferror1002OOB,1,mean,na.rm=T)
dferror1002OOB$se<-apply(dferror1002OOB,1,sd,na.rm=T)

save(dferror1002, file ="./data/data_antprofiler/missForest/ManualError_dferror1002NRMSE.RData")
save(dferror1002OOB, file ="./data/data_antprofiler/missForest/ManualError_dferror1002OOB.RData")

# plot error from OOB and NRMSE
require(reshape2)
load("./data/data_antprofiler/missForest/ManualError_dferror1002OOB.RData")
load("./data/data_antprofiler/missForest/ManualError_dferror1002NRMSE.RData")
dferror1002NRMSE<-dferror1002

    # prepare datasets for ggplot
dferror1002OOB$trait <- rownames(dferror1002OOB)
oobmelt <- cbind(melt(dferror1002OOB[,c("trait","moy")]),melt(dferror1002OOB[,c("trait","se")]))
oobmelt$variable<-oobmelt$variable<-oobmelt$trait<-NULL
names(oobmelt)<-c("mean","trait","Sd")
oobmelt$Error_type <- "OOB"

dferror1002NRMSE$trait <- rownames(dferror1002NRMSE)
nrmsemelt <- cbind(melt(dferror1002NRMSE[,c("trait","moy")]),melt(dferror1002NRMSE[,c("trait","se")]))
nrmsemelt$variable<-nrmsemelt$variable<-nrmsemelt$trait<-NULL
names(nrmsemelt)<-c("mean","trait","Sd")
nrmsemelt$Error_type <- "NRMSE"

traiterrors<-rbind(oobmelt,nrmsemelt)
    # plot
ggplot(traiterrors,aes(x=trait,y=mean,color=Error_type)) + 
geom_point(position=position_dodge(width=0.5),size = 2.5) + 
geom_errorbar(aes(ymin=mean - Sd, ymax=mean + Sd),position=position_dodge(width=0.5),width = 0.1) +
  xlab("") +
  ylab("Mean error +/- SD") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))



### Final imputation using optimal parameters per traits
load("./data/data_antprofiler/df1002.RData")
load("./data/data_antprofiler/dfall.RData")
load("./data/data_antprofiler/index1002.RData")

  # plot missing values
pdf(file= "./data/data_antprofiler/glmulti/m1904/missmap.pdf", onefile=T, paper='A4', width=10, height=25)
par(mar = c(14,2,2,2))
missmap(dfall[rownames(dfall) %in% sp1002, 1:11], x.cex = 1.5, y.cex = 0.5, col = c("beige", "forestgreen" ), legend = F)
dev.off()


#best parameters are ntree=100 and mtry=6

load("./data/data_antprofiler/missForest/error_bestk_selection.RData") #optimal number of eigenvectors
bestk #1 or 25
names(bestk) <- c("IndependentFoundation","SuperColonial","Ubiquitous", "ColonyFoundation", "DisturbanceSpecialist","NestingType")

#select columns
k=1
df <- dfall[rownames(dfall) %in% sp1002, 
c(1:24, 31:170)]

#100 imputations with different number of eigenvectors per trait
errors <- list()
counter=1
for (i in (1:100)){ # i = imputation ID
  for(k in c(1,2,20,21)){ # k = number of eigenvectors
    o <- missForest(df[,c(1:(24+k))], maxiter = 15, ntree = 100, mtry = 6, variablewise=T)
    dfimputed <- o$ximp
    dfimputed$species <- rownames(dfimputed)
    save(dfimputed, file = paste("./data/data_antprofiler/missForest/dfimputed708_i", i, "_k", k, ".RData", sep = ""))
    errorsdf <- data.frame(nrep=i,neigen=k,trait=names(df[,c(1:(24+k))]),error=o$OOB)
    errors[[counter]] <- errorsdf
    print(paste(i,k))
    save(errors,file = paste("./data/data_antprofiler/missForest/dfimputed708_k", k, ".RData", sep = ""))
    counter=counter+1
    }
}

#aassemble traits with optimal numbers of eigenvectors
dirpath<-"data/data_antprofiler/missForest/"

load("./data/data_antprofiler/missForest/error_bestk_selection.RData") #eigenvectors
names(bestk) <- c("IndependentFoundation","SuperColonial","Ubiquitous", "ColonyFoundation", "DisturbanceSpecialist","NestingType")
bestk

#open all files
bestk

#bestk=1 supercolonialist, disturbance specialist, ubiquitous
all.files<-intersect(list.files(dirpath,full.names=T,pattern="dfimputed"),
                     list.files(dirpath,full.names=T,pattern="k1"))
mylist1<- lapply(all.files, function(x) {
  load(file = x)
  get(ls()[ls()!= "filename"])
})

mylist1 <- lapply(mylist1, function(x) x[(names(x) %in% c("SuperColonial", "Ubiquitous","DisturbanceSpecialist"))])

#bestk=25 IndependentFoundation
all.files<-intersect(list.files(dirpath,full.names=T,pattern="dfimputed"),
                      list.files(dirpath,full.names=T,pattern="k25.R"))
mylist2<- lapply(all.files, function(x) {
  load(file = x)
  get(ls()[ls()!= "filename"])
})

mylist2 <- lapply(mylist2, function(x) x[(names(x) %in% c("IndependentFoundation"))])


#####list of datasets with all columns
all.list <- list()
for (i in 1:100){
  all.list[[i]] <- cbind(mylist1[[i]],mylist2[[i]])
}

save(all.list, file = "./data/data_antprofiler/all_imputed_datasets.RData")



### Evaluate the mean and sd values of OOB error per trait
load("./data/data_antprofiler/missForest/error_bestk_selection.RData") #eigenvectors
names(bestk) <- c("IndependentFoundation","SuperColonial","Ubiquitous", "ColonyFoundation", "DisturbanceSpecialist","NestingType")
bestk

dirpath<-"data/data_antprofiler/missForest/" #where the datasets were saved

# load all files and put them in a list
all.files<-list.files(dirpath,full.names=T,pattern="dfimputed708_k")
mylist1<- lapply(all.files, function(x) {
  load(file = x)
  get(ls()[ls()!= "filename"])
})

mylist1e<-mylist1[[1]]
mylist1e<-rbindlist(mylist1e)
mylist1e<-mylist1e[trait %in% c("SuperColonial", "Ubiquitous","DisturbanceSpecialist")]

mylist2e<-mylist1[[2]]
mylist2e<-rbindlist(mylist2e)
mylist2e<-mylist2e[trait %in% c("IndependentFoundation")]

# calculate mean and sd OOB error for each trait
allerr<-rbindlist(list(mylist1e,mylist2e))
allerr$nrep<-allerr$neigen<-NULL
allerr<-allerr[,error.sd:=sd(error),by="trait"]
allerr<-allerr[,error:=mean(error),by="trait"]
allerr<-unique(allerr)

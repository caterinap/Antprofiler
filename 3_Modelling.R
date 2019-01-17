#### Third part of the analysis: model the invasiveness and extract predictions ####
library(ggplot2)
library(pscl)
library(MASS)
library(GGally)
library(car)
library(glmulti)
library(MuMIn)
library(visreg)
library(boot)
library(pROC)

## load 100 imputed dataset
load("./data/data_antprofiler/dfall.RData")
load("./data/data_antprofiler/all_imputed_datasets.RData")
dfall$species<-rownames(dfall)

#add invasive information from whole dataset (this information was removed for the imputation)
all.list2 <- list()
for (i in 1:100){
  all.list[[i]]$species <- rownames(all.list[[i]])
  all.list2[[i]] <- merge(all.list[[i]],dfall[,c("species","Invasive")],by="species")

  all.list2[[i]]$SuperColonial<-as.factor(all.list2[[i]]$SuperColonial)
  all.list2[[i]]$IndependantFundation<-as.factor(all.list2[[i]]$IndependantFundation)
  all.list2[[i]]$Ubiquitous<-as.factor(all.list2[[i]]$Ubiquitous)
  all.list2[[i]]$DisturbanceSpecialist<-as.factor(all.list2[[i]]$DisturbanceSpecialist)
  all.list2[[i]]$ColonyFoundation<-as.factor(all.list2[[i]]$ColonyFoundation)
  all.list2[[i]]$Invasive<-as.factor(all.list2[[i]]$Invasive)
}


#Full model on one dataset
df <- all.list2[[1]] #select only one dataset to check the model

m1 <- glm(Invasive ~  SuperColonial + Ubiquitous + ColonyFoundation + DisturbanceSpecialist, 
                    data = df,family=binomial(logit))
summary(m1)
vif(m1)
visreg(m1)

summary(dispmod::glm.binomial.disp(m1))

#model assumptions
E2 <- resid(m1, type="pearson")
F2 <- fitted(m1, type="response")
plot(x=F2, y=E2, xlab="fitted values", ylab="Pearson residuals")
abline(h=0, lty=2)

# Cook's distance
plot(cooks.distance(m1), ylim=c(0,1), ylab="Cook distance values", type="h")

plot(m1)

#find best model based on aicc
m2 <- glmulti(m1, level = 1, crit="aicc", family = binomial(logit))

summary(m2) #best model is full model
weightable(m2) #all models

#some models are comparable (AICc<2) -> average them
mm1 <- glm(Invasive ~ 1 + SuperColonial + ColonyFoundation + DisturbanceSpecialist, data = df,family=binomial)
mm2 <- glm(Invasive ~ 1 + SuperColonial + Ubiquitous + ColonyFoundation + DisturbanceSpecialist, data = df,family=binomial)

#average the models
mm.avg <- model.avg(mm1,mm2)
summary(mm.avg)

summary(m1) #averaged and full model are very similar

#cross-validation - leave-one-out
cost <- function(r, pi = 0) mean(abs(r - pi) > 0.5) 
m1.cv <- cv.glm(data = df, m1, cost, K = 1002)  # leave-one-out cross validation

m1.cv$delta #error rate [1] 0.007982040

#cross-validation - groups of 19 
crossvalerr<-vector()
for (i in 1:100){
m1.cv <- cv.glm(data = df, m1, cost, K = 19)
print(m1.cv$delta)
crossvalerr[i] <- m1.cv$delta
}
mean(crossvalerr) #0.007804391
sd(crossvalerr) #0.0004991028


#area under the curve
m.roc <- roc(df$Invasive, predict(m1, backtransform = TRUE))
plot(m.roc)
m.roc #AUC= 0.951

#pseudo R-squared
r.squaredGLMM(m1) #0.5314497

#chisq test
anova(m1,test='Chisq')

# predictions
preds<-predict(m1,type="response",se.fit=T)
df$predictedInva<-preds$fit
df$predictedInva_SE<-preds$se.fit

#Super-invasives: 
invasives <- df[df$Invasive == 1,]
quant05 <- quantile(invasives$predictedInva, probs = seq(0, 1, 0.05))[2]
potentialinvasives <- df[df$Invasive == 0 & df$predictedInva >= quant05,]

write.csv(potentialinvasives,"potentialInvaders.csv",row.names=F)

##########repeat the same glm on 100 datasets and identify the invaders############

alldf<-list()
potentialinv<-list()

for (i in 1:100){
df <- all.list2[[i]]
m1 <- glm(Invasive ~  SuperColonial + Ubiquitous + ColonyFoundation + DisturbanceSpecialist, 
          data = df,family=binomial(logit))
preds<-predict(m1,type="response",se.fit=T)
df$predictedInva<-preds$fit
df$predictedInva_SE<-preds$se.fit
qq <- quantile(df[df$Invasive == 1,]$predictedInva, probs = seq(0, 1, 0.05))[2]
potentialinv[[i]] <- data.frame(n=i,invs=df[df$Invasive == 0 & df$predictedInva >= qq,])
alldf[[i]] <- df
}

all.inv<-do.call(rbind,potentialinv)
unique(all.inv$invs.species)
sapply(all.inv, function(x) length(unique(x)))
#in how many models are these species considered as invasives?
summary(as.factor(all.inv$invs.species))

#summary table
require(data.table)
all.inv <- data.table(all.inv[,.(invs.species,invs.predictedInva,invs.predictedInva_SE)])
all.inv$count<-1
all.inv<-all.inv[, list(count=sum(count), invs.predictedInva=mean(invs.predictedInva)),
   by=c("invs.species")]
write.csv(all.inv,"potentialInvaders_100models.csv",row.names=F)
save(alldf,file = "./data/data_antprofiler/all.predictions100GLm.RData")


#leave-one-out invasive at a time and predict it

alldf<-list()
potentialinv<-list()
df<-all.list2[[1]]
invasp<-unique(df[df$Invasive==1,]$species)
counter=1

for (j in invasp){
for (i in 1:100){
  tryCatch({
  df <- all.list2[[i]]
  df <- df[df$species != j,] #remove one invasive
  m1 <- glm(Invasive ~  SuperColonial + Ubiquitous + ColonyFoundation + DisturbanceSpecialist, 
            data = df,family=binomial(logit))
  preds<-NA
  preds<-predict(m1,type="response",se.fit=T)
  df$predictedInva<-preds$fit
  df$predictedInva_SE<-preds$se.fit
  qq <- quantile(df[df$Invasive == 1,]$predictedInva, probs = seq(0, 1, 0.05))[2]
  potentialinv[[counter]] <- data.frame(n=i,invs=df[df$Invasive == 0 & df$predictedInva >= qq,],removedsp=j)
  df$removedsp<-j
  alldf[[counter]] <- df
  }, error=function(e){})
  counter=counter+1
  print(i);print(j);print(counter)
}
 }


all.inv<-do.call(rbind,potentialinv)
unique(all.inv$invs.species)
sapply(all.inv, function(x) length(unique(x)))
#in how many models are these species considered as invasives?
summary(as.factor(all.inv$invs.species))

#summary table
require(data.table)
all.inv<-data.table(all.inv)
all.inv <- all.inv[,.(invs.species,invs.predictedInva,invs.predictedInva_SE)]
all.inv$count<-1
all.inv<-all.inv[, list(count=sum(count), invs.predictedInva=mean(invs.predictedInva)),
                 by=c("invs.species")]
all.inv[,percmod:=count/1900*100]
write.csv(all.inv,"potentialInvaders_100models_LOO.csv",row.names=F)



#leave-one-out invasive at a time and predict it #########

alldf<-list()
potentialinv<-list()
df<-all.list2[[1]]
invasp<-unique(df[df$Invasive==1,]$species)
counter=1

for (j in invasp){
  for (i in 1:100){
    tryCatch({
      df <- all.list2[[i]]
      df[df$species==j,]$Invasive<-0 #transform invasive into non invasive
      m1 <- glm(Invasive ~  SuperColonial + Ubiquitous + ColonyFoundation + DisturbanceSpecialist, 
                data = df,family=binomial(logit))
      preds<-NA
      preds<-predict(m1,type="response",se.fit=T)
      df$predictedInva<-preds$fit
      df$predictedInva_SE<-preds$se.fit
      qq <- quantile(df[df$Invasive == 1,]$predictedInva, probs = seq(0, 1, 0.05))[2]
      potentialinv[[counter]] <- data.frame(n=i,invs=df[df$Invasive == 0 & df$predictedInva >= qq,],removedsp=j)
      df$removedsp<-j
      alldf[[counter]] <- df
    }, error=function(e){})
    counter=counter+1
    print(i);print(j);print(counter)
  }
}

all.inv<-do.call(rbind,potentialinv)
intersect(unique(all.inv$invs.species),invasp)
sapply(all.inv, function(x) length(unique(x)))
#in how many models are these species considered as invasives?
summary(as.factor(all.inv$invs.species))

#summary table
require(data.table)
all.inv<-data.table(all.inv)
all.inv <- all.inv[,.(invs.species,invs.predictedInva,invs.predictedInva_SE)]
all.inv$count<-1
all.inv<-all.inv[, list(count=sum(count), invs.predictedInva=mean(invs.predictedInva)),
                 by=c("invs.species")]
all.inv[,percmod:=count/1900*100]
write.csv(all.inv,"potentialInvaders_100models_INvaZeros.csv",row.names=F)

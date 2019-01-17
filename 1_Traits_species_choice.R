#### First part of the analysis: select the traits and the species ####
require(data.table)

##load data
load("./data/data_antprofiler/dfall.RData")
  #how many NAs per trait?
apply(dfall,2, function(x) sum(is.na(x))/nrow(dfall)*100)
      
## Step1: check the correlations between traits and invasive status
ind<-combn(colnames(dfall)[-1],2)
listres<-list()
for (i in 1:NCOL(ind)){
  tt<-chisq.test(dfall[,ind[1,i]],dfall[,ind[2,i]])
  dfres<-data.frame(tr1=ind[1,i],tr2=ind[2,i],chi=tt$statistic,
                    pval=tt$p.value)
  listres[[i]]<-dfres
}

listres<-rbindlist(listres)
listres<-listres[pval>0.05 | tr1=="Invasive" | tr2=="Invasive"]
listres[pval<0.06 & tr2=="Invasive"]

unique(unique(listres$tr1),unique(listres$tr2))
         

##Step2: select species in order to have less than 60% NAs in target traits
df <- data.table(dfall)
  #sort by increasing number of NA in columns
y = sort(df[,lapply(.SD, function(x) sum(is.na(x)))])
setcolorder(df, names(y))
  #sort by increasing number of NA in rows
df[, idx := rowSums(is.na(df))] # count nr of NA in rows
df = df[order(idx),] # sort
df$species <- row.names(dfall)
  #select species with less than, e.g., 4 NAs per row                 
df<-df[idx<4]
apply(df,2, function(x) sum(is.na(x))/nrow(df)*100) #how many NAs per trait?
  #if more than 60% of NAs or too few species -> repeat with more/less NAs per row (i.e. NAs per species for all traits)

  #make sure to add back invasive species (in case they were not selected because of too many NAs)
ivs <- dfall[dfall$Invasive==1,]
ivs$species<-rownames(ivs)
setcolorder(ivs, names(df))

df2<-data.table(rbind(df,ivs))
df2<-unique(df2,by="species")      

save(df2, file = "./data/data_antprofiler/df1002.RData")    

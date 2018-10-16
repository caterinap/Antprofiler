#### First part of the analysis: explore the dataset, select the traits and the species ####
require(data.table)

#load data
load("./data/data_antprofiler/dfall.RData")
apply(dfall,2, function(x) sum(is.na(x))/nrow(dfall)*100)

#example to select species in order to have less than 60% NAs in target traits
df <- data.table(dfall) #dfall is the table with all traits

y = sort(df[,lapply(.SD, function(x) sum(is.na(x)))]) # nr of NA in columns, increasing
setcolorder(df, names(y))

df[, idx := rowSums(is.na(df))] # count nr of NA in rows
df = df[order(idx),] # sort by nr of NA in rows
df$species <- row.names(dfall)
df<-df[idx<4]
apply(df,2, function(x) sum(is.na(x))/nrow(df)*100)
df[, idx := NULL] # idx not needed anymore

#add back invasive species (in case they were not selected)
ivs <- dfall[dfall$Invasive==1,]
ivs$species<-rownames(ivs)
setcolorder(ivs, names(df))

df2<-data.table(rbind(df,ivs))
df2<-unique(df2,by="species")

#correlations between traits
listres<-list()
for (i in 1:NCOL(ind)){
  tt<-chisq.test(df1002[,ind[1,i]],df1002[,ind[2,i]])
  dfres<-data.frame(tr1=ind[1,i],tr2=ind[2,i],chi=tt$statistic,
                    pval=tt$p.value)
  listres[[i]]<-dfres
}


listres<-rbindlist(listres)
listres<-listres[pval>0.05 | tr1=="Invasive" | tr2=="Invasive"]
listres[pval<0.06 & tr2=="Invasive"]

unique(unique(listres$tr1),unique(listres$tr2))
      
      

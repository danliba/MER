rm(list = ls())
setwd('D:/Maestria/MER/SOTON/clases/deepSeaEco/assigment2')


dir()
library(readxl)
library(vegan)
library(ade4)

data0<-read_excel("ROVtransect_DATA.xlsx", sheet = "abundance")
View(data0)
data1<-data0[,5:16]
data1[is.na(data1)]<-0
View(data1)


data_bio_t<-t(data1)
View(data_bio_t)

par(mfrow=c(1,1))
boxplot(data_bio_t)
#Primero con todos los datos fisicos del 2003 al 2015 promedio anual
#row.names(data_c_bio)<-c(2003:2015)

# here we reduce the distance
datalt_bio_t<-log(data_bio_t+1)
#we reduce the distance between points
#distort the distances
View(datalt_bio_t)
plot(datalt_bio_t)  

#step 4
#Calculate the matrix of association using the coefficient
#matdist<-as.dist(datalt)

matdist_bio_t<-vegdist(datalt_bio_t,method='bray')

#Step 5: Apply clustering method and genearte the dendrogram
# group average clustering method
#View(matdist)

LS<-hclust(matdist_bio_t,method = 'single')
LC<-hclust(matdist_bio_t,method='complete')
GA<-hclust(matdist_bio_t,method='average')

par(mfrow=c(1,3))
plot(LS,ylab='Bray Curtis Method',
     xlab='stations',main='Single linkage')
abline(b=0,a=0.25,col='red')

plot(LC,ylab='Bray Curtis Maximum distance',
     xlab='stations',main='Clusters Analysis - Complete linkage')
abline(b=0,a=0.98,col='red')

plot(as.dendrogram(LC), horiz = T)

plot(GA,ylab='Bray Curtis Method',
     xlab='stations',main='Group Average linkage')
abline(b=0,a=0.25,col='red')

groupe<-cutree(LC,4);groupe

## most representative species
library(labdsv)
data7<-data1[,colSums(data1)!=0] #delete the columns with 0 values
indval(data7,groupe)->IV;IV

IV$relfrq#faithfulness 
IV$relabu #specificity
IV$indval #value of indval from 0 to 1

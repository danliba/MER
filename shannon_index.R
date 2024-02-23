rm(list = ls())
setwd('D:/Maestria/MER/SOTON/clases/deepSeaEco/assigment2')


dir()
library(readxl)
library(vegan)
library(ade4)

data0<-read_excel("ROVtransect_DATA.xlsx", sheet = "abundance")
View(data0)
data1<-data0[,5:17]
data1[is.na(data1)]<-0
View(data1)

shannon<-0 #cummulative shannnon 
for (ii in 1:dim(data1)[1]){
  shannon[ii]<-diversity(data1[ii,],'shannon')
}

##PLOT
windows(width = 12, height = 6)  # Set the width and height of the figure
par(mfrow = c(1, 1))  # Set the layout to a single plot

plot(data0$`Distance along transect (m)`,shannon,pch=15,col='black',xlab='years',ylab='',
     xlim=c(1,303),axes=FALSE,type='b',main='Diversidad Funcional Iquitos')
axis(2, ylim=c(0,1),col="black",las=1)  
axis(side=1,at=data0$time,labels=data0$time)

box()
mtext("Shannon Index",side=2,line=2.5)
grid(lty = "solid")  # Major grid lines
#grids
minor_ticks <- seq(1, 303, by = 1)
abline(v = minor_ticks, col = "gray", lty = "dotted", lwd = 0.5)

par(new = TRUE)
#plot(data0$time,shannon/max((shannon)), pch=16,  xlab="", ylab="", ylim=c(0.7,1), 
#     axes=FALSE, type="b", col="red")
plot(data0$`Distance along transect (m)`,shannon/log(dim(data1)[2]), pch=16,  xlab="", ylab="", ylim=c(0.1,0.8), 
     axes=FALSE, type="b", col="red")
mtext("Shannon Evennes",side=4,col="red",line=4) 
axis(4, col="red",col.axis="red",las=0)

#legend
legend("bottomright", legend = c("Shannon Index", "Shannon Evenness"), 
       pch = c(15, 16), col = c("black", "red"), lty = c(1, 1))


shannon_eveness<-shannon/log(dim(data1)[2])

shannon_iquitos<-data.frame(cbind(shannon,shannon_eveness))
library("writexl")
write_xlsx(shannon_iquitos,"Shannon.xlsx")

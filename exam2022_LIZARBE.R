##############################################################
####################### Exam 2022 ############################
##############################################################

rm(list = ls())

setwd('D:/Maestria/MER/Bordeaux/clases/Environmental_data_treatment_and_modelling/exam_practice/exam2022/')
dir()

data1<-read.table("DATASET.txt",header=TRUE,row.names=1);data1
View(data1)

library(vegan)
library(ade4)
#library(labdsv)
#Question 1: Using the group average clustering method and a distance of 3.8, 
#how many groups of stations can you identify from physico-chemical variables? 
#You will ensure that you follow the different steps to perform correctly your classification.
head(data1)

#first we separte the physico-chemical data 
data_pc<-data1[,8:16]

View(data_pc)
boxplot(data_pc)

# here we reduce the distance
datalt<-log(data_pc+1)
#we reduce the distance between points
#distort the distances
View(datalt)
plot(datalt)

#step 4
#Calculate the matrix of association using the coefficient
library(vegan)
matdist<-vegdist(datalt,method='euclidian')

#Step 5: Apply clustering method and genearte the dendrogram
# group average clustering method

GA<-hclust(matdist,method='average')

plot(GA,ylab='Euclidian distance',
     xlab='stations',main='Group Average linkage')
abline(b=0,a=3,col='red')
groupe<-cutree(GA,3.8);groupe

#Conclusion: For the physico-chemical point of view we can observe that there are 3 clusters 
#The first grouping A_a A_b A_c A_d A_e A_f B_d C_a C_b C_c C_d, Second cluster: B_c
#Third cluster: D_a D_b D_c D_d D_e 

############################################################################
# Question 2: Which marshes are the groups associated with? Do the NO3 values 
#depend on the groups at the level of significance alfa = 5%? 
#############################################################################
View(data1)
#Group 1: is associated with Marsh A and C 
#Group 2: is associated wit Marsh  B
#Group 3: is associated with Marsh D

#Do the NO3 values depend on the groups at the level of significance alfa = 5%? 
#Now we will test if NO3 is similar in all the groups

#We will prove normality first
NO3<-data1[,13]

#we extract NO3 and do somple plotting
data1[data1$March=='A',13]->A; shapiro.test(A) #
data1[data1$March=='B',13]->B; shapiro.test(B) #p<0.05
data1[data1$March=='C',13]->C; shapiro.test(C) #
data1[data1$March=='D',13]->D; shapiro.test(D) #p<0.05

boxplot(A,B,C,D)
shapiro.test(NO3) #p <0.05 is not normal

#we will perform a non-parametric one way anova --> Kurskall wallis

kruskal.test(data1$NO3~data1$March,data=data1)
#Ho: NO3 is constant and doesnt varie along the sites
#Ha: NO3 varies along the sites
#p<0.05 , we reject the null hypothesis and accept the alternative
#At 0.05 significance level, we conclude that the mean NO3 along the 4 March is not the same
#This means that NO3 values depend on the group

########################################################################################
#Question 3: Perform a principal component analysis from physico-chemical variables. 
#By basing your analysis on the 2 first factorial axes, can you identify the variables 
#that explain these groups of stations?
########################################################################################
aa<-scale(data_pc[,1])
#First we rearrange the matrix for the PCA
View(data_pc)

#Now we do the previous steps for PCA
###1: We probe the multinormality of the data ############

windows()
par(mfrow=c(3,3))

#plot(density.default(aa),main="depth");shapiro.test(data_pc[,1]) #pval >0.05
plot(density.default(scale(data_pc[,1])),main="depth");shapiro.test(data_pc[,1]) #pval >0.05
plot(density.default(scale(data_pc[,2])),main="Irradiance");shapiro.test(data_pc[,2]) #pval >0.05
plot(density.default(scale(data_pc[,3])),main="Optical depth");shapiro.test(data_pc[,3]) #pval >0.05
plot(density.default(scale(data_pc[,4])),main="Temp");shapiro.test(data_pc[,4]) #pval >0.05
plot(density.default(scale(data_pc[,5])),main="NO2");shapiro.test(data_pc[,5]) #pval >0.05
plot(density.default(scale(data_pc[,6])),main="NO3");shapiro.test(data_pc[,6]) #pval < 0.05
plot(density.default(scale(data_pc[,7])),main="PO4");shapiro.test(data_pc[,7]) #pval < 0.05
plot(density.default(scale(data_pc[,8])),main="NP");shapiro.test(data_pc[,8]) #pval < 0.05
plot(density.default(scale(data_pc[,9])),main="Turb");shapiro.test(data_pc[,9]) #pval < 0.05

#They are not normal but we will continue anyways because the teacher said this is 
#for educational purposes

####Step 2: Let's prove monotonic relation
# Assumption 2 : linear or monotonic relationship between variables
plot(data_pc)
# there is no linear relationship between all of the variables
# we will still apply the PCA

##Now we do the STEPS for the PCA 
######### STEP 1 - Transformation of DATA #########
# Check the range of variation of parameters
windows()
par(mfrow=c(1,2)) 
boxplot(data_pc, horizontal=F)
# as the variables have high diversity ranges we will now standardize the data
data5<-scale(data_pc)
boxplot(data5, horizontal=F)

#all scaled

######### STEP 2&3 - Computation of the covariance matrix and spectral decomposition#########
library(ade4)
# pcaBA<-dudi.pca(BA3) 
# This option requires to provide interactively an answer
# about the number of kept axes
# Another option :
pcaMER<-dudi.pca(data5,scannf=FALSE,nf=3)
pcaMER$eig

eigpercent<-pcaMER$eig/sum(pcaMER$eig)*100;eigpercent
barplot(eigpercent,names.arg=c("F1","F2","F3","F4","F5","F6","F7","F8","F9"))
# The three first factorial axes explain approximately 72% of the 
# total variance associated with the dataset. 
# This means that we can reduce the analysis to this 3 factorial axes

#graphical representations of the PCA
par(mfrow=c(1,2)) 
s.corcircle(pcaMER$co, xax=1, yax=2) # Variable space in the F1-F2 plan
s.label(pcaMER$li, xax=1, yax=2) # Observation space in the F1-F2 plan

## Can you identify the variables that explain these groups of stations?
#### Conclusion ########
#We observe that NO3, NP (nitrate phosphate ratio), Turb and positively correlated and PO4 
# negatively correlated with Factor 1. 
#Factor 1 represents the turbidity concentration and the gradient of Phosphate-No3, due to 
#agricultural activities in the area
#Factor 2: is related with optical depth and temperature
#Stations in group 1 ( Marsh A) are positively related to Factor 2
#Stations in group 3 (Marsh D) are positively related to factor 1, this means that that
#stations Da and Dc could be explained by the proximity of urban and agricultural spaces

###############################################################################
#Question 4: Which species are specific to these groups for an indicator value of
#0.8? Do we obtain the same result from a classification carried out on the abundance data? 

## Phytoplankton data
phyto<-data1[,17:50]

View(phyto)
boxplot(phyto)

# here we reduce the distance
datalt<-log(phyto+1)
#we reduce the distance between points
#distort the distances
View(datalt)

#step 4
#Calculate the matrix of association using the coefficient
library(vegan)
matdist<-vegdist(datalt,method='bray')

#Step 5: Apply clustering method and genearte the dendrogram
# group average clustering method

LS<-hclust(matdist,method = 'single')
LC<-hclust(matdist,method='complete')
GA<-hclust(matdist,method='average')

par(mfrow=c(1,3))
plot(LS,ylab='Dissimilarity of Bray-Curtis',
     xlab='stations',main='HCA-Single linkage')
abline(b=0,a=0.8,col='red')

plot(LC,ylab='Dissimilarity of Bray-Curtis',
     xlab='stations',main='HCA-Complete linkage')
abline(b=0,a=0.8,col='red')

plot(GA,ylab='Dissimilarity of Bray-Curtis',
     xlab='stations',main='CHA-Group Average linkage')
abline(b=0,a=0.8,col='red')

groupe<-cutree(GA,3);groupe

#We will use complete linkage

library(labdsv)

data7<-phyto[,colSums(phyto)!=0] #delete the columns with 0 values
#We need to remove species absent for any stations!
indval(data7,groupe)->IV;IV
IV$relfrq#faithfulness 
IV$relabu #specificity
IV$indval #value of indval from 0 to 1


##Group 1:cycl
#Group 2:
#Group 3: Navi

#are representative
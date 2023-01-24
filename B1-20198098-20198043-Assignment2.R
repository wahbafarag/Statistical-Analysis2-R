# Wahby Mustafa Farag      20198098
# Reem Ali Abdelrahman     20198043

#install.packages("corrplot")
#BiocManager::install("preprocessCore")

library(antiProfilesData)
library(corrplot)
library(preprocessCore)

dt <- antiProfilesData::apColonData

fdata<-fData(dt)
pdata<-pData(dt)
edata<-as.data.frame(exprs(dt))

#edata
#pdata
#fdata

#1.a Show the type of each column
sapply(pdata, class)

sapply(edata, class)
for (i in colnames(edata)){
  print(class(edata[,i]))
}


# 1.b Show column names and rows name

colnames(pdata)
rownames(pdata)
# colnames(edata) same as rownames(pdata)
rownames(edata)


# 1.c Calculate summary of each column

sapply(pdata, summary)
sapply(edata, summary)


# 1.d Show frequency of categorical data, taking into the consideration, NA values frequency if any. 

table(pdata$Tissue, useNA = "ifany")
table(pdata$SubType, useNA = "ifany")
table(pdata$Status, useNA = "ifany")
#table(pdata$ClinicalGroup, useNA = "ifany")
table(length(pdata$ClinicalGroup[is.na(pdata$ClinicalGroup)]))

# 1.e Calculate the correlation and covariance between the first 10 columns only of our data set and draw full correlation matrix

data <- edata[,1:10]
cov(data,y=NULL ,use="everything")
corl <- cor(data,y=NULL,use="everything")
corrplot(corl, type = "upper",tl.col = "black",order = "hclust", tl.srt = 45)
colPalette<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = corl, col = colPalette, symm = TRUE)



# 1.f For both genes: GSM95478,GSM95473 show the plot with a line of their relation.

x <- edata[,"GSM95478"]
y <- edata[,"GSM95473"]
plot(x, y,pch=1 ,main = "Line and Relation Plot",xlab = "GSM95478",ylab = "GSM95473")
abline(lm(y ~ x), col = "purple",lwd=4)

#-----------------------------------------------------------------------------

#Question Two :

#Using PCA and SVD, Prove by plotting and values that both can return the same result by
#suitable normalization.


nor = normalize.quantiles(as.matrix(edata))
pc1 = prcomp(nor)
cntEdata = t(t(nor) - colMeans(nor))
svd_1 = svd(cntEdata)

plot(pc1$rotation,svd_1$v,col="black")
plot(pc1$rotation[,1],svd1$v[,1],col="black")

par(mfrow = c(1,2))
plot(pc1$rotation,col= "black" ,main = "PCA Values")
plot(svd_1$v, col = "black",main = "SVD Values")
par(mfrow = c(1,1))


#-----------------------------------------------------------------------------

#Question 3:
#  256 visual artists were surveyed to find out their zodiac sign. The results were: Aries (29),
#Taurus (24), Gemini (22), Cancer (19), Leo (21), Virgo (18), Libra (19), Scorpio (20), Sagittarius
#(23), Capricorn (18), Aquarius (20), Pisces (23).


#3.1) Test the hypothesis that zodiac signs are evenly distributed across visual artists:


res <- as.factor(c(rep("Aries",29), rep("Taurus", 24), rep("Gemini", 22),rep("Cancer", 19), rep("Leo", 21), rep("Virgo",18),rep("Libra", 19),rep("Scorpio", 20), rep("Sagittarius", 23),rep("Capricorn", 18),rep("Aquarius", 20), rep("Pisces",23)))
hyp=c(rep(1,12))
hyp=hyp/sum(hyp)
chisq.test(table(res),p=hyp)
 

#3.2) Explicitly mention your H1 and Ho assumption.

# p-value is greater than 0.05 so we accept it
# h0 is unif dist over zodiac
# h1 is null uniformly dist over zodiac (1-h0)


#-----------------------------------------------------------------------------

# Question 4:
# Plot hierarchical clusters on our first 10 columns of edata and apply the kmeans to all the edata
#columns and show the centroid of the result.

normalizeData = log2(edata + 1) 
omitEdata <- na.omit(normalizeData)

ecDist = dist(t(omitEdata[,1:10]))
hiClust = hclust(ecDist)
plot(hiClust,hang = -1,main = "Hierarchical Clusters",ylab = "Height",xlab = "Distance")


clust = kmeans(omitEdata,centers=3)
#names(clust)
#dim(clust$centers)
clust1Cent = clust$centers[,1]
clust2Cent = clust$centers[,2]
clust3Cent = clust$centers[,2]
clust4Cent = clust$centers[,4]
table(clust$cluster)

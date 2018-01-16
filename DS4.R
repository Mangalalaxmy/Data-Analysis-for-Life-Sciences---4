# Exercises for DS4
# Distance Exercises
#1
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
str(e)
head(tissue)
#How many biological replicates for hippocampus?
table(tissue)

#2
#What is the distance between samples 3 and 45?
x = e[,3]
y = e[,45]
sqrt(crossprod(x-y))

#3
#What is the distance between gene 210486_at and 200805_at?
x = e["210486_at",]
y = e["200805_at",]
sqrt(sum((x-y)^2))

#4
#If I run the command (don't run it!)
d = as.matrix(dist( e))
#How many cells (number of rows times number of columns) would this matrix have?
nrow(e)^2

#5
#Compute the distance between all pairs of samples:
d = dist(t(e))
#Read the help file for dist.
#How many distances are stored in d? (Hint: What is the length of d)?
length(d)

#6
#Why is the answer above not ncol(e)^2?
#Because we take advantage of symmetry: only lower triangular matrix is stored thus only ncol(e)*(ncol(e)-1)/2 values. correct

#Projection Exercises
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)

#1
#Suppose you want to make an MA plot of the first two samples y = geneExpression[,1:2]. Which of the following projections of y gives us new coordinates such that column 2 versus column 1 is an MA plot?
y = geneExpression[,1:2]
plot((y[,1]+y[,2])/2, y[,1]-y[,2])
#y(1,1 : 1,-1)

#2
#Say Y is M×N, in the SVD Y=UDV??? which of the following is NOT correct?
# D are the coordinates of the projection U???Y

#SVD Exercises
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
s = svd(e)
signflips = sample(c(-1,1),ncol(e),replace=TRUE)
signflips
newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
all.equal( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

#1
#Compute the SVD of e
s = svd(e)
#Now compute the mean of each row:
m = rowMeans(e)
#What is the correlation between the first column of U and m?
cor(s$u[,1], m)

#2
#In the above question we saw how the first column relates to the mean of the rows of e. Note that if we change these means, the distances between columns do not change. Here some R code showing how changing the means does not change the distances:
newmeans = rnorm(nrow(e)) ##random values we will add to create new means
newe = e+newmeans ##we change the means
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45]))
#So we might as well make the mean of each row 0 since it does not help us approximate the column distances. We will define y as the detrended e and recompute the SVD:
y = e - rowMeans(e)
s = svd(y)
#We showed that UDV??? is equal to y up to numerical error
resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))
#The above can be made more efficient in two ways. First, using the crossprod and second not creating a diagonal matrix. Note that in R we can multiply a matrix x by vector a. The result is a matrix with row i equal to x[i,]*a[i]. Here is an example to illustrate this.
x=matrix(rep(c(1,2),each=5),5,2)
x
x*c(1:5)
#Note that the above code is actually equivalent to:
sweep(x,1,1:5,"*")
#This means that we don't have to convert s$d into a matrix to obtain DV???.
#Which of the following gives us the same as diag(s$d)%*%t(s$v)?
#s$d * t(s$v)

#3
#If we define vd = t(s$d * t(s$v)) then which of the following is not the same UDV??? :
s$u %*% s$d * t(s$v)

#4
#Let z = s$d * t(s$v). We showed derivation demonstrating that because U is orthogonal the distance between e[,3] and e[,45] is the same as the distance between y[,3] and y[,45] which is the same as z[,3] and z[,45]
z = s$d * t(s$v)
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))
#Note that the columns z have 189 entries, compared to 22,215 for e.
#What is the difference (in absolute value) between the actual distance sqrt(crossprod(e[,3]-e[,45])) and the approximation using only two dimension of z
actdist = sqrt(crossprod(e[,3]-e[,45]))
approxdist = sqrt(crossprod(z[1:2,3]-z[1:2,45]))
abs(actdist - approxdist)

#5
#What is the minium number of dimensions we need to use for the approximation in SVD Exercises #4 to be within 10% or less?
ds = 1:189
actdist = sqrt(crossprod(e[,3]-e[,45]))
approxdist = sapply(ds,function(d){
  sqrt(crossprod(z[1:d,3,drop=FALSE]-z[1:d,45,drop=FALSE] )) 
})
percentdiff = 100 * abs((actdist - approxdist)/actdist)
min(ds[which(percentdiff < 10)])
plot(ds,percentdiff)
abline(h=10,v=7)

#6
#Compute distances between sample 3 and all other samples:
distances = sqrt(apply(e[,-3]-e[,3],2,crossprod))
#Recompute this distance using the 2 dimensional approximation.
#What is the Spearman correlation between this approximate distance and the actual distance?
approxdistances = sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
cor(distances, approxdistances, method = "spearman")
plot(distances,approxdistances)

#MDS Exercises
library(tissuesGeneExpression)
data(tissuesGeneExpression)
#In these exercise we will demonstrate the relantionship between the SVD and the output of mdscale, the function in R that performs MDS.
#Using the z we computed in SVD Exercises #4
y = e - rowMeans(e)
s = svd(y)
z = s$d * t(s$v)
#we can make an mds plot
library(rafalib)
ftissue = factor(tissue)
mypar(1,1)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
#Now run the function cmdscale on the original data
d = dist(t(e))
mds = cmdscale(d)

#1
#What is the correlation between the first row of z and the first column in mds?
cor(z[1,], mds[,1])

#2
#What is the correlation between the second row of z and the second column od mds?
cor(z[2,], mds[,2])

#3
#Note that the mds plot is not the same:
library(rafalib)
ftissue = factor(tissue)
mypar(1,2)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
plot(mds[,1],mds[,2],col=as.numeric(ftissue))
#Given the answer to MDS Exercises #1, what do we have to do to z[1,] and z[2,] to get a practically identical plot?
#multiply z[1,] and z[2,] by -1

#4
#Load the following dataset
library(GSE5859Subset)
data(GSE5859Subset)
#Compute the svd and compute z
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)
#Which dimension of z most correlates with the outcome sampleInfo$group?
which.max(cor(sampleInfo$g,t(z)))
plot(z[1,])
plot(sampleInfo$group)

#5
#Load the following dataset
library(GSE5859Subset)
data(GSE5859Subset)
#Compute the svd and compute z
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)
#What is this max correlation?
max(cor(sampleInfo$g,t(z)))

#6
#Load the following dataset
library(GSE5859Subset)
data(GSE5859Subset)
#Compute the svd and compute z
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)
#Which dimension of z has the second highest correlates with the outcome sampleInfo$group?
which.max(cor(sampleInfo$g,t(z))[-1]) + 1
sort(cor(sampleInfo$g,t(z)))
tail(sort(which(a<1)))

#7
#Note these measurements were made during two months:
sampleInfo$date
#We can extract the month this way:
month = format( sampleInfo$date, "%m")
month = factor( month)
#Which dimension of z has the highest correlates with the outcome month
which.max(cor(as.numeric(month), t(z)))
#What is this correlation?
max(cor(as.numeric(month), t(z)))

#8
#In MDS Exercises #7 we saw that that one of the dimensions was highly correlated to the sampleInfo$group. Now take the 5th column of U and stratify by the gene chromosome. Remove chrUn and make a boxplot of the values of U6 stratified by chromosome.
#Which chromosome looks different from the rest? Copy and paste the name as it appears in geneAnnotation
mypar()
U = s$u
strata = split(U[,6],geneAnnotation$CHR)
strata = strata[which(names(strata)!="chrUn")]
boxplot(strata)
boxplot(strata,range=0)
boxplot(strata,range=0,ylim=c(-0.025,0.025))
medians = sapply(strata,median)
names(strata)[ which.max(abs(medians)) ]

#Hierarchical Clustering Exercises
#1
#Create a random matrix with no correlation in the following way:
set.seed(1)
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n
#Run hierarchical clustering on this data with the hclust function with default parameters to cluster the columns. Create a dendrogram.
#From the dendrogram which pairs of samples are the furthest away from each other?
d = dist(t(x))
ho = hclust(d)
plot(ho)

#2
#Set the seed at 1, set.seed(1) and replicate the creation of this matrix 100 times
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
#then perform hierarchical clustering as in the solution to question 2.4.1 and find the number of clusters if you use cutree at height 143. Note that this number is a random variable.
#Based on the Monte Carlo simulation, what is the standard error of this random variable?
set.seed(1)
hc = replicate(100, {x = matrix(rnorm(m*n),m,n)
d=dist(t(x))
ho=hclust(d)
length(unique(cutree(ho, h=143)))})
plot(table(hc))
popsd(hc)

#K-Means Exercises
#1
#Run kmeans with 4 centers for the blood RNA data:
library(GSE5859Subset)
data(GSE5859Subset)
#Set the seed to 10, set.seed(10) right before running kmeans with 5 centers.
#Explore the relationship of clusters and information in sampleInfo. Which of th following best described what you find:
mds=cmdscale(dist(t(geneExpression)))
set.seed(10)
result=kmeans(t(geneExpression),5)
mypar(1,1)
plot(mds,bg=result$cl,pch=21)
table(sampleInfo$group,result$cluster)
table(sampleInfo$date,result$cluster)
##looks better if we re-order:
table(sampleInfo$date,result$cluster)[,c(4,1,5,3,2)]

#Heat map Exercises
#1
#Load the data:
library(GSE5859Subset)
data(GSE5859Subset)
#Pick the 25 genes with the highest across sample variance. This function might help
install.packages("matrixStats")
library(matrixStats)
?rowMads ##we use mads due to a outlier sample
#While a heatmap function is included in R, we recommend the heatmap.2 function from the gplots package on CRAN because it is a bit more customized. For example, it stretches to fill the window.
library(gplots)
#Use heatmap.2 to make a heatmap showing the sampleInfo$group with color, the date as labels, the rows labelled with chromosome, and scaling the rows.
#What do we learn from this heatmap?
library(RColorBrewer)
library(rafalib)
variances = rowMads(geneExpression)
selected = order(variances, decreasing = TRUE)[1:25]
hmcol = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
gcols = brewer.pal(3,"Dark2")[sampleInfo$g+1]
labcol = gsub("2005-","",sampleInfo$date)
heatmap.2(geneExpression[selected,], labCol=labcol, trace="none", scale="row", labRow=geneAnnotation$CHR[selected], ColSideColors=gcols, Col=hmcol)
 #OR
##load libraries
library(rafalib)
library(gplots)
library(matrixStats)
library(RColorBrewer)
##make colors
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(25)
gcol=brewer.pal(3,"Dark2")
gcol=gcol[sampleInfo$g+1]

##make lables: remove 2005 since it's common to all
labcol= gsub("2005-","",sampleInfo$date)  

##pick highly variable genes:
sds =rowMads(geneExpression)
ind = order(sds,decreasing=TRUE)[1:25]

## make heatmap
heatmap.2(geneExpression[ind,],
          col=cols,
          trace="none",
          scale="row",
          labRow=geneAnnotation$CHR[ind],
          labCol=labcol,
          ColSideColors=gcol,
          key=FALSE)
#A group of chrY genes are higher in group 0 and appear to drive the clustering. Within those clusters there appears to be clustering by month.

#2
#Create a large data set of random data that is completely independent of sampleInfo$group like this:
set.seed(17)
m = nrow(geneExpression)
n = ncol(geneExpression)
x = matrix(rnorm(m*n),m,n)
g = factor(sampleInfo$g )
#Create two heatmaps with these data. Show the group g either with labels or colors.
#1. Taking the 50 genes with smallest p-values obtained with rowttests
#2. Taking the 50 genes with largest standard deviations.
#Which of the following statements is true:
mypar(1,2)
t = rowttests(x,g)
ps = order(t$p.value)[1:50]
cols = colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
heatmap.2(x[ps,], labCol = g, col=cols, trace="none")
s = rowSds(x)
sds = order(-s,decreasing = TRUE)[1:50]
heatmap.2(x[sds,], labCol = g, col=cols, trace="none")

#Conditional Expectation Exercises
#1
#Generate some random data to imitate heights for men (0) and women (1):
n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
#Treating the data generated above as the population, if we know someone is 176 cm tall, what it the probability that this person is a woman: 
#Pr(Y=1|X=176)=E(Y|X=176)?
mean(y[x==176])

#2
#Now make a plot of E(Y|X=x) for x=seq(160,178) using the data generated in Conditional Expectation Exercises #1.
#Suppose for each height x you predict 1 (female) if Pr(Y|X=x)>0.5 and 0 (male) otherwise. What is the largest height for which you predict female ?
xops = seq(160,178)
p = sapply(xops,function(x0) mean(y[x==x0]))
plot(xops, p)
abline(h=0.5)
abline(v=168)

#Smoothing Exercises
n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]
#Set the seed at 5, set.seed(5) and take a random sample of 250 individuals from the population like this:
set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]

#1
#Use loess to estimate f(x)=E(Y|X=x) using the default parameters. What is the predicted f(168)?
fit=loess(Y~X)
predict(fit,newdata=data.frame(X=168))
##Here is a plot
xs = seq(160,178)
Pr =sapply(xs,function(x0) mean(Y[X==x0]))
plot(xs,Pr)
fitted=predict(fit,newdata=data.frame(X=xs))
lines(xs,fitted)

#2
#The loess estimate above is a random variable thus we should compute its standard error. Use Monte Carlo simulation to compute the standard error of your estimate of f(168).
#Set the seed to 5, set.seed(5) and perform 1000 simulation of the computations performed in question 2.7.1. Report the the SE of the loess based estimate.
set.seed(5)
B=1000
reps = replicate(B,{N = 250
   ind = sample(length(y),N)
   Y = y[ind]
   X = x[ind]
   fit=loess(Y~X)
   est = predict(fit,newdata=data.frame(X=168))
   return (est)
   })
popsd(reps)

#KNN and CrossValidation exercises
library(GSE5859Subset)
data(GSE5859Subset)
#And define the outcome and predictors. To make the problem more difficult we will only consider autosomal genes:
y = factor(sampleInfo$group)
X = t(geneExpression)
out = which(geneAnnotation$CHR%in%c("chrX","chrY"))
X = X[,-out]
#Note, you will also need to load the following:
library(caret)

#1
#Set the seed to 1 set.seed(1) then use the createFolds function in the caret package to create 10 folds of y.
#What is the 2nd entry in the fold 3?
set.seed(1)
ids = createFolds(y, k=10)
ids[[3]][2]
sapply(ids, function(i) table(y[i]))

#2
#For the following questions we are going to use kNN. We are going to consider a smaller set of predictors by filtering genes using t-tests. Specifically, we will perform a t-test and select the m genes with the smallest p-values.
#Let m=8 and k=5 and train kNN by leaving out the second fold idx[[2]]
#How many mistakes do we make on the test set? Remember it is indispensable that you perform the ttest on the training data.
library(class)
library(genefilter)
m=8
k=5
ind = idx[[2]]
pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
ind2 = order(pvals)[1:m]
predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
sum(predict!=y[ind])

#3
#Now run the code for kNN and Cross Validation Exercises #2 for all 10 folds and keep track of the errors. What is our error rate (number of errors divided by number of predictions) ?
library(class)
library(genefilter)
m=8
k=5
result = sapply(idx,function(ind){
  pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
  ind2 = order(pvals)[1:m]
  predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
  sum(predict!=y[ind])
})
sum(result)/length(y)

#4
#Now we are going to select the best values of k and m. Use the expand grid function to try out the following values:
ms=2^c(1:11)
ks=seq(1,9,2)
params = expand.grid(k=ks,m=ms)
#Now use apply or a loop to obtain error rates for each of these pairs of parameters. Which pair of parameters minimizes the error rate?
#k=
#m=
errors = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})
params[which.min(errors),]
##make a plot and confirm its just one min:
errors = matrix(errors,5,11)
library(rafalib)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks))

#5
#Repeat question kNN and Cross Validation Exercises #4 but now perform the t-test filtering before the cross validation. Note how this biases the entire result and gives us much lower estimated error rates.
#What is the minimum error rate?
pvals = rowttests(t(X),factor(y))$p.val
errors = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})
min(errors)
##make a plot and compare to previous question
errors = matrix(errors,5,11)
library(rafalib)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks))

#6
#Repeat the cross-validation we performed in question kNN and Cross Validation Exercises #4 but now instead of defining y as sampleInfo$group use:
y = factor(as.numeric(format( sampleInfo$date, "%m")=="06"))
#What is the minimum error rate now?
errors = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})
min(errors)
##make a plot and confirm its just one min:
errors = matrix(errors,5,11)
library(rafalib)
mypar(1,1)
matplot(ms,t(errors),type="l",log="x")
legend("topright",as.character(ks),lty=seq_along(ks),col=seq_along(ks))

#Confounding Exercises
#Install the latest version of the dagdata packge from the genomicsclass github repository. Load the admissions data from the dagdata package:
library(dagdata)
data(admissions)
#Familiarize yourself with this table:
print( admissions )
#You can also obtain this data directly from here.
#To install the library(GSE5859):
library(devtools)
install_github("genomicsclass/GSE5859")

#1
#Let's compute the proportion of men who were accepted:
Mindex = which(admissions$Gender==1)
Maccepted= sum(admissions$Number[Mindex] * admissions$Percent[Mindex]/100)
Mapplied = sum(admissions$Number[Mindex])
Maccepted/Mapplied
#What is the proportion of women that were accepted?
Windex = which(admissions$Gender==0)
Waccepted= sum(admissions$Number[Windex] * admissions$Percent[Windex]/100)
Wapplied = sum(admissions$Number[Windex])
Waccepted/Wapplied

#2
#Now that we have observed different acceptance rates between genders, test for the significance of this result.
#If you perform an independence test, what is the p-value?
#Hint: create a table that has the totals for accepted and not-accepted by gender then use chisq.test
#This difference actually led to a lawsuit.
#Now notice that looking at the data by major, the differences disappear.
index = admissions$Gender==1
men = admissions[index,]
women = admissions[!index,]
print( data.frame( major=admissions[1:6,1],men=men[,3], women=women[,3]) )
#How can this be? This is referred to as Simpson's Paradox.
#In the following questions we will try to decipher why this is happening.
Ms_accepted = admissions$Number[Mindex] * admissions$Percent[Mindex]/100
Ws_accepted = admissions$Number[Windex] * admissions$Percent[Windex]/100
accepted = cbind(Ms_accepted, Ws_accepted)
pval = chisq.test(accepted)$p.value

#3
#We can quantify how "hard" a major is using the percent of students that were accepted. Compute the percent that were accepted (regardless of gender) to each major and call this vector H
#Which is the hardest major? (enter a letter)
Femp = subset(admissions$Percent, admissions$Gender==0)
Memp = subset(admissions$Percent, admissions$Gender==1)
majorp = data.frame(major=admissions[1:6,1], men=men[,3], women=women[,3])
majorp$H = (majorp$men+majorp$women)/2
table(majorp$major,majorp$H)

major = admissions[1:6,1]
men = admissions[1:6,]
women =admissions[7:12,]
H = (men$Number*men$Percent/100 + women$Number*women$Percent/100) / (men$Number+women$Number)
major[which.min(H)]

#4
#What proportion gets in for the major from Confounding Exercises #3?
#0.065

#5
#For men, what is the correlation between the number of applications across majors and H
majorp$Mn=subset(admissions$Number, admissions$Gender==1)
majorp$Wn=subset(admissions$Number, admissions$Gender==0)
majorp$Totaln = majorp$Mn + majorp$Wn
cor(majorp$Mn,H)

#6
#For women, what is the correlation between the number of applications across majors and H
cor(majorp$Wn,H)

#7
#Given the answers to Confounding Exercises #5 and #6, which best explains the differences in admission percentages when we combine majors?
#There is confounding between gender and preference for "hard" majors: females are more likely to apply to harder majors.

#Confounding Exercises in Genomics
#Load the data for this gene expression dataset:
library(Biobase)
library(GSE5859)
data(GSE5859)
#Note that this is the original dataset from which we selected the subset used in GSE5859Subset.  You can obtain it from the genomicsclass GitHub repository.
#We can extract the gene expression data and sample information table using the Bio conductor functions exprs and pData like this:
geneExpression = exprs(e)
sampleInfo = pData(e)

#1
#Familiarize yourself with the sampleInfo table. Note that some samples were processed at different times. This is an extraneous variable and should not affect the values in geneExpression. However, as we have seen in previous analyses it does appear to have an effect so we will explore this here.
#You can extract the year from each date like this:
year = format(sampleInfo$date,"%y")
#Note there are
length( unique(year) )
#unique years for which we have data.
#For how many of these years do we have more than one ethnicity represented?
tab = table(year,sampleInfo$ethnicity)
print(tab)
x = rowSums( tab != 0)
sum( x >= 2)

#2
#Repeat the above exercise but now instead of year consider the month as well. Specifically, instead of the year variable defined above use:
month.year = format(sampleInfo$date,"%m%y")
#For what proportion of these month.year values do we have more than one ethnicity represented?
tab = table(month.year,sampleInfo$ethnicity)
print(tab)
x = rowSums( tab != 0)
mean( x >= 2)
1/24

#3
#Perform a t-test (use rowttests) comparing CEU samples processed in 2002 to those processed in 2003. Then use the qvalue package to obtain q-values for each gene.
#How many genes have q-values < 0.05?
library(genefilter)
library(qvalue)
year = factor( format(sampleInfo$date,"%y") )
index = which(year%in% c("02","03") & sampleInfo$ethnicity=="CEU")
year = droplevels(year[index])
pval = rowttests(geneExpression[ ,index], year)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)
print( qvalue(pval)$pi0 )

#4
#Now perform a t-test (use rowttests) comparing CEU samples processed in 2003 to CEU samples processed in 2004. Then use the qvalue package to obtain q-values for each gene.
#How many genes have q-values < 0.05?
year = factor( format(sampleInfo$date,"%y") )
index = which(year%in% c("03","04") & sampleInfo$ethnicity=="CEU")
year = droplevels(year[index])
pval = rowttests(geneExpression[ ,index], year)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

#5
#Now we are going to compare ethnicities as was done in the original publication in which these data were first presented. Use the rowttests function to compare the ASN population to the CEU population. Once again, use the qvalue function to obtain q-values.
#How many genes have q-values < 0.05?
g <- factor(sampleInfo$ethnicity) 
index = which(sampleInfo$ethnicity%in% c("ASN", "CEU"))
g = droplevels((g[index]))
pval = rowttests(geneExpression[ ,index], g)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

#6
#Note that over 80% of genes are called differentially expressed between ethnic groups. However, due to the confounding with processing date, we need to confirm these differences are actually due to ethnicity. This will not be easy due to the almost perfect confounding. However, above we noted that two groups were represented in 2005. Just like we stratified by majors to remove the "major effect" in our admissions example, here we can stratify by year and perform a t-test comparing ASN and CEU, but only for samples processed in 2005.
#How many genes have q-values < 0.05?
year = factor( format(sampleInfo$date,"%y") )
g <- factor(sampleInfo$ethnicity) 
index = which(sampleInfo$ethnicity%in% c("ASN", "CEU") & year == "05")
g = droplevels((g[index]))
pval = rowttests(geneExpression[ ,index], g)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

#7
#To provide a more balanced comparison we repeat the analysis but now taking 3 random CEU samples from 2002. Repeat the analysis above but comparing the ASN from 2005 to three random CEU samples from 2002. Set the seed at 3, set.seed(3)
#How many genes have q-values < 0.05?
library(qvalue)
library(genefilter)
year = factor( format(sampleInfo$date,"%y") )
index1 = which(sampleInfo$ethnicity=="ASN" & year=="05")
set.seed(3)
index2 = sample( which(sampleInfo$ethnicity == "CEU" & year=="02"), 3)
index = c( index1, index2)
g = droplevels(sampleInfo$ethnicity[index])
pval = rowttests(geneExpression[ ,index], g)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

#Modeling Batch effects Exercises
#For the dataset we have been working with, models do not help due to the almost perfect confounding. This is one reason we created the subset dataset:
library(GSE5859Subset)
data(GSE5859Subset)
#Here we purposely confounded month and group (sex) but not completely:
sex = sampleInfo$group
month = factor( format(sampleInfo$date,"%m"))
table( sampleInfo$group, month)
#1
#Using the functions rowttests and qvalue compare the two groups, in this case males and females so coded in sex. Because this is a smaller dataset, which decreases our power, we will use a more lenient FDR cut-off of 10%.
#How many gene have q-values less than 0.1?
library(qvalue)
library(genefilter)
pvals = rowttests(geneExpression,factor(sampleInfo$g))$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals<0.1)

#2
#Note that sampleInfo$group here represents males and females. Thus we expect differences to be on chrY and, for genes that escape inactivation, chrX. Note that we do not expect many autosomal genes to be different between males and females. This gives us an opportunity to evaluate false and true positives with experimental data. For example, we evaluate results using the proportion genes of the list that are on chrX or chrY.
#For the list of genes with q<0.1 calculated in Modeling Batch Effects Exercises #1, what proportion of genes are on chrX or chrY?
XYprop = which((qvals<0.1)%in%geneAnnotation$CHR == "chrX" | geneAnnotation$CHR == "chrY")
length(XYprop)/length(impvals)
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

#3
#Now for the autosomal genes (not on chrX and chrY) for which q-value < 0.1 perform a t-test comparing samples processed in June to those processed in October.
#What proportion of these have p-values < 0.05?
index = which(qvals<0.1 & !geneAnnotation$CHR%in%c("chrX","chrY"))
month = factor( format(sampleInfo$date,"%m"))
pvals = rowttests(geneExpression[index,],month)$p.value
mean(pvals<0.05)

#4
#The above result shows that the great majority of the autosomal genes show differences due to processing data. This provides further evidence that confounding is resulting in false positives. So we are going to try to model the month effect to better estimate the sex effect. We are going to use a linear model:
#Which of the following creates the appropriate design matrix?
X = model.matrix(~sex+month)

#5
#Now use the X defined above to fit a regression model using lm for each gene. Note that you can obtain p-values for estimated parameters using summary. Here is an example:
X = model.matrix(~sex+month)
i = 234
y = geneExpression[i,]
fit = lm(y~X-1)
summary(fit)$coef
#How many of the q-values for the group comparison are <0.1 now?
pvals <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
} ) )
qvals <- qvalue(pvals)$qvalue
ct <- which(qvals<0.1)
length(ct)

#6
#With this new list, what proportion of these are chrX and chrY?
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

#7
#Now, from the linear model in Modeling Batch Effects Exercises #6, extract the p-values related to the coefficient representing the October versus June differences using the same linear model.
#How many of the q-values for the month comparison are < 0.1 now?
pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X)
  summary(fit)$coef[3,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1)

#Factor Analysis Exercises
#We will continue to use this dataset:
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)
#and define
y = geneExpression - rowMeans(geneExpression)
#Compute and plot an image of the correlation for each sample. Make two image plots of these correlations. In the first one, plot the correlation as image. In the second, order the samples by date and then plot the an image of the correlation. The only difference in these plots is the order in which the samples are plotted.
#1
#Based on these plots, which of the following you would say is true:
library(rafalib)
sex = sampleInfo$group
mypar(1,2)
cors = cor(y)
image(cors)
o = order(sampleInfo$date)
image(cors[o,o])
#The fact that in the plot ordered by month we see two groups mainly driven by month and within these, we see subgroups driven by date seems to suggest date more than month per se are the hidden factors.

#2
#Based on the correlation plots above, we could argue that there are at least two hidden factors. Using PCA estimate these two factors. Specifically, apply the svd to y and use the first two PCs as estimates.
#Which command gives us these estimates?
pcs = svd(y)$v[,1:2]

#3
#Plot each of the estimated factor ordered by date. Use color to denote month. The first factor is clearly related to date.
#Which of the following appear to be most different according to this factor?
pcs = svd(y)$v[,1:2]
o = order(sampleInfo$date)
cols = as.numeric(month)[o]
mypar(2,1)
for(i in 1:2){
  plot(pcs[o,i],col=cols,xaxt="n",xlab="")
  label = gsub("2005-","",sampleInfo$date[o])
  axis(1,1:ncol(y),label,las=2)
}
#June 23 and June 27 

#4
#Use the svd function to obtain the principal components (PCs) for our detrended gene expression data y:
#How many principal components (PCs) explain more than 10% each of the variability?
s = svd(y)
varexplained = s$d^2/ sum(s$d^2)
plot(varexplained)
sum(varexplained>0.10)

#5
#Which PC most correlates (negative or positive correlation) with month?
month = factor( format(sampleInfo$date,"%m"))
cors = cor( as.numeric(month),s$v)
plot(t(cors))
which.max(abs(cors))
#What is this correlation (in absolute value)?
max(abs(cors))

#6
#Which PC most correlates (negative or positive correlation) with sex?
sex = factor( format(sampleInfo$group))
cors = cor( as.numeric(sex),s$v)
plot(t(cors))
which.max(abs(cors))
#What is this correlation (in absolute value)?
max(abs(cors))


#7
#Now instead of using month, which we have shown does not quite describe the batch, add the two estimated factors in Factor Analysis Exercises #6 to the linear model we used in previous exercises.
X <- model.matrix(~sex+s$v[,1:2])
#Apply this model to each gene, and compute q-values for the sex difference.
#How many q-values are <0.1 for the sex comparison?
pvals <- t( sapply(1:nrow(geneExpression),function(j){
  y <- geneExpression[j,]
  fit <- lm(y~X-1)
  summary(fit)$coef[2,4]
} ) )
qvals <- qvalue(pvals)$qvalue
ct <- which(qvals<0.1)
length(ct)
#What proportion of the genes are on chrX and chrY?
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

#SVA Exercises
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library(sva)
library(Biobase)
library(GSE5859Subset)
data(GSE5859Subset)
#1
#In the previous section we estimated factors using PCA. But we noted that the first factor was correlated with our outcome of interest:
s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])
#As in the previous questions we are interested in finding genes that are differentially expressed between the two groups (males and females in this case). Here we learn to use SVA to estimate these effects while using a factor analysis approach to account for batch effects.
#The svafit function estimates factors, but downweighting the genes that appear to correlate with the outcome of interest. It also tries to estimate the number of factors and returns the estimated factors like this:
sex = sampleInfo$group
mod = model.matrix(~sex)
svafit = sva(geneExpression,mod)
head(svafit$sv)
#Note that the resulting estimated factors are not that different from the PCs
for(i in 1:ncol(svafit$sv)){
  print( cor(s$v[,i],svafit$sv[,i]) )
}
#Now fit a linear model to estimate the difference between males and females for each gene but that instead accounting for batch effects using month it includes the factors estimated by sva in the model. Use the qvalue function to estimate q-values.
#How many genes have q-value < 0.1?
library(qvalue)
library(sva)
X= model.matrix(~sex+svafit$sv)
pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1)

#2
#What proportion of the genes from SVA Exercises #1 are from chrY or chrX?
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

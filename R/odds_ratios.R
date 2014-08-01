# Simple, internal function to check input datasets for compatability with later functions in library(sppairs)
or.check<-function(dataset)
{
if(dim(dataset)[2]!=2){stop("input does not have two columns")}		# check size
test<-as.character(apply(dataset, 2, class))		# check contains numbers
for(i in 1:2){if(any(c("numeric", "integer")==test[i])==F)
	{stop(paste("column", i, "is not numeric", sep=" "))}}
	if(any(dataset>1)){dataset<-make.binary(dataset)}		# check if any values >1; if so, convert to binary
invisible(dataset)		# return the (possibly corrected) dataset for use in later functions.
}


# Calculate asymmetric odds ratio from a contingency table
or.contingency<-function(dataset)	
{
dataset<-or.check(dataset)		# will either correct the dataset, or stop this function with an error
cont.table <-as.matrix(table(dataset))[2:1, 2:1]
c<-cont.table[1, 1]; d<-cont.table[1, 2]; e<-cont.table[2, 1]; f<-cont.table[2, 2]
odds.ratio<-(c/d) / ( (c+e) / (d+f) )
return(odds.ratio)
}


# Calculate odds ratio from results of glm or glmer. Internal to or.glmer and or.glm
or.regression<-function(b, z0, z1)
{
library(boot)	# for logit() & inv.logit
# regardless of the method used, calculate the odds ratio
odds.ratio<-exp( z1- logit(
	(b * inv.logit(z1)) +
	((1-b) * inv.logit(z0))
)	# end logit
)	# end exp
return(odds.ratio)
}


# Odds ratio calculation using lme4
or.glmer<-function(dataset, random.effect)
{
dataset<-or.check(dataset)		# will either correct the dataset, or stop this function with an error
b<-(1/dim(dataset)[1])*sum(dataset[, 2])		# proportion of rows at which sp. B occurred
library(lme4)
model0<-glmer(dataset[, 1]~1 + (1|random.effect), family=binomial(link="logit"))
model1<-glmer(dataset[, 1]~dataset[, 2] + (1|random.effect), family=binomial(link="logit"))
z0<-as.numeric(fixef(model0))
z1<-sum(as.numeric(fixef(model1)))
odds.ratio<-or.regression(b, z0, z1)
return(odds.ratio)
}


# Odds ratio calculation using glm
or.glm<-function(dataset)
{
dataset<-or.check(dataset)		# will either correct the dataset, or stop this function with an error
b<-(1/dim(dataset)[1])*sum(dataset[, 2])		# proportion of rows at which sp. B occurred
model0<-glm(dataset[, 1]~1, family=binomial(link="logit"))
model1<-glm(dataset[, 1]~dataset[, 2], family=binomial(link="logit"))
z0<-as.numeric(coef(model0)[1])	# intercept; occurrence of sp. A in the absence of sp. B
z1<-sum(as.numeric(coef(model1)))	# covariate; occurrence of sp. A in the presence of sp. B
odds.ratio<-or.regression(b, z0, z1)
return(odds.ratio)
}
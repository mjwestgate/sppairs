# Functions for calculating odds ratios

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
for(i in 1:2){dataset[, i]<-factor(dataset[, i], levels=c(0, 1), labels=c("0", "1"))}	 # avoids errors with 100% zeros or ones
count.table <-as.matrix(table(dataset))[2:1, 2:1]
pc.table<-(100/sum(count.table))*count.table
c<-pc.table[1, 1]; d<-pc.table[1, 2]; e<-pc.table[2, 1]; f<-pc.table[2, 2]
# odds.ratio<-(c/d) / ( (c+e) / (d+f) )
odds.ratio<-(c / e) / ( (c + d) / (e + f) )
return(odds.ratio)
}


# Calculate symmetric odds ratio
or.symmetric<-function(dataset)
{
dataset<-or.check(dataset)
cont.table <-as.matrix(table(dataset))[2:1, 2:1]
pc.table<-(100/sum(cont.table))*cont.table
c<-pc.table[1, 1]; d<-pc.table[1, 2]; e<-pc.table[2, 1]; f<-pc.table[2, 2]
odds.ratio<-(c/e) / (d/f) # == (c/d)/(e/f)
return(odds.ratio)
}


# Calculate odds ratio from results of glm or glmer. Internal to or.glmer and or.glm
or.regression<-function(b, z0, z1)
{
val1<-(b * inv.logit(z1)) 
val2<-((1-b) * inv.logit(z0))
odds.ratio<-exp(z1 - logit(val1 + val2))
return(odds.ratio)
}


# Odds ratio calculation using lme4
or.glmer<-function(dataset, random.effect, complex=FALSE)
{
b<-(1/nrow(dataset))* sum(dataset[, 2])		# proportion of rows at which sp. B occurred
model<-glmer(dataset[, 1]~dataset[, 2] + (1|random.effect), family=binomial(link="logit"),
	control=glmerControl(optimizer="bobyqa")) # Note: optimizer doesn't make much difference in most cases
# work out if there is a convergence error
if(length(model@optinfo$conv$lme4$messages)>0){converge.warning<-TRUE
	}else{converge.warning<-FALSE}
z0<-as.numeric(fixef(model))[1] # coef of a model with fixed effects, but ignoring their effects
z1<-sum(as.numeric(fixef(model)))
odds.ratio<-or.regression(b, z0, z1)
if(complex){
	return(c(b=b, 
		intercept=z0, slope=as.numeric(fixef(model))[2], 
		converge.warning= converge.warning, odds=odds.ratio))
}else{return(odds.ratio)}
}


# Odds ratio calculation using glm
or.glm<-function(dataset)
{
dataset<-or.check(dataset)		# will either correct the dataset, or stop this function with an error
b<-(1/dim(dataset)[1])*sum(dataset[, 2])		# proportion of rows at which sp. B occurred
model<-glm(dataset[, 1]~dataset[, 2], family=binomial(link="logit"))
z0<-as.numeric(coef(model)[1])	# intercept; occurrence of sp. A in the absence of sp. B
z1<-sum(as.numeric(coef(model)))	# intercept + slope; occurrence of sp. A in the presence of sp. B
odds.ratio<-or.regression(b, z0, z1)
if(complex){
	return(c(b=b, intercept=z0, slope=as.numeric(coef(model)[2]), odds=odds.ratio))
}else{return(odds.ratio)}
}
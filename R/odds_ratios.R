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


# Pairwise odds ratio calculation
spaa<-function(dataset, method, rarity.cutoff, random.effect)
{

# set some default behaviour
if(missing(method)){
	if(missing(random.effect)){method<-"contingency"}else{method<-"glmer"}
}else{	# i.e. is method is specified
	if(method=="glmer" & missing(random.effect)){
	method <-"contingency"
	cat("Warning: method 'glmer' not available without specifying 'random.effect'. Switched to method 'contingency'")}
}
if(missing(rarity.cutoff))rarity.cutoff<-0.1
if(any(dataset>1)){dataset<-make.binary(dataset)}		# check if any values >1; if so, convert to binary
  
# apply rarity cutoff 
occu.result<-prop.occupied(dataset)
dataset<-dataset[, which(occu.result>rarity.cutoff)]
  
# create an object to export frequency information
frequency.result<-occu.result[which(occu.result>rarity.cutoff)]
frequency.result<-data.frame(
	name=names(frequency.result),
	frequency=as.numeric(frequency.result),
	stringsAsFactors=FALSE)
#frequency.result$species<-as.character(frequency.result$species)
  
# create combn of remaining species
n.species<-dim(dataset)[2]
combinations<-t(combn(c(1:n.species), 2))
combinations<-rbind(combinations, combinations[, c(2, 1)])
combinations.final<-data.frame(
	col1=combinations[, 1],
	col2=combinations[, 2],
	sp1=colnames(dataset)[combinations[, 1]],
	sp2=colnames(dataset)[combinations[, 2]],
	odds=rep(NA, dim(combinations)[1]),
	stringsAsFactors=FALSE)
combinations.final$odds<-as.numeric(combinations.final$odds)
n.rows<-dim(combinations.final)[1]

# run analysis in loop
for(i in 1: n.rows){
	dataset.thisrun<-dataset[, c(combinations.final$col1[i], combinations.final$col2[i])]
	or.simple<-or.contingency(dataset.thisrun)	# run as a test: req. as glmer will fail if no overlap in presences.
	switch(method,
		"contingency"={combinations.final$odds[i]<-or.simple},
		"glm"={if(any(c(0, Inf)==or.simple)){combinations.final$odds[i]<-or.simple
			}else{combinations.final$odds[i]<-or.glm(dataset.thisrun)}},
		"glmer"={if(any(c(0, Inf)==or.simple)){combinations.final$odds[i]<-or.simple
			}else{combinations.final$odds[i]<-or.glmer(dataset.thisrun, random.effect)}})
	}
  
# create and return objects necessary for plotting
output<-list(
	values=list(
		method=method,
		observations=dim(dataset)[1],
		rarity.cutoff= rarity.cutoff,
		species= n.species,
		combinations= n.rows),
	species=frequency.result,
	combinations=combinations.final)

return(output)
}


# Internal function to convert the 'long' output from pairwise.odds.ratios() into a matrix
pairwise.or.matrix<-function(or.object)
{
# work out species info
n.species<-max(or.object$col1)
species.names<-rep("none", n.species)
for(i in 1:n.species){species.names[i]<-as.character(or.object$sp1[which(or.object$col1==i)[1]])}
  
# 	create & fill output
output.matrix<-matrix(data=NA, n.species, n.species)
colnames(output.matrix)<-species.names
rownames(output.matrix)<-species.names
for(i in 1:dim(or.object)[1]){output.matrix[or.object$col1[i], or.object$col2[i]]<-or.object$odds[i]}
  
return(output.matrix)
}
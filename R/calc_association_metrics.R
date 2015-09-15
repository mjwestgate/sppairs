# Functions for calculating the association between a pair of vectors (mainly from odds ratios)

# Calculate asymmetric odds ratio from a contingency table
or.asymmetric<-function(dataset, run.check=FALSE)	
{
if(run.check)dataset<-or.check(dataset)		# will either correct the dataset, or stop this function with an error
for(i in 1:2){dataset[, i]<-factor(dataset[, i], levels=c(0, 1), labels=c("0", "1"))}	 # avoids errors with 100% zeros or ones
count.table <-as.matrix(table(dataset))
pc.table<-(100/sum(count.table))*count.table
# values for calculation
both.pres<-pc.table[2, 2]
ApresBabs<-pc.table[2, 1]
Bpres<-sum(pc.table[, 2])
Babs<-sum(pc.table[, 1])
# calculate and return
odds.ratio<-(both.pres/ApresBabs) / (Bpres / Babs)
return(odds.ratio)
}


# Calculate symmetric odds ratio
or.symmetric<-function(dataset, run.check= FALSE)
{
if(run.check)dataset<-or.check(dataset)	
cont.table <-as.matrix(table(dataset))
pc.table<-(100/sum(cont.table))*cont.table
# inputs
both.pres<-pc.table[2, 2]
both.abs<-pc.table[1, 1]
ApresBabs<-pc.table[2, 1]
AabsBpres<-pc.table[1, 2]
# calculate and return
odds.ratio<-(both.pres / AabsBpres) / (ApresBabs / both.abs)
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
or.glmer<-function(dataset, random.effect, complex=FALSE, run.check= FALSE)
{
if(run.check)dataset<-or.check(dataset)	
b<-(1/nrow(dataset))* sum(dataset[, 1])		# proportion of rows at which sp. B occurred
model<-glmer(dataset[, 2]~dataset[, 1] + (1|random.effect), family=binomial(link="logit"),
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
		converge.warning= converge.warning, or.contingency=or.contingency(dataset), value=odds.ratio))
}else{return(odds.ratio)}
}


# Odds ratio calculation using glm
or.glm<-function(dataset, complex=FALSE, run.check= FALSE)
{
if(run.check)dataset<-or.check(dataset)	# will either correct the dataset, or stop this function with an error
b<-(1/dim(dataset)[1])*sum(dataset[, 1])		# proportion of rows at which sp. B occurred
model<-glm(dataset[, 2]~dataset[, 1], family=binomial(link="logit"))
z0<-as.numeric(coef(model)[1])	# intercept; occurrence of sp. A in the absence of sp. B
z1<-sum(as.numeric(coef(model)))	# intercept + slope; occurrence of sp. A in the presence of sp. B
odds.ratio<-or.regression(b, z0, z1)
if(complex){return(c(b=b, intercept=z0, slope=as.numeric(coef(model)[2]), value=odds.ratio))
}else{return(odds.ratio)}
}


# mutual information - based on entropy::mi.plugin
mutual.information<-function(dataset, run.check= FALSE){
	if(run.check)dataset<-or.check(dataset)
	for(i in 1:2){dataset[, i]<-factor(dataset[, i], levels=c(0, 1), labels=c("0", "1"))}	 # avoids errors with 100% zeros or ones
	count.table <-as.matrix(table(dataset))
	pc.table<-(1/sum(count.table))*count.table
	entropy<-function(x){
		y<-x/sum(x)
		keep<-which(y>0)
		H<-(-sum(y[keep]*log(y[keep], base=2)))
		return(H)}
	H1<-entropy(rowSums(pc.table))
	H2<-entropy(colSums(pc.table))
	H12<-entropy(pc.table)
	result<-H1+H2-H12
	return(result)
	}
# Note: equvalent code using library(entropy) would be: mi.plugin(pc.table, unit="log2")
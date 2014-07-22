# Convert an abundance dataset to presence/absence.
make.binary<-function(dataset, threshold)
{
dataset<-as.data.frame(dataset)
if(missing(threshold)){	# for converting abundance to presence/absence data, this works well
	for(i in 1: dim(dataset)[2]){if(any(dataset[, i]>0))dataset[which(dataset[, i]>0), i]<-1}
}else{	# otherwise, cut by threshold
	for(i in 1: dim(dataset)[2]){dataset[, i]<-as.numeric(
		cut(dataset[, i], breaks=c(-Inf, threshold, Inf), include.lowest=TRUE))-1}}
return(dataset)
}


# Calculate the proportion of rows occupied (where rows represent either visits or sites)
prop.occupied<-function(dataset, random.effect)
{
dataset<-as.data.frame(dataset)	# forces dim() to function for ncol=1
if(missing(random.effect)){
	occupied<-(1/dim(dataset)[1])*apply(dataset, 2, sum)
}else{	# i.e. if 'site' is given as a random effect
	n.sites<-length(levels(random.effect))
	n.occupied<-apply(dataset, 2, function(x)
	{length(which(xtabs(x~random.effect)>0))})
	occupied<-(1/n.sites)*n.occupied}
return(occupied)
}
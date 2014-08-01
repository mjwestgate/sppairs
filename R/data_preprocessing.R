# Convert an abundance dataset to presence/absence.
make.binary<-function(dataset, threshold)
{
dataset<-as.data.frame(dataset)	 # avoid errors if single col. entered
if(missing(threshold)){	# for converting abundance to presence/absence data, this works well
	dataset<-apply(dataset, 2, function(x){
		if(any(x>0)){x[which(x>0)]<-1}
		return(x)})
}else{	# otherwise, cut by threshold
	dataset<-apply(dataset, 2, function(x){
		x<-as.numeric(cut(x, breaks=c(-Inf, threshold, Inf), include.lowest=TRUE))-1
		return(x)})
	}
dataset<-as.data.frame(dataset)	# as apply converts dataset to matrix
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
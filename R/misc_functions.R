# Miscellaneous useful functions


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
prop.occupied<-function(dataset)
{
dataset<-as.data.frame(dataset)	# forces dim() to function for ncol=1
occupied<-(1/dim(dataset)[1])*apply(dataset, 2, sum)
return(occupied)
}


# function to make a square (asymmetric) distance matrix of odds ratios
make.or.matrix<-function(
	spaa.object	# result from spaa()
	){
	if(class(spaa.object)!="spaa"){stop("input not of class spaa")}

	spp.names<-spaa.object$species$species
	n.spp<-length(spp.names)

	# work out what to do for an asymmetric matrix
	result<-matrix(data=NA, nrow= n.spp, ncol= n.spp)
	colnames(result)<-spp.names
	rownames(result)<-spp.names
	for(i in 1:nrow(spaa.object$combinations)){
		sp1<-spaa.object$combinations$sp1[i]
		sp2<-spaa.object$combinations$sp2[i]
		row.i<-which(spp.names==sp1)
		col.i<-which(spp.names==sp2)
		if(spaa.object$asymmetric){
			result[row.i, col.i]<-spaa.object$combinations$odds[i]
		}else{
			result[row.i, col.i]<-spaa.object$combinations$odds[i]
			result[col.i, row.i]<-spaa.object$combinations$odds[i]}
		}

	return(result)
	}
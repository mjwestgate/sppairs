# Core functions

# SPAA
# Pairwise association calculation
spaa<-function(x, method="or.symmetric", asymmetric=FALSE, ...)#. further information passed to method (e.g. random.effect, complex)
{

if(ncol(x)<2){stop("Error: No species are sufficiently common to analyse")}

# create combn of remaining species
n.species<-dim(x)[2]
combinations<-t(combn(c(1:n.species), 2))
if(asymmetric){	combinations<-rbind(combinations, combinations[, c(2, 1)])}
combinations.final<-data.frame(
	col1=combinations[, 1],
	col2=combinations[, 2],
	sp1=colnames(x)[combinations[, 1]],
	sp2=colnames(x)[combinations[, 2]],
	stringsAsFactors=FALSE)

# calculate pairwise association using the specified method
result<-apply(combinations.final[, 1:2], 1, FUN=function(y, source, method, ...){
	y<-as.numeric(y)
	dataset.thisrun<-source[, y]
	value<-do.call(method, list(dataset.thisrun, ...))
	return(value)
	}, source=x, method=method)
combinations.final$value<-as.numeric(result)
return(combinations.final[, -c(1:2)])

}


## FUNCTIONS TO PROCESS INPUTS ##
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
# note: this defaults to FALSE, but can be called by spaa by setting run.check=TRUE


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


# function to ensure dataset is in a useful format for binary analysis
# This function is not called by any others; but may be useful for preprocessing; e.g. spaa(clean.dataset(x))
# most useful purpose is to cut very rare or very common species.
clean.dataset<-function(x, make.binary=TRUE, cutoff.min=0.1, cutoff.max=1, ...){ #... allows threshold to be passed to make.binary
	if(make.binary){x<-make.binary(x, ...)}		# check if any values >1; if so, convert to binary
	# apply rarity cutoff 
	occu.result<-prop.occupied(x)
	keep.entries<-which(c(occu.result>=cutoff.min & occu.result<=cutoff.max)==TRUE)
	if(length(keep.entries)>0){
		x<-x[, keep.entries]
		return(x)
	}else{stop("Error: no columns have frequencies within the specified limits")}
}




## SUMMARY STATISTICS ON INPUTS
# Calculate the proportion of rows occupied (where rows represent either visits or sites)
prop.occupied<-function(dataset)
{
dataset<-as.data.frame(dataset)	# forces dim() to function for ncol=1
occupied<-(1/dim(dataset)[1])*apply(dataset, 2, sum)
return(occupied)
}


# binary entropy function
# returns a vector of the same length as the input, giving the binary entropy (in Shannons)
# may be useful in combinations with prop.occupied
binary.entropy<-function(x){ # input is a vector of proportions
	entropy<-function(x){if(x==0){return(0)}else{return(-x*(log(x, base=2)))}}
	sapply(x, FUN=function(x){sum(c(entropy(x), entropy(1-x)))})
	}



## FUNCTIONS TO PROCESS OUTPUTS
# function to make a square (asymmetric) distance matrix of odds ratios
make.or.matrix<-function(
	input	# result from spaa()
	){
	# work out properties of the input
	spp.names<-unique(c(input[, 1], input[, 2]))
	n.spp<-length(spp.names)
	if(nrow(input)==choose(n.spp, 2)){asymmetric<-FALSE}else{asymmetric<-TRUE}
	# create a matrix
	result<-matrix(data=NA, nrow= n.spp, ncol= n.spp)
	colnames(result)<-spp.names
	rownames(result)<-spp.names
	# fill with a loop
	for(i in 1:nrow(input)){
		sp1<-input$sp1[i]
		sp2<-input$sp2[i]
		row.i<-which(spp.names==sp1)
		col.i<-which(spp.names==sp2)
		if(asymmetric){
			result[row.i, col.i]<-input$value[i]
		}else{
			result[row.i, col.i]<-input$value[i]
			result[col.i, row.i]<-input$value[i]}
	}
	return(result)
	}
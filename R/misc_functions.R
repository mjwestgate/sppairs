# Miscellaneous useful functions - should probably make these available to users


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


# function to make a square (asymmetric) distance matrix of odds ratios
make.or.matrix<-function(
	spaa.object	# result from spaa()
	){
	if(class(spaa.object)!="spaa"){stop("input not of class spaa")}

	# get the data we are interested in; split into two separate dataset (upper and lower triangles)
	dataset<-spaa.object$combinations
	halfway<-dim(dataset)[1]*0.5
	dset1<-dataset[c(1:halfway), ]
	dset2<-dataset[c((halfway+1):dim(dataset)[1]), ]

	# create a fake distance matrix with the appropriate information
	species<-sort(unique(dataset$sp1))
	n.species<-length(species)
	test.dist<-as.dist(matrix(data=0, nrow=n.species, ncol=n.species), upper=FALSE)
	attr(test.dist, "Labels")<-species

	# overwrite this with correct data
	dist1<-test.dist; dist2<-test.dist
	dist1[1:dim(combn(n.species, 2))[2]]<-dset1$odds
	dist2[1:dim(combn(n.species, 2))[2]]<-dset2$odds

	# convert back to matrix
	dist1<-as.matrix(dist1); dist2<-as.matrix(dist2)
	# any(c(dist1==t(dist1))==FALSE)	# check symmetry

	# merge, such that dist1 provides lower; dist2 provides upper; separated by NAs
	final.matrix<-matrix(data=NA, nrow=n.species, ncol=n.species)
	colnames(final.matrix)<-species
	rownames(final.matrix)<-species
	for(i in 1:n.species){	
		if(i>1){
			rows<-which(c(1:n.species)<i)
			final.matrix[rows, i]<-dist2[rows, i]}
		if(i<n.species){
			rows<-which(c(1:n.species)>i)
			final.matrix[rows, i]<-dist1[rows, i]}
		}
	# which(c(final.matrix ==t(final.matrix))==TRUE) # some values identical; but very few

	return(final.matrix)
	}
# Pairwise odds ratio calculation
spaa<-function(dataset, method, random.effect, complex=FALSE, clean=TRUE)
{

# SET DEFAULTS
if(missing(method)){
	if(missing(random.effect)){method<-"contingency"}else{method<-"glmer"}
}else{	# i.e. is method is specified
	if(method=="glmer" & missing(random.effect)){
	method <-"contingency"
	cat("Warning: method 'glmer' not available without specifying 'random.effect'. Switched to method 'contingency'")}
}
# if method not recognised, set to or.contingency
method.list<-c("contingency", "glm", "glmer", "pearsons", "spearmans", "or.symmetric", "mutual.info")
if(any(method.list==method)==F){	
	cat(paste("Warning: method '", method, "' not recognised. Switching to method 'contingency'", sep=""))
	method<-"contingency"
	}
# set symmetry behavior (avoids duplicating symmetric calculations)
if(any(c("or.symmetric", "pearson", "spearman", "mutual.info")==method)){asymmetric<-FALSE}else{asymmetric<-TRUE}
# If requested, clean the dataset
if(length(clean)==1 & is.logical(clean[1])){
	if(clean(dataset))dataset<-clean.dataset(dataset) # run using defaults
}else{
if(is.list(clean)){
	dataset<-do.call("clean.dataset", append(list(dataset=dataset), clean))}}
if(ncol(dataset)==0){stop("Error: No species are sufficiently common to analyse")}

# create combn of remaining species
n.species<-dim(dataset)[2]
combinations<-t(combn(c(1:n.species), 2))
if(asymmetric){	combinations<-rbind(combinations, combinations[, c(2, 1)])}
combinations.final<-data.frame(
	col1=combinations[, 1],
	col2=combinations[, 2],
	sp1=colnames(dataset)[combinations[, 1]],
	sp2=colnames(dataset)[combinations[, 2]],
	stringsAsFactors=FALSE)

# RUN ANALYSIS
# note: do.call() would be more efficient than switch() for method selection
result<-apply(combinations.final[, 1:2], 1, FUN=function(x, source, method){
	x<-as.numeric(x)
	dataset.thisrun<-source[, x]
	result<-switch(method, 
		"contingency"={or.contingency(dataset.thisrun, run.check=FALSE)},
		"glm"={or.glm(dataset.thisrun, complex, run.check=FALSE)},
		"glmer"={or.glmer(dataset.thisrun, random.effect, complex, run.check=FALSE)},
		"or.symmetric"={or.symmetric(dataset.thisrun, run.check=FALSE)},
		"pearson"={cor(x=dataset.thisrun[, 1], y= dataset.thisrun[, 2], method="pearson")},
		"spearman"={cor(x=dataset.thisrun[, 1], y= dataset.thisrun[, 2], method="spearman")},
		"mutual.info"={mutual.information(dataset.thisrun, run.check=FALSE)}
		)
	return(result)
	}, source=dataset, method=method)

if(complex){
	result<-as.matrix(as.data.frame(result))
	result<-t(result)
	rownames(result)<-NULL
	result<-as.data.frame(result)
	combinations.final<-cbind(combinations.final, result)
}else{combinations.final$value<-as.numeric(result)}

# return only a data.frame with the requisite calculations
return(combinations.final[, -c(1:2)])

}
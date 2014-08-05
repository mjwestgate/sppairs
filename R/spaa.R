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
if(any(c("contingency", "glm", "glmer")==method)==F){	# if method not recognised, set to 
	cat(paste("Warning: method '", method, "' not recognised. Switching to method 'contingency'", sep=""))
	method<-"contingency"
	}
if(missing(rarity.cutoff))rarity.cutoff<-0.1
if(any(dataset>1)){dataset<-make.binary(dataset)}		# check if any values >1; if so, convert to binary
  
# apply a cutoff to exclude spp. present at all (or none) of the sites
dataset<-dataset[, which(c(apply(dataset, 2, sum)==dim(dataset)[1])==FALSE)]	# remove 100% occupied
dataset<-dataset[, which(c(apply(dataset, 2, function(x){length(which(x==0))})==dim(dataset)[1])==FALSE)]	# 0%
# apply rarity cutoff 
occu.result<-prop.occupied(dataset)
dataset<-dataset[, which(occu.result>rarity.cutoff)]

# create an object to export frequency information
frequency.result<-occu.result[which(occu.result>rarity.cutoff)]
frequency.result<-data.frame(
	species=names(frequency.result),
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

class(output)<-"spaa"

return(output)
}
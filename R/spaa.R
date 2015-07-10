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
# if method not recognised, set to or.contingency
if(any(c("contingency", "glm", "glmer", "pearsons", "spearmans", "or.symmetric")==method)==F){	
	cat(paste("Warning: method '", method, "' not recognised. Switching to method 'contingency'", sep=""))
	method<-"contingency"
	}
if(missing(rarity.cutoff))rarity.cutoff<-0.1
if(any(dataset>1)){dataset<-make.binary(dataset)}		# check if any values >1; if so, convert to binary
if(any(c("or.symmetric", "pearsons", "spearmans")==method)){asymmetric<-FALSE}else{asymmetric<-TRUE}

# apply a cutoff to exclude spp. present at all (or none) of the sites
dataset<-dataset[, which(c(apply(dataset, 2, sum)==dim(dataset)[1])==FALSE)]	# remove 100% occupied
dataset<-dataset[, which(c(apply(dataset, 2, function(x){
	length(which(x==0))})==dim(dataset)[1])==FALSE)]	# 0%
# apply rarity cutoff 
occu.result<-prop.occupied(dataset)
dataset<-dataset[, which(occu.result>rarity.cutoff)]
if(ncol(dataset)==0){stop("Error: No species are sufficiently common to analyse")}

# create an object to export frequency information
frequency.result<-occu.result[which(occu.result>rarity.cutoff)]
frequency.result<-data.frame(
	species=names(frequency.result),
	frequency=as.numeric(frequency.result),
	stringsAsFactors=FALSE)
  
# create combn of remaining species
n.species<-dim(dataset)[2]
combinations<-t(combn(c(1:n.species), 2))
if(asymmetric){
	combinations<-rbind(combinations, combinations[, c(2, 1)])}
combinations.final<-data.frame(
	col1=combinations[, 1],
	col2=combinations[, 2],
	sp1=colnames(dataset)[combinations[, 1]],
	sp2=colnames(dataset)[combinations[, 2]],
	odds=rep(NA, dim(combinations)[1]),
	stringsAsFactors=FALSE)
combinations.final$odds<-as.numeric(combinations.final$odds)
n.rows<-dim(combinations.final)[1]

# run analysis using apply
result<-apply(combinations.final[, 1:2], 1, FUN=function(x, source, method){
	x<-as.numeric(x)
	dataset.thisrun<-source[, x]
	or.simple<-or.contingency(dataset.thisrun)	
	result<-switch(method, 
		"contingency"={or.simple},
		"glm"={
			if(any(c(0, Inf)==or.simple)){or.simple
			}else{or.glm(dataset.thisrun)}},
		"glmer"={if(any(c(0, Inf)==or.simple)){or.simple
			}else{or.glmer(dataset.thisrun, random.effect)}},
		"or.symmetric"={or.symmetric(dataset.thisrun)},
		"pearsons"={cor(x=dataset.thisrun[, 1], y= dataset.thisrun[, 2], method="pearson")},
		"spearmans"={cor(x=dataset.thisrun[, 1], y= dataset.thisrun[, 2], method="spearman")})
	return(result)
	}, source=dataset, method=method)
combinations.final$odds<-as.numeric(result)

# run analysis in loop
# for(i in 1: n.rows){
	# dataset.thisrun<-dataset[, c(combinations.final$col1[i], combinations.final$col2[i])]
	# or.simple<-or.contingency(dataset.thisrun)	# run as a test: req. as glmer will fail if no overlap in presences.
	# switch(method,
		# "contingency"={combinations.final$odds[i]<-or.simple},
		# "glm"={if(any(c(0, Inf)==or.simple)){combinations.final$odds[i]<-or.simple
			# }else{combinations.final$odds[i]<-or.glm(dataset.thisrun)}},
		# "glmer"={if(any(c(0, Inf)==or.simple)){combinations.final$odds[i]<-or.simple
			# }else{combinations.final$odds[i]<-or.glmer(dataset.thisrun, random.effect)}},
		# "or.symmetric"={combinations.final$odds[i]<-or.symmetric(dataset.thisrun)},
		# "pearsons"={combinations.final$odds[i]<-cor(
			# x=dataset.thisrun[, 1], y= dataset.thisrun[, 2], method="pearson")},
		# "spearmans"={combinations.final$odds[i]<-cor(
			# x=dataset.thisrun[, 1], y= dataset.thisrun[, 2], method="spearman")}
		# )
	# }
  
# create and return objects necessary for plotting
output<-list(
	values=list(
		method=method,
		observations=dim(dataset)[1],
		rarity.cutoff= rarity.cutoff,
		species= n.species,
		combinations= n.rows),
	species=frequency.result,
	combinations=combinations.final,
	asymmetric=asymmetric)

class(output)<-"spaa"

return(output)
}


# summary
summary.spaa<-function(object){
	print.spaa(object)
	frequencies<-as.numeric(xtabs(rep(1, dim(object$combination)[1])~
		cut(object$combination$odds, breaks=c(0, 0.000001, 1/9, 1/6, 1/3, 3, 6, 9, 10^8, Inf))))
	print(data.frame(odds.ratio=c("0", "<1/9", "<1/6", "<1/3", ">3", ">6", ">9", "Inf"),
		n.pairs=frequencies[c(1:4, 6:9)]))
}


# print
print.spaa<-function(object){
	cat(paste(object$values$combinations, "combinations of", object$values$species, "species, ", sep=" "))
	cat(paste("calculated using method '", object$values$method, 
		"' with cutoff = ", object$values$rarity.cutoff, sep=""))
	}



# Internal function to convert the 'long' output from spaa() into a matrix
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

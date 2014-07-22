#' Odds Ratio Test 
#' 
# Simple, internal function to check input datasets for compatability with later functions in library(naa). Checks the number of columns; whether they are numeric; and whether they are binary.
#' @param dataset a two-column matrix or data.frame containing presence/absence data for two species
#' @return Gives either an error message, or invisibly returns 'dataset' (possibly corrected by 'make.binary')
odds.ratio.test<-function(dataset)
{
if(dim(dataset)[2]!=2){stop("input does not have two columns")}		# check size
test<-as.character(apply(dataset, 2, class))		# check contains numbers
for(i in 1:2){if(any(c("numeric", "integer")==test[i])==F)
	{stop(paste("column", i, "is not numeric", sep=" "))}}
	if(any(dataset>1)){dataset<-make.binary(dataset)}		# check if any values >1; if so, convert to binary
invisible(dataset)		# return the (possibly corrected) dataset for use in later functions.
}

#' Odds ratio calculation using a contingency table
#' 
#' Calculate asymmetric odds ratio from a contingency table
#' @param dataset Two column data.frame, with presence/absence for two species
#' @return Value of the odds ratio between those species
or.contingency<-function(dataset)	
{
dataset<-odds.ratio.test(dataset)		# will either correct the dataset, or stop this function with an error
cont.table <-as.matrix(table(dataset))[2:1, 2:1]
c<-cont.table[1, 1]; d<-cont.table[1, 2]; e<-cont.table[2, 1]; f<-cont.table[2, 2]
odds.ratio<-(c/d) / ( (c+e) / (d+f) )
return(odds.ratio)
}

#' Odds ratio calculation using regression
#' 
#' Calculate asymmetric odds ratio using logistic regression. Random effects can be specified to account for multiple visits to the same site (using glmer in lme4).
#' @param dataset Two column data.frame containing species occurrence. Note that first column is species A and second species B, where A is the species subject to investigation (i.e. does p/a of B affect A?)
#' @param random.effect # if specified, should be a factor to be used as random effect in glmer(lme4). If missing, glm() is used to calculate the odds ratio
#' @return value of the odds ratio between those species
or.regression<-function(dataset, random.effect)
{
dataset<-odds.ratio.test(dataset)		# will either correct the dataset, or stop this function with an error
b<-(1/dim(dataset)[1])*sum(dataset[, 2])		# proportion of sites at which sp. B occurred
require(boot)	# for logit() & inv.logit
  
# To call or.regression without random effects (note this should give similar results to or.contingency()):
if(missing(random.effect))	
	{
	model0<-glm(dataset[, 1]~1, family=binomial(link="logit"))
	model1<-glm(dataset[, 1]~dataset[, 2], family=binomial(link="logit"))
	z0<-as.numeric(coef(model0)[1])	# intercept; occurrence of sp. A in the absence of sp. B
	z1<-sum(as.numeric(coef(model1)))	# covariate; occurrence of sp. A in the presence of sp. B
# Alternatively (i.e. if random effects are specified), call glmer(lme4)
}else{
	require(lme4)
	model0<-glmer(dataset[, 1]~1 + (1|random.effect), family=binomial(link="logit"))
	model1<-glmer(dataset[, 1]~dataset[, 2] + (1|random.effect), family=binomial(link="logit"))
	z0<-as.numeric(fixef(model0))
	z1<-sum(as.numeric(fixef(model1)))
}	# end else
  
# regardless of the method used, calculate the odds ratio
odds.ratio<-exp( z1- logit(
	(b * inv.logit(z1)) +
	((1-b) * inv.logit(z0))
)	# end logit
)	# end exp
  
return(odds.ratio)
} 


#' Pairwise odds ratio calculation
#'
#' Calculates a number of odds ratios for all (sufficiently common) species. Basically a wrapper function for or.contingency() and or.regression().
#' @param dataset data.frame where each column is a species. Rows can be sites or visits; if visits, a grouping factor is a good idea.
#' @param random.effect Grouping variable for visits, if given. Passed to or.regression()
#' @param rarity.cutoff Minimum proportion of occupied sites/visits for a species (column) to be included. Defaults to 0.1.
#' @param quiet If TRUE, code does not give a progress report. Defaults to FALSE.
#' @return A list containing the components: [[1]] the frequency at which species occur at sites; [[2]] pairwise odds ratios in long format; [[3]] pairwise odds ratios in wide format
pairwise.odds.ratios<-function(dataset, random.effect, rarity.cutoff, quiet)
{
if(missing(quiet))quiet<-FALSE
if(missing(rarity.cutoff))rarity.cutoff<-0.1
if(any(dataset>1)){dataset<-make.binary(dataset)}		# check if any values >1; if so, convert to binary
  
# apply rarity cutoff 
occu.result<-prop.occupied(dataset)
dataset<-dataset[, which(occu.result>rarity.cutoff)]
  
# create an object to export frequency information
frequency.result<-occu.result[which(occu.result>rarity.cutoff)]
frequency.result<-data.frame(
	species=names(frequency.result),
	value=as.numeric(frequency.result))
frequency.result$species<-as.character(frequency.result$species)
  
# create combn of remaining species
n.species<-dim(dataset)[2]
combinations<-t(combn(c(1:n.species), 2))
combinations<-rbind(combinations, combinations[, c(2, 1)])
combinations.final<-data.frame(
	col1=combinations[, 1],
	col2=combinations[, 2],
	sp1=colnames(dataset)[combinations[, 1]],
	sp2=colnames(dataset)[combinations[, 2]],
	odds=rep(NA, dim(combinations)[1]))
combinations.final$odds<-as.numeric(combinations.final$odds)
  
# run odds.ratio calculation for each pair
n.rows<-dim(combinations.final)[1]
cat(paste("calculating", dim(combinations.final)[1], "combinations of", n.species, "species.", sep=" ")) 
if(quiet==FALSE){
	pings<-seq(20, n.rows, 20)
	cat(" n(rows) calculated = ")}

# run loop
for(i in 1:dim(combinations.final)[1]){
	if(quiet==FALSE){	# create an output meter
  		if(any(pings==i)){	
    		thisrun<-which(pings==i)
    		cat(paste(pings[thisrun], " ... ", sep=""))}}

	# then run analysis
	dataset.thisrun<-dataset[, c(combinations.final$col1[i], combinations.final$col2[i])]
	or.simple<-or.contingency(dataset.thisrun)	# run as a test: req. as glmer will fail if no overlap in presences.
	if(missing(random.effect)){combinations.final$odds[i]<-or.simple
	}else{
		if(any(c(0, Inf)==or.simple)){combinations.final$odds[i]<-or.simple
		}else{combinations.final$odds[i]<-or.regression(dataset.thisrun, random.effect)}
	}	# end if(random effect specified)
}	# end loop (i.e. all combinations of spp.)
  
# create and return objects necessary for plotting
output<-list(
	frequency=frequency.result,
	result.long=combinations.final,
	result.wide=pairwise.or.matrix(combinations.final))

return(output)
}


#' Pairwise odds ratio matrix
#'
#' Internal function to convert the 'long' output from pairwise.odds.ratios() into a matrix
#' @param or.object data.frame containing columns listing the two species being compared, and the odds ratio between them
#' @return A matrix
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
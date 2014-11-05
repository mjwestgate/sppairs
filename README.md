# sppairs
# Species Pairwise Association Analysis (SPAA)

SPAA uses odds ratios to calculate and visualise pairwise associations between species from presence/absence data. It is a simple method that can be applied to site by species matrices in a range of contexts, and provides an alternative to dissimilarity-based approaches. This package gives a simple implementation of the SPAA approach, which you can read a  full description of in this article:

    Lane, P.W., D.B. Lindenmayer, P.S. Barton, W. Blanchard & M.J. Westgate (2014) 
    "Visualisation of species pairwise associations: a case study of surrogacy in 
    bird assemblages" Ecology and Evolution. DOI: 10.1002/ece3.1182
    
Article available here (OA): http://onlinelibrary.wiley.com/doi/10.1002/ece3.1182/full   
Package details:   
Website: http://martinwestgate.com/code/sppairs/   
Authors: Martin Westgate (<martinjwestgate@gmail.com>) and Peter Lane   

# Example

```
# Load required packages
library(devtools)
install_github('sppairs', 'mjwestgate')
library(sppairs)

# get some data
library(cooccur)
data(beetles)

# Test odds ratio calculation for a single set of species
# First, try a simple approach using contingency tables
or.contingency(beetles[, c(1, 10)])
# Alternatively, the odds ratios can be calculated from logistic regression:
or.glm(beetles[, c(1, 10)])	
# using mixed models to account for repeat visits (requires lme4)
groups<-as.factor(rep(c(1:10), length.out=dim(beetles)[1])) # create grouping variable for example purposes ONLY
or.glmer(beetles[, c(1, 10)], random.effect=groups)

# Calculate odds ratios for all pairs of species.
or.test<-spaa(beetles, rarity.cutoff=0.1)
# NOTE: 
	# 1. This can take a long time; particularly for many spp.
	# 2. glmer in lme4 v1+ often gives error messages. 


# Some basic summaries:
# Histogram of % sites occupied
hist(or.test$species$frequency, 
	las=1, xlim=c(0, 1),
	xlab="proportion", ylab="number of species", main="proportion of sites occupied by each spp.")
# Histogram of odds ratios
hist(or.test$combinations$odds,
	las=1,
	breaks=c(c(0,xlim=c(0, 4),  0.333, 0.666), 0, seq(1, 4, 0.5), Inf), 
	xlab="odds ratio", main="Results of pairwise.odds.ratio() for all sites and spp.")
# how many spp. in each indicator category?
length(which(or.test$combinations$odds>3))	#  n pairs of significant indicators
length(which(or.test$combinations$odds<0.33))	# n pairs of contra-indicators
or.test$combinations[which(or.test$combinations$odds>3), ] # show pairs of species that meet +ve indicator criteria


# Draw a simple plot
plot(or.test)
plot(or.test, draw.frequencies="none", add.key="none") # Simple version, no key

# Alternatively, calculate point and line values, then plot
point.coords<-spaa.points(or.test)	# avoids issue whereby igraph gives a different arrangement each time
line.values<-spaa.lines(or.test)		# sets line properties in a sensible way.
plot.spaa(list(point.coords, line.values))	# default behaviour (as above)

# Customize plot appearance
key.matrix<-matrix(data=c(
0, 0.85, 0, 1,
0.85, 1, 0, 1,
0, 0.2, 0.02, 0.32,	# 1. bottom left	
0.18, 0.40, 0.02, 0.22),
nrow=4, ncol=4, byrow=TRUE)
rownames(key.matrix)<-c("network", "species", "points", "lines")

plot.spaa(list(point.coords, line.values), plot.control=list(
key.placement=key.matrix,
line.cols=brewer.pal(5, "Blues")[2:5],
line.breaks=c(3, 5, 7, 9, Inf),
line.widths=rep(2, 4),	# note: length(line.widths) must = length(line.cols)
point.breaks=c(0, 0.2, 0.3, 0.4, 0.5, 1))
)
```

SPAA uses odds ratios to calculate and visualise pairwise associations between species from presence/absence data. It is a simple method that can be applied to site by species matrices in a range of contexts, and provides an alternative to dissimilarity-based approaches.

This package gives a simple implementation of Lane's approach. *Improvements will be made over time, but you should treat this as in beta for the time being.* You can read a full description of SPAA in this article:

    Lane, P.W., D.B. Lindenmayer, P.S. Barton, W. Blanchard & M.J. Westgate (2014) 
    "Visualisation of species pairwise associations: a case study of surrogacy in 
    bird assemblages" Ecology and Evolution. DOI: 10.1002/ece3.1182
    
Available here: http://onlinelibrary.wiley.com/doi/10.1002/ece3.1182/full

Author: Martin Westgate <martinjwestgate@gmail.com> & Peter Lane

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
# Alternatively, the symmetric odds ratio can be calculated from regression:
# using glm(); should give similar results to above
or.regression(beetles[, c(1, 10)])	
# using mixed models to account for repeat visits (requires lme4)
groups<-as.factor(rep(c(1:10), length.out=dim(beetles)[1])) # create grouping variable for example purposes ONLY
or.regression(beetles[, c(1, 10)], random.effect=groups)

# Calculate odds ratios for all pairs of species.
or.test<-pairwise.odds.ratios(beetles, rarity.cutoff=0.1)
# NOTE: 
	# 1. This can take a long time; particularly for many spp.
	# 2. glmer in lme4 v1+ often gives error messages. 

# Some basic summaries:
# Histogram of % sites occupied
hist(or.test$frequency[, 2], las=1, xlim=c(0, 1),
	xlab="proportion", ylab="number of species", main="proportion of sites occupied by each spp.")
# Histogram of odds ratios
hist(or.test$result.long$odds, breaks=c(c(0, 0.333, 0.666), 0, seq(1, 4, 0.5), Inf), xlim=c(0, 4), las=1,
	xlab="odds ratio", main="Results of pairwise.odds.ratio() for all sites and spp.")
# how many spp. in each indicator category?
length(which(or.test$result.long$odds>3))	#  n pairs of significant indicators
length(which(or.test$result.long$odds<0.33))	# n pairs of contra-indicators
or.test$result.long[which(or.test$result.long$odds>3), ] # show pairs of species that meet +ve indicator criteria

# Calculate point and line values for plot
point.coords<-or.points(or.test)	# note: igraph gives a different arrangement each time
line.values<-or.lines(or.test)		# sets line properties in a sensible way.

# Draw
plot.pairs(points=point.coords, lines=line.values, draw.frequencies="none", add.key="none") # Simple version, no key
plot.pairs(points=point.coords, lines=line.values)	# default behaviour (point, line behaviour as in text).
plot.pairs(points=point.coords, lines=line.values, draw.frequencies="both", add.key="both") # with a simple line key as well
```

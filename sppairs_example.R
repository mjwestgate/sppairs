# This script uses naa to generate some of the results from the paper describing that technique.

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
# create grouping variable for example purposes ONLY
groups<-as.factor(rep(c(1:10), length.out=dim(beetles)[1]))
or.regression(beetles[, c(1, 10)], random.effect=groups)

# Calculate odds ratios for all pairs of species.
or.test<-pairwise.odds.ratios(beetles, rarity.cutoff=0.1)
# NOTES: 
	# 1. This can take a long time; particularly for many spp.
	# 2. glmer in lme4 v1+ often gives error messages. 

# Some basic summaries:
# Histogram of % sites occupied
hist(or.test$frequency[, 2], las=1, xlim=c(0, 1),
	xlab="proportion", ylab="number of species", main="proportion of sites occupied by each spp.")
# Histogram of odds ratios
hist(or.test$result.long$odds, breaks=c(c(0, 0.333, 0.666), seq(1, 4, 0.5), Inf), xlim=c(0, 4), las=1,
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
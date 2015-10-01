# sppairs
# Species Pairwise Association Analysis (SPAA)

SPAA was designed to calculate the degree of spatial association/co-occurrence between species from presence/absence data. The current implementation allows users to call any function that might plausibly be used to measure association between each pair of vectors. This allows users to call the original odds ratio -based methods; to call other functions such as cor(); or to build and call their own functions. A description of the original method is available in this article:

    Lane, P.W., D.B. Lindenmayer, P.S. Barton, W. Blanchard & M.J. Westgate (2014) 
    "Visualisation of species pairwise associations: a case study of surrogacy in 
    bird assemblages" Ecology and Evolution. DOI: 10.1002/ece3.1182
    
Article available here (OA): http://onlinelibrary.wiley.com/doi/10.1002/ece3.1182/full   
Package details:   
Website: http://martinwestgate.com/code/sppairs/   
Authors: Martin Westgate (<martinjwestgate@gmail.com>) and Peter Lane   

Note: Automated plotting of pairwise association networks is retained in this version for backwards-compatability, but no longer retains any functionality beyond that shipped with igraph. Users should consider calling igraph directly for cleaner results.

# Example

```
# Load required packages
library(devtools)
install_github('mjwestgate/sppairs')
library(sppairs)

# get some data
library(cooccur)
data(beetles)

# Test odds ratio calculation for a single set of species
# First, try a simple approach using contingency tables
or.symmetric(beetles[, c(1, 10)])
# Alternatively, the odds ratios can be calculated from logistic regression:
or.glm(beetles[, c(1, 10)])	
# using mixed models to account for repeat visits (requires lme4)
groups<-as.factor(rep(c(1:10), length.out=dim(beetles)[1])) # create grouping variable for example purposes ONLY
or.glmer(beetles[, c(1, 10)], random.effect=groups)

# Calculate odds ratios for all pairs of species.
or.test<-spaa(clean.dataset(beetles, min.cutoff=0.1))
# NOTE: 
	# 1. This can take a long time; particularly for many spp.
	# 2. glmer in lme4 v1+ often gives error messages. 

```

# sppairs
# Species Pairwise Association Analysis (SPAA)

NOTE: I DO NOT RECOMMEND THAT YOU USE THIS (BETA) VERSION OF SPAA.
This version is:
1. greatly modified from earlier versions - do not expect earlier code to work at all, or if it does, to give similar results
2. in active development - do not expect stable (or even sensible) results at present

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

```

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
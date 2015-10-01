# function to create a quick plot of results
plot.spaa<-function(input, ...){
# set some defaults
method<-attr(input, "association.function")
or.methods<-c("or.symmetric", "or.asymmetric", "or.glm", "or.glmer")
if(is.null(method)==FALSE & any(or.methods== method)){
	line.breaks<-c(0, 1/6, 1/3, 3, 6, Inf)
	line.colors<-brewer.pal(5, "RdBu")[5:1] 
	line.colors[3]<-NA
}else{ # i.e. if you are just provided data and it doesn't state it is OR-related
	span.zero<-min(input[, 3], na.rm=TRUE)<0 & max(input[, 3], na.rm=TRUE)>0
	if(span.zero){ # when there are +ve and -ve values only
		max.val<-sqrt(max(input[which(input[, 3]!=Inf), 3]^2, na.rm=TRUE))
		seq.vals<-seq(0, max.val, length.out=3)
		line.breaks<-c(-seq.vals[3:2], seq.vals)
		line.colors<-brewer.pal(5, "RdBu")[5:1] 
		if(any(input[, 3]==Inf))line.breaks[5]<-Inf
		if(any(input[, 3]==-Inf))line.breaks[1]<-(-Inf)
	}else{	 # when only +ve or -ve values
		line.breaks<-seq(
			min(input[which(input[, 3]!=-Inf), 3], na.rm=TRUE),
			max(input[which(input[, 3]!=Inf), 3], na.rm=TRUE),
			length.out=4)
		if(any(input[, 3]==Inf))line.breaks[4]<-Inf
		if(any(input[, 3]==-Inf))line.breaks[1]<-(-Inf)
		line.colors<-brewer.pal(3, "Purples") 
	}
}
input$color<-line.colors[cut(input[, 3], breaks=line.breaks, labels=FALSE)]
graph.data<-graph_from_data_frame(input, directed=TRUE, vertices=NULL)
plot(graph.data, ...)
}


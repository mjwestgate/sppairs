# Plot results of pairwise.odds.ratios, using results from or.points and or.lines
plot.spaa<-function(object, plot.control, graph.function, ...)
{
# set behaviour for different inputs
if(class(object)=="spaa"){
	if(missing(plot.control)){
		threshold<-3
	}else{
		if(length(plot.control$threshold)==0){threshold<-3
		}else{threshold<-plot.control$threshold}
	}
	input<-list(
		points=spaa.points(object, threshold, graph.function, ...),
		lines=spaa.lines(object, threshold))
}else{if(class(object)=="list"){
	input<-object
	if(length(names(input))==0){names(input)<-c("points", "lines")}
	}}

# set error messages
if(dim(input$points)[1]==0)stop("Error: no strong inter-species associations identified")

# set plot attributes
input<-set.plot.attributes(input, plot.control) # set plot attributes/defaults
draw.sppairs(input)
invisible(input)
}	# end plot.spaa



# a key to species labels
draw.spp.key<-function(input, labels){
	label.size<-max(nchar(input$points$species))
	if(label.size<10){x.min<-0.8}else{x.min<-0.5}

	# work out y values
	y.key<-rep(0, dim(input$points)[1])
	for(i in 2:dim(input$points)[1]){
		y.key[i]<-y.key[(i-1)]+sum(c(input$points$radius[c((i-1), i)], input$plot.control$point.unit*0.5))}
	y.key<-max(y.key)-y.key
	y.min<-0-input$points$radius[dim(input$points)[1]]
	y.max<-max(y.key)+input$points$radius[1]

	# draw plot
	par(mar=rep(0, 4), cex=0.7)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(x.min, 1.1), ylim=c(y.min, y.max), asp=1)
	mtext("Species", side=3, font=2, cex=input$plot.control$text.size*0.7, line=-1, adj=0.5)

	# add points
	for(i in 1:dim(input$points)[1]){
		draw.circle(x=1, y=y.key[i], 
			r=input$points$radius[i], 
			bg=input$points$col[i], 
			col=input$plot.control$point.label[1])} 
	# and labels for those points
	text(x=1, y=y.key,
		labels=c(1:dim(input$points)[1]), 
		cex=input$plot.control$text.size, 
		col=input$plot.control$point.label[2])

	# add descriptions
	if(missing(labels)){species.labels<-input$points$species
	}else{species.labels<-labels}

	text(x=(1-max(input$points$radius)), y=y.key, 
		labels=species.labels, 
		cex=input$plot.control$text.size,
		pos=2)	
	}



# Function to add point key to plot
draw.point.key<-function(input){

	# set labels
	# dataset<-input$plot.control
	nrow<-length(input$plot.control$point.breaks)
	values<-format(input$plot.control$point.breaks)
	point.labels<-paste(values[1:(nrow-1)], values[2:nrow], sep=" - ")
	npoints<-nrow-1

	# work out y values
	y.key<-rep(0, npoints)
	point.radii<-input$plot.control$point.size*input$plot.control$point.unit
	for(i in 2:npoints){
		y.key[i]<-y.key[(i-1)]+sum(c(point.radii[c((i-1), i)], input$plot.control$point.unit*0.5))}
	y.key<-max(y.key)-y.key
	y.min<-0-point.radii[length(point.radii)]
	y.max<-max(y.key)+point.radii[1]
	x.max<-1+max(point.radii)

	# set up plot
	par(mar=c(0, 0, 1, 0), cex=0.7)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(0.5, x.max), ylim=c(y.min, y.max), asp=1)

	# add points
	for(i in 1:npoints){
		draw.circle(x=1, y=y.key[i], 
			r= point.radii[i], 
			bg=input$plot.control$point.cols[i], 
			col=input$plot.control$point.label[1])} 

	text(x=(1-max(point.radii)), y= y.key, pos=2, 
		cex=input$plot.control$text.size, labels=point.labels, col=input$plot.control$point.label[2])
	mtext("Occupancy", side=3, font=2, cex=input$plot.control$text.size*0.7, adj=0.5, line=0)
	}



# Function to draw key to lines
draw.line.key<-function(input, labels){
	dataset<-input$plot.control
	nrow<-length(dataset$line.cols)

	# work out how to display odds ratios
	values<-dataset$line.breaks

	# ID unplotted values (between 1/threshold - threshold
	# threshold.loc<-which(values==dataset$threshold)
	inv.threshold.loc<-which(values==(1/dataset$threshold))

	# format for printing	
	values2<-c(1/values[values<1], values[values>=1])
	low.vals<-which(values2<100)
	high.vals<-which(values2>=100)
	low.c<-as.character(as.integer(values2[low.vals]))
	high.c<-as.character(values2[high.vals])
	simple.values<-c(high.c, low.c)[order(c(high.vals, low.vals))]	
	# add 1/ to low values
	neg.vals<-which(values<1)
	simple.values[neg.vals]<-paste("1/", simple.values[neg.vals], sep="")
	if(simple.values[1]=="Inf/1")simple.values[1]<-0

	if(missing(labels)){
		line.labels<-paste(simple.values[1:nrow], simple.values[2:(nrow+1)], sep=" - ")[-inv.threshold.loc]
	}else{line.labels<-labels}
	plot.cols<-dataset$line.cols[-inv.threshold.loc]
	plot.widths<-dataset$line.widths[-inv.threshold.loc]

	# add plot
	par(mar=c(0, 3, 1, 0), cex=0.7)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(0, 1), ylim=c(0, (nrow-0.7)))
	for(i in 1: nrow){
		lines(x=c(0.1, 1), y=rep(c(1: nrow)[i], 2), 
			col=plot.cols[i], lwd=plot.widths[i], lend="butt")}
	axis(2, at=c(1:(nrow-1)), labels=line.labels, lwd=0, line=-1, 
		cex.axis=input$plot.control$text.size, las=1)
	mtext("Odds ratio", side=3, font=2, line=0, cex=input$plot.control$text.size*0.7, adj=0.5)
	}


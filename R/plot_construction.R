# function to set plot defaults, and overwrite if new data is provided (altered from circleplot)
	# NOTE: the above does not allow arbitrary point colours as in circleplot.
set.plot.attributes<-function(
	input,
	plot.control,
	draw.frequencies
	)
	{

	# set defaults depending on whether sets use.frequencies or not
	if(draw.frequencies){
		point.breaks<-seq(0, 1, 0.2) #c(0, 0.1, 0.25, 0.5, 0.75, 1)
		point.cols<-brewer.pal(5, "Spectral")[5:1]
		point.size<-seq(1, 2, length.out=5)
	}else{
		point.breaks<-c(0, 1); point.cols<-"white"; point.size<-1}

	# work out point sizes
	range.values<-c(max(input$points$x)-min(input$points$x), 
		max(input$points$y)-min(input$points$y))
	point.unit<-max(range.values)*0.01

	# set some defaults
	plot.defaults<-list(
		threshold=3,
		point.label="grey30",	# use for boundary and label
		point.cols=point.cols,
		point.breaks=point.breaks,
		point.size=point.size,
		point.unit= point.unit,
		line.breaks=c(0, 0.000001, 1/9, 1/6, 1/3, 3, 6, 9, 10^8, Inf),
		line.cols=c("black", brewer.pal(7, "RdBu")[7:1], "magenta"),
		line.widths=rep(2, 9),
		text.size=0.7,
		key.placement=matrix(data=NA, nrow=1, ncol=4)
		)

	# overwrite these values where others are provided
	if(missing(plot.control)==FALSE){
		names.provided<-names(plot.control)
		for(i in 1:length(plot.defaults)){
			if(any(names.provided==names(plot.defaults)[i])){
			entry.thisrun<-which(names.provided==names(plot.defaults)[i])
			plot.defaults[[i]]<-plot.control[[entry.thisrun]]
			}}}

	# set point colours and sizes according to breaks given in plot.defaults
	# note that these vary categorically, not continuously as before
	point.categories<-cut(input$points$freq, breaks=plot.defaults$point.breaks, 
		include.lowest=TRUE, labels=FALSE)
	input$points$col<-plot.defaults$point.cols[point.categories]
	#if(length(plot.defaults$point.size)==1){plot.defaults$point.size<-rep(plot.defaults$point.size, 2)}
	#point.sizes<-seq(plot.defaults$point.size[1], plot.defaults$point.size[2], 
	#	length.out=length(plot.defaults$point.cols))	# change to make line.width and point.size behaviour the same
	input$points$cex<-plot.defaults$point.size[point.categories]	
	input$points$radius<-input$points$cex*point.unit

	# set line values
	line.classes<-cut(input$lines$odds, breaks= plot.defaults$line.breaks, 
		include.lowest=TRUE, labels=FALSE)
	input$lines$col<-plot.defaults$line.cols[line.classes]
	input$lines$lwd<-plot.defaults$line.widths[line.classes]

	# set line order
	input$lines$order<-input$lines$odds
	low.rows<-which(input$lines$odds<1)
	input$lines$order[low.rows]<-1/input$lines$order[low.rows]
	input$lines<-input$lines[order(input$lines$order, decreasing=FALSE), ]

	# add plot.control to input
	input$plot.control<-plot.defaults

	return(input)
	}
	


# Plot results of pairwise.odds.ratios, using results from or.points and or.lines
plot.spaa<-function(object, plot.control, draw.frequencies)
{
# set defaults
if(missing(draw.frequencies))draw.frequencies<-TRUE

# set behaviour for different inputs
if(class(object)=="spaa"){
	if(missing(plot.control)){threshold<-3
	}else{
		if(length(plot.control$threshold==0)){threshold<-3
		}else{threshold<-plot.control$threshold}
	}
	input<-list(
		points=spaa.points(object, threshold),
		lines=spaa.lines(object, threshold))
} else if(class(object)=="list"){
	input<-object
	if(length(names(input))==0){names(object)<-c("points", "lines")}
	}

# set error messages
if(dim(input$points)[1]==0)stop("Error: no strong inter-species associations identified")

# set plot attributes
input<-set.plot.attributes(input, plot.control, draw.frequencies) # set plot attributes/defaults

# set draw behaviour depending on key requests
if(any(is.na(input$plot.control$key.placement))){draw.sppairs(input)
}else{
	place.matrix<-input$plot.control$key.placement
	command.list<-rownames(place.matrix)
	split.screen(place.matrix)
	for(i in 1:dim(place.matrix)[1]){
		screen(i)
		switch(command.list[i],
			network={draw.sppairs(input)},
			species={draw.spp.key(input)},
			points={draw.point.key(input)},
			lines={draw.line.key(input)})
		}
	close.screen(all.screens=TRUE)
	}
}	# end plot.pairs




# draw a plot of the netowrk given by sppairs
draw.sppairs<-function(input)
{
# how much should we adjust point sizes etc.?
interpoint.dist<-mean(c((max(input$point$x)-min(input$points$x)), 
	(max(input$point$y)-min(input$point$y))))

# set plot region
par(mar=rep(0.5, 4), cex=0.7)
plot(x=input$points$x, y=input$points$y, type="n", ann=F, axes=F, asp=1)

# run loop to draw lines showing significant effects
for(i in 1:dim(input$lines)[1])
	{
	# which species are connected by this line?
	sp1<-which(input$points$species==input$lines$sp1[i])
	sp2<-which(input$points$species==input$lines$sp2[i])
	dataset<-input$points[c(sp1, sp2), ]
	arrow.points<-find.new.points(dataset)	
	arrows.default(arrow.points, col=input$lines$col[i], 
		lwd=input$lines$lwd[i], code=input$lines$arrow.code[i])
	}
  
# add points
for(i in 1:dim(input$points)[1]){
	draw.circle(input$points$x[i], input$points$y[i], 
		r=input$points$radius[i], 
		bg=input$points$col[i], 
		col=input$plot.control$point.label)} 
# and labels for those points
text(input$points$x, input$points$y, 
	labels=c(1: dim(input$points)[1]), 	# NOTE: consider changing this for consistency between plots
	cex=input$plot.control$text.size, col=input$plot.control$point.label)

}	# end draw.sppairs()



# a key to species labels
draw.spp.key<-function(input){
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
			col=input$plot.control$point.label)} 
	# and labels for those points
	text(x=1, y=y.key,
		labels=c(1:dim(input$points)[1]), 
		cex=input$plot.control$text.size, 
		col=input$plot.control$point.label)

	# add descriptions
	text(x=(1-max(input$points$radius)), y=y.key, 
		labels=input$points$species, 
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
			col=input$plot.control$point.label)} 

	text(x=(1-max(point.radii)), y= y.key, pos=2, 
		cex=input$plot.control$text.size, labels=point.labels)
	mtext("Occupancy", side=3, font=2, cex=input$plot.control$text.size*0.7, adj=0.5, line=0)
	}



# Function to draw key to lines
draw.line.key<-function(input){
	dataset<-input$plot.control
	nrow<-length(dataset$line.cols)

	# work out how to display odds ratios
	values<-dataset$line.breaks
	values2<-c(1/values[values<1], values[values>=1])
	low.vals<-which(values2<100)
	high.vals<-which(values2>=100)
	low.c<-as.character(as.integer(values2[low.vals]))
	high.c<-as.character(values2[high.vals])
	simple.values<-c(high.c, low.c)[order(c(high.vals, low.vals))]	
	# add 1/ to low values
	neg.vals<-which(values<1)
	simple.values[neg.vals]<-paste(simple.values[neg.vals], "/1", sep="")
	if(simple.values[1]=="Inf/1")simple.values[1]<-0

	line.labels<-paste(simple.values[1:nrow], simple.values[2:(nrow+1)], sep=" - ")

	# add plot
	par(mar=c(0, 3, 1, 0), cex=0.7)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(0, 1), ylim=c(0, nrow))
	for(i in 1: nrow){
		lines(x=c(0.1, 1), y=rep(c(1: nrow)[i], 2), 
			col= dataset$line.cols[i], lwd= dataset$line.widths[i], lend="butt")}
	axis(2, at=c(1:nrow), labels=line.labels, lwd=0, line=-1, 
		cex.axis=input$plot.control$text.size, las=1)
	mtext("Odds ratio", side=3, font=2, line=0, cex=input$plot.control$text.size*0.7, adj=0.5)
	}


# function to set plot defaults, and overwrite if new data is provided (altered from circleplot)
	# NOTE: the above does not allow arbitrary point colours as in circleplot.
set.plot.attributes<-function(
	input,
	plot.control
	)
	{
	# work out point sizes
	range.values<-c(max(input$points$x)-min(input$points$x), 
		max(input$points$y)-min(input$points$y))
	point.unit<-max(range.values)*0.01

	# set some defaults
	plot.defaults<-list(
		threshold=3,
		point.label=rep("grey30", 2),	# point border and label colour, in that order; or a single value
		point.cols=brewer.pal(5, "Spectral")[5:1],
		point.breaks=seq(0, 1, 0.2),
		point.size=seq(1, 2, length.out=5),
		line.breaks=c(0, 0.000001, 1/9, 1/6, 1/3, 3, 6, 9, 10^8, Inf),
		line.cols=c("black", brewer.pal(7, "RdBu")[7:1], "magenta"),
		line.widths=rep(1.5, 9),
		line.curvature=1.5,
		arrow.size=1.5,
		arrow.angle=20,
		text.size=0.7
		)

	# USER CONTROL
	# overwrite these values where others are provided
	if(missing(plot.control)==FALSE){
		names.provided<-names(plot.control)
		for(i in 1:length(plot.defaults)){
			if(any(names.provided==names(plot.defaults)[i])){
			entry.thisrun<-which(names.provided==names(plot.defaults)[i])
			if(length(entry.thisrun)>1){
				stop(paste("multiple values supplied for plot.control$", names(plot.defaults)[i], sep=""))
			}else{plot.defaults[[i]]<-plot.control[[entry.thisrun]]}
			}}
		}

	# ERROR CATCHING
	# allow 'global' setting of point and line attributes,
	# via single values as inputs to point.size, point.cols, point.label, line.width
	if(length(plot.defaults$point.label)==1){plot.defaults$point.label<-rep(plot.defaults$point.label, 2)}
	if(length(plot.defaults$point.size)==1){
		plot.defaults$point.size<-rep(plot.defaults$point.size, length(plot.defaults$point.breaks)-1)}
	if(length(plot.defaults$point.cols)==1){
		plot.defaults$point.cols<-rep(plot.defaults$point.cols, length(plot.defaults$point.breaks)-1)}
	if(length(plot.defaults$line.widths)==1){
		plot.defaults$line.widths<-rep(plot.defaults$line.widths, length(plot.defaults$line.breaks)-1)}

	# cut point values as needed
	point.categories<-cut(input$points$freq, breaks=plot.defaults$point.breaks, 
		include.lowest=TRUE, labels=FALSE)
	input$points$col<-plot.defaults$point.cols[point.categories]
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

	# append unit size to plot.control
	plot.defaults<-append(plot.defaults, list(point.unit= point.unit)) # should not be available to users

	# use point.unit to scale arrow size
	plot.defaults$arrow.size<-plot.defaults$arrow.size*plot.defaults$point.unit

	# EXPORT
	# add plot.control to input
	input$plot.control<-plot.defaults

	return(input)
	}
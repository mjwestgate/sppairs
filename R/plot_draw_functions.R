# create a function to draw points as polygons given location and attribute data
draw.circle<-function(x, y, r, bg, col){
	degrees<-seq(0, 360, length.out=100)
	radians<-degrees*(pi/180)
	x.new<-x + (r * cos(radians))
	y.new<-y + (r * sin(radians))
	polygon(x.new, y.new, col=bg, border=col)
	}


# create a function that draws an arrowhead with the specified direction
arrow.coords<-function(
	direction, # left or right
	angle,	# how 'sharp' is the arrowhead?
	length,
	location,
	rotation,
	reverse
	)
	{
	# set some defaults
	if(missing(angle))angle<-20		# note this gives a total angle of 40 deg. (2*angle)
	if(missing(length))length<-1
	if(missing(location))location<-c(0, 0)
	if(missing(rotation))rotation<-0	# assume units in radians
	centre<-c(0, 0)
	angle<-angle*(pi/180)	# assume units are in degrees, and convert to radians

	# calculate x and y coordinates
	xlim<-c(centre[1]-(length*0.5), centre[1]+(length*0.5))	
	height<-length*tan(angle)
	ylim<-c(centre[2]-height, centre[2]+height)	

	# arrange for a left-facing arrow
	arrow<-data.frame(
		x=c(xlim[1], xlim[2], xlim[2], xlim[1]),
		y=c(centre[2], ylim[1], ylim[2], centre[2]))
	if(reverse){arrow$x<-arrow$x*-1}	# flip to face right if required

	# rotate
	arrow$new.x<-(arrow$x * cos(rotation))-(arrow$y * sin(rotation))
	arrow$new.y<-(arrow$y * cos(rotation))+(arrow$x * sin(rotation))

	# adjust position to match location
	arrow$new.x<-arrow$new.x + as.numeric(location[1])
	arrow$new.y<-arrow$new.y + as.numeric(location[2])

	# ensure only one set of x and y are given, return
	arrow<-arrow[, 3:4]; colnames(arrow)<-c("x", "y")
	return(arrow)
	}
#plot(1~1, type="n", bty="l", ann=FALSE, axes=FALSE, xlim=c(-1, 1), ylim=c(-1, 1))
#polygon(arrow.coords(direction=2, angle=30*(pi/180)), col="blue", border=NA)


# function to return a curved line, if a straight line would overlap other points
line.adjust<-function(
	all.data,
	curvature
	)
	{
	#(1) find the slope of the line between two points
	vertex.rows<-which(all.data$match==1)
	model<-lm(y~x, data= all.data[vertex.rows, ])
	midpoint.x<-mean(all.data$x[vertex.rows])

	# set remaining behaviour, allowing different rules if points are directly above or below one another
	#if(is.na(coef(model)[2])){coef.thisrun<-0
	#}else{
	coef.thisrun<-coef(model)[2]  #}
	midpoint.y<-as.numeric(coef(model)[1]+(coef.thisrun*midpoint.x))

	# 2) set the midpoint of that line as the origin
	all.data$origin.x<-all.data$x-midpoint.x
	all.data$origin.y<-all.data$y-midpoint.y

	# 3) rotate all points such that this line is horizontal
	theta<-as.numeric(-atan(coef.thisrun))
	all.data$new.x<-(all.data$origin.x * cos(theta))-(all.data$origin.y * sin(theta))
	all.data$new.y<-(all.data$origin.y * cos(theta))+(all.data$origin.x * sin(theta))
	# note: new.x is used instead of x because your 2-line approach otherwise causes errors

	# 3a. exclude all points outside x range of (rotated) input data
	exclude.test1<-all.data$new.x > max(all.data$new.x[all.data$match==1])
	if(any(exclude.test1)){all.data<-all.data[-which(exclude.test1==TRUE), ]}
	exclude.test2<-all.data$new.x < min(all.data$new.x[all.data$match==1])
	if(any(exclude.test2)){all.data<-all.data[-which(exclude.test2==TRUE), ]}

	# 4) show boundaries of points, continue if there is a need for curved lines
	size.test<-dim(all.data)[1]>2 
	all.data$ysq<-sqrt(all.data$new.y^2)-(all.data$radius*curvature)		
	# note: 1.5 is a proximity ratio - closer than half of radius, overlap =TRUE
	overlap.test<-any(c(all.data$match==0 & all.data$ysq<0))

	if(c(size.test & overlap.test)){
	# 5) if a line passes too close to one or more points, curve it to maximize distance from all points
		# accomplish this by passing the smallest distance possible from the closest point
		col<-which(colnames(all.data)=="new.x")
		match.points<-all.data[all.data$match==1, c(col, (col+1))]
		included.rows<-which(all.data$match==0)
		row<-c(1:dim(all.data)[1])[included.rows][which.min(all.data$ysq[included.rows])]
		clearance<-all.data$radius[row]*curvature
		new.point<-c(all.data$new.x[row], all.data$new.y[row]+(sign(all.data$new.y[row])* -1 * clearance))
		match.points<-rbind(match.points, new.point)
		colnames(match.points)<-c("x", "y")

		# model using quadratic curve
		model2<-lm(y ~ x + I(x^2), data=match.points)
		x.vals<-seq(match.points$x[1], match.points$x[2], length.out=101)
		y.vals<-coef(model2)[1] + (coef(model2)[2] * x.vals) + (coef(model2)[3] * (x.vals^2)) 
		parabola.initial<-data.frame(x=x.vals, y=y.vals)

		# rotate back to original orientation
		parabola.final<-parabola.initial
		parabola.final$x<-(parabola.initial$x * cos(-theta))-(parabola.initial$y * sin(-theta))
		parabola.final$y<-(parabola.initial$y * cos(-theta))+(parabola.initial$x * sin(-theta))

		# relocate back to initial origin
		parabola.final$x<-parabola.final$x+ midpoint.x
		parabola.final$y<-parabola.final$y+ midpoint.y	
		result<-parabola.final

	}else{ # if there is no overlap, return a straight line - no transformation required.
		new.vertices<-which(all.data$match==1) # as number of rows has changed
		x.final<-seq(min(all.data$x[new.vertices]), max(all.data$x[new.vertices]), length.out=101)
		#if(is.na(coef(model)[2])){
		#	y.final<-seq(min(all.data$y[new.vertices]), max(all.data$y[new.vertices]), length.out=101)
		#}else{
			y.final<-predict.lm(model, newdata=data.frame(x=x.final), se.fit=FALSE)  # }
		result<-data.frame(x=x.final, y=as.numeric(y.final))
		}

	# return a list
	result<-list(angle=-theta, line=result)

	return(result)
	}	# end points.adjust




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
	dataset<-input$points
	dataset$match<-0
	dataset$match[c(sp1, sp2)]<-1
	result<-line.adjust(dataset, input$plot.control$line.curvature)
		lines(result$line$x, result$line$y, col=input$lines$col[i], lwd=input$lines$lwd[i])

	# calculate arrows
	# is input$lines$sp1[i] to the left of sp2? If not, angle should be reversed
	reverse<-input$points$x[sp2] > input$points$x[sp1]
	# angle given by line.adjust goes from sp1 > sp2
	#}else{ 	# the reverse is true

	# add an arrow
	if(input$lines$arrow.code[i]==1){
		arrow.thisrun<-arrow.coords(
			angle=input$plot.control$arrow.angle,
			direction=input$lines$arrow.code[i],
			length=input$plot.control$arrow.size,
			location=as.numeric(result$line[51, ]),
			rotation=result$angle,
			reverse=reverse)
		polygon(arrow.thisrun, col=input$lines$col[i], border=NA)
		}
	}	# end lines
  
# add points
for(i in 1:dim(input$points)[1]){
	draw.circle(input$points$x[i], input$points$y[i], 
		r=input$points$radius[i], 
		bg=input$points$col[i], 
		col=input$plot.control$point.label[1])} 
# try different approach?
#do.call(draw.circle, input$points...

# and labels for those points
if(is.na(input$plot.control$point.label[2])==FALSE){
	text(input$points$x, input$points$y, 
		labels=c(1: dim(input$points)[1]), 	# NOTE: consider changing this for consistency between plots
		cex=input$plot.control$text.size, col=input$plot.control$point.label[2])}

}	# end draw.sppairs()



# Function used in draw.sppairs to find the edges of a circle that intersect a line
# OBSOLETE
find.new.points<-function(dataset){
	# use lm to find the slope of the line joining these points
	model<-lm(y~x, data=dataset)
	theta<-atan(coef(model)[2])
		
	# create a data.frame to export the new points
	edge.values<-as.data.frame(matrix(data=NA, nrow=4, ncol=3))
	colnames(edge.values)<-c("point", "x", "y")
	edge.values$point<-rep(c(1, 2), each=2)

	# calculate points on the circumference bisected by this line
	for(i in 1:2){
		rows<-list(c(1, 2), c(3, 4))[[i]]
		delta.x<-dataset$radius[i]*cos(theta)
		delta.y<-dataset$radius[i]*sin(theta)
		edge.values$x[rows]<-c(dataset$x[i]-delta.x, dataset$x[i]+delta.x)
		edge.values$y[rows]<-c(dataset$y[i]-delta.y, dataset$y[i]+delta.y)
		}

	# work out which set of points are closest together
	comparison.matrix<-as.matrix(dist(edge.values[, 2:3]))[1:2, 3:4]
	selection<-which.min(comparison.matrix)
	point1<-rep(c(1, 2), 2)[selection]; point2<-rep(c(3, 4), each=2)[selection]

	coordinates<-c(edge.values$x[point1], edge.values$x[point2], 
		edge.values$y[point1], edge.values$y[point2])
	return(coordinates)
}
# Internal function to get sensible line widths
lwd.scale<-function(vector){
  vector<-vector-min(vector)
  vector<-vector/max(vector)	# 0-1 range
  vector<-(vector*2)+1}	# 1-3 range



# Simple internal function to make basic arrows appear in a standard manner
arrows.default<-function(coordinates, ...){
  arrows(x0= coordinates[1], x1= coordinates[2], y0= coordinates[3], y1= coordinates[4],
         length=0.05, angle=30, ljoin=1, lend="butt", ...)}



# create a function to draw points as polygons given location and attribute data
draw.circle<-function(x, y, r, bg, col){
	degrees<-seq(0, 360, length.out=100)
	radians<-degrees*(pi/180)
	x.new<-x + (r * cos(radians))
	y.new<-y + (r * sin(radians))
	polygon(x.new, y.new, col=bg, border=col)
	}



# Function used in draw.sppairs to find the edges of a circle that intersect a line
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
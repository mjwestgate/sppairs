# Internal function to get sensible line widths
lwd.scale<-function(vector){
  vector<-vector-min(vector)
  vector<-vector/max(vector)	# 0-1 range
  vector<-(vector*2)+1}	# 1-3 range


# Simple internal function to make basic arrows appear in a standard manner
arrows.default<-function(coordinates, ...){
  arrows(x0= coordinates[1], x1= coordinates[2], y0= coordinates[3], y1= coordinates[4],
         length=0.05, angle=30, ...)}


# Function to shorten arrow lengths, so that arrowheads are always visible
shorten.line<-function(coordinates, reduction)	
{
if(missing(reduction))reduction<-0.05

# calculate length of the hypotenuse
delta.x<-(coordinates[2]-coordinates[1])
delta.y<-(coordinates[4]-coordinates[3])
hypotenuse.new<-sqrt(delta.x^2 + delta.y^2)-reduction

# calculate the angle between the hypotenuse and the adjacent
theta<-atan2(delta.y, delta.x)
x.new<-(hypotenuse.new*cos(theta))+coordinates[1]
y.new<-(hypotenuse.new*sin(theta))+coordinates[3]

# return new coordinates
return(c(
	x0=coordinates[1], x1=x.new,
	y0=coordinates[3],y1=y.new))
}


# Function to perpendicularly offset lines, for cases where there are arrows in both directions. Not yet implemented.
offset.line<-function(coordinates, offset)
{
if(missing(offset))offset<-0.05

# calulate slope of the line
delta.x<-(coordinates[2]-coordinates[1])
delta.y<-(coordinates[4]-coordinates[3])
theta<-atan2(delta.y, delta.x)

# calculate a perpendicular offset
new.theta<-theta+(pi*0.5)
x.offset<-offset*cos(new.theta)
y.offset<-offset*sin(new.theta)
coordinates.new<-c(
	coordinates[1:2]+x.offset,
	coordinates[3:4]+y.offset)

return(coordinates.new)			
}
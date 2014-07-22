# Use igraph to calculate sensible point locations
or.points<-function(dataset, threshold)
{
require(igraph)
if(missing(threshold))threshold<-3

# select input data
species.results<-dataset$frequency	# what are the species we are looking at?
or.matrix<-dataset$result.wide	# extract wide-format odds ratios from input data

# ID interacting spp.
result.positive<-make.binary(or.matrix, threshold=threshold)	# positively-interacting species
# and for contra-indicators
lower.threshold<-(1/threshold)
or.matrix.neg<-(lower.threshold-or.matrix)
result.negative<-make.binary(or.matrix.neg, threshold=0)

# combine -ve and +ve results into one matrix; ID most strongly interacting spp.
result<-as.matrix(result.positive)+ as.matrix(result.negative)
result.freq<-apply(
	cbind(colSums(result, na.rm=TRUE), rowSums(result, na.rm=TRUE)), 1, sum)
rows<-which(result.freq==0)
result<-result[-rows, -rows]
species.results<-species.results[-rows, ]	

# get plot coordinates using igraph
net<-graph.adjacency(result, mode="directed", weighted=TRUE, diag=FALSE)
layout<-layout.fruchterman.reingold(net)
points<-as.data.frame(cbind(species.results, layout))
colnames(points)[2:4]<-c("freq", "x", "y")

return(points)
}


# Function for setting the attributes of lines, to be drawn between points given by or.points()
or.lines<-function(dataset, threshold, reduce)
{
if(missing(threshold))threshold<-3
if(missing(reduce))reduce<-TRUE

# make some changes to result.long for plotting lines
line.data<-dataset$result.long[, 3:5]
line.data$sp1<-as.character(line.data$sp1)
line.data$sp2<-as.character(line.data$sp2)
# remove values outside of threshold range
line.data<-rbind(
	line.data[line.data$odds<=(1/threshold) & line.data$odds>=0, ],
	line.data[line.data$odds>=threshold & line.data$odds<=Inf, ])

# add new columns to line data to allow prettier plots
line.data<-cbind(line.data, as.data.frame(matrix(data=0, nrow=dim(line.data)[1], ncol=5)))
colnames(line.data)[4:8]<-c("difference", "arrow.code", "offset", "colour", "width")
line.data$arrow.code<-2

# run loop to get properties of lines representing strong effects
for(i in 1:dim(line.data)[1])
{
test<-c(line.data$sp2==line.data$sp1[i] & line.data$sp1==line.data$sp2[i])
if(any(test==TRUE)){
	other.row<-which(test==TRUE)
	if(other.row>i){		# selects only later rows
		# behaviour now depends on 'reduce'
		if(reduce){
			values<-line.data$odds[c(i, other.row)]
			line.data$arrow.code[i]<-3
			line.data$arrow.code[other.row]<-99
			line.data$odds[i]<-mean(values)
			line.data$difference[i]<-max(values)-min(values)   
		}else{	# i.e. if reduce==FALSE
			line.data$offset[c(i, other.row)]<-1
	}}}}

# if required, get rid of excess rows
if(reduce){if(any(line.data$arrow.code==99)){line.data<-line.data[-which(line.data$arrow.code==99), ]}}

# set colours & widths
rows.pos<-which(c(line.data$odds>threshold & line.data$odds!=Inf)==TRUE)
line.data$colour[rows.pos]<-"red"
line.data$width[rows.pos]<-lwd.scale(log10(line.data$odds[rows.pos]-threshold+1))
rows.neg<-which(c(line.data$odds<threshold & line.data$odds!=0)==TRUE)
line.data$colour[rows.neg]<-"blue"
line.data$width[rows.neg]<-lwd.scale(log10(1/(line.data$odds[rows.neg]+0.01)))

# set colours & widths for special cases
if(any(line.data$odds==Inf)){
	rows.inf<-which(line.data$odds==Inf)
	line.data$colour[rows.inf]<-"magenta"
	line.data$width[rows.inf]<-2}		
if(any(line.data$odds==0)){
	rows.inf<-which(line.data$odds==0)
	line.data$colour[rows.inf]<-"black"
	line.data$width[rows.inf]<-2}		

return(line.data)
}


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



# Plot results of pairwise.odds.ratios, using results from or.points and or.lines
plot.pairs<-function(points, lines, draw.frequencies, add.key)
{
# set defaults
if(missing(draw.frequencies))draw.frequencies<-"both"
if(missing(add.key))add.key<-"species"

# how much should the arrows be adjusted?
interpoint.dist<-mean(c((max(points[, 3])-min(points[, 3])), (max(points[, 4])-min(points[, 4]))))

# set plot region
if(any(c("objects", "species")==add.key)){
	split.screen(matrix(c(0, 0.7, 0, 1, 0.7, 1, 0, 1), nrow=2, ncol=4, byrow=T))
	screen(1)	
}else{if(add.key=="both"){
	split.screen(matrix(c(0, 0.7, 0, 1, 0.7, 1, 0.2, 1, 0.7, 1, 0, 0.2), nrow=3, ncol=4, byrow=T))
	screen(1)	
}}

par(mar=c(1,1,1,1), cex=0.7)
plot(x=points$x, y=points$y, type="n", ann=F, axes=F, asp=1)

# run loop to draw lines showing significant effects
for(i in 1:dim(lines)[1])
{
# which species are connected by this line?
sp1<-which(points$species==lines$sp1[i])
sp2<-which(points$species==lines$sp2[i])

# get line coordinates, with reduction in line length proportional to terminal point szie
if(any(c("size", "both")==draw.frequencies)){
	reduction.thisrun<-(points$freq[sp2]+interpoint.dist)*0.02
}else{reduction.thisrun<-0.02*interpoint.dist}

line.thisrun<-shorten.line(c(x0=points$x[sp1], x1=points$x[sp2],
	y0=points$y[sp1], y1=points$y[sp2]), 
	reduction=reduction.thisrun)
if(lines$arrow.code[i]==3){	# if arrows on both sides, repeat length reduction
	line.thisrun<-shorten.line(c(x0=line.thisrun[2], x1=line.thisrun[1],
		y0=line.thisrun[4], y1=line.thisrun[3]), 
		reduction=reduction.thisrun)}

# if unidirectional arrows occur in both directions, offset them. Note: Does not occur by default.
if(lines$offset[i]==1){line.thisrun<-offset.line(line.thisrun, offset=interpoint.dist*0.01)}

# draw arrows
arrows.default(line.thisrun, col=lines$colour[i], 
	lwd=lines$width[i], code=lines$arrow.code[i])
}
  
# set point sizes/colours
point.cols<-as.character(cut(points$freq, breaks=c(0, 0.036, 0.1, 0.25, 0.5, 0.75, 1),
	labels=c("blue", "green", "yellow", "brown1", "orange", "red")))
point.edge.cols<-as.character(cut(points$freq, breaks=c(0, 0.036, 0.1, 0.25, 0.5, 0.75, 1),
	labels=c("darkblue", "darkgreen", "darkgoldenrod1", "brown3", "darkorange", "darkred")))
point.cex<-2+(4*points$freq)

# overwrite these calculations if not required
if(draw.frequencies=="none"){point.cols<-"white"; point.edge.cols<-"black"; point.cex<-3}
if(draw.frequencies=="size"){point.cols<-"white"; point.edge.cols<-"black"}
if(draw.frequencies=="colour"){point.cex<-3}

# add numbered points
points(points$x, points$y, pch=21, bg=point.cols, col=point.edge.cols, cex= point.cex)
text(points$x, points$y, labels=c(1: dim(points)[1]), cex=0.7, col="black")


# ADD KEYS (according to behaviour set by add.keys)
# add species list first
if(any(c("species", "both")==add.key))
{
screen(2)
par(mar=rep(0, 4), cex=0.6)
plot(1~1, ann=F, axes=F, type="n", xlim=c(0, 1), ylim=c(0, 1.05))
text(x=c(0, 0.1), y=rep(1.05, 2), labels=c("#", "Species"), font=2, pos=4)
text(x=0, y=seq(1, 0, length.out=dim(points)[1]), labels=c(1:dim(points)[1]), pos=4)
text(x=0.1, y=seq(1, 0, length.out=dim(points)[1]), labels=points$species, pos=4)		
}

if(any(c("objects", "both")==add.key))
	{
	if(add.key=="objects"){screen(2)}else{screen(3)}

	# calculate line widths
	pos.vals<-lines[which(lines$colour=="red"), ]
	pos.labels<-round(seq(min(pos.vals$odds), max(pos.vals$odds), length.out=4), 2)
	pos.lwd<-seq(1, 3, length.out=4)
	neg.vals<-lines[which(lines$colour=="blue"), ]
	neg.labels<-round(seq(max(neg.vals$odds), min(neg.vals$odds), length.out=4), 2)
	neg.lwd<-seq(1, 3, length.out=4)

	# add plot
	par(mar=c(2, 0, 0, 0), cex=0.5)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(0, 9), ylim=c(0, 1.1))
	for(i in 1:4){lines(x=rep(c(4:1)[i], 2), y=c(0, 0.7), col="blue", lwd= neg.lwd[i], lend="butt")}
	for(i in 1:4){lines(x=rep(c(5:8)[i], 2), y=c(0, 0.7), col="red", lwd= pos.lwd[i], lend="butt")}
	axis(1, at=c(1:8), labels=c(neg.labels[4:1], pos.labels), lwd=0, line=-1, cex.axis=1)
	arrows(x0=4, x1=1, y0=0.8, y1=0.8, length=0.05)
	arrows(x0=5, x1=8, y0=0.8, y1=0.8, length=0.05)
	text(x=c(2.5, 6.5), rep(0.8, 2), labels=c("negative", "positive"), pos=3)
	mtext("odds ratio", side=2, line=-1.3, cex=0.5, adj=0.3)
}

close.screen(all.screens=TRUE)
}
# Use igraph to calculate sensible point locations
or.points<-function(dataset, threshold)
{
library(igraph)
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
if(length(rows)>0)
	{result<-result[-rows, -rows]
	species.results<-species.results[-rows, ]	}

# get plot coordinates using igraph
net<-graph.adjacency(result, mode="directed", weighted=TRUE, diag=FALSE)
layout<-layout.fruchterman.reingold(net)
points<-as.data.frame(cbind(species.results, layout))
colnames(points)[2:4]<-c("freq", "x", "y")

# set some colours
bg.list<-brewer.pal(8, "Set2")[c(3, 5, 6, 2, 4)] #c("deepskyblue", "aquamarine", "")
col.list<-c("blue", "green", "yellow", "orange", "red")
break.values<-c(0, 0.1, 0.25, 0.5, 0.75, 1)

points$bg<-as.character(cut(points$freq, 
	breaks= break.values,
	labels= bg.list))
points$col<-as.character(cut(points$freq, breaks= break.values,
	labels= col.list))
points$cex<-2+(4*points$freq)

return(points)
}


# Function for setting the attributes of lines, to be drawn between points given by or.points()
or.lines<-function(dataset, threshold, reduce)
{
library(RColorBrewer)
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
line.breaks<-c(0, 0.000001, 1/9, 1/6, 1/3, 3, 6, 9, 10^8, Inf)
line.cols<-c("black", brewer.pal(7, "RdBu")[7:1], "magenta")
line.widths<-c(2, 3, 2, 1, 0, 1, 2, 3, 2)
line.data$colour<-as.character(cut(line.data$odds,
 	breaks=line.breaks, labels=line.cols,
	include.lowest=TRUE))
line.data$width<-line.widths[as.numeric(as.character(cut(test.lines$odds,
 	breaks=line.breaks, labels=c(1:9), include.lowest=TRUE)))]

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
if(missing(draw.frequencies))draw.frequencies<-TRUE
if(missing(add.key))add.key<-"species"

# set error messages
if(dim(points)[1]==0)stop("Error: no strong inter-species associations identified")

switch(add.key,
	"none"={draw.sppairs(points, lines, draw.frequencies)},
	"species"={
		split.screen(matrix(c(0, 0.7, 0, 1, 0.7, 1, 0, 1), nrow=2, ncol=4, byrow=T))	
		screen(1); draw.sppairs(points, lines, draw.frequencies)
		screen(2); draw.spp.key(points, draw.frequencies)
		close.screen(all.screens=TRUE)},
	"object"={
		split.screen(matrix(c(0, 0.8, 0, 1, 0.8, 1, 0, 1), nrow=2, ncol=4, byrow=T))
		screen(1); draw.sppairs(points, lines, draw.frequencies)
		screen(2); draw.object.key(draw.frequencies)
		close.screen(all.screens=TRUE)},
	"both"={
		split.screen(matrix(c(0, 0.6, 0, 1, 0.6, 0.8, 0, 1, 0.8, 1, 0, 1), nrow=3, ncol=4, byrow=T))
		screen(1); draw.sppairs(points, lines, draw.frequencies)
		screen(2); draw.spp.key(points, draw.frequencies)
		screen(3); draw.object.key(draw.frequencies)
		close.screen(all.screens=TRUE)}
	)
}	# end plot.pairs


draw.object.key<-function(draw.frequencies)
	{
	if(draw.frequencies){
		screens<-split.screen(matrix(c(0, 1, 0, 0.5, 0, 1, 0.5, 1), nrow=2, ncol=4, byrow=T))
		screen(screens[1]); draw.line.key()
		screen(screens[2]); draw.point.key()
	}else{draw.line.key()}
	}



draw.sppairs<-function(points, lines, draw.frequencies)
{
# how much should we adjust point sizes etc.?
interpoint.dist<-mean(c((max(points[, 3])-min(points[, 3])), (max(points[, 4])-min(points[, 4]))))

# set plot region
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
  
# add numbered points
if(draw.frequencies==TRUE){
	points(points$x, points$y, pch=21, bg=points$bg, col=points$col, cex=points$cex)
}else{
	points(points$x, points$y, pch=21, bg="white", col="black", cex=3)}
text(points$x, points$y, labels=c(1: dim(points)[1]), cex=0.7, col="black")

}	# end draw.sppairs()



draw.spp.key<-function(points, draw.frequencies){
	par(mar=rep(0, 4), cex=0.6)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(0, 1), ylim=c(0, 1.05))
	text(x=0.1, y=1.05, labels="Species", font=2, pos=4)
	if(draw.frequencies==TRUE){
		points(rep(0, dim(points)[1]), seq(1, 0, length.out=dim(points)[1]), 
		pch=21, bg=points$bg, col=points$col, cex=points$cex)
	}else{
		points(rep(0, dim(points)[1]), seq(1, 0, length.out=dim(points)[1]), 
		pch=21, bg="white", col="black", cex=3)}
	text(x=0, y=seq(1, 0, length.out=dim(points)[1]), labels=c(1:dim(points)[1]))
	text(x=0.1, y=seq(1, 0, length.out=dim(points)[1]), labels=points$species, pos=4)	
	}


draw.point.key<-function(){
	bg.list<-brewer.pal(8, "Set2")[c(3, 5, 6, 2, 4)] 
	col.list<-c("blue", "green", "yellow", "orange", "red")
	break.values<-c(0, 0.1, 0.25, 0.5, 0.75, 1)
	point.labels<-paste(break.values[1:5], break.values[2:6], sep=" - ")
	point.cex<-2+(4*break.values[2:6])

	par(mar=c(0, 0, 1, 0), cex=0.5)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(0.2, 1), ylim=c(0.5, 5.5))
	points(x=rep(0.9, 5), y=c(1:5), pch=21,
		bg=bg.list, col=col.list, cex=point.cex)
	text(x=rep(0.8, 5), y=c(1:5), pos=2, cex=1,
		labels= point.labels)
	mtext("spp. occurrence", side=3, cex=0.5, adj=0.5)
	}


draw.line.key<-function(){
	library(RColorBrewer)
	line.breaks<-c(0, 0.000001, 1/9, 1/6, 1/3, 3, 6, 9, 10^8, Inf)
	line.labels<-c("0", "<1/9", "<1/6", "<1/3", ">3", ">6", ">9", "Inf")
	line.cols<-c("black", brewer.pal(7, "RdBu")[c(7:5, 3:1)], "magenta")
	line.widths<-c(2, 3, 2, 1, 1, 2, 3, 2)
	# add plot
	par(mar=c(0, 3, 0, 0), cex=0.5)
	plot(1~1, ann=F, axes=F, type="n", xlim=c(0, 1), ylim=c(0, 9))
	for(i in 1:8){
		lines(x=c(0.1, 1), y=rep(c(1:8)[i], 2), 
			col= line.cols[i], lwd= line.widths[i], lend="butt")}
	axis(2, at=c(1:8), labels=line.labels, lwd=0, line=-1, cex.axis=1, las=1)
	mtext("odds ratio", side=3, line=-1.3, cex=0.5, adj=0.5)
	}


# Internal function to convert the 'long' output from spaa() into a matrix
pairwise.or.matrix<-function(or.object)
{
# work out species info
n.species<-max(or.object$col1)
species.names<-rep("none", n.species)
for(i in 1:n.species){species.names[i]<-as.character(or.object$sp1[which(or.object$col1==i)[1]])}
  
# 	create & fill output
output.matrix<-matrix(data=NA, n.species, n.species)
colnames(output.matrix)<-species.names
rownames(output.matrix)<-species.names
for(i in 1:dim(or.object)[1]){output.matrix[or.object$col1[i], or.object$col2[i]]<-or.object$odds[i]}
  
return(output.matrix)
}


# Use igraph to calculate sensible point locations
or.points<-function(dataset, threshold)
{
library(igraph)
if(missing(threshold))threshold<-3

# select input data
species.results<-dataset$species	# what are the species we are looking at?
or.matrix<-pairwise.or.matrix(dataset$combinations)	# extract wide-format odds ratios from input data

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
line.data$width<-line.widths[as.numeric(as.character(cut(line.data$odds,
 	breaks=line.breaks, labels=c(1:9), include.lowest=TRUE)))]

return(line.data)
}
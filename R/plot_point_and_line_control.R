# AVAILABLE TO USERS:

# Use igraph to calculate sensible point locations for one or more spaa objects
spaa.points<-function(dataset, 
	threshold, 
	graph.function, 	# what function should be called to generate point locations?
	return.matrix,		# logical; if TRUE, returns an adjacency matrix; if FALSE (the default) returns a data.frame of point coordinates
	... # optional information passed to igraph::layout.fruchterman.reingold
	)
{
if(missing(threshold))threshold<-3
if(missing(graph.function))graph.function<-"layout.fruchterman.reingold"
if(missing(return.matrix))return.matrix<-FALSE

# run code for a single spaa result
if(class(dataset)=="spaa"){
	species.results<-dataset$species	# what are the species we are looking at?
	result<-adjacency.test(dataset, threshold)
}else{if(class(dataset)=="list"){
	
	n.inputs<-length(dataset)

	# get species list & mean frequency
	species.data<-list()
	for(i in 1: n.inputs){species.data[[i]]<-dataset[[i]]$species}
	species.list<-list()
	for(i in 1: n.inputs){species.list[[i]]<-species.data[[i]]$species}	
	species.list<-sort(unique(unlist(species.list)))
	n.species<-length(species.list)
	species.results<-data.frame(species=species.list, frequency=0, stringsAsFactors=FALSE)
	for(i in 1:n.inputs){
		colnames(species.data[[i]])[2]<-paste("x", i, sep="")
		species.results<-merge(species.results, species.data[[i]], by="species", all=TRUE)}
	species.results$frequency<-apply(species.results[, c(3:dim(species.results)[2])], 1, 
		FUN=function(x){mean(x, na.rm=TRUE)})
	species.results<-species.results[, 1:2]

	# calculate adjacency
	distance.array<-array(data=NA, dim=c(n.species, n.species, n.inputs), 
		dimnames=list(species.list, species.list, c(1:n.inputs)))
	for(i in 1:n.inputs){distance.array[, , i]<-adjacency.test(dataset[[i]], threshold, species.list)}
	result<-matrix(data=NA, nrow=n.species, ncol=n.species)
		colnames(result)<-species.list
	for(i in 1:n.species){result[, i]<-apply(distance.array[, i, ], 1, FUN=function(x){sum(x, na.rm=TRUE)})}
}}	# end list code

# remove uninformative spp
result.freq<-apply(
	cbind(colSums(result, na.rm=TRUE), rowSums(result, na.rm=TRUE)), 1, sum)
rows<-which(result.freq==0)
if(length(rows)>0)
	{result<-result[-rows, -rows]
	species.results<-species.results[-rows, ]	}

if(return.matrix){return(result)
}else{
	# get plot coordinates 
	if(graph.function=="graph.grid"){
		layout<-graph.grid(result)
	}else{
		net<-graph.adjacency(result, mode="max", weighted=NULL, diag=FALSE)
		layout<-do.call(graph.function, list(net, ...)) #layout.fruchterman.reingold(net) # layout.kamada.kawai(net)
	}
	points<-as.data.frame(cbind(species.results, layout))
	colnames(points)[2:4]<-c("freq", "x", "y")
	
	return(points)
	}

}



# Function for setting the attributes of lines, to be drawn between points given by or.points()
spaa.lines<-function(dataset, threshold)
{
if(missing(threshold))threshold<-3

# make some changes to result.long for plotting lines
line.data<-dataset$combinations[, 3:5]	# remove spp numbers
# line.data$sp1<-as.character(line.data$sp1)
# line.data$sp2<-as.character(line.data$sp2)
# remove values outside of threshold range
line.data<-rbind(
	line.data[which(line.data$odds<=(1/threshold)), ], # & line.data$odds>=0, ],
	line.data[which(line.data$odds>=threshold), ]) # & line.data$odds<=Inf, ])

# add new columns to line data to allow prettier plots
line.data$difference<-0; line.data$arrow.code<-1

# run loop to get properties of lines representing strong effects
for(i in 1:dim(line.data)[1])
{
repeat.test<-c(line.data$sp2==line.data$sp1[i] & line.data$sp1==line.data$sp2[i])
if(any(repeat.test==TRUE)){
	other.row<-which(repeat.test==TRUE)
	if(other.row>i){		# selects only later rows
			# is line uni- or bi-directional?
			values<-line.data$odds[c(i, other.row)]
			line.data$arrow.code[i]<-2
			line.data$arrow.code[other.row]<-NA
			line.data$odds[i]<-mean(values)
			line.data$difference[i]<-max(values)-min(values)   
	}}}

if(any(is.na(line.data$arrow.code))){
	line.data<-line.data[-which(is.na(line.data$arrow.code)), ]}
return(line.data)
}



## INTERNAL FUNCTIONS ##

# create a function to run inside spaa.points that calculates point adjacency
adjacency.test<-function(dataset,	# must be a single spaa() result
	threshold,
	species.list  # if given, a list of species to determine output matrix dimensions an colnames
	)
	{
	# select input data
	or.matrix<-pairwise.or.matrix(dataset$combinations)	# extract wide-format odds ratios from input data
	if(missing(species.list)==FALSE){
		corrected.matrix<-matrix(data=NA, nrow=length(species.list), ncol=length(species.list))
		rownames(corrected.matrix)<-species.list
		colnames(corrected.matrix)<-species.list
		for(i in 1:length(species.list)){
			if(any(colnames(or.matrix)==species.list[i])){
				col<-which(colnames(or.matrix)==species.list[i])
				data.thisrun<-data.frame(species=rownames(or.matrix), value=or.matrix[, col], stringsAsFactors=FALSE)
				corrected.matrix[, i]<-as.numeric(
					merge(data.frame(species=species.list, stringsAsFactors=FALSE), 
					data.thisrun, by="species", all=TRUE)$value)
			}}
		or.matrix<-corrected.matrix
	}	# end else

	# ID interacting spp.
	result.positive<-make.binary(or.matrix, threshold=threshold)	# positively-interacting species
	# and for contra-indicators
	lower.threshold<-(1/threshold)
	or.matrix.neg<-(lower.threshold-or.matrix)
	result.negative<-make.binary(or.matrix.neg, threshold=0)
	
	# combine -ve and +ve results into one matrix; ID most strongly interacting spp.
	result<-as.matrix(result.positive)+ as.matrix(result.negative)
	return(result)
	}



# create a wrapper function to detect outlying points in spaa.points and return a more contrained result.
# this is a messy function, and should be used with caution; hence not able to be called yet.
graph.kmeans.reduced<-function(input)
	{
	values<-layout.fruchterman.reingold(input)
	distances<-dist(values, upper=TRUE)
	mean.distances<-kmeans(distances, 2) # apply(distances, 1, function(x){mean(x, na.rm=TRUE)})

	# identify rows that belong to well- and poorly- represented clusters
	frequencies<-as.matrix(xtabs(rep(1, length(mean.distances$cluster))~mean.distances$cluster))
	rarest.value<-as.numeric(rownames(frequencies)[which.min(frequencies)])
	distant.points<-as.numeric(names(which(mean.distances$cluster==rarest.value)))
	close.points<-c(1:dim(values)[1])[-distant.points]

	# find central point of largest cluster
	centre.point<-apply(values[close.points, ], 2, mean)
	# find mean interpoint distance within this cluster
	cluster.distances<-as.matrix(dist(values[close.points, ], upper=TRUE, diag=FALSE))
		for(i in 1:dim(cluster.distances)[1]){cluster.distances[i, i]<-NA}
	mean.dist<-mean(apply(cluster.distances, 1, FUN=function(x){min(x, na.rm=TRUE)}))
	# find max distance from the centre
	max.dist<-quantile(sqrt((values[close.points, 1]-centre.point[1])^2 + (values[close.points, 2]-centre.point[2])^2), 0.9)
	radius<-max.dist + mean.dist	
	# find centre of smaller cluster
	outer.point<-apply(values[distant.points, ], 2, mean)

	# use this information to place the outlying points equidistantly along a circle surrounding the central cluster
	for(i in 1:length(distant.points)){
		centres<-as.data.frame(rbind(centre.point, values[distant.points[i], ]))
			colnames(centres)<-c("x", "y")
		model<-lm(y~x, data=centres)
		theta<-as.numeric(atan(coef(model)[2]))
		y.diff<-radius*cos(theta)
		x.diff<-radius*sin(theta)
		values[distant.points[i], ]<-c(
			centre.point[1]-x.diff, centre.point[2]-y.diff)	
		}
	return(values)	
	}



# function to put points on a grid
graph.grid<-function(
	input	,	# association matrix
	x.expansion		# as plot(asp=1), we may need to change the coordinates of the actual points
	){
	# set defaults
	if(missing(x.expansion))x.expansion<-1.2

	# work out the dimensions of the grid
	input.size<-dim(input)[1]

	# try an approach on incrementally increasing sized parabolas
	# first determine how many you need
	size.ID<-data.frame(n.parabolas=c(1:28), n.entries=cumsum(c(3:30))*2)
	row<-min(which(size.ID$n.entries>=input.size))
	n<-size.ID$n.parabolas[row]

	# then work out coordinates for that many points
	heights<-seq(0, 1, length.out=n+1)[c(2: (n+1))]
	widths<-(heights*0.5) + 0.5
	result.list<-list()
	for(i in 1:n){
		data.thisrun<-data.frame(x=c(-widths[i], 0, widths[i]), y=c(0, heights[i], 0))
		model<-lm(y~x+ I(x**2), data=data.thisrun)
		xvals<-seq(-widths[i], widths[i], length.out=(3+i))
		x.interval<-xvals[2]-xvals[1]
		xvals<-xvals[c(1:(2+i))] + (x.interval*0.5)
		yvals<-as.numeric(predict(model, newdata=data.frame(x=xvals)))
		result<-data.frame(x=xvals, y=yvals)
		result.list[[i]]<-result
		}
	point.dframe<-do.call(rbind, result.list)
	point.dframe<-rbind(point.dframe, 
		data.frame(x=point.dframe$x, y=-1*point.dframe$y))

	# rotate all points by a fixed amount
	rotation<-0.123
	point.dframe$new.x<-(point.dframe$x * cos(rotation))-(point.dframe$y * sin(rotation))
	point.dframe$new.y<-(point.dframe$y * cos(rotation))+(point.dframe$x * sin(rotation))
	point.dframe<-point.dframe[, 3:4]
	colnames(point.dframe)<-c("x", "y")

	# order points by increasing distance from origin
	point.dframe$dist<-sqrt(point.dframe$x^2 + point.dframe$y^2)
	point.dframe<-point.dframe[order(point.dframe$dist, decreasing=FALSE), ]

	# order species by decreasing number of connections
	spp.dist<-cbind(
		apply(input, 1, FUN=function(x){sum(x, na.rm=TRUE)}),
		apply(input, 2, FUN=function(x){sum(x, na.rm=TRUE)}))
	spp.freq<-data.frame(
		species=rownames(spp.dist),
		n=apply(spp.dist, 1, sum),
		stringsAsFactors=FALSE)
	rownames(spp.freq)<-NULL
	spp.freq<-spp.freq[order(spp.freq$n, decreasing=TRUE), ]

	# merge point and species data
	point.dframe<-point.dframe[c(1:dim(spp.freq)[1]), ]	# ensure same size
	result<-cbind(point.dframe, spp.freq)

	# reorder such that spp are in initial order
	result<-result[order(result$species), 1:2]
	# result<-jitter(result)	# line drawing doesn't work if points are directly above/below one another 
		# as coef(model))[2]==NA
	return(result)
	}

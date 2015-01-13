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



# create a function to run inside spaa.points that calculates point adjacency
adjacency.test<-function(dataset,	# must be a single spaa() result
	threshold,
	species.list  # if given, a list of species to determine output matrix dimensions an colnames
	)
	{
	# select input data
	or.matrix<-pairwise.or.matrix(dataset$combinations)	# extract wide-format odds ratios from input data
	
	# add extra rows/columns if necessary
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
	}
					
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



# Use igraph to calculate sensible point locations for one or more spaa objects
spaa.points<-function(dataset, threshold)
{
if(missing(threshold))threshold<-3

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

# get plot coordinates using igraph
net<-graph.adjacency(result, mode="directed", weighted=TRUE, diag=FALSE)
layout<-layout.kamada.kawai(net)# layout.fruchterman.reingold(net)
points<-as.data.frame(cbind(species.results, layout))
colnames(points)[2:4]<-c("freq", "x", "y")

return(points)
}


# Function for setting the attributes of lines, to be drawn between points given by or.points()
spaa.lines<-function(dataset, threshold, reduce)
{
if(missing(threshold))threshold<-3
if(missing(reduce))reduce<-TRUE

# make some changes to result.long for plotting lines
line.data<-dataset$combinations[, 3:5]
line.data$sp1<-as.character(line.data$sp1)
line.data$sp2<-as.character(line.data$sp2)
# remove values outside of threshold range
line.data<-rbind(
	line.data[line.data$odds<=(1/threshold) & line.data$odds>=0, ],
	line.data[line.data$odds>=threshold & line.data$odds<=Inf, ])

# add new columns to line data to allow prettier plots
line.data<-cbind(line.data, as.data.frame(matrix(data=0, nrow=dim(line.data)[1], ncol=2)))
colnames(line.data)[4:5]<-c("difference", "arrow.code") #, "offset", "colour", "width")
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

return(line.data)
}
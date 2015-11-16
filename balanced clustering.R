library(plyr)
set.seed(4)
setwd("C:/Users/ading001c/Documents/Geocode/CAR Module Realignment/")
#data <- read.csv(file.choose(),header=TRUE)
Geocode = read.csv("Geocode - Copy.csv", header=T)

n.all=nrow(Geocode)
n.cluster=length(unique(Geocode$Cluster))

###Input###
## Sampling ##
sample.size=30
sample.tmp=lapply(seq(n.cluster),function(k) Geocode[sample(Geocode[Geocode$Cluster==k,][,1], sample.size),])
sample=arrange(do.call(rbind,sample.tmp),ID)
non_sample=Geocode[!(Geocode$ID  %in% sample[,1]),]
non_sample$Cluster<-NA

##poly-stable##
n.row=nrow(non_sample)

# cluster mean#
centers <- aggregate(cbind(sample$Latitude,sample$Longitude)~sample$Cluster, FUN=mean)

# required minimum cluster size m
m=100

### Method ###

#mark properties in non-sample as free if no cluster assigned
free=function(prop.ID){ 
  free=is.na(non_sample[non_sample$ID==prop.ID,4])
    return (free)
}


## calculate distance of all the points from each of the clusters
    # Convert to radian function
    radians = function(theta=0){return(theta * pi / 180)}
    G <- radians(non_sample[,c(2,3)])
    C <- radians(centers[,c(2,3)])    
    R <- 3958.756 # Earth mean radius [miles]
    # Define the Grear Circle Distance function
    dist <- function(xind, yind)
      acos(sin(G[xind, 1]) * sin(C[yind, 1]) + 
             cos(G[xind, 1]) * cos(C[yind, 1]) * cos(C[yind, 2] - G[xind, 2])) * R

    dist_mat <- outer(seq.int(n.row), seq.int(n.cluster), dist)
    dist_mat1 <- cbind(dist_mat,non_sample$ID)
    colnames(dist_mat1) <-c(seq(30),"ID")
    dist_mat1=as.data.frame(dist_mat1)

#sort all properties based on distance
    cluster.pref.ls<-as.data.frame(lapply(seq(n.cluster), 
                            function(h) arrange(as.data.frame(non_sample$ID),dist_mat[,h])))
    prop.pref.ls<-as.data.frame(lapply(seq(n.row), function(k) 
                            arrange(as.data.frame(seq(30)),dist_mat[k,])))
#calculate number of points to be assigned for each cluster
    open.slots=m-table(sample$Cluster)
#index in the sorted list
    index<- rep(1,30)

points.remain=sum(open.slots)

iteration=0
while (points.remain>0){
    for (h in 1:n.cluster){ 
      open.slots.old=open.slots
        if(open.slots[[h]]>0){
            for (i in index[h]:(index[h]+open.slots[[h]]-1)){#clusters propose
                proposed<- cluster.pref.ls[i,h]
                #points dispose
                if (free(proposed)==TRUE){
                  #add item i to cluster h
                  non_sample[non_sample$ID==proposed,4]<-h
                  open.slots[[h]]=open.slots[[h]]-1}#end if
                else if((free(proposed)==FALSE)&&
                          (dist_mat1[dist_mat1$ID==proposed,h] < 
                             dist_mat1[dist_mat1$ID==proposed,non_sample[non_sample$ID==proposed,4]])){
                  #item i is engaged to cluster h*, but i is closer to cluster h
                  #remove item i from cluster h*
                  old.cluster=non_sample[non_sample$ID==proposed,4]
                  open.slots[[old.cluster]]=open.slots[[old.cluster]]+1
                  #move item i to cluster h
                  non_sample[non_sample$ID==proposed,4]<-h
                  open.slots[[h]]=open.slots[[h]]-1}#end else             
            }#end for
            index[h] <- index[h]+open.slots.old[[h]] #Update index
        }#end if
    }#end for
    points.remain <- sum(open.slots)
    iteration=iteration+1
}#end while


#After poly-stable terminates, the remaining points are greedily 
  #assigned to their nearest clusters
remain.points<-subset(non_sample,is.na(non_sample$Cluster)) #subset ponits with NA values
remain.points$Cluster<-apply(dist_mat[dist_mat1$ID  %in% remain.points[,1],],MARGIN=1,which.min)


ret.val=arrange(rbind(sample,subset(non_sample,!is.na(non_sample$Cluster)),remain.points),ID)

#Update cluster centers
centers <- aggregate(cbind(ret.val$Latitude,ret.val$Longitude)~ret.val$Cluster, FUN=mean)

#Plot data

  plot(ret.val$Longitude, ret.val$Latitude,col=ret.val$Cluster, pch=16, cex=.6,
       xlab="Longitude",ylab="Latitude")
  library(maps)
  map("state", add=T)
  map("county", add=T)
  text(centers[,c(3,2)],label=c(1:28,31,32),cex=.8,col="red",font=2)

###refine###

#max weighted edge function
maxweight.edge = function(V1,V2){
  tmp.edge <- edge.list[ (edge.list$from==V1) & (edge.list$to==V2), ]
  max <- max(unlist(tmp.edge$weight))
  maxweight.edge <- tmp.edge [ tmp.edge$weight == max, ]
  return (maxweight.edge)
}

#search cycle function
 find.cycle=function(x=1){
   order<-c()
   visited<-c()
   arr <- tail(adj.mat[x,],ncol(adj.mat)-x)
   while(sum(arr)>0){
     #if edge exits
     col<- which(arr > 0)[[1]]+x                   #get col index
     order[1]<- rownames(adj.mat)[1]               #Initiate node order 
     visited[1]<-1
     visited<-c(visited,col)                       #mark as visited                      
     order<- c(order, rownames(adj.mat)[col])      #Extend node order
     x=col                                         #next row to visit
     arr <- tail(adj.mat[x,],ncol(adj.mat)-x)
   }

   cycle.formed= adj.mat[tail(order,1),head(order,1)] > 0
   
   while(cycle.formed == FALSE){
     order<-head(order,-1)
     cycle.formed = adj.mat[tail(order,1),head(order,1)] > 0
   }

   ret.val = list("order" = order, "visited" = visited) 
   return(ret.val)
 }


  max.iter=100
  theta=0
  Iteration = 0
  while(theta>0){
    options(error=recover)
    theta=0
    #Update cluster centers
    centers <- aggregate(cbind(ret.val$Latitude,ret.val$Longitude)~ret.val$Cluster, FUN=mean)
    
    #Update dist_mat_rf
    G <- radians(ret.val[,c(2,3)])
    C <- radians(centers[,c(2,3)])    
    dist_mat_rf <- outer(seq.int(n.all), seq.int(n.cluster), dist)
    size=table(ret.val$Cluster)
    
    ##Individual Reassginments
    opt.cluster <- apply(dist_mat_rf,MARGIN=1,which.min)               #Calculate optimal cluster
    swap.candidates<-list()
    
    for (i in 1:n.all){
      if (opt.cluster[[i]] != ret.val$Cluster[[i]]){                    #if x has a better option
        if(size[ret.val$Cluster[[i]]] > m){                               #cluster can afford to give
          ret.val$Cluster[[i]] <- opt.cluster[[i]]                        #move x to new cluster
          theta=theta+1                                                 #termination threshold
          size=table(ret.val$Cluster)
        }
        else{
          swap.candidates[[i]] <- ret.val[i,]                           #add x to hold_mat
        }
      }
    }
    swap.candidates <- do.call(rbind, swap.candidates)
    
    ## Group Reassignments
    # Construct potential assignment graph G=(V,E) from {H(x)} for all x in X, such that
    # 1.Every cluster becomes a node.
    # 2.An edge (A,B) is added for every point in A that is closer to 
    #  the center of B (so that there are potentially multiple edges between the same pair of nodes)
    #  its weight should be proportional to the benefit gained by moving it.
    
    #Create edge list("ID","from","to","benefits")
    sub.dist.mat <- as.data.frame(dist_mat_rf[cbind(dist_mat_rf,ret.val$ID)[,31] %in% swap.candidates$ID,])    #subset distance matrix
    
    edge.list<-list()
    for(i in 1:nrow(sub.dist.mat)){
      Cl_col<-swap.candidates[i,4]
      tmp<-sub.dist.mat[i,Cl_col]
      for(k in 1:n.cluster){
        benefits<-tmp-sub.dist.mat[i,k]
        if (benefits>0){
          edge.list[[length(edge.list)+1]]<-c(swap.candidates[i,c(1,4)],k,benefits)
        }   
      }  
    }
    
    edge.list<- do.call(rbind, edge.list)
    colnames(edge.list) <-c("ID","from","to","weight")
    edge.list <- as.data.frame(edge.list)
    
    #Create igraph object from edgelist
    library(igraph)
    V <- unique(c(edge.list$from, edge.list$to))
    graph<- graph.data.frame(edge.list[,c(2:4)],directed=TRUE)
    E(graph)$weight=as.numeric(edge.list$weight)
    E(graph)$ID=edge.list$ID
    
    #Get cycles
    SCC<-clusters(graph,mode="strong")
    cycles<-list()
    for (i in 1:SCC$no){
      if (length(which(SCC$membership==i))>1){
        cycles[[length(cycles)+1]]<- V(graph)$name[SCC$membership==i]
      }
    }
       
    #find cycle by visiting the adj.mat -- not optimal
    #investigate in Hamiltonian Cycle and Eulerian Cycle for optimization
    
    cycles <- lapply(1:SCC$no, function(x) {
      sg <- induced.subgraph(graph, SCC$membership == x)
      n <- sum(SCC$membership == x)
      col <- rep(c(1, 0), c(1, n-1))  # Used to grab isomorphism starting at 1
      sg.idx <- graph.get.subisomorphisms.vf2(sg, graph.ring(n, directed=TRUE), col, col)[[1]]
      which(SCC$membership == x)[sg.idx]
    })
    
    
    node.order<-list()
    for (i in 1:length(cycles)){
      vids<-cycles[[i]]
      subgraph <- induced.subgraph(graph, vids)
      adj.mat<- get.adjacency(subgraph, edges=F, names=TRUE,sparse=F)
      fc <- find.cycle()
      #return result
      node.order[[i]] <- fc$order
    }
    node.order[[1]] <- as.character(c(10,9,1,12,4,3,13))
    node.order[[1]] <- as.character(c(20,34))
    #sub edgelist
    subEL<-list()
    for (i in 1:length(node.order)){
      data<- node.order[[i]]
      subEL[[i]]<-as.data.frame(cbind(data,append(data[-1], head(data,1), after=length(data[-1]))))
      colnames(subEL[[i]])=c("from","to")
    }

    #points to move
    points.holder<-list()
    for( i in 1:length(subEL)){
      n=nrow(subEL[[i]])
      for (k in 1:n){
        points.holder[[length(points.holder)+1]]<-maxweight.edge(subEL[[i]][k,1], subEL[[i]][k,2])
      }
    }    
    points.to.move<-data.frame(matrix(unlist(points.holder), ncol=4, byrow=T))
    colnames(points.to.move)=c("ID","from","to","benefit")
    points.to.move<-arrange(points.to.move,ID)

    #move points
    ret.val.cp=ret.val
    for (i in 1:nrow(points.to.move)){
      ret.val$Cluster[ret.val$ID == points.to.move$ID[i]]<-points.to.move$to[i] 
    }
    
    
    plot(ret.val$Longitude, ret.val$Latitude,col=ret.val$Cluster,main="Result After", pch=16, cex=.6,
         xlab="Longitude",ylab="Latitude")
    library(maps)
    map("state", add=T)
    map("county", add=T)
    text(centers[,c(3,2)],label=c(1:28,31,32),cex=.8,col="red",font=2)
    
    #configuration
    theta=theta+nrow(points.to.move)
    Iteration=Iteration+1
    cat ("Iteration:", Iteration ,"points moved =", theta, "gain =" , sum (points.to.move$benefit))
    if(Iteration>max.iter) break

  }#end outer while loop


library(cluster)
clusplot(cbind(ret.val$Longitude, ret.val$Latitude), ret.val$Cluster, color=T,shade=TRUE, 
         labels=0, lines=0,col.p=ret.val$Cluster,
         xlab="Longitude",ylab="Latitude",cex=1)

plot(ret.val.pl$Longitude, ret.val.pl$Latitude,col=ret.val.pl$Cluster, pch=16, cex=.6,
     xlab="Longitude",ylab="Latitude")
library(maps)
map("state", add=T)
map("county", add=T)
text(centers[,c(3,2)],label=c(1:28,31,32),cex=.8,col="red",font=2)

text(ret.val.pl[,c(3,2)],label=ret.val.pl$ID,cex=.1,col="white",font=2)

#plotGoogleMap
coordinates(ret.val)<-~Longitude+Latitude # convert to SPDF
proj4string(ret.val) <- CRS("+init=epsg:4326") 
ic <- iconlabels(ret.val$Cluster,height=10,icon=TRUE,colPalette=rainbow(30))
m <- plotGoogleMaps(ret.val,zcol="Cluster",filename='myMap.htm',iconMarker=ic, 
                    colPalette=rainbow(30),mapTypeId='ROADMAP',layerName='Clusters')

#Plor graph
V(graph)$color <- rainbow(SCC$no)[SCC$membership]
plot.igraph(graph,vertex.label=V(graph)$name,layout=layout.kamada.kawai, 
            vertex.label.color="black",edge.color="black",mark.groups=split(1:vcount(graph), SCC$membership),
            edge.arrow.size=0.15,edge.curved=TRUE)

write.table(ret.val, paste("ret.val.csv", sep=""), sep=",", row.names=F)
write.table(centers, paste("centers.csv", sep=""), sep=",", row.names=F)

plot(non_sample$Longitude, non_sample$Latitude,col=non_sample$Cluster, pch=16, cex=.6,
     xlab="Longitude",ylab="Latitude")

test <- subset(non_sample, 
                !(is.na(Cluster) )
)
clusplot(cbind(test$Longitude, test$Latitude), test$Cluster, color=T,shade=TRUE, 
         labels=0, lines=0,col.p=test$Cluster,
         xlab="Longitude",ylab="Latitude",cex=1)

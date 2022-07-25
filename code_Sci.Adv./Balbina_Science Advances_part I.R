# Palmeirim et al. 2022, Science Advances

# Summary of the script:

# 1 - Network draws (Figure 3)
# 2 - Landscape-level properties (Figure 4)
    # 2.1 Connectance
    # 2.2 Modularity
    # 2.3 Nestedness
    # 2.4 Robustness
        # 2.2.1 Largest-to-smallest
        # 2.2.2 Smallest-to-largest
        # 2.2.3 Random
    # 2.5 Figure 4
        # 2.5.1 Simulated values        
        # 2.5.2 Observed values
    # 2.6 Figure 5
# 3 - Site- and species-level properties (node properties, Figure 6)
    # 3.1 Site and species normalised degree
        # 3.1.1 Site normalised degree
        # 3.1.2 Species normalised degree
    # 3.2 Site and species nestedness contribution

# Figure 7

# Figures from Supplementary Materials: S1 to S4



# ************
  
# 1 - network draws (Figure 3)

setwd("~/Desktop/code_Sci.Adv.")

# loading the packages needed
library(bipartite); library (networkD3); library(reshape); library(reshape2)
library(igraph); library(magrittr);  library(plyr); library(ggplot2); library(vegan)
library(MuMIn); options(na.action = "na.fail"); options(error=recover)
library(fitdistrplus); library(nlme)

# 1.1 All taxa

# loading the data
df<-read.csv("all_taxa.csv", header=T, row.names = 1)
traits<-read.csv("all_taxa_trait_rede.csv", header=T)

#Transform the interaction matrix into an igraph object
df1<-  df[,colSums(df) > 0] # make sure all columns have interactions with lines
str(df1) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1),df1)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1]))
netlist[,2]=as.character(paste(netlist[,2]))
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
#plot(gad)
#Check the main properties of the igraph object
V(g)
E(g)
#ordering data frames so that the order of traits matches the order of species in graph g
guild <- traits[order(match(traits$area, V(g)$name)),] 
str(guild)
# check whether igraph object and traits df have same length
setdiff(guild$area, V(g)$name)
setdiff(V(g)$name, guild$area)
#removing species not present in the habitat matrix
del<-(setdiff(traits$area, V(g)$name)) 
guild <- traits[!(traits[,1] %in% c(del)),]
# check again, both objects must have same lenght
setdiff(guild$area, V(g)$name)
setdiff(V(g)$name, guild$area)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="sm"]="circle"
V(g)$shape[V(g)$guild=="liz"]="circle"
V(g)$shape[V(g)$guild=="mid.large"]="circle"
V(g)$shape[V(g)$guild=="birds"]="circle"
V(g)$shape[V(g)$guild=="frogs"]="circle"
V(g)$shape[V(g)$guild=="beetles"]="circle"
V(g)$shape[V(g)$guild=="bees"]="circle"
V(g)$shape[V(g)$guild=="trees"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="sm"]="tomato"
V(g)$color[V(g)$guild=="liz"]="light green"
V(g)$color[V(g)$guild=="mid.large"]="dark red"
V(g)$color[V(g)$guild=="birds"]="orange"
V(g)$color[V(g)$guild=="frogs"]="light blue"
V(g)$color[V(g)$guild=="beetles"]="blue"
V(g)$color[V(g)$guild=="bees"]="gold"
V(g)$color[V(g)$guild=="trees"]="yellowgreen"
plot(g,
     # Set the drawing mode
     layout=layout.fruchterman.reingold,
     # Set title
     main='Balbina habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape,
     vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=2000,
     vertex.label.cex=0.5,
     #   vertex.size=log(igraph::degree(g)),
     vertex.size=log(igraph::degree(g)),)

# 1.2 

# set directory
# data on 8 taxonomic groups across 22 islands and 3 CF sites
# sps--site matrices (binary data)
df.sm<-read.csv("small_mammals.csv", header = T, row.names = 1) 
df.liz<-read.csv("lizards.csv", header = T, row.names = 1) 
df.mid.large<-read.csv("mid.large.csv", header = T, row.names = 1) 
df.birds<-read.csv("birds.csv", header = T, row.names = 1) 
df.frogs<-read.csv("frogs.csv", header = T, row.names = 1) 
df.beetles<-read.csv("beetles2.csv", header = T, row.names = 1) 
df.bees<-read.csv("bees.csv", header = T, row.names = 1) 
df.trees<-read.csv("trees.csv", header = T, row.names = 1) 

# traits matrices
traits.sm<-read.csv("small_mammals_trait_rede.csv", header = T) 
traits.liz<-read.csv("lizards_trait_rede.csv", header = T) 
traits.mid.large<-read.csv("mid.large_trait_rede.csv", header = T) 
traits.birds<-read.csv("birds_trait_rede.csv", header = T) 
traits.frogs<-read.csv("frogs_trait_rede.csv", header = T) 
traits.beetles<-read.csv("beetles_trait_rede.csv", header = T) 
traits.bees<-read.csv("bees_trait_rede.csv", header = T) 
traits.trees<-read.csv("trees_trait_rede.csv", header = T) 

# Small mammals

df1.sm<-df.sm[,colSums(df.sm) > 0]
str(df1.sm) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.sm),df1.sm)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 44 = 25 sites + 19 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F); plot(g)
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.sm[order(match(traits.sm$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.sm$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.sm$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.sm$area, V(g)$name)) 
guild <- traits.sm[!(traits.sm[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="tomato"
E(g)$width <- 0.8
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=1000, vertex.label.cex=0.5, vertex.size=25*(igraph::degree(g, normalized = TRUE)),)

# Lizards

df1.liz<-df.liz[,colSums(df.liz) > 0]
str(df1.liz) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.liz),df1.liz)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 44 = 25 sites + 14 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.liz[order(match(traits.liz$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.liz$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.liz$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.liz$area, V(g)$name)) 
guild <- traits.liz[!(traits.liz[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="purple"
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=100, vertex.label.cex=0.6, vertex.size=25*(igraph::degree(g, normalized = TRUE)),)

# Mid-large mammals

df1.mid.large<-df.mid.large[,colSums(df.mid.large) > 0]
str(df1.mid.large) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.mid.large),df1.mid.large)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 53 = 25 sites + 28 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F); plot(g)
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.mid.large[order(match(traits.mid.large$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.mid.large$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.mid.large$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.mid.large$area, V(g)$name)) 
guild <- traits.mid.large[!(traits.mid.large[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="dark red"
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=100, vertex.label.cex=0.5, vertex.size=25*(igraph::degree(g, normalized = TRUE)),)

# Understorey birds

df1.birds<-df.birds[,colSums(df.birds) > 0]
str(df1.birds) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.birds),df1.birds)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 144 = 25 sites + 114 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
plot(g)
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.birds[order(match(traits.birds$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.birds$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.birds$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.birds$area, V(g)$name)) 
guild <- traits.birds[!(traits.birds[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="orange"
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=100, vertex.label.cex=0.5, vertex.size=25*(igraph::degree(g, normalized = TRUE)),)

# Frogs

df1.frogs<-df.frogs[,colSums(df.frogs) > 0]
str(df1.frogs) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.frogs),df1.frogs)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 144 = 25 sites + 114 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
plot(g)
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.frogs[order(match(traits.frogs$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.frogs$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.frogs$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.frogs$area, V(g)$name)) 
guild <- traits.frogs[!(traits.frogs[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="lightsteelblue2"
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=100, vertex.label.cex=0.65, vertex.size=25*(igraph::degree(g, normalized = TRUE)),)

# Dung beetles

df1.beetles<-df.beetles[,colSums(df.beetles) > 0]
str(df1.beetles) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.beetles),df1.beetles)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 144 = 25 sites + 114 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
plot(g)
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.beetles[order(match(traits.beetles$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.beetles$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.beetles$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.beetles$area, V(g)$name)) 
guild <- traits.beetles[!(traits.beetles[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="blue"
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=100, vertex.label.cex=0.8, 
     vertex.size=100 *(igraph::degree(g, normalized = TRUE)),)

# Orchid bees

df1.bees<-df.bees[,colSums(df.bees) > 0]
str(df1.bees) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.bees),df1.bees)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 144 = 25 sites + 114 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
plot(g)
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.bees[order(match(traits.bees$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.bees$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.bees$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.bees$area, V(g)$name)) 
guild <- traits.bees[!(traits.bees[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="gold"
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=100, vertex.label.cex=0.5, vertex.size=25*(igraph::degree(g, normalized = TRUE)),)

# Trees

df1.trees<-df.trees[,colSums(df.trees) > 0]
str(df1.trees) # check data
## preparing igraph object to igraph
net = cbind.data.frame(reference=row.names(df1.trees),df1.trees)
netlist = melt(net, na.rm = T) # transform to edge list
colnames(netlist) = c("site", "species", "weight")
netlist[,1]=as.character(paste(netlist[,1])) # sites
netlist[,2]=as.character(paste(netlist[,2])) # species
netlist2 <- subset(netlist, weight > 0) # remove non-existent interactions
nodes <- unique(data.frame(nodes = c(netlist2[,1], netlist2[,2])))
head(nodes)# check again
tail(nodes) # and again
nrow(nodes) # 144 = 25 sites + 114 sps
links = netlist2
## transform to igraph object
g <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
plot(g)
#Check the main properties of the igraph object
V(g); E(g)
#ordering dataframes so that the order of traits matches the order of species in graph g
guild <- traits.trees[order(match(traits.trees2$area, V(g)$name)),]; str(guild)
# check whether igraph object and traits df have same length
setdiff(traits.trees$area, V(g)$name) # character(0)
setdiff(V(g)$name, traits.trees$area) # character(0)
#removing species not present in the habitat matrix
del<-(setdiff(traits.trees$area, V(g)$name)) 
guild <- traits.trees[!(traits.trees[,1] %in% c(del)),]
# check again, both objects must have same length
setdiff(guild$area, V(g)$name) # character(0)
setdiff(V(g)$name, guild$area) # character(0)
## create a vector for traits to be graphed
V(g)$guild=as.character(guild$type)
V(g)$shape=V(g)$guild
V(g)$shape[V(g)$guild=="site"]="square"
V(g)$shape[V(g)$guild=="species"]="circle"
## to create a vector colors
V(g)$color=V(g)$guild
V(g)$color[V(g)$guild=="site"]="dark green" ## blue
V(g)$color[V(g)$guild=="species"]="yellowgreen"
plot(g,layout=layout.fruchterman.reingold,
     # main='Small mammals habitat-network',
     # Set node attributes
     vertex.shape=V(g)$shape, vertex.color = V(g)$color,
     # Set link colors
     edge.color = "light gray",
     # Set link curvature from 0 to 1
     edge.curved=0.3,
     # Set nodes parameters
     vertex.label.dist=1000, vertex.label.cex=0.5,
     vertex.size=25*(igraph::degree(g, normalized = TRUE)),)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# 2 - landscape-level properties (Figure 4)

# loading the packages needed
library(bipartite); library(vegan); library(reshape2); library(igraph)
library(networkD3); library(reshape2)

# 2.1 Connectance

# small mammals 

# observed value
connect.sm<-networklevel(df1.sm, level="both", weighted=FALSE, index="connectance")
# simulated value given the null model "r1"
null.sm.r1<-oecosimu(df1.sm,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.sm.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.sm.r1$oecosimu$simulated

# lizards 

# observed value
connect.liz<-networklevel(df1.liz, level="both", weighted=FALSE, index="connectance")
# simulated value given the null model "r1"
null.liz.r1<-oecosimu(df1.liz,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.liz.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.liz.r1$oecosimu$simulated

# mid-large mammals

# observed value
connect.mid.large<-networklevel(df1.mid.large, level="both", weighted=FALSE, index="connectance")
# simulated value given the null model "r1"
null.mid.large.r1<-oecosimu(df1.mid.large,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.mid.large.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.mid.large.r1$oecosimu$simulated

# understorey birds

# observed value
connect.birds<-networklevel(df1.birds, level="both", weighted=FALSE, index="connectance")
# simulated value given the null model "r1"
null.birds.r1<-oecosimu(df1.birds,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.birds.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.birds.r1$oecosimu$simulated

# frogs

# observed value
connect.frogs<-networklevel(df1.frogs, level="both", weighted=FALSE, index="connectance") 
# simulated value given the null model "r1"
null.frogs.r1<-oecosimu(df1.frogs,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.frogs.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.frogs.r1$oecosimu$simulated

# dung beetles

# observed value
connect.beetles<-networklevel(df1.beetles, level="both", weighted=FALSE, index="connectance") 
# simulated value given the null model "r1"
null.beetles.r1<-oecosimu(df1.beetles,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.beetles.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.beetles.r1$oecosimu$simulated

# bees

# observed value
connect.bees<-networklevel(df1.bees, level="both", weighted=FALSE, index="connectance") 
# simulated value given the null model "r1"
null.bees.r1<-oecosimu(df1.bees,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.bees.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.bees.r1$oecosimu$simulated

# trees

# observed value
connect.trees<-networklevel(df1.trees, level="both", weighted=FALSE, index="connectance") 
# simulated value given the null model "r1"
null.trees.r1<-oecosimu(df1.trees,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.trees.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.trees.r1$oecosimu$simulated

# all taxa

# observed value
connect.all.taxa<-networklevel(df, level="both", weighted=FALSE, index="connectance") 
# simulated value given the null model "r1"
null.all.taxa.r1<-oecosimu(df,networklevel, "r1", index="connectance",statistic = "c-score", nsim=100)
null.all.taxa.r1$oecosimu$z 
# obtaining the simulated values needed for the plot
null.all.taxa.r1$oecosimu$simulated


# 2.2 Modularity (using Beckett 2016' algorithm)

# Small mammals 

# observed value
mod.sm<-computeModules(df1.sm)
listModuleInformation(mod.sm)
printoutModuleInformation(mod.sm)

# simulated value given the null model "r1"
nulls.r1.sm<- (simulate(vegan::nullmodel(df1.sm, method="r1"), nsim = 100))
mod.nulls.r1.sm <- apply(nulls.r1.sm, 3, computeModules)
like.nulls.r1.sm <- sapply(mod.nulls.r1.sm, function(x) x@likelihood)
z_nulls.r1.sm <- (mod.sm@likelihood - mean(like.nulls.r1.sm))/sd(like.nulls.r1.sm)

# Lizards 

# observed value
mod.liz<-computeModules(df1.liz)
listModuleInformation(mod.liz)
printoutModuleInformation(mod.liz)

# simulated value given the null model "r1"
nulls.r1.liz<- (simulate(vegan::nullmodel(df1.liz, method="r1"), nsim = 100))
mod.nulls.r1.liz <- apply(nulls.r1.liz, 3, computeModules)
like.nulls.r1.liz <- sapply(mod.nulls.r1.liz, function(x) x@likelihood)
z_nulls.r1.liz <- (mod.liz@likelihood - mean(like.nulls.r1.liz))/sd(like.nulls.r1.liz)

# Mid-large mammals 

# observed value
mod.mid.large<-computeModules(df1.mid.large)
listModuleInformation(mod.mid.large)
printoutModuleInformation(mod.mid.large)

# simulated value given the null model "r1"
nulls.r1.mid.large<- (simulate(vegan::nullmodel(df1.mid.large, method="r1"), nsim = 100))
mod.nulls.r1.mid.large <- apply(nulls.r1.mid.large, 3, computeModules)
like.nulls.r1.mid.large <- sapply(mod.nulls.r1.mid.large, function(x) x@likelihood)
z_nulls.r1.mid.large <- (mod.mid.large@likelihood - mean(like.nulls.r1.mid.large))/sd(like.nulls.r1.mid.large)

# Birds

# observed value
mod.birds<-computeModules(df1.birds)
listModuleInformation(mod.birds)
printoutModuleInformation(mod.birds)

# simulated value given the null model "r1"
nulls.r1.birds<- (simulate(vegan::nullmodel(df1.birds, method="r1"), nsim = 100))
mod.nulls.r1.birds <- apply(nulls.r1.birds, 3, computeModules)
like.nulls.r1.birds <- sapply(mod.nulls.r1.birds, function(x) x@likelihood)
z_nulls.r1.birds <- (mod.birds@likelihood - mean(like.nulls.r1.birds))/sd(like.nulls.r1.birds)

# Frogs

# observed value
mod.frogs<-computeModules(df1.frogs)
listModuleInformation(mod.frogs)
printoutModuleInformation(mod.frogs)

# simulated value given the null model "r1"
nulls.r1.frogs<- (simulate(vegan::nullmodel(df1.frogs, method="r1"), nsim = 100))
mod.nulls.r1.frogs <- apply(nulls.r1.frogs, 3, computeModules)
like.nulls.r1.frogs <- sapply(mod.nulls.r1.frogs, function(x) x@likelihood)
z_nulls.r1.frogs <- (mod.frogs@likelihood - mean(like.nulls.r1.frogs))/sd(like.nulls.r1.frogs)


# Dung beetles

# observed value
mod.beetles<-computeModules(df1.beetles)
listModuleInformation(mod.beetles)
printoutModuleInformation(mod.beetles)

# simulated value given the null model "r1"
nulls.r1.beetles<- (simulate(vegan::nullmodel(df1.beetles, method="r1"), nsim = 100))
mod.nulls.r1.beetles <- apply(nulls.r1.beetles, 3, computeModules)
like.nulls.r1.beetles <- sapply(mod.nulls.r1.beetles, function(x) x@likelihood)
z_nulls.r1.beetles <- (mod.beetles@likelihood - mean(like.nulls.r1.beetles))/sd(like.nulls.r1.beetles)


# Orchid bees

# observed value
mod.bees<-computeModules(df1.bees)
listModuleInformation(mod.bees)
printoutModuleInformation(mod.bees)

# simulated value given the null model "r1"
nulls.r1.bees<- (simulate(vegan::nullmodel(df1.bees, method="r1"), nsim = 100))
mod.nulls.r1.bees <- apply(nulls.r1.bees, 3, computeModules)
like.nulls.r1.bees <- sapply(mod.nulls.r1.bees, function(x) x@likelihood)
z_nulls.r1.bees <- (mod.bees@likelihood - mean(like.nulls.r1.bees))/sd(like.nulls.r1.bees)

# Trees

# observed value
mod.trees<-computeModules(df1.trees)
listModuleInformation(mod.trees)
printoutModuleInformation(mod.bees)

# simulated value given the null model "r1"
nulls.r1.trees<- (simulate(vegan::nullmodel(df1.trees, method="r1"), nsim = 100))
mod.nulls.r1.trees <- apply(nulls.r1.trees, 3, computeModules)
like.nulls.r1.trees <- sapply(mod.nulls.r1.trees, function(x) x@likelihood)
z_nulls.r1.trees <- (mod.trees@likelihood - mean(like.nulls.r1.trees))/sd(like.nulls.r1.trees)

# All taxa

# observed value
mod.all.taxa<-computeModules(df)
listModuleInformation(mod.all.taxa)
printoutModuleInformation(mod.all.taxa)

# simulated value given the null model "r1"
nulls.r1.all.taxa<- (simulate(vegan::nullmodel(df1.all.taxa, method="r1"), nsim = 100))
mod.nulls.r1.all.taxa <- apply(nulls.r1.all.taxa, 3, computeModules)
like.nulls.r1.all.taxa <- sapply(mod.nulls.r1.all.taxa, function(x) x@likelihood)
z_nulls.r1.all.taxa <- (mod.all.taxa@likelihood - mean(like.nulls.r1.all.taxa))/sd(like.nulls.r1.all.taxa)

# 2.3 Nestedness (NODF)

# Small mammals

# observed value
nested_sm<-networklevel(df.sm, index="NODF")

# simulated value given the null model "r1"
null_r1_sm<-oecosimu(df.sm, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_sm$oecosimu$z #
# simulated values needed for the plot
null_r1.sm<-null_r1_sm$oecosimu$simulated
mean(null_r1.sm)

# Lizards

# observed value
nested_liz<-networklevel(df.liz, index="NODF")

# simulated value given the null model "r1"
null_r1_liz<-oecosimu(df.liz, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_liz$oecosimu$z #
# simulated values needed for the plot
null_r1.liz<-null_r1_liz$oecosimu$simulated
mean(null_r1.liz)


# Mid-large mammals

# observed value
nested_mid.large<-networklevel(df.mid.large, index="NODF")

# simulated value given the null model "r1"
null_r1_mid.large<-oecosimu(df.mid.large, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_mid.large$oecosimu$z #
# simulated values needed for the plot
null_r1.mid.large<-null_r1_mid.large$oecosimu$simulated
mean(null_r1.mid.large)

# Birds

# observed value
nested_birds<-networklevel(df.birds, index="NODF")

# simulated value given the null model "r1"
null_r1_birds<-oecosimu(df.birds, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_birds$oecosimu$z #
# simulated values needed for the plot
null_r1.birds<-null_r1_birds$oecosimu$simulated
mean(null_r1.birds)

# Frogs

# observed value
nested_frogs<-networklevel(df.frogs, index="NODF")

# simulated value given the null model "r1"
null_r1_frogs<-oecosimu(df.frogs, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_frogs$oecosimu$z #
# simulated values needed for the plot
null_r1.frogs<-null_r1_frogs$oecosimu$simulated
mean(null_r1.frogs)

# Dung beetles

# observed value
nested_beetles<-networklevel(df.beetles, index="NODF")

# simulated value given the null model "r1"
null_r1_beetles<-oecosimu(df.beetles, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_beetles$oecosimu$z #
# simulated values needed for the plot
null_r1.beetles<-null_r1_beetles$oecosimu$simulated
mean(null_r1.beetles)

# Orchid bees

# observed value
nested_bees<-networklevel(df.bees, index="NODF")

# simulated value given the null model "r1"
null_r1_bees<-oecosimu(df.bees, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_bees$oecosimu$z #
# simulated values needed for the plot
null_r1.bees<-null_r1_bees$oecosimu$simulated
mean(null_r1.bees)

# Trees

# observed value
nested_trees<-networklevel(df.trees, index="NODF")

# simulated value given the null model "r1"
null_r1_trees<-oecosimu(df.trees, networklevel,"r1",index= "NODF", statistic = "C.score", nsimul = 100)
null_r1_trees$oecosimu$z #
# simulated values needed for the plot
null_r1.trees<-null_r1_trees$oecosimu$simulated
mean(null_r1.trees)

# 2.4 Robustness

# 2.4.1 Largest-to-smallest

# extinction order
largest_to_smallest<-c(25,24,23,18,9,11,4,22,5,8,17,12,20,19,3,13,14,6,1,10,16,21,15,2,7)

# small mammals

# pattern of secondary extinctions
sm.sec.ext.lower<-second.extinct(df.sm, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(sm.sec.ext.lower, pch=19, cex=0.5)
# observed value
sm.rob.lower<-robustness(sm.sec.ext.lower)
# calculate random matrices
null_matrix_sm<-(simulate(vegan::nullmodel(df.sm, method="r1"), nsim = 100))
# transform array to a list
list_sm <- alply(null_matrix_sm,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_sm,second.extinct, participant="lower",
                  details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.sm.lower<-mean(null.rob_r1)
sd.sm.lower<-sd(null.rob_r1)
z_score_sm.lower<-(sm.rob.lower-mean.sm.lower)/sd.sm.lower

# Lizards

# pattern of secondary extinctions
liz.sec.ext.lower<-second.extinct(df.liz, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(liz.sec.ext.lower, pch=19, cex=0.5)
# observed value
liz.rob.lower<-robustness(liz.sec.ext.lower)
# calculate random matrices
null_matrix_liz<-(simulate(vegan::nullmodel(df.liz, method="r1"), nsim = 100))
# transform array to a list
list_liz = alply(null_matrix_liz,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_liz,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.liz.lower<-mean(null.rob_r1)
sd.liz.lower<-sd(null.rob_r1)
z_score_liz.lower<-(liz.rob.lower-mean.liz.lower)/sd.liz.lower

# Mid-large mammals

# pattern of secondary extinctions
mid.large.sec.ext.lower<-second.extinct(df.mid.large, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(mid.large.sec.ext.lower, pch=19, cex=0.5)
# observed value
mid.large.rob.lower<-robustness(mid.large.sec.ext.lower)
# calculate random matrices
null_matrix_mid.large<-(simulate(vegan::nullmodel(df.mid.large, method="r1"), nsim = 100))
# transform array to a list
list_mid.large = alply(null_matrix_mid.large,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_mid.large,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.mid.large.lower<-mean(null.rob_r1)
sd.mid.large.lower<-sd(null.rob_r1)
z_score_mid.large.lower<-(mid.large.rob.lower-mean.mid.large.lower)/sd.mid.large.lower

# Birds

# pattern of secondary extinctions
birds.sec.ext.lower<-second.extinct(df.birds, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(mid.large.sec.ext.lower, pch=19, cex=0.5)
# observed value
birds.rob.lower<-robustness(birds.sec.ext.lower)
# calculate random matrices
null_matrix_birds<-(simulate(vegan::nullmodel(df.birds, method="r1"), nsim = 100))
# transform array to a list
list_birds = alply(null_matrix_birds,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_birds,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.birds.lower<-mean(null.rob_r1)
sd.birds.lower<-sd(null.rob_r1)
z_score_birds.lower<-(birds.rob.lower-mean.birds.lower)/sd.birds.lower


# Frogs

# pattern of secondary extinctions
frogs.sec.ext.lower<-second.extinct(df.frogs, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(frogs.sec.ext.lower, pch=19, cex=0.5)
# observed value
frogs.rob.lower<-robustness(frogs.sec.ext.lower)
# calculate random matrices
null_matrix_frogs<-(simulate(vegan::nullmodel(df.frogs, method="r1"), nsim = 100))
# transform array to a list
list_frogs = alply(null_matrix_frogs,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_frogs,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.frogs.lower<-mean(null.rob_r1)
sd.frogs.lower<-sd(null.rob_r1)
z_score_frogs.lower<-(frogs.rob.lower-mean.frogs.lower)/sd.frogs.lower


# Dung-beetles - here we used the dung-beetles data excluding sites with zero dung-beetles records

df.beetles2<-read.csv("beetles_no_zeros.csv", header = T, row.names = 1) 
largest_to_smallest_beetles<-c(17,18,19,13,6,7,2,16,3,5,12,8,14,1,9,10,4,11,15)

# pattern of secondary extinctions
beetles.sec.ext.lower<-second.extinct(df.beetles2, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest_beetles)
slope.bipartite(beetles.sec.ext.lower, pch=19, cex=0.5)
# observed value
beetles.rob.lower<-robustness(beetles.sec.ext.lower)
# calculate random matrices
null_matrix_beetles<-(simulate(vegan::nullmodel(df.beetles2, method="r1"), nsim = 100))
# transform array to a list
list_beetles = alply(null_matrix_beetles,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_beetles,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest_beetles)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.beetles.lower<-mean(null.rob_r1)
sd.beetles.lower<-sd(null.rob_r1)
z_score_beetles.lower<-(beetles.rob.lower-mean.beetles.lower)/sd.beetles.lower


# Orchid bees

# pattern of secondary extinctions
bees.sec.ext.lower<-second.extinct(df.bees, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(bees.sec.ext.lower, pch=19, cex=0.5)
# observed value
bees.rob.lower<-robustness(bees.sec.ext.lower)
# calculate random matrices
null_matrix_bees<-(simulate(vegan::nullmodel(df.bees, method="r1"), nsim = 100))
# transform array to a list
list_bees <- alply(null_matrix_bees,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_bees,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.bees.lower<-mean(null.rob_r1)
sd.bees.lower<-sd(null.rob_r1)
z_score_bees.lower<-(bees.rob.lower-mean.bees.lower)/sd.bees.lower


# Trees

# pattern of secondary extinctions
trees.sec.ext.lower<-second.extinct(df.trees, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(trees.sec.ext.lower, pch=19, cex=0.5)
# observed value
trees.rob.lower<-robustness(trees.sec.ext.lower)
# calculate random matrices
null_matrix_trees<-(simulate(vegan::nullmodel(df.trees, method="r1"), nsim = 100))
# transform array to a list
list_trees <- alply(null_matrix_bees,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_trees,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.trees.lower<-mean(null.rob_r1)
sd.trees.lower<-sd(null.rob_r1)
z_score_trees.lower<-(trees.rob.lower-mean.trees.lower)/sd.trees.lower

# All

# pattern of secondary extinctions
all.sec.ext.lower<-second.extinct(df, participant="lower",details="TRUE",method="external", ext.row = largest_to_smallest)
slope.bipartite(all.sec.ext.lower, pch=19, cex=0.5)
# observed value
all.rob.lower<-robustness(all.sec.ext.lower)
# calculate random matrices
null_matrix_all<-(simulate(vegan::nullmodel(df, method="r1"), nsim = 100))
# transform array to a list
list_all <- alply(null_matrix_all,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_all,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = largest_to_smallest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.all.lower<-mean(null.rob_r1)
sd.all.lower<-sd(null.rob_r1)
z_score_all.lower<-(all.rob.lower-mean.all.lower)/sd.all.lower


# 2.4.2 Smallest-to-largest

# extinction order
smallest_to_largest<-c(7,2,15,21,16,10,1,6,14,13,3,19,20,12,17,8,5,22,4,11,9,18,23,24,25)

# Small mammals

# pattern of secondary extinctions
sm.sec.ext.lower<-second.extinct(df.sm, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(sm.sec.ext.lower, pch=19, cex=0.5)
# observed value
sm.rob.lower<-robustness(sm.sec.ext.lower)
# calculate random matrices
null_matrix_sm<-(simulate(vegan::nullmodel(df.sm, method="r1"), nsim = 100))
# transform array to a list
list_sm <- alply(null_matrix_sm,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_sm,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.sm.lower<-mean(null.rob_r1)
sd.sm.lower<-sd(null.rob_r1)
z_score_sm.lower<-(sm.rob.lower-mean.sm.lower)/sd.sm.lower

# Lizards

# pattern of secondary extinctions
liz.sec.ext.lower<-second.extinct(df.liz, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(liz.sec.ext.lower, pch=19, cex=0.5)
# observed value
liz.rob.lower<-robustness(liz.sec.ext.lower)
# calculate random matrices
null_matrix_liz<-(simulate(vegan::nullmodel(df.liz, method="r1"), nsim = 100))
# transform array to a list
list_liz <- alply(null_matrix_liz,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_liz,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.liz.lower<-mean(null.rob_r1)
sd.liz.lower<-sd(null.rob_r1)
z_score_liz.lower<-(liz.rob.lower-mean.liz.lower)/sd.liz.lower

# Mid-large mammals

# pattern of secondary extinctions
mid.large.sec.ext.lower<-second.extinct(df.mid.large, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(mid.large.sec.ext.lower, pch=19, cex=0.5)
# observed value
mid.large.rob.lower<-robustness(mid.large.sec.ext.lower)
# calculate random matrices
null_matrix_mid.large<-(simulate(vegan::nullmodel(df.mid.large, method="r1"), nsim = 100))
# transform array to a list
list_mid.large <- alply(null_matrix_mid.large,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_mid.large,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.mid.large.lower<-mean(null.rob_r1)
sd.mid.large.lower<-sd(null.rob_r1)
z_score_mid.large.lower<-(mid.large.rob.lower-mean.mid.large.lower)/sd.mid.large.lower

# Birds

# pattern of secondary extinctions
birds.sec.ext.lower<-second.extinct(df.birds, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(birds.sec.ext.lower, pch=19, cex=0.5)
# observed value
birds.rob.lower<-robustness(birds.sec.ext.lower)
# calculate random matrices
null_matrix_birds<-(simulate(vegan::nullmodel(df.birds, method="r1"), nsim = 100))
# transform array to a list
list_birds <- alply(null_matrix_birds,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_birds,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.birds.lower<-mean(null.rob_r1)
sd.birds.lower<-sd(null.rob_r1)
z_score_birds.lower<-(birds.rob.lower-mean.birds.lower)/sd.birds.lower

# Frogs

# pattern of secondary extinctions
frogs.sec.ext.lower<-second.extinct(df.frogs, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(frogs.sec.ext.lower, pch=19, cex=0.5)
# observed value
frogs.rob.lower<-robustness(frogs.sec.ext.lower)
# calculate random matrices
null_matrix_frogs<-(simulate(vegan::nullmodel(df.frogs, method="r1"), nsim = 100))
# transform array to a list
list_frogs <- alply(null_matrix_frogs,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_frogs,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.frogs.lower<-mean(null.rob_r1)
sd.frogs.lower<-sd(null.rob_r1)
z_score_frogs.lower<-(frogs.rob.lower-mean.frogs.lower)/sd.frogs.lower

# Dung-beetles

smallest_to_largest_beetles<-c(15,11,4,10,9,1,14,8,12,5,3,16,2,7,6,13,19,18,17)

# pattern of secondary extinctions
beetles.sec.ext.lower<-second.extinct(df.beetles2, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest_beetles)
slope.bipartite(beetles.sec.ext.lower, pch=19, cex=0.5)
# observed value
beetles.rob.lower<-robustness(beetles.sec.ext.lower)
# calculate random matrices
null_matrix_beetles<-(simulate(vegan::nullmodel(df.beetles2, method="r1"), nsim = 100))
# transform array to a list
list_beetles <- alply(null_matrix_beetles,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_beetles,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest_beetles)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.beetles.lower<-mean(null.rob_r1)
sd.beetles.lower<-sd(null.rob_r1)
z_score_beetles.lower<-(beetles.rob.lower-mean.beetles.lower)/sd.beetles.lower

# Orchid bees

# pattern of secondary extinctions
bees.sec.ext.lower<-second.extinct(df.bees, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(bees.sec.ext.lower, pch=19, cex=0.5)
# observed value
bees.rob.lower<-robustness(bees.sec.ext.lower)
# calculate random matrices
null_matrix_bees<-(simulate(vegan::nullmodel(df.bees, method="r1"), nsim = 100))
# transform array to a list
list_bees <- alply(null_matrix_bees,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_bees,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.bees.lower<-mean(null.rob_r1)
sd.bees.lower<-sd(null.rob_r1)
z_score_bees.lower<-(bees.rob.lower-mean.bees.lower)/sd.bees.lower

# Trees

# pattern of secondary extinctions
trees.sec.ext.lower<-second.extinct(df.trees, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(trees.sec.ext.lower, pch=19, cex=0.5)
# observed value
trees.rob.lower<-robustness(trees.sec.ext.lower)
# calculate random matrices
null_matrix_trees<-(simulate(vegan::nullmodel(df.trees, method="r1"), nsim = 100))
# transform array to a list
list_trees <- alply(null_matrix_trees,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_trees,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.trees.lower<-mean(null.rob_r1)
sd.trees.lower<-sd(null.rob_r1)
z_score_trees.lower<-(trees.rob.lower-mean.trees.lower)/sd.trees.lower

# All

# pattern of secondary extinctions
all.sec.ext.lower<-second.extinct(df, participant="lower",details="TRUE",method="external", ext.row = smallest_to_largest)
slope.bipartite(all.sec.ext.lower, pch=19, cex=0.5)
# observed value
all.rob.lower<-robustness(all.sec.ext.lower)
# calculate random matrices
null_matrix_all<-(simulate(vegan::nullmodel(df.all, method="r1"), nsim = 100))
# transform array to a list
list_all <- alply(null_matrix_all,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_all,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = smallest_to_largest)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.all.lower<-mean(null.rob_r1)
sd.all.lower<-sd(null.rob_r1)
z_score_all.lower<-(all.rob.lower-mean.all.lower)/sd.all.lower


# 2.4.3 Random

# Extinction order
random<-c(14,15,2,11,5,12,23,17,16,10,13,21,7,9,24,1,22,20,3,4,18,6,19,25,8)

# Small mammals

# pattern of secondary extinctions
sm.sec.ext.lower<-second.extinct(df.sm, participant="lower", details="TRUE", method="external", ext.row = random)
slope.bipartite(sm.sec.ext.lower, pch=19, cex=0.5)
# observed value
sm.rob.lower<-robustness(sm.sec.ext.lower)
# calculate random matrices
null_matrix_sm<-(simulate(vegan::nullmodel(df.sm, method="r1"), nsim = 100))
# transform array to a list
list_sm <- alply(null_matrix_sm,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_sm,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.sm.lower<-mean(null.rob_r1)
sd.sm.lower<-sd(null.rob_r1)
z_score_sm.lower<-(sm.rob.lower-mean.sm.lower)/sd.sm.lower

# Lizards

# pattern of secondary extinctions
liz.sec.ext.lower<-second.extinct(df.liz, participant="lower", details="TRUE", method="external", ext.row = random)
slope.bipartite(liz.sec.ext.lower, pch=19, cex=0.5)
# observed value
liz.rob.lower<-robustness(liz.sec.ext.lower)
# calculate random matrices
null_matrix_liz<-(simulate(vegan::nullmodel(df.liz, method="r1"), nsim = 100))
# transform array to a list
list_liz <- alply(null_matrix_liz,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_liz,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.liz.lower<-mean(null.rob_r1)
sd.liz.lower<-sd(null.rob_r1)
z_score_liz.lower<-(liz.rob.lower-mean.liz.lower)/sd.liz.lower

# Mid-large mammals

# pattern of secondary extinctions
mid.large.sec.ext.lower<-second.extinct(df.mid.large, participant="lower", details="TRUE", method="external", ext.row = random)
slope.bipartite(mid.large.sec.ext.lower, pch=19, cex=0.5)
# observed value
mid.large.rob.lower<-robustness(mid.large.sec.ext.lower)
# calculate random matrices
null_matrix_mid.large<-(simulate(vegan::nullmodel(df.mid.large, method="r1"), nsim = 100))
# transform array to a list
list_mid.large <- alply(null_matrix_mid.large,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_mid.large,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.mid.large.lower<-mean(null.rob_r1)
sd.mid.large.lower<-sd(null.rob_r1)
z_score_mid.large.lower<-(mid.large.rob.lower-mean.mid.large.lower)/sd.mid.large.lower

# Birds

# pattern of secondary extinctions
birds.sec.ext.lower<-second.extinct(df.birds, participant="lower", details="TRUE", method="external", ext.row = random)
slope.bipartite(birds.sec.ext.lower, pch=19, cex=0.5)
# observed value
birds.rob.lower<-robustness(birds.sec.ext.lower)
# calculate random matrices
null_matrix_birds<-(simulate(vegan::nullmodel(df.birds, method="r1"), nsim = 100))
# transform array to a list
list_birds <- alply(null_matrix_birds,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_birds,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.birds.lower<-mean(null.rob_r1)
sd.birds.lower<-sd(null.rob_r1)
z_score_birds.lower<-(birds.rob.lower-mean.birds.lower)/sd.birds.lower

# Frogs

# pattern of secondary extinctions
frogs.sec.ext.lower<-second.extinct(df.frogs, participant="lower", details="TRUE", method="external", ext.row = random)
slope.bipartite(frogs.sec.ext.lower, pch=19, cex=0.5)
# observed value
frogs.rob.lower<-robustness(frogs.sec.ext.lower)
# calculate random matrices
null_matrix_frogs<-(simulate(vegan::nullmodel(df.frogs, method="r1"), nsim = 100))
# transform array to a list
list_frogs <- alply(null_matrix_frogs,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_frogs,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.frogs.lower<-mean(null.rob_r1)
sd.frogs.lower<-sd(null.rob_r1)
z_score_frogs.lower<-(frogs.rob.lower-mean.frogs.lower)/sd.frogs.lower

# Dung-beetles

random_beetles<-c(14,15,2,11,5,12,17,16,10,13,7,9,1,3,4,18,6,19,8)

# pattern of secondary extinctions
beetles.sec.ext.lower<-second.extinct(df.beetles2, participant="lower", details="TRUE", method="external", ext.row = random_beetles)
slope.bipartite(beetles.sec.ext.lower, pch=19, cex=0.5)
# observed value
beetles.rob.lower<-robustness(beetles.sec.ext.lower)
# calculate random matrices
null_matrix_beetles<-(simulate(vegan::nullmodel(df.beetles, method="r1"), nsim = 100))
# transform array to a list
list_beetles <- alply(null_matrix_beetles,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_beetles,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random_beetles)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.beetles.lower<-mean(null.rob_r1)
sd.beetles.lower<-sd(null.rob_r1)
z_score_beetles.lower<-(beetles.rob.lower-mean.beetles.lower)/sd.beetles.lower

# Orchid bees

# pattern of secondary extinctions
bees.sec.ext.lower<-second.extinct(df.bees, participant="lower", details="TRUE", method="external", ext.row = random)
slope.bipartite(bees.sec.ext.lower, pch=19, cex=0.5)
# observed value
bees.rob.lower<-robustness(bees.sec.ext.lower)
# calculate random matrices
null_matrix_bees<-(simulate(vegan::nullmodel(df.bees, method="r1"), nsim = 100))
# transform array to a list
list_bees <- alply(null_matrix_bees,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_bees,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.bees.lower<-mean(null.rob_r1)
sd.bees.lower<-sd(null.rob_r1)
z_score_bees.lower<-(bees.rob.lower-mean.bees.lower)/sd.bees.lower

# Trees

# pattern of secondary extinctions
trees.sec.ext.lower<-second.extinct(df.trees, participant="lower", details="TRUE", method="external", ext.row = random)
slope.bipartite(trees.sec.ext.lower, pch=19, cex=0.5)
# observed value
trees.rob.lower<-robustness(trees.sec.ext.lower)
# calculate random matrices
null_matrix_trees<-(simulate(vegan::nullmodel(df.trees, method="r1"), nsim = 100))
# transform array to a list
list_bees <- alply(null_matrix_trees,3) 
# use lapply to run second.extinct and robustness in a list
list_sec <- lapply(list_trees,second.extinct, participant="lower",
                   details="TRUE", method="external", ext.row = random)
# simulated values
rob <- lapply(list_sec, robustness)
# calculate z-score
null.rob_r1<-unlist(rob) # simulated values 
mean.trees.lower<-mean(null.rob_r1)
sd.trees.lower<-sd(null.rob_r1)
z_score_trees.lower<-(trees.rob.lower-mean.trees.lower)/sd.trees.lower

# 2.5 Figure 4

# 2.5.1 Simulated values (obtained using the null model r1)
null.models<-read.table("null_models.txt", header = T) 

# Connectance
null.models$taxa<-factor(null.models$taxa, levels =unique(null.models$taxa[order(null.models$connect_r1 )]))
p.fill2<-c("blue","orange","light blue","black","yellowgreen","gold","tomato","purple", "dark red")
p2 <- ggplot(null.models, aes(null.models[,1], null.models[,2], colour= null.models[,1]))
p2 + geom_boxplot(colour = p.fill2) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Connectance") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') 

# Modularity
p.fill5<-c("dark red", "yellowgreen","black","purple", "tomato", "gold","light blue", "orange","blue")
null.models$taxa<-factor(null.models$taxa, levels =unique(null.models$taxa[order(null.models$mod_r1 )]))
p5 <- ggplot(null.models, aes(null.models[,1], null.models[,3], colour= null.models[,1]))
p5 + geom_boxplot(colour= p.fill5) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Modularity") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') +
  scale_y_continuous(limits = c(0.1,0.4)) 

# Nestedness
p.fill7<-c("blue","orange","light blue","yellowgreen", "black","gold", "tomato", "purple", "dark red")
null.models$taxa<-factor(null.models$taxa, levels =unique(null.models$taxa[order(null.models$nest_r1)]))

p7 <- ggplot(null.models, aes(null.models[,1], null.models[,4], colour= null.models[,1]))
p7 + geom_boxplot(colour = p.fill7) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Nestedness") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') +
  scale_y_continuous(limits = c(30,90))

# Robustness (sites) - ordered from largest to smallest
p.fill8<-c("blue","orange","light blue", "yellowgreen","black", "tomato", "purple", "gold", "dark red")
null.models$taxa<-factor(null.models$taxa, levels =unique(null.models$taxa[order(null.models$robustness_LL_r1_large_to_small)]))
p8 <- ggplot(null.models, aes(null.models[,1], null.models[,5], colour= null.models[,1]))
p8 + geom_boxplot(colour = p.fill8) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Robustness") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') 
  

# 2.5.2 Observed values
obs<-read.table("observed_landscape_metrics.txt", header = T) 

# Connectance
obs2<-as.data.frame(obs[1:9,])
obs2$taxa<-factor(obs2$taxa, levels =unique(obs2$taxa[order(obs2$value)]))
p3 <- ggplot(obs2, aes(obs2[,1], obs2[,2]))
p3 + geom_point(size = 3.5, stroke = 0.65, color = "black",  alpha = 1, shape = 16) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Connectance") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') +
  scale_y_continuous(limits = c(0.15,0.5)) 

# Modularity
obs3<-as.data.frame(obs[10:18,])
obs3$taxa<-factor(obs3$taxa, levels =unique(obs3$taxa[order(obs3$value)]))
p6 <- ggplot(obs3, aes(obs3[,1], obs3[,2]))
p6 + geom_point(size = 3.5, stroke = 0.65, color = "black",  alpha = 1, shape = 16) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Modularity") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') +
  scale_y_continuous(limits = c(0.1,0.4)) 

# Nestedness
obs4<-as.data.frame(obs[19:27,])
obs4$taxa<-factor(obs4$taxa, levels =unique(obs4$taxa[order(obs4$value)]))
p8 <- ggplot(obs4, aes(obs4[,1], obs4[,2]))
p8 + geom_point(size = 3.5, stroke = 0.65, color = "black",  alpha = 1, shape = 16) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Nestedness") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') +
  scale_y_continuous(limits = c(30,90)) 

# Robustness
obs5<-as.data.frame(obs[28:36,])
obs5$taxa<-factor(obs5$taxa, levels =unique(obs5$taxa[order(obs5$value)]))

p8 <- ggplot(obs5, aes(obs5[,1], obs5[,2]))
p8 + geom_point(size = 3.5, stroke = 0.65, color = "black",  alpha = 1, shape = 16) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Robustness") +
  theme(axis.title.x = element_text(size=0,color='black'),
        axis.title.y = element_text(size=18,color='black'), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) +
  theme(legend.position='none') 

# 2.6 Figure 5

# Largest-to-smallest
rob<-read.table("robustness_plot_largest_to_smallest.txt", header = T) 
robustness<-as.data.frame(rob)
p.fill2<-c("orange","light blue","yellowgreen","gold","tomato","purple", "dark red","blue")
p1 <- ggplot(robustness, aes(x = -log10(robustness[,1]), y = robustness[,2], group = robustness[,3], color = robustness[,3]))
p1 + geom_point(size = 3, alpha = 0.5) +  
  geom_smooth(span = 1.5, se = F) +
  scale_colour_manual(values=c("gold","blue","orange","light blue","purple","dark red","tomato","yellowgreen")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y=expression(paste("% remaining species")),
       x=expression(paste("Remaining forest area (ha)"))) +
  theme(legend.position="none") +
  theme(plot.title = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=18,color='black', face="bold"),
        axis.title.y = element_text(size=18,color='black', face="bold"), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15))

# Smallest-to-largest
rob2<-read.table("robustness_plot_smallest_to_largest.txt", header = T) 
robustness<-as.data.frame(rob2)
p.fill2<-c("orange","light blue","yellowgreen","gold","tomato","purple", "dark red","blue")
p1 <- ggplot(robustness, aes(x = log10(robustness[,2]), y = robustness[,3], group = robustness[,4], color = robustness[,4]))
p1 + geom_point(size = 3, alpha = 0.5) +  
  geom_smooth(span = 1.5, se = F) +
  scale_colour_manual(values=c("gold","blue","orange","light blue","purple","dark red","tomato","yellowgreen")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y=expression(paste("% remaining species")),
       x=expression(paste("Remaining forest area (ha)"))) +
  theme(legend.position="none") +
  theme(plot.title = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=18,color='black', face="bold"),
        axis.title.y = element_text(size=18,color='black', face="bold"), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15))

# 3 - Site and species-level properties (node properties)

# 3.1 Normalised degree

# 3.1.1. Site-level
# Small mammals
k_lower.sm<-specieslevel(df.sm, index="normalised degree",level="lower") 
# Lizard
k_lower.liz<-specieslevel(df.liz, index="normalised degree",level="lower") 
# Mid-large mammals
k_lower.mid.large<-specieslevel(df.mid.large, index="normalised degree",level="lower")
# Understorey birds
k_lower.birds<-specieslevel(df.birds, index="normalised degree",level="lower") 
# Frogs
k_lower.frogs<-specieslevel(df.frogs, index="normalised degree",level="lower") 
# Dung-beetles
k_lower.beetles<-specieslevel(df.beetles, index="normalised degree",level="lower") 
# Orchid bees
k_lower.bees<-specieslevel(df.bees, index="degree",level="lower") 
# Trees
k_lower.trees<-specieslevel(df.trees, index="normalised degree",level="lower") 
# All.taxa
k_lower.all.taxa<-specieslevel(df, index="degree",level="lower") 

# 3.1.1. Species-level
# Small mammals
k_higher.sm<-specieslevel(df.sm, index="normalised degree",level="higher") 
# Lizards
k_higher.liz<-specieslevel(df.liz, index="normalised degree",level="higher") 
# Mid-large mammals
k_higher.mid.large<-specieslevel(df.mid.large, index="normalised degree",level="higher")
# Understorey birds
k_higher.birds<-specieslevel(df.birds, index="normalised degree",level="higher") 
# Frogs
k_higher.frogs<-specieslevel(df.frogs, index="normalised degree",level="higher") 
# Dung-beetles
k_higher.beetles<-specieslevel(df.beetles, index="normalised degree",level="higher") 
# Orchid bees
k_higher.bees<-specieslevel(df.bees, index="degree",level="higher") 
# Trees
k_higher.trees<-specieslevel(df.trees, index="normalised degree",level="higher") 
# All taxa
k_higher.all.taxa<-specieslevel(df, index="degree",level="higher") 

# 3.2 Site and species nestedness contribution

# Small mammals
nest.cont.sm<-nestedcontribution(df.sm, nsimul = 100)
# Lizards
nest.cont.liz<-nestedcontribution(df.liz, nsimul = 100)
# Mid-large mammals
nest.cont.mid.large<-nestedcontribution(df.mid.large, nsimul = 100)
# Understorey birds
nest.cont.birds<-nestedcontribution(df.birds, nsimul = 100)
# Frogs
nest.cont.frogs<-nestedcontribution(df.frogs, nsimul = 100)
# Dung-beetles
nest.cont.beetles<-nestedcontribution(df.beetles, nsimul = 100)
# Orchid bees
nest.cont.bees<-nestedcontribution(df.bees, nsimul = 100)
# Trees
nest.cont.trees<-nestedcontribution(df.trees, nsimul = 100)
# All taxa
nest.cont.all<-nestedcontribution(df, nsimul = 100)


###############################################################################


# 4.3  Figure 5
#pannel A
nest.fig<-read.table("nest_cont_sites_figure2.txt", header = T) 

nest.plot<-ggballoonplot(nest.fig, y="variables",
                         font.label = list(type=2,size = 44, color = "black"),
                         color = "grey", fill = "grey", size.range = c(3, 20)) 

#pannel B
nest.fig<-read.table("nest_cont_sites_figure2.txt", header = T) 

nest.plot<-ggballoonplot(nest.fig, y="variables",
                         font.label = list(type=2,size = 44, color = "black"),
                         color = "grey", fill = "grey", size.range = c(3, 20)) 


# Figure 6

# 6 - Site-level modelling

# Figure 7


##################################################################################
# Figures from Supplementary Materials

node_sites<-read.table("node_sites_metrics.txt", header = T) 

# Figure S1 - Histograms showing normalised degree at the site level, color-coded by forest area

# Small mammals
sites.sm<-node_sites[node_sites$Taxonomic_group == "sm",]
p1 <- ggplot(sites.sm, aes(reorder(sites.sm$site_code, sites.sm$normalised_degree), sites.sm$normalised_degree, fill = log10(area)))
p1 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95)) 

# Lizards
sites.liz<-node_sites[node_sites$Taxonomic_group == "liz",]
p2 <- ggplot(sites.liz, aes(reorder(sites.liz$site_code, sites.liz$normalised_degree), sites.liz$normalised_degree, fill = log10(area)))
p2 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95)) 

# Mid-large mammals
sites.mid.large<-node_sites[node_sites$Taxonomic_group == "mid.large",]
p3 <- ggplot(sites.mid.large, aes(reorder(sites.mid.large$site_code, sites.mid.large$normalised_degree), sites.mid.large$normalised_degree, fill = log10(area)))
p3 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10)) +
  scale_y_continuous(limits = c(0,0.95)) 

# Understorey birds
sites.birds<-node_sites[node_sites$Taxonomic_group == "birds",]
p4 <- ggplot(sites.birds, aes(reorder(sites.birds$site_code, sites.birds$normalised_degree), sites.birds$normalised_degree, fill = log10(area)))
p4 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95)) 

# Frogs
sites.frogs<-node_sites[node_sites$Taxonomic_group == "frogs",]
p5 <- ggplot(sites.frogs, aes(reorder(sites.frogs$site_code, sites.frogs$normalised_degree), sites.frogs$normalised_degree, fill = log10(area)))
p5 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95)) 

# Dung beetles
sites.beetles<-node_sites[node_sites$Taxonomic_group == "beetles",]
p6 <- ggplot(sites.beetles, aes(reorder(sites.beetles$site_code, sites.beetles$normalised_degree), sites.beetles$normalised_degree, fill = log10(area)))
p6 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95)) 

# Orchid bees
sites.bees<-node_sites[node_sites$Taxonomic_group == "bees",]
p7 <- ggplot(sites.bees, aes(reorder(sites.bees$site_code, sites.bees$normalised_degree), sites.bees$normalised_degree, fill = log10(area)))
p7 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95)) 

# Trees
sites.trees<-node_sites[node_sites$Taxonomic_group == "trees",]
p8 <- ggplot(sites.trees, aes(reorder(sites.trees$site_code, sites.trees$normalised_degree), sites.trees$normalised_degree, fill = log10(area)))
p8 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95)) 

# All taxa
sites.all<-node_sites[node_sites$Taxonomic_group == "all_taxa",]
p8 <- ggplot(sites.all, aes(reorder(sites.all$site_code, sites.all$normalised_degree), sites.all$normalised_degree, fill = log10(area)))
p8 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(0,0.95))

# Figure S2 - Histograms showing nestedness contribution at the site level, color-coded by forest area

# Small mammals
sites.sm<-node_sites[node_sites$Taxonomic_group == "sm",]
p1 <- ggplot(sites.sm, aes(reorder(site_code, nested_contribution), nested_contribution, fill = log10(area)))
p1 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(-1,10)) 

# Lizards
sites.liz<-node_sites[node_sites$Taxonomic_group == "liz",]
p2 <- ggplot(sites.liz, aes(reorder(site_code, nested_contribution),nested_contribution, fill = log10(area)))
p2 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(-1,10)) 

# Mid-large mammals
sites.mid.large<-node_sites[node_sites$Taxonomic_group == "mid.large",]
p3 <- ggplot(sites.mid.large, aes(reorder(site_code, nested_contribution), nested_contribution, fill = log10(area)))
p3 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10)) +
  scale_y_continuous(limits = c(-1,10)) 

# Understorey birds
sites.birds<-node_sites[node_sites$Taxonomic_group == "birds",]
p4 <- ggplot(sites.birds, aes(reorder(sites.birds$site_code, sites.birds$nested_contribution), sites.birds$nested_contribution, fill = log10(area)))
p4 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10)) +
  scale_y_continuous(limits = c(-1,10)) 

# Frogs
sites.frogs<-node_sites[node_sites$Taxonomic_group == "frogs",]
p5 <- ggplot(sites.frogs, aes(reorder(sites.frogs$site_code, sites.frogs$nested_contribution), sites.frogs$nested_contribution, fill = log10(area)))
p5 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(-1,10)) 

# Dung beetles
sites.beetles<-node_sites[node_sites$Taxonomic_group == "beetles",]
p6 <- ggplot(sites.beetles, aes(reorder(sites.beetles$site_code, sites.beetles$nested_contribution), sites.beetles$nested_contribution, fill = log10(area)))
p6 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(-1,10)) 

# Orchid bees
sites.bees<-node_sites[node_sites$Taxonomic_group == "bees",]
p7 <- ggplot(sites.bees, aes(reorder(sites.bees$site_code, sites.bees$nested_contribution), sites.bees$nested_contribution, fill = log10(area)))
p7 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10))+
  scale_y_continuous(limits = c(-1,10)) 

# Trees
sites.trees<-node_sites[node_sites$Taxonomic_group == "trees",]
p8 <- ggplot(sites.trees, aes(reorder(sites.trees$site_code, sites.trees$nested_contribution), sites.trees$nested_contribution, fill = log10(area)))
p8 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10)) +
  scale_y_continuous(limits = c(-1,10)) 

# All taxa
sites.all<-node_sites[node_sites$Taxonomic_group == "all_taxa",]
p8 <- ggplot(sites.all, aes(reorder(sites.all$site_code, sites.all$nested_contribution), sites.all$nested_contribution, fill = log10(area)))
p8 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="darkgreen", mid = "yellow", low="red", midpoint = 2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 10),
    axis.text.x=element_text(size=10)) +
  scale_y_continuous(limits = c(-1,10)) 


# Figure S3 - Histograms showing normalised degree at the species level, color-coded by species body mass

node_species<-read.table("node_species_metrics.txt", header = T) 

# Small mammals
species.sm<-node_species[node_species$taxa == "sm",]
p1 <- ggplot(species.sm, aes(reorder(species.sm$species, species.sm$normalized.degree),species.sm$normalized.degree, fill = log10(body_size)))
p1 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.8) +
  #theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12))+
  scale_y_continuous(limits = c(0,1)) 

# Lizards
species.liz<-node_species[node_species$taxa == "liz",]
p2 <- ggplot(species.liz, aes(reorder(species.liz$species, species.liz$normalized.degree),species.liz$normalized.degree, fill =  log10(svl)))
p2 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.8) +
  #  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12)) +
  scale_y_continuous(limits = c(0,1)) 

# Mid-large mammals
species.mid.large<-node_species[node_species$taxa == "mid.large",]
p3 <- ggplot(species.mid.large, aes(reorder(species.mid.large$species, species.mid.large$normalized.degree),species.mid.large$normalized.degree, fill = log10(body_mass)))
p3 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 0.7) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12)) +
  scale_y_continuous(limits = c(0,1)) 

# Understorey birds
species.birds<-node_species[node_species$taxa == "birds",]
p4 <- ggplot(species.birds, aes(reorder(species.birds$species, species.birds$normalized.degree),species.birds$normalized.degree, fill=log10(body_size)))
p4 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.2) +
  #theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 5),
    axis.text.x=element_text(size=12)) +
  scale_y_continuous(limits = c(0,1)) 

# Frogs
species.frogs<-read.table("node_species_metrics_frogs.txt", header = T) 
p5 <- ggplot(species.frogs, aes(reorder(species.frogs$species, species.frogs$normalized.degree),species.frogs$normalized.degree, fill=log10(body_size_mm)))
p5 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.7) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12)) +
  scale_y_continuous(limits = c(0,1)) 

# Dung beetles
species.beetles<-node_species[node_species$taxa == "beetles",]
p6 <- ggplot(species.beetles, aes(reorder(species.beetles$species, species.beetles$normalized.degree),species.beetles$normalized.degree,fill=log10(body_size)))
p6 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = -2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12)) +
  scale_y_continuous(limits = c(0,1)) 

# Orchid bees
species.bees<-node_species[node_species$taxa == "bees",]
p7 <- ggplot(species.bees, aes(reorder(species.bees$species, species.bees$normalized.degree),species.bees$normalized.degree,fill=body_size))
p7 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 15) +
  #theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12)) +
  scale_y_continuous(limits = c(0,1)) 

# Trees
species.trees<-node_species[node_species$taxa == "trees",]
p8 <- ggplot(species.trees, aes(reorder(species.trees$species, species.trees$normalized.degree),species.trees$normalized.degree, fill=body_size))
p8 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 3.5) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=1)) +
  scale_y_continuous(limits = c(0,1)) 

# Figure S4 - Histograms showing nestedness contribution at the species level, color-coded by species body mass

# Small mammals
species.sm<-node_species[node_species$taxa == "sm",]
p1 <- ggplot(species.sm, aes(reorder(species.sm$species, species.sm$nested_contribution),species.sm$nested_contribution, fill = log10(body_size)))
p1 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.8) +
  #theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12))

# Lizards
species.liz<-node_species[node_species$taxa == "liz",]
p2 <- ggplot(species.liz, aes(reorder(species.liz$species, species.liz$nested_contribution),species.liz$nested_contribution, fill =  log10(body_size)))
p2 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.8) +
  #  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12)) 

# Mid-large mammals
species.mid.large<-node_species[node_species$taxa == "mid.large",]
p3 <- ggplot(species.mid.large, aes(reorder(species.mid.large$species, species.mid.large$nested_contribution),species.mid.large$nested_contribution, fill = log10(body_size)))
p3 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 0.7) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12))

# Understorey birds
species.birds<-node_species[node_species$taxa == "birds",]
p4 <- ggplot(species.birds, aes(reorder(species.birds$species, species.birds$nested_contribution),species.birds$nested_contribution, fill=log10(body_size)))
p4 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.2) +
  #theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 5),
    axis.text.x=element_text(size=12))

# Frogs
species.frogs<-read.table("node_species_metrics_frogs.txt", header = T) 
p5 <- ggplot(species.frogs, aes(reorder(species.frogs$species, species.frogs$nested_contribution),species.frogs$nested_contribution, fill=log10(body_size_mm)))
p5 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 1.7) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12))

# Dung beetles
species.beetles<-node_species[node_species$taxa == "beetles",]
p6 <- ggplot(species.beetles, aes(reorder(species.beetles$species, species.beetles$nested_contribution),species.beetles$nested_contribution,fill=log10(body_size)))
p6 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = -2) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12)) 

# Orchid bees
species.bees<-node_species[node_species$taxa == "bees",]
p7 <- ggplot(species.bees, aes(reorder(species.bees$species, species.bees$nested_contribution),species.bees$nested_contribution,fill=body_size))
p7 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 15) +
  #theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=12))

# Trees
species.trees<-node_species[node_species$taxa == "trees",]
p8 <- ggplot(species.trees, aes(reorder(species.trees$species, species.trees$nested_contribution),species.trees$nested_contribution, fill=body_size))
p8 + geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  scale_fill_gradient2(high="blue", mid = "yellow", low="brown", midpoint = 3.5) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(#title = element_text(size =15),
    axis.title.y = element_text(size  =0),
    axis.title.x = element_text(size =0),
    axis.text.y=element_text(size = 12),
    axis.text.x=element_text(size=1))


# Figure S5 - Histograms showing nestedness contribution at the species level, color-coded by species body mass

# pannel A
p1 <- ggplot(node_species, aes(x = body_size_st, y = normalized.degree, group = node_species[,1], color = node_species[,1])) 
p1 + geom_point(size = 3, alpha = 0.5) +  
  geom_smooth(data = node_species, aes(x = body_size_st, y = normalized.degree), 
              method = "lm", formula = y ~ x , se = F) +
  #geom_smooth(span = 5, se = F) +
  scale_colour_manual(values=c("gold","blue","orange","light blue","purple","dark red","tomato","yellowgreen")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y="Normalised degree",
       x= "Relativized species size (%)") +
  theme(legend.position="none") +
  theme(plot.title = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=18,color='black', face="bold"),
        axis.title.y = element_text(size=18,color='black', face="bold"), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15)) 

# pannel B
p2 <- ggplot(node_species, aes(x = body_size_st, y = nested_contribution, group = node_species[,1], color = node_species[,1])) 
p2 + geom_point(size = 3, alpha = 0.5) +  
  geom_smooth(data = node_species, aes(x = body_size_st, y = nested_contribution), 
              method = "lm", formula = y ~ x , se = F) +
  #  geom_smooth(span = 0.5, se = F) +
  scale_colour_manual(values=c("gold","blue","orange","light blue","purple","dark red","tomato","yellowgreen")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y="Nestedness contribution",
       x= "Relativized species size (%)") +
  # scale_y_continuous(limits = c(-0.1,1.001)) +
  theme(legend.position="none") +
  theme(plot.title = element_text(size=15, face="bold"),
        axis.title.x = element_text(size=18,color='black', face="bold"),
        axis.title.y = element_text(size=18,color='black', face="bold"), 
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15))


### Lines 1-10:Opening phylogenetic packages that will be needed for analyses####
library(phytools)
library(ape)
library(geiger)
library(nlme)
library(picante)
library(caper)
library(Rmpfr)
library(OUwie)
library(qpcR)

###Lines 14-39: Loading functions to be able to open phylogenetic trees####

### function to be able to read trees
read.newick<-function(file="",text){
  # check to see if reading from file
  if(file!="") text<-scan(file,sep="\n",what="character")
  if(length(text)>1){
    tree<-lapply(text,newick)
    class(tree)<-"multiPhylo"
  } else tree<-newick(text)
  return(tree)
}

# main Newick string function
# written by Liam J. Revell 2013
newick<-function(text){
  text<-unlist(strsplit(text, NULL))
  tip.label<-vector(mode="character")
  node.label<-vector(mode="character") 
  edge<-matrix(c(1,NA),1,2) 
  edge.length<-vector()
  currnode<-1
  Nnode<-currnode
  i<-j<-k<-1
  while(text[i]!=";"){
    if(text[i]=="("){
      if(j>nrow(edge)) edge<-rbind(edge,c(NA,NA))
      edge[j,1]<-currnode
      i<-i+1
      # is the next element a label?
      if(is.na(match(text[i],c("(",")",",",":",";")))){
        temp<-getLabel(text,i)
        tip.label[k]<-temp$label
        i<-temp$end
        edge[j,2]<--k
        k<-k+1
        # is there a branch length?
        if(text[i]==":"){
          temp<-getEdgeLength(text,i)
          edge.length[j]<-temp$edge.length
          i<-temp$end
        }	
      } else if(text[i]=="("){
        Nnode<-Nnode+1 # creating a new internal node
        currnode<-Nnode
        edge[j,2]<-currnode # move to new internal node
      }
      j<-j+1
    } else if(text[i]==")"){
      i<-i+1
      # is the next element a label?
      if(is.na(match(text[i],c("(",")",",",":",";")))){
        temp<-getLabel(text,i)
        node.label[currnode]<-temp$label
        i<-temp$end
      }
      # is there a branch length?
      if(text[i]==":"){
        temp<-getEdgeLength(text,i)
        if(currnode>1){ 
          ii<-match(currnode,edge[,2])
          edge.length[ii]<-temp$edge.length
        } else root.edge<-temp$edge.length
        i<-temp$end
      }	
      if(currnode>1) currnode<-edge[match(currnode,edge[,2]),1] # move down the tree
    } else if(text[i]==","){
      if(j>nrow(edge)) edge<-rbind(edge,c(NA,NA))
      edge[j,1]<-currnode
      i<-i+1
      # is the next element a label?
      if(is.na(match(text[i],c("(",")",",",":",";")))){
        temp<-getLabel(text,i)
        tip.label[k]<-temp$label
        i<-temp$end
        edge[j,2]<--k
        k<-k+1
        # is there a branch length?
        if(text[i]==":"){
          temp<-getEdgeLength(text,i)
          edge.length[j]<-temp$edge.length
          i<-temp$end
        }
      } else if(text[i]=="("){
        Nnode<-Nnode+1 # creating a new internal node
        currnode<-Nnode
        edge[j,2]<-currnode # move to internal node
      }
      j<-j+1
    }
  }
  Ntip<-k-1
  edge[edge>0]<-edge[edge>0]+Ntip
  edge[edge<0]<--edge[edge<0]
  edge.length[is.na(edge.length)]<-0
  node.label[is.na(node.label)]<-""
  if(length(node.label)==0) node.label<-NULL
  # assemble into "phylo" object
  tree<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,edge.length=edge.length,node.label=node.label)
  class(tree)<-"phylo"
  return(tree)
}

# function gets label
# written by Liam J. Revell 2011-2013
getLabel<-function(text,start,stop.char=c(",",":",")")){
  i<-0
  label<-vector()
  while(is.na(match(text[i+start],stop.char))){
    label[i+1]<-text[i+start]
    i<-i+1
  }
  return(list(label=paste(label,collapse=""),end=i+start))
}

# function gets branch length
# written by Liam J. Revell 2011-2013
getEdgeLength<-function(text,start){
  i<-start+1; m<-1
  temp<-vector()
  while(is.na(match(text[i],c(",",")",";")))){
    temp[m]<-text[i]
    i<-i+1
    m<-m+1
  }
  return(list(edge.length=as.numeric(paste(temp,collapse="")),end=i))
}

####Lines 141-382: OUwie and model selection for wrasses####

###
### LOAD TREES
### 

#labridae
##Load tree
labridtree<-read.newick(file="Baliga_MCC.nwk")
labridtree<-collapse.singles(labridtree)

plot(labridtree)
is.rooted(labridtree) ## Yes, tree is rooted
labridtree$tip.label
is.ultrametric(labridtree)
## TRUE


#full labrid data
labrid_data<-read.csv("wrassedata.csv", header=1)
head(labrid_data)
rownames(labrid_data)<-labrid_data$Species

### PRUNE THE TREE to just include your species data
labridspecies<-row.names(labrid_data)
labridspecies ### 81 species
pruned_labrid<-drop.tip(labridtree, setdiff(labridtree$tip.label, labridspecies))
setdiff(labridtree$tip.label, labridspecies)
plotTree(pruned_labrid)
pruned_labrid$tip.label

### Order Data
ordered_labrid_data <-labrid_data[pruned_labrid$tip.label,] # this puts your data into the same order as the tips on the tree
ordered_labrid_data

head(ordered_labrid_data)

#### Dataset
names(ordered_labrid_data)

loginput<-ordered_labrid_data[,8]
names(loginput)<-row.names(ordered_labrid_data)

logoutput<-ordered_labrid_data[,9]
names(logoutput)<-row.names(ordered_labrid_data)

logcoupler<-ordered_labrid_data[,10]
names(logcoupler)<-row.names(ordered_labrid_data)

logkt<- ordered_labrid_data[,11]
names(logkt)<-row.names(ordered_labrid_data)

diet<-ordered_labrid_data[,6]
names(diet)<-row.names(ordered_labrid_data)

type<-ordered_labrid_data[,7] ### identify discrete character for regime differences
names(type)<-row.names(ordered_labrid_data)

type


## make simmap
#### Running OUwie

##### You need to make a diet type matrix for all analyses
type <- as.matrix(ordered_labrid_data)[, "type"]
type ## should list your species, and whether they are shallow or deep

### need to make stochastic character maps of habitat type using the make.simmap function
trees_type<-make.simmap(pruned_labrid,type, nsim=500) ### Make 500 stochastic character maps for diet category
trees_type ### should say "500 phylogenetic trees with mapped discrete characters"



#### make separate data frames for all your traits

### 1. input
loginput_OUwie <- data.frame(Genus_species = rownames(ordered_labrid_data), Reg = type, X = as.numeric(loginput))
loginput_OUwie

##### 2. output
logoutput_OUwie <- data.frame(Genus_species = rownames(ordered_labrid_data), Reg = type, X = as.numeric(logoutput))
logoutput_OUwie


##### 3. coupler
logcoupler_OUwie <- data.frame(Genus_species = rownames(ordered_labrid_data), Reg = type, X = as.numeric(logcoupler))
logcoupler_OUwie


#### kt
logkt_OUwie <- data.frame(Genus_species = rownames(ordered_labrid_data), Reg = type, X = as.numeric(logkt))
logkt_OUwie




#### Run your models

#input
BM1_loginput<-lapply(trees_type,OUwie,data= loginput_OUwie,model="BM1",simmap.tree=TRUE)
BMS_loginput<-lapply(trees_type,OUwie,data= loginput_OUwie,model="BMS",simmap.tree=TRUE)
OU1_loginput<-lapply(trees_type,OUwie,data= loginput_OUwie,model="OU1",simmap.tree=TRUE)
OUM_loginput<-lapply(trees_type,OUwie,data= loginput_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_loginput<-lapply(trees_type,OUwie,data= loginput_OUwie,model="OUMV",simmap.tree=TRUE)



#output
BM1_logoutput<-lapply(trees_type,OUwie,data= logoutput_OUwie,model="BM1",simmap.tree=TRUE)
BMS_logoutput<-lapply(trees_type,OUwie,data= logoutput_OUwie,model="BMS",simmap.tree=TRUE)
OU1_logoutput<-lapply(trees_type,OUwie,data= logoutput_OUwie,model="OU1",simmap.tree=TRUE)
OUM_logoutput<-lapply(trees_type,OUwie,data= logoutput_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_logoutput<-lapply(trees_type,OUwie,data= logoutput_OUwie,model="OUMV",simmap.tree=TRUE)



#coupler
BM1_logcoupler<-lapply(trees_type,OUwie,data= logcoupler_OUwie,model="BM1",simmap.tree=TRUE)
BMS_logcoupler<-lapply(trees_type,OUwie,data= logcoupler_OUwie,model="BMS",simmap.tree=TRUE)
OU1_logcoupler<-lapply(trees_type,OUwie,data= logcoupler_OUwie,model="OU1",simmap.tree=TRUE)
OUM_logcoupler<-lapply(trees_type,OUwie,data= logcoupler_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_logcoupler<-lapply(trees_type,OUwie,data= logcoupler_OUwie,model="OUMV",simmap.tree=TRUE)


#kt
BM1_logkt<-lapply(trees_type,OUwie,data= logkt_OUwie,model="BM1",simmap.tree=TRUE)
BMS_logkt<-lapply(trees_type,OUwie,data= logkt_OUwie,model="BMS",simmap.tree=TRUE)
OU1_logkt<-lapply(trees_type,OUwie,data= logkt_OUwie,model="OU1",simmap.tree=TRUE)
OUM_logkt<-lapply(trees_type,OUwie,data= logkt_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_logkt<-lapply(trees_type,OUwie,data= logkt_OUwie,model="OUMV",simmap.tree=TRUE)



#### Get mean AIC scores for each model

#kt
meanaicc_BM1_logkt<-t(sapply(BM1_logkt,function(x) x$AICc))
meanaicc_BM1_logkt
mean(meanaicc_BM1_logkt) # -138.2987

meanaicc_BMS_logkt<-t(sapply(BMS_logkt,function(x) x$AICc))
meanaicc_BMS_logkt
mean(meanaicc_BMS_logkt) # -136.9242

meanaicc_OU1_logkt<-t(sapply(OU1_logkt,function(x) x$AICc))
meanaicc_OU1_logkt
mean(meanaicc_OU1_logkt) # -147.9847

meanaicc_OUM_logkt<-t(sapply(OUM_logkt,function(x) x$AICc))
meanaicc_OUM_logkt
mean(meanaicc_OUM_logkt) # -151.4389

meanaicc_OUMV_logkt<-t(sapply(OUMV_logkt,function(x) x$AICc))
meanaicc_OUMV_logkt
mean(meanaicc_OUMV_logkt) # -151.1844

#input
meanaicc_BM1_loginput<-t(sapply(BM1_loginput,function(x) x$AICc))
meanaicc_BM1_loginput
mean(meanaicc_BM1_loginput) # -212.9015

meanaicc_BMS_loginput<-t(sapply(BMS_loginput,function(x) x$AICc))
meanaicc_BMS_loginput
mean(meanaicc_BMS_loginput) # -211.7291

meanaicc_OU1_loginput<-t(sapply(OU1_loginput,function(x) x$AICc))
meanaicc_OU1_loginput
mean(meanaicc_OU1_loginput) # -224.4719

meanaicc_OUM_loginput<-t(sapply(OUM_loginput,function(x) x$AICc))
meanaicc_OUM_loginput
mean(meanaicc_OUM_loginput) # -226.2652

meanaicc_OUMV_loginput<-t(sapply(OUMV_loginput,function(x) x$AICc))
meanaicc_OUMV_loginput
mean(meanaicc_OUMV_loginput) # -226.1387


#output
meanaicc_BM1_logoutput<-t(sapply(BM1_logoutput,function(x) x$AICc))
meanaicc_BM1_logoutput
mean(meanaicc_BM1_logoutput) # -246.7342

meanaicc_BMS_logoutput<-t(sapply(BMS_logoutput,function(x) x$AICc))
meanaicc_BMS_logoutput
mean(meanaicc_BMS_logoutput) # -245.1012

meanaicc_OU1_logoutput<-t(sapply(OU1_logoutput,function(x) x$AICc))
meanaicc_OU1_logoutput
mean(meanaicc_OU1_logoutput) # -253.443

meanaicc_OUM_logoutput<-t(sapply(OUM_logoutput,function(x) x$AICc))
meanaicc_OUM_logoutput
mean(meanaicc_OUM_logoutput) # -254.8644

meanaicc_OUMV_logoutput<-t(sapply(OUMV_logoutput,function(x) x$AICc))
meanaicc_OUMV_logoutput
mean(meanaicc_OUMV_logoutput) # -253.7929


#coupler
meanaicc_BM1_logcoupler<-t(sapply(BM1_logcoupler,function(x) x$AICc))
meanaicc_BM1_logcoupler
mean(meanaicc_BM1_logcoupler) # -261.8418

meanaicc_BMS_logcoupler<-t(sapply(BMS_logcoupler,function(x) x$AICc))
meanaicc_BMS_logcoupler
mean(meanaicc_BMS_logcoupler) # -261.6427

meanaicc_OU1_logcoupler<-t(sapply(OU1_logcoupler,function(x) x$AICc))
meanaicc_OU1_logcoupler
mean(meanaicc_OU1_logcoupler) # -259.6839

meanaicc_OUM_logcoupler<-t(sapply(OUM_logcoupler,function(x) x$AICc))
meanaicc_OUM_logcoupler
mean(meanaicc_OUM_logcoupler) # -260.5624

meanaicc_OUMV_logcoupler<-t(sapply(OUMV_logcoupler,function(x) x$AICc))
meanaicc_OUMV_logcoupler
mean(meanaicc_OUMV_logcoupler) # -259.551


### Calculate AIC weights for your AIC scores.
### paste your AIC values BM1, BMS, OU1, OUM, OUMV <- in that order
#### IF THERE IS A NEGATIVE SIGN MAKE SURE TO INCLUDE THAT. Be very, very careful.

library(qpcR)
logkt_aic<-c(-138.2987,-136.9242,-147.9847,-151.4389,-151.1844)

loginput_aic<- c(-212.9015,-211.7291,-224.4719,-226.2652,-226.1387)

logoutput_aic<- c(-246.7342,-245.1012,-253.443,-254.8644,-253.7929)

logcoupler_aic<- c(-261.8418,-261.6427,-259.6839,-260.5624,-259.551) 

akaike.weights(logkt_aic) ### OUM but OUMV works well too...
akaike.weights(loginput_aic) #OUM fits best but OUMV works well too as does OU1...
akaike.weights(logoutput_aic) ##OUM fits best but OUMV and OU1 fit well too...
akaike.weights(logcoupler_aic) #BM1 but all models well supported???
##### This will give you delta AIC (0=lowest AIC; best-fit model) and weights (all weights add up to 1)


#####Lines 383-616: Ouwie and model selection for stomatopods#####

###
### LOAD TREES
### 

#stomatopoda
##Load tree
stomatopodtree<-read.newick(file="stomatopod.nwk")
stomatopodtree<-collapse.singles(stomatopodtree)
plot(stomatopodtree)
is.rooted(stomatopodtree)
## Yes, tree is rooted
stomatopodtree$tip.label
is.ultrametric(stomatopodtree)
## TRUE


#full stomatopod data
stomatopod_data<-read.csv("stomatopodata.csv", header=1)
head(stomatopod_data)
rownames(stomatopod_data)<-stomatopod_data$Species

### PRUNE THE TREE to just include your species data
stomatopodspecies<-row.names(stomatopod_data)
stomatopodspecies ### 34 species
pruned_stomatopod<-drop.tip(stomatopodtree, setdiff(stomatopodtree$tip.label, stomatopodspecies))
setdiff(pruned_stomatopod$tip.label, stomatopodspecies)
plotTree(pruned_stomatopod)
pruned_stomatopod$tip.label

### Order Data
ordered_stomatopod_data <-stomatopod_data[pruned_stomatopod$tip.label,] # this puts your data into the same order as the tips on the tree
ordered_stomatopod_data

head(ordered_stomatopod_data)

#### Dataset
names(ordered_stomatopod_data)

loginput_ms<-ordered_stomatopod_data[,6]
names(loginput_ms)<-row.names(ordered_stomatopod_data)

logoutput_ms<-ordered_stomatopod_data[,7]
names(logoutput_ms)<-row.names(ordered_stomatopod_data)

logcoupler_ms<-ordered_stomatopod_data[,8]
names(logcoupler_ms)<-row.names(ordered_stomatopod_data)

logkt_ms<- ordered_stomatopod_data[,9]
names(logkt_ms)<-row.names(ordered_stomatopod_data)

type_ms<-ordered_stomatopod_data[,10] ### identify the column where your data are
names(type)<-row.names(ordered_stomatopod_data)

type_ms

## make simmap
#### Running OUwie

##### You need to make a habitat type matrix for all analyses
type_ms <- as.matrix(ordered_stomatopod_data)[, "type"]
type_ms ## should list your species, and whether they are shallow or deep

### need to make stochastic character maps of habitat type using the make.simmap function
trees_type_ms<-make.simmap(pruned_stomatopod,type_ms, nsim=500) ### Make 500 stochastic character maps for habitat depth
trees_type_ms ### should say "500 phylogenetic trees with mapped discrete characters"


#### make separate data frames for all your traits

### 1. input
loginput_ms_OUwie <- data.frame(Genus_species = rownames(ordered_stomatopod_data), Reg = type_ms, X = as.numeric(loginput_ms))
loginput_ms_OUwie

##### 2. output
logoutput_ms_OUwie <- data.frame(Genus_species = rownames(ordered_stomatopod_data), Reg = type_ms, X = as.numeric(logoutput_ms))
logoutput_ms_OUwie


##### 3. coupler
logcoupler_ms_OUwie <- data.frame(Genus_species = rownames(ordered_stomatopod_data), Reg = type_ms, X = as.numeric(logcoupler_ms))
logcoupler_ms_OUwie


#### kt
logkt_ms_OUwie <- data.frame(Genus_species = rownames(ordered_stomatopod_data), Reg = type_ms, X = as.numeric(logkt_ms))
logkt_ms_OUwie




#### Run your models

#input
BM1_loginput_ms<-lapply(trees_type_ms,OUwie,data= loginput_ms_OUwie,model="BM1",simmap.tree=TRUE)
BMS_loginput_ms<-lapply(trees_type_ms,OUwie,data= loginput_ms_OUwie,model="BMS",simmap.tree=TRUE)
OU1_loginput_ms<-lapply(trees_type_ms,OUwie,data= loginput_ms_OUwie,model="OU1",simmap.tree=TRUE)
OUM_loginput_ms<-lapply(trees_type_ms,OUwie,data= loginput_ms_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_loginput_ms<-lapply(trees_type_ms,OUwie,data= loginput_ms_OUwie,model="OUMV",simmap.tree=TRUE)

#output
BM1_logoutput_ms<-lapply(trees_type_ms,OUwie,data= logoutput_ms_OUwie,model="BM1",simmap.tree=TRUE)
BMS_logoutput_ms<-lapply(trees_type_ms,OUwie,data= logoutput_ms_OUwie,model="BMS",simmap.tree=TRUE)
OU1_logoutput_ms<-lapply(trees_type_ms,OUwie,data= logoutput_ms_OUwie,model="OU1",simmap.tree=TRUE)
OUM_logoutput_ms<-lapply(trees_type_ms,OUwie,data= logoutput_ms_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_logoutput_ms<-lapply(trees_type_ms,OUwie,data= logoutput_ms_OUwie,model="OUMV",simmap.tree=TRUE)


#coupler
BM1_logcoupler_ms<-lapply(trees_type_ms,OUwie,data= logcoupler_ms_OUwie,model="BM1",simmap.tree=TRUE)
BMS_logcoupler_ms<-lapply(trees_type_ms,OUwie,data= logcoupler_ms_OUwie,model="BMS",simmap.tree=TRUE)
OU1_logcoupler_ms<-lapply(trees_type_ms,OUwie,data= logcoupler_ms_OUwie,model="OU1",simmap.tree=TRUE)
OUM_logcoupler_ms<-lapply(trees_type_ms,OUwie,data= logcoupler_ms_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_logcoupler_ms<-lapply(trees_type_ms,OUwie,data= logcoupler_ms_OUwie,model="OUMV",simmap.tree=TRUE)

#kt
BM1_logkt_ms<-lapply(trees_type_ms,OUwie,data= logkt_ms_OUwie,model="BM1",simmap.tree=TRUE)
BMS_logkt_ms<-lapply(trees_type_ms,OUwie,data= logkt_ms_OUwie,model="BMS",simmap.tree=TRUE)
OU1_logkt_ms<-lapply(trees_type_ms,OUwie,data= logkt_ms_OUwie,model="OU1",simmap.tree=TRUE)
OUM_logkt_ms<-lapply(trees_type_ms,OUwie,data= logkt_ms_OUwie,model="OUM",simmap.tree=TRUE)
OUMV_logkt_ms<-lapply(trees_type_ms,OUwie,data= logkt_ms_OUwie,model="OUMV",simmap.tree=TRUE)


#### Get mean AIC scores for each model

#kt
meanaicc_BM1_logkt_ms<-t(sapply(BM1_logkt_ms,function(x) x$AICc))
meanaicc_BM1_logkt_ms
mean(meanaicc_BM1_logkt_ms) # -82.21707

meanaicc_BMS_logkt_ms<-t(sapply(BMS_logkt_ms,function(x) x$AICc))
meanaicc_BMS_logkt_ms
mean(meanaicc_BMS_logkt_ms) # -80.75224

meanaicc_OU1_logkt_ms<-t(sapply(OU1_logkt_ms,function(x) x$AICc))
meanaicc_OU1_logkt_ms
mean(meanaicc_OU1_logkt_ms) # -80.10237

meanaicc_OUM_logkt_ms<-t(sapply(OUM_logkt_ms,function(x) x$AICc))
meanaicc_OUM_logkt_ms
mean(meanaicc_OUM_logkt_ms) # -86.1127

meanaicc_OUMV_logkt_ms<-t(sapply(OUMV_logkt_ms,function(x) x$AICc)) 
meanaicc_OUMV_logkt_ms
mean(meanaicc_OUMV_logkt_ms) # -88.70393


#input
meanaicc_BM1_loginput_ms<-t(sapply(BM1_loginput_ms,function(x) x$AICc))
meanaicc_BM1_loginput_ms
mean(meanaicc_BM1_loginput_ms) # -121.4827

meanaicc_BMS_loginput_ms<-t(sapply(BMS_loginput_ms,function(x) x$AICc))
meanaicc_BMS_loginput_ms
mean(meanaicc_BMS_loginput_ms) # -119.7558

meanaicc_OU1_loginput_ms<-t(sapply(OU1_loginput_ms,function(x) x$AICc))
meanaicc_OU1_loginput_ms
mean(meanaicc_OU1_loginput_ms) # -121.4161

meanaicc_OUM_loginput_ms<-t(sapply(OUM_loginput_ms,function(x) x$AICc))
meanaicc_OUM_loginput_ms
mean(meanaicc_OUM_loginput_ms) # -125.5831

meanaicc_OUMV_loginput_ms<-t(sapply(OUMV_loginput_ms,function(x) x$AICc))
meanaicc_OUMV_loginput_ms
mean(meanaicc_OUMV_loginput_ms) # -122.8416


#output
meanaicc_BM1_logoutput_ms<-t(sapply(BM1_logoutput_ms,function(x) x$AICc))
meanaicc_BM1_logoutput_ms
mean(meanaicc_BM1_logoutput_ms) # -72.82617

meanaicc_BMS_logoutput_ms<-t(sapply(BMS_logoutput_ms,function(x) x$AICc))
meanaicc_BMS_logoutput_ms
mean(meanaicc_BMS_logoutput_ms) # -71.81533

meanaicc_OU1_logoutput_ms<-t(sapply(OU1_logoutput_ms,function(x) x$AICc))
meanaicc_OU1_logoutput_ms
mean(meanaicc_OU1_logoutput_ms) # -76.15634

meanaicc_OUM_logoutput_ms<-t(sapply(OUM_logoutput_ms,function(x) x$AICc))
meanaicc_OUM_logoutput_ms
mean(meanaicc_OUM_logoutput_ms) # -81.02534

meanaicc_OUMV_logoutput_ms<-t(sapply(OUMV_logoutput_ms,function(x) x$AICc))
meanaicc_OUMV_logoutput_ms
mean(meanaicc_OUMV_logoutput_ms) # -81.69081



#coupler
meanaicc_BM1_logcoupler_ms<-t(sapply(BM1_logcoupler_ms,function(x) x$AICc))
meanaicc_BM1_logcoupler_ms
mean(meanaicc_BM1_logcoupler_ms) # -142.4579

meanaicc_BMS_logcoupler_ms<-t(sapply(BMS_logcoupler_ms,function(x) x$AICc))
meanaicc_BMS_logcoupler_ms
mean(meanaicc_BMS_logcoupler_ms) # -140.0665

meanaicc_OU1_logcoupler_ms<-t(sapply(OU1_logcoupler_ms,function(x) x$AICc))
meanaicc_OU1_logcoupler_ms
mean(meanaicc_OU1_logcoupler_ms) # -144.5267

meanaicc_OUM_logcoupler_ms<-t(sapply(OUM_logcoupler_ms,function(x) x$AICc))
meanaicc_OUM_logcoupler_ms
mean(meanaicc_OUM_logcoupler_ms) # -142.038

meanaicc_OUMV_logcoupler_ms<-t(sapply(OUMV_logcoupler_ms,function(x) x$AICc))
meanaicc_OUMV_logcoupler_ms
mean(meanaicc_OUMV_logcoupler_ms) # -143.6003



### Calculate AIC weights for your AIC scores.
### paste your AIC values BM1, BMS, OU1, OUM, OUMV <- in that order
#### IF THERE IS A NEGATIVE SIGN MAKE SURE TO INCLUDE THAT. Be very, very careful.
library(qpcR)
logkt_ms_aic<-c(-82.21707,-80.75224,-80.10237,-86.1127,-88.70393)

loginput_ms_aic<- c(-121.4827,-119.7558,-121.4161,-125.5831,-122.8416)

logoutput_ms_aic<- c(-72.82617,-71.81533,-76.15634,-81.02534,-81.69081) 

logcoupler_ms_aic<- c(-142.4579,-140.0665,-144.5267,-142.038,-143.6003) 

akaike.weights(logkt_ms_aic) #OUMV fits best
akaike.weights(loginput_ms_aic) #OUM fits best
akaike.weights(logoutput_ms_aic) #OUMV fits best
akaike.weights(logcoupler_ms_aic) #OU1 fits best


#Lines 618-643: Creating functions of models/parameters for best fit models for wrasses#####

#kt
solution_OUMV_logkt<-t(sapply(OUMV_logkt,function(x) x$solution))
theta_OUMV_logkt<-t(sapply(OUMV_logkt,function(x) x$theta))

solution_OUM_logkt<-t(sapply(OUM_logkt,function(x) x$solution)) #best
theta_OUM_logkt<-t(sapply(OUM_logkt,function(x) x$theta))

#input
solution_OUMV_loginput<-t(sapply(OUMV_loginput,function(x) x$solution))
theta_OUMV_loginput<-t(sapply(OUMV_loginput,function(x) x$theta))

solution_OUM_loginput<-t(sapply(OUM_loginput,function(x) x$solution)) #best
theta_OUM_loginput<-t(sapply(OUM_loginput,function(x) x$theta))

#output
solution_OUMV_logoutput<-t(sapply(OUMV_logoutput,function(x) x$solution))
theta_OUMV_logoutput<-t(sapply(OUMV_logoutput,function(x) x$theta))

solution_OUM_logoutput<-t(sapply(OUM_logoutput,function(x) x$solution)) #best
theta_OUM_logoutput<-t(sapply(OUM_logoutput,function(x) x$theta))

#coupler
solution_BM1_logcoupler<-t(sapply(BM1_logcoupler,function(x) x$solution))
theta_BM1_logcoupler<-t(sapply(BM1_logcoupler,function(x) x$theta))


# Lines 646-662: Creating functions of models/parameters for best fit models in stomatopods ####
solution_OUMV_logkt_ms<-t(sapply(OUMV_logkt_ms,function(x) x$solution))
theta_OUMV_logkt_ms<-t(sapply(OUMV_logkt_ms,function(x) x$theta))

solution_OUM_loginput_ms<-t(sapply(OUM_loginput_ms,function(x) x$solution))
theta_OUM_loginput_ms<-t(sapply(OUM_loginput_ms,function(x) x$theta))

solution_OUMV_logoutput_ms<-t(sapply(OUMV_logoutput_ms,function(x) x$solution))
theta_OUMV_logoutput_ms<-t(sapply(OUMV_logoutput_ms,function(x) x$theta))

solution_OUMV_logcoupler_ms<-t(sapply(OUMV_logcoupler_ms,function(x) x$solution))
theta_OUMV_logcoupler_ms<-t(sapply(OUMV_logcoupler_ms,function(x) x$theta))

solution_OU1_logcoupler_ms<-t(sapply(OU1_logcoupler_ms,function(x) x$solution)) #best
theta_OU1_logcoupler_ms<-t(sapply(OU1_logcoupler_ms,function(x) x$theta)) #best



#Lines 664-957: Using simulated data to see which models we can fit for wrasses#####


regime_labrid <- as.matrix(ordered_labrid_data)[, "type"]
regime_labrid ## should list your species, and whether they are shallow or deep

trees_regime<-make.simmap(pruned_labrid,regime_labrid, nsim=100) 
trees_regime ### should say "100 phylogenetic trees with mapped discrete characters"

class(trees_regime)


#Simulate an Ornstein-Uhlenbeck model with different state means
#and sigma^2 per selective regime, but same alpha
alpha_bm1=c(1e-10, 1e-10)
sigma.sq_bm1=c(0.45,0.45)
theta0=1.0
theta_bm1=c(0,0)

alpha_bms=c(1e-10, 1e-10)
sigma.sq_bms=c(0.45,0.9)
theta0=1.0
theta_bms=c(0,0)

alpha_ou1=c(0.1,0.1)
sigma.sq_ou1=c(0.9,0.9)
theta0=1.0
theta_ou1=c(1.0,1.0)

alpha_oum=c(1.0, 1.0)
sigma.sq_oum=c(0.45,0.45)
theta0=1.0
theta_oum=c(1.0,2.0)

alpha_oumv=c(1.0, 1.0)
sigma.sq_oumv=c(0.45,0.9)
theta0=1.0
theta_oumv=c(1.0,2.0)

alpha_oumva=c(1.0, .5)
sigma.sq_oumva=c(0.45,0.9)
theta0=1.0
theta_oumva=c(1.0,2.0)


#making simulations
#BM1
ktdata_sim_bm1<- OUwie.sim(trees_regime[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_bm1, sigma.sq=sigma.sq_bm1, theta0=theta0, theta=theta_bm1)

#BMS
ktdata_sim_bms<- OUwie.sim(trees_regime[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_bms, sigma.sq=sigma.sq_bms, theta0=theta0, theta=theta_bms)

#OU1
ktdata_sim_ou1<- OUwie.sim(trees_regime[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_ou1, sigma.sq=sigma.sq_ou1, theta0=theta0, theta=theta_ou1)

#OUM
ktdata_sim_oum<- OUwie.sim(trees_regime[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_oum, sigma.sq=sigma.sq_oum, theta0=theta0, theta=theta_oum)

#OUMV
ktdata_sim_oumv<- OUwie.sim(trees_regime[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_oumv, sigma.sq=sigma.sq_oumv, theta0=theta0, theta=theta_oumv)

#OUMVA
ktdata_sim_oumva<- OUwie.sim(trees_regime[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_oumva, sigma.sq=sigma.sq_oumva, theta0=theta0, theta=theta_oumva)


#making new data frames with simulated data
#bm1
simkt_bm1<- data.frame(ordered_labrid_data$Species,ordered_labrid_data$type, ktdata_sim_bm1$X)

#bms
simkt_bms<- data.frame(ordered_labrid_data$Species,ordered_labrid_data$type, ktdata_sim_bms$X)

#ou1
simkt_ou1<- data.frame(ordered_labrid_data$Species,ordered_labrid_data$type, ktdata_sim_ou1$X)

#oum
simkt_oum<- data.frame(ordered_labrid_data$Species,ordered_labrid_data$type, ktdata_sim_oum$X)


#oumv
simkt_oumv<- data.frame(ordered_labrid_data$Species,ordered_labrid_data$type, ktdata_sim_oumv$X)

#oumva
simkt_oumva<- data.frame(ordered_labrid_data$Species,ordered_labrid_data$type, ktdata_sim_oumva$X)


#running OUwie with simulated data
#bm1 
bm1_BM1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bm1,model=c("BM1"),simmap.tree=TRUE))
bm1_BMS_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bm1,model=c("BMS"),simmap.tree=TRUE))
bm1_OU1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bm1,model=c("OU1"),simmap.tree=TRUE))
bm1_OUM_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bm1,model=c("OUM"),simmap.tree=TRUE))
bm1_OUMV_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bm1,model=c("OUMV"),simmap.tree=TRUE))
bm1_OUMVA_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bm1,model=c("OUMVA"),simmap.tree=TRUE))

#bms
bms_BM1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bms,model=c("BM1"),simmap.tree=TRUE))
bms_BMS_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bms,model=c("BMS"),simmap.tree=TRUE))
bms_OU1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bms,model=c("OU1"),simmap.tree=TRUE))
bms_OUM_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bms,model=c("OUM"),simmap.tree=TRUE))
bms_OUMV_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bms,model=c("OUMV"),simmap.tree=TRUE))
bms_OUMVA_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_bms,model=c("OUMVA"),simmap.tree=TRUE))

#ou1
ou1_BM1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("BM1"),simmap.tree=TRUE))
ou1_BMS_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("BMS"),simmap.tree=TRUE))
ou1_OU1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OU1"),simmap.tree=TRUE))
ou1_OUM_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OUM"),simmap.tree=TRUE))
ou1_OUMV_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OUMV"),simmap.tree=TRUE))
ou1_OUMVA_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OUMVA"),simmap.tree=TRUE))

#oum
oum_BM1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oum,model=c("BM1"),simmap.tree=TRUE))
oum_BMS_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oum,model=c("BMS"),simmap.tree=TRUE))
oum_OU1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oum,model=c("OU1"),simmap.tree=TRUE))
oum_OUM_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oum,model=c("OUM"),simmap.tree=TRUE))
oum_OUMV_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oum,model=c("OUMV"),simmap.tree=TRUE))
oum_OUMVA_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oum,model=c("OUMVA"),simmap.tree=TRUE))

#oumv
oumv_BM1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumv,model=c("BM1"),simmap.tree=TRUE))
oumv_BMS_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumv,model=c("BMS"),simmap.tree=TRUE))
oumv_OU1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumv,model=c("OU1"),simmap.tree=TRUE))
oumv_OUM_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumv,model=c("OUM"),simmap.tree=TRUE))
oumv_OUMV_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumv,model=c("OUMV"),simmap.tree=TRUE))
oumv_OUMVA_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumv,model=c("OUMVA"),simmap.tree=TRUE))

#oumva
oumva_BM1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumva,model=c("BM1"),simmap.tree=TRUE))
oumva_BMS_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumva,model=c("BMS"),simmap.tree=TRUE))
oumva_OU1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumva,model=c("OU1"),simmap.tree=TRUE))
oumva_OUM_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumva,model=c("OUM"),simmap.tree=TRUE))
oumva_OUMV_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumva,model=c("OUMV"),simmap.tree=TRUE))
oumva_OUMVA_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_oumva,model=c("OUMVA"),simmap.tree=TRUE))

#aic scores
#bm1
meanaicc_bm1_BM1_sim<-t(sapply(bm1_BM1_sim,function(x) x$AICc))
meanaicc_bm1_BM1_sim
mean(meanaicc_bm1_BM1_sim) #406.7591 #yay

meanaicc_bm1_BMS_sim<-t(sapply(bm1_BMS_sim,function(x) x$AICc))
meanaicc_bm1_BMS_sim
mean(meanaicc_bm1_BMS_sim) #408.822

meanaicc_bm1_OU1_sim<-t(sapply(bm1_OU1_sim,function(x) x$AICc))
meanaicc_bm1_OU1_sim
mean(meanaicc_bm1_OU1_sim) #408.9191

meanaicc_bm1_OUM_sim<-t(sapply(bm1_OUM_sim,function(x) x$AICc))
meanaicc_bm1_OUM_sim
mean(meanaicc_bm1_OUM_sim) #410.5854

meanaicc_bm1_OUMV_sim<-t(sapply(bm1_OUMV_sim,function(x) x$AICc))
meanaicc_bm1_OUMV_sim
mean(meanaicc_bm1_OUMV_sim) #412.7324

meanaicc_bm1_OUMVA_sim<-t(sapply(bm1_OUMVA_sim,function(x) x$AICc))
meanaicc_bm1_OUMVA_sim
mean(meanaicc_bm1_OUMVA_sim) #394.6426 .....

sim_bm1<- c(419.8318,421.3853,420.4748,420.549,422.3636,423.0157)
akaike.weights(sim_bm1)

#bms
meanaicc_bms_BM1_sim<-t(sapply(bms_BM1_sim,function(x) x$AICc))
meanaicc_bms_BM1_sim
mean(meanaicc_bms_BM1_sim) #447.631

meanaicc_bms_BMS_sim<-t(sapply(bms_BMS_sim,function(x) x$AICc))
meanaicc_bms_BMS_sim
mean(meanaicc_bms_BMS_sim) #445.3701 #yay

meanaicc_bms_OU1_sim<-t(sapply(bms_OU1_sim,function(x) x$AICc))
meanaicc_bms_OU1_sim
mean(meanaicc_bms_OU1_sim) #449.6103

meanaicc_bms_OUM_sim<-t(sapply(bms_OUM_sim,function(x) x$AICc))
meanaicc_bms_OUM_sim
mean(meanaicc_bms_OUM_sim) #449.5676

meanaicc_bms_OUMV_sim<-t(sapply(bms_OUMV_sim,function(x) x$AICc))
meanaicc_bms_OUMV_sim
mean(meanaicc_bms_OUMV_sim) #449.5676

meanaicc_bms_OUMVA_sim<-t(sapply(bms_OUMVA_sim,function(x) x$AICc))
meanaicc_bms_OUMVA_sim
mean(meanaicc_bms_OUMVA_sim) #297.6972 but why????

sim_bms<- c(458.6416,454.6759,460.7994,462.9802,459.1082,185.6256)
akaike.weights(sim_bms)

#ou1
meanaicc_ou1_BM1_sim<-t(sapply(ou1_BM1_sim,function(x) x$AICc))
meanaicc_ou1_BM1_sim
mean(meanaicc_ou1_BM1_sim) #389.9667

meanaicc_ou1_BMS_sim<-t(sapply(ou1_BMS_sim,function(x) x$AICc))
meanaicc_ou1_BMS_sim
mean(meanaicc_ou1_BMS_sim) #385.6754 

meanaicc_ou1_OU1_sim<-t(sapply(ou1_OU1_sim,function(x) x$AICc))
meanaicc_ou1_OU1_sim
mean(meanaicc_ou1_OU1_sim) #367.171#yay

meanaicc_ou1_OUM_sim<-t(sapply(ou1_OUM_sim,function(x) x$AICc))
meanaicc_ou1_OUM_sim
mean(meanaicc_ou1_OUM_sim) #369.2455

meanaicc_ou1_OUMV_sim<-t(sapply(ou1_OUMV_sim,function(x) x$AICc))
meanaicc_ou1_OUMV_sim
mean(meanaicc_ou1_OUMV_sim) #365.2558 ####but why....

meanaicc_ou1_OUMVA_sim<-t(sapply(ou1_OUMVA_sim,function(x) x$AICc))
meanaicc_ou1_OUMVA_sim
mean(meanaicc_ou1_OUMVA_sim) #-1320069 nope...

sim_ou1<- c(363.6772,365.7432,343.7365,345.9011,348.1397,-1.524019e+13)
akaike.weights(sim_ou1)

#oum
meanaicc_oum_BM1_sim<-t(sapply(oum_BM1_sim,function(x) x$AICc))
meanaicc_oum_BM1_sim
mean(meanaicc_oum_BM1_sim) #148.4085

meanaicc_oum_BMS_sim<-t(sapply(oum_BMS_sim,function(x) x$AICc))
meanaicc_oum_BMS_sim
mean(meanaicc_oum_BMS_sim) #147.5537

meanaicc_oum_OU1_sim<-t(sapply(oum_OU1_sim,function(x) x$AICc))
meanaicc_oum_OU1_sim
mean(meanaicc_oum_OU1_sim) #124.8072 

meanaicc_oum_OUM_sim<-t(sapply(oum_OUM_sim,function(x) x$AICc))
meanaicc_oum_OUM_sim
mean(meanaicc_oum_OUM_sim) #107.1089 yay

meanaicc_oum_OUMV_sim<-t(sapply(oum_OUMV_sim,function(x) x$AICc))
meanaicc_oum_OUMV_sim
mean(meanaicc_oum_OUMV_sim) #109.1171

meanaicc_oum_OUMVA_sim<-t(sapply(oum_OUMVA_sim,function(x) x$AICc))
meanaicc_oum_OUMVA_sim
mean(meanaicc_oum_OUMVA_sim) #-50108944476 nope...

#oumv
meanaicc_oumv_BM1_sim<-t(sapply(oumv_BM1_sim,function(x) x$AICc))
meanaicc_oumv_BM1_sim
mean(meanaicc_oumv_BM1_sim) #209.2446

meanaicc_oumv_BMS_sim<-t(sapply(oumv_BMS_sim,function(x) x$AICc))
meanaicc_oumv_BMS_sim
mean(meanaicc_oumv_BMS_sim) # 205.5904

meanaicc_oumv_OU1_sim<-t(sapply(oumv_OU1_sim,function(x) x$AICc))
meanaicc_oumv_OU1_sim
mean(meanaicc_oumv_OU1_sim) #161.7155

meanaicc_oumv_OUM_sim<-t(sapply(oumv_OUM_sim,function(x) x$AICc))
meanaicc_oumv_OUM_sim
mean(meanaicc_oumv_OUM_sim) #139.2369

meanaicc_oumv_OUMV_sim<-t(sapply(oumv_OUMV_sim,function(x) x$AICc))
meanaicc_oumv_OUMV_sim
mean(meanaicc_oumv_OUMV_sim) #137.2267 #yay

meanaicc_oumv_OUMVA_sim<-t(sapply(oumv_OUMVA_sim,function(x) x$AICc))
meanaicc_oumv_OUMVA_sim
mean(meanaicc_oumv_OUMVA_sim) # -2.923225e+12 nope...

#oumva
meanaicc_oumva_BM1_sim<-t(sapply(oumva_BM1_sim,function(x) x$AICc))
meanaicc_oumva_BM1_sim
mean(meanaicc_oumva_BM1_sim) #291.0797

meanaicc_oumva_BMS_sim<-t(sapply(oumva_BMS_sim,function(x) x$AICc))
meanaicc_oumva_BMS_sim
mean(meanaicc_oumva_BMS_sim) #286.1763

meanaicc_oumva_OU1_sim<-t(sapply(oumva_OU1_sim,function(x) x$AICc))
meanaicc_oumva_OU1_sim
mean(meanaicc_oumva_OU1_sim) #232.366

meanaicc_oumva_OUM_sim<-t(sapply(oumva_OUM_sim,function(x) x$AICc))
meanaicc_oumva_OUM_sim
mean(meanaicc_oumva_OUM_sim) #228.6443

meanaicc_oumva_OUMV_sim<-t(sapply(oumva_OUMV_sim,function(x) x$AICc))
meanaicc_oumva_OUMV_sim
mean(meanaicc_oumva_OUMV_sim) # 226.2742

meanaicc_oumva_OUMVA_sim<-t(sapply(oumva_OUMVA_sim,function(x) x$AICc))
meanaicc_oumva_OUMVA_sim
mean(meanaicc_oumva_OUMVA_sim) #-50070509555

#Lines 958-end: Using simulated data to see which models we can fit for wrasses####


regime_stomatopod <- as.matrix(ordered_stomatopod_data)[, "type"]
regime_stomatopod ## should list your species, and whether they are shallow or deep

### need to make stochastic character maps of habitat type using the make.simmap function
trees_regime_ms<-make.simmap(pruned_stomatopod,regime_stomatopod, nsim=100) ### Make 500 stochastic character maps for habitat depth
trees_regime_ms ### should say "500 phylogenetic trees with mapped discrete characters"

class(trees_regime_ms)

#making simulated data
#BM1
ktdata_sim_bm1_ms<- OUwie.sim(trees_regime_ms[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_bm1, sigma.sq=sigma.sq_bm1, theta0=theta0, theta=theta_bm1)

#BMS
ktdata_sim_bms_ms<- OUwie.sim(trees_regime_ms[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_bms, sigma.sq=sigma.sq_bms, theta0=theta0, theta=theta_bms)

#OU1
ktdata_sim_ou1_ms<- OUwie.sim(trees_regime_ms[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_ou1, sigma.sq=sigma.sq_ou1, theta0=theta0, theta=theta_ou1)

#OUM
ktdata_sim_oum_ms<- OUwie.sim(trees_regime_ms[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_oum, sigma.sq=sigma.sq_oum, theta0=theta0, theta=theta_oum)

#OUMV
ktdata_sim_oumv_ms<- OUwie.sim(trees_regime_ms[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_oumv, sigma.sq=sigma.sq_oumv, theta0=theta0, theta=theta_oumv)

#OUMVA
ktdata_sim_oumva_ms<- OUwie.sim(trees_regime_ms[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_oumva, sigma.sq=sigma.sq_oumva, theta0=theta0, theta=theta_oumva)


#making new data frames with simulated data
#bm1
ktdata_sim_bm1_ms<- data.frame(ordered_stomatopod_data$Species,ordered_stomatopod_data$type, ktdata_sim_bm1_ms$X)
ktdata_sim_bms_ms<- data.frame(ordered_stomatopod_data$Species,ordered_stomatopod_data$type, ktdata_sim_bms_ms$X)
ktdata_sim_ou1_ms<- data.frame(ordered_stomatopod_data$Species,ordered_stomatopod_data$type, ktdata_sim_ou1_ms$X)
ktdata_sim_oum_ms<- data.frame(ordered_stomatopod_data$Species,ordered_stomatopod_data$type, ktdata_sim_oum_ms$X)
ktdata_sim_oumv_ms<- data.frame(ordered_stomatopod_data$Species,ordered_stomatopod_data$type, ktdata_sim_oumv_ms$X)
ktdata_sim_oumva_ms<- data.frame(ordered_stomatopod_data$Species,ordered_stomatopod_data$type, ktdata_sim_oumva_ms$X)

#running OUwie with simulated data
#KT 
#bm1 
bm1_BM1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bm1_ms,model=c("BM1"),simmap.tree=TRUE))
bm1_BMS_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bm1_ms,model=c("BMS"),simmap.tree=TRUE))
bm1_OU1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bm1_ms,model=c("OU1"),simmap.tree=TRUE))
bm1_OUM_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bm1_ms,model=c("OUM"),simmap.tree=TRUE))
bm1_OUMV_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bm1_ms,model=c("OUMV"),simmap.tree=TRUE))
bm1_OUMVA_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bm1_ms,model=c("OUMVA"),simmap.tree=TRUE))

#bms
bms_BM1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bms_ms,model=c("BM1"),simmap.tree=TRUE))
bms_BMS_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bms_ms,model=c("BMS"),simmap.tree=TRUE))
bms_OU1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bms_ms,model=c("OU1"),simmap.tree=TRUE))
bms_OUM_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bms_ms,model=c("OUM"),simmap.tree=TRUE))
bms_OUMV_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bms_ms,model=c("OUMV"),simmap.tree=TRUE))
bms_OUMVA_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_bms_ms,model=c("OUMVA"),simmap.tree=TRUE))

#ou1
ou1_BM1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_ou1_ms,model=c("BM1"),simmap.tree=TRUE))
ou1_BMS_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_ou1_ms,model=c("BMS"),simmap.tree=TRUE))
ou1_OU1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_ou1_ms,model=c("OU1"),simmap.tree=TRUE))
ou1_OUM_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_ou1_ms,model=c("OUM"),simmap.tree=TRUE))
ou1_OUMV_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_ou1_ms,model=c("OUMV"),simmap.tree=TRUE))
ou1_OUMVA_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_ou1_ms,model=c("OUMVA"),simmap.tree=TRUE))

#oum
oum_BM1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("BM1"),simmap.tree=TRUE))
oum_BMS_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("BMS"),simmap.tree=TRUE))
oum_OU1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OU1"),simmap.tree=TRUE))
oum_OUM_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OUM"),simmap.tree=TRUE))
oum_OUMV_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OUMV"),simmap.tree=TRUE))
oum_OUMVA_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OUMVA"),simmap.tree=TRUE))

#oumv
oumv_BM1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumv_ms,model=c("BM1"),simmap.tree=TRUE))
oumv_BMS_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumv_ms,model=c("BMS"),simmap.tree=TRUE))
oumv_OU1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumv_ms,model=c("OU1"),simmap.tree=TRUE))
oumv_OUM_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumv_ms,model=c("OUM"),simmap.tree=TRUE))
oumv_OUMV_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumv_ms,model=c("OUMV"),simmap.tree=TRUE))
oumv_OUMVA_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumv_ms,model=c("OUMVA"),simmap.tree=TRUE))

#oumva
oumva_BM1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumva_ms,model=c("BM1"),simmap.tree=TRUE))
oumva_BMS_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumva_ms,model=c("BMS"),simmap.tree=TRUE))
oumva_OU1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumva_ms,model=c("OU1"),simmap.tree=TRUE))
oumva_OUM_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumva_ms,model=c("OUM"),simmap.tree=TRUE))
oumva_OUMV_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumva_ms,model=c("OUMV"),simmap.tree=TRUE))
oumva_OUMVA_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oumva_ms,model=c("OUMVA"),simmap.tree=TRUE))

#comparison of aic to see best fitting model. model under which data was simulated be the best fitting.

#bm1
meanaicc_bm1_BM1_sim_ms<-t(sapply(bm1_BM1_sim_ms,function(x) x$AICc))
meanaicc_bm1_BM1_sim_ms
mean(meanaicc_bm1_BM1_sim_ms) #183.36

meanaicc_bm1_BMS_sim_ms<-t(sapply(bm1_BMS_sim_ms,function(x) x$AICc))
meanaicc_bm1_BMS_sim_ms
mean(meanaicc_bm1_BMS_sim_ms) #185.3559

meanaicc_bm1_OU1_sim_ms<-t(sapply(bm1_OU1_sim_ms,function(x) x$AICc))
meanaicc_bm1_OU1_sim_ms
mean(meanaicc_bm1_OU1_sim_ms) #185.7876

meanaicc_bm1_OUM_sim_ms<-t(sapply(bm1_OUM_sim_ms,function(x) x$AICc))
meanaicc_bm1_OUM_sim_ms
mean(meanaicc_bm1_OUM_sim_ms) #188.2847

meanaicc_bm1_OUMV_sim_ms<-t(sapply(bm1_OUMV_sim_ms,function(x) x$AICc))
meanaicc_bm1_OUMV_sim_ms
mean(meanaicc_bm1_OUMV_sim_ms) #190.647

meanaicc_bm1_OUMVA_sim_ms<-t(sapply(bm1_OUMVA_sim_ms,function(x) x$AICc))
meanaicc_bm1_OUMVA_sim_ms
mean(meanaicc_bm1_OUMVA_sim_ms) #190.4572

sim_bm1_ms<- c(419.8318,421.3853,420.4748,420.549,422.3636,423.0157)
akaike.weights(sim_bm1_ms)

#bms
meanaicc_bms_BM1_sim_ms<-t(sapply(bms_BM1_sim_ms,function(x) x$AICc))
meanaicc_bms_BM1_sim_ms
mean(meanaicc_bms_BM1_sim_ms) #204.3592

meanaicc_bms_BMS_sim_ms<-t(sapply(bms_BMS_sim_ms,function(x) x$AICc))
meanaicc_bms_BMS_sim_ms
mean(meanaicc_bms_BMS_sim_ms) #204.333

meanaicc_bms_OU1_sim_ms<-t(sapply(bms_OU1_sim_ms,function(x) x$AICc))
meanaicc_bms_OU1_sim_ms
mean(meanaicc_bms_OU1_sim_ms) #206.7868

meanaicc_bms_OUM_sim_ms<-t(sapply(bms_OUM_sim_ms,function(x) x$AICc))
meanaicc_bms_OUM_sim_ms
mean(meanaicc_bms_OUM_sim_ms) #209.1473

meanaicc_bms_OUMV_sim_ms<-t(sapply(bms_OUMV_sim_ms,function(x) x$AICc))
meanaicc_bms_OUMV_sim_ms
mean(meanaicc_bms_OUMV_sim_ms) # 209.484

meanaicc_bms_OUMVA_sim_ms<-t(sapply(bms_OUMVA_sim_ms,function(x) x$AICc))
meanaicc_bms_OUMVA_sim_ms
mean(meanaicc_bms_OUMVA_sim_ms) #211.9325

sim_bms_ms<- c(458.6416,454.6759,460.7994,462.9802,459.1082,185.6256)
akaike.weights(sim_bms)

#ou1
meanaicc_ou1_BM1_sim_ms<-t(sapply(ou1_BM1_sim_ms,function(x) x$AICc))
meanaicc_ou1_BM1_sim_ms
mean(meanaicc_ou1_BM1_sim_ms) #150.3567

meanaicc_ou1_BMS_sim_ms<-t(sapply(ou1_BMS_sim_ms,function(x) x$AICc))
meanaicc_ou1_BMS_sim_ms
mean(meanaicc_ou1_BMS_sim_ms) #151.8195

meanaicc_ou1_OU1_sim_ms<-t(sapply(ou1_OU1_sim_ms,function(x) x$AICc))
meanaicc_ou1_OU1_sim_ms
mean(meanaicc_ou1_OU1_sim_ms) #146.8343

meanaicc_ou1_OUM_sim_ms<-t(sapply(ou1_OUM_sim_ms,function(x) x$AICc))
meanaicc_ou1_OUM_sim_ms
mean(meanaicc_ou1_OUM_sim_ms) #149.3897

meanaicc_ou1_OUMV_sim_ms<-t(sapply(ou1_OUMV_sim_ms,function(x) x$AICc))
meanaicc_ou1_OUMV_sim_ms
mean(meanaicc_ou1_OUMV_sim_ms) #150.6897

meanaicc_ou1_OUMVA_sim_ms<-t(sapply(ou1_OUMVA_sim_ms,function(x) x$AICc))
meanaicc_ou1_OUMVA_sim_ms
mean(meanaicc_ou1_OUMVA_sim_ms) #150.3465

sim_ou1<- c(363.6772,365.7432,343.7365,345.9011,348.1397,-1.524019e+13)
akaike.weights(sim_ou1)

#oum
meanaicc_oum_BM1_sim_ms<-t(sapply(oum_BM1_sim_ms,function(x) x$AICc))
meanaicc_oum_BM1_sim_ms
mean(meanaicc_oum_BM1_sim_ms) #64.6109

meanaicc_oum_BMS_sim_ms<-t(sapply(oum_BMS_sim_ms,function(x) x$AICc))
meanaicc_oum_BMS_sim_ms
mean(meanaicc_oum_BMS_sim_ms) #66.09303

meanaicc_oum_OU1_sim_ms<-t(sapply(oum_OU1_sim_ms,function(x) x$AICc))
meanaicc_oum_OU1_sim_ms
mean(meanaicc_oum_OU1_sim_ms) #65.39562

meanaicc_oum_OUM_sim_ms<-t(sapply(oum_OUM_sim_ms,function(x) x$AICc))
meanaicc_oum_OUM_sim_ms
mean(meanaicc_oum_OUM_sim_ms) #57.34175

meanaicc_oum_OUMV_sim_ms<-t(sapply(oum_OUMV_sim_ms,function(x) x$AICc))
meanaicc_oum_OUMV_sim_ms
mean(meanaicc_oum_OUMV_sim_ms) #60.10721

meanaicc_oum_OUMVA_sim_ms<-t(sapply(oum_OUMVA_sim_ms,function(x) x$AICc))
meanaicc_oum_OUMVA_sim_ms
mean(meanaicc_oum_OUMVA_sim_ms) #35.59733 ###wtf???

#oumv
meanaicc_oumv_BM1_sim_ms<-t(sapply(oumv_BM1_sim_ms,function(x) x$AICc))
meanaicc_oumv_BM1_sim_ms
mean(meanaicc_oumv_BM1_sim_ms) #91.30437

meanaicc_oumv_BMS_sim_ms<-t(sapply(oumv_BMS_sim_ms,function(x) x$AICc))
meanaicc_oumv_BMS_sim_ms
mean(meanaicc_oumv_BMS_sim_ms) # 90.42694

meanaicc_oumv_OU1_sim_ms<-t(sapply(oumv_OU1_sim_ms,function(x) x$AICc))
meanaicc_oumv_OU1_sim_ms
mean(meanaicc_oumv_OU1_sim_ms) #86.43928

meanaicc_oumv_OUM_sim_ms<-t(sapply(oumv_OUM_sim_ms,function(x) x$AICc))
meanaicc_oumv_OUM_sim_ms
mean(meanaicc_oumv_OUM_sim_ms) #76.43124 ##this should not be the case

meanaicc_oumv_OUMV_sim_ms<-t(sapply(oumv_OUMV_sim_ms,function(x) x$AICc))
meanaicc_oumv_OUMV_sim_ms
mean(meanaicc_oumv_OUMV_sim_ms) #79.13219  ####oum does better????

meanaicc_oumv_OUMVA_sim_ms<-t(sapply(oumv_OUMVA_sim_ms,function(x) x$AICc))
meanaicc_oumv_OUMVA_sim_ms
mean(meanaicc_oumv_OUMVA_sim_ms) #-884.3681????

#oumva
meanaicc_oumva_BM1_sim_ms<-t(sapply(oumva_BM1_sim_ms,function(x) x$AICc))
meanaicc_oumva_BM1_sim_ms
mean(meanaicc_oumva_BM1_sim_ms) #102.8744

meanaicc_oumva_BMS_sim_ms<-t(sapply(oumva_BMS_sim_ms,function(x) x$AICc))
meanaicc_oumva_BMS_sim_ms
mean(meanaicc_oumva_BMS_sim_ms) #98.23815

meanaicc_oumva_OU1_sim_ms<-t(sapply(oumva_OU1_sim_ms,function(x) x$AICc))
meanaicc_oumva_OU1_sim_ms
mean(meanaicc_oumva_OU1_sim_ms) #100.2095

meanaicc_oumva_OUM_sim_ms<-t(sapply(oumva_OUM_sim_ms,function(x) x$AICc))
meanaicc_oumva_OUM_sim_ms
mean(meanaicc_oumva_OUM_sim_ms) #89.38316

meanaicc_oumva_OUMV_sim_ms<-t(sapply(oumva_OUMV_sim_ms,function(x) x$AICc))
meanaicc_oumva_OUMV_sim_ms
mean(meanaicc_oumva_OUMV_sim_ms) #76.61978

meanaicc_oumva_OUMVA_sim_ms<-t(sapply(oumva_OUMVA_sim_ms,function(x) x$AICc))
meanaicc_oumva_OUMVA_sim_ms
mean(meanaicc_oumva_OUMVA_sim_ms) #-2905.722

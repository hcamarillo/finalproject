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

#simulatons for mantis shrimp####
regime_stomatopod <- as.matrix(ordered_stomatopod_data)[, "type"]
regime_stomatopod ## should list your species, and whether they are shallow or deep

### need to make stochastic character maps of habitat type using the make.simmap function
trees_regime_ms<-make.simmap(pruned_stomatopod,regime_stomatopod, nsim=100) ### Make 500 stochastic character maps for habitat depth
trees_regime_ms ### should say "500 phylogenetic trees with mapped discrete characters"

class(trees_regime_ms)

#Simulate an Ornstein-Uhlenbeck model with different state means
#and sigma^2 per selective regime, but same alpha

alpha_oum=c(1.0, 1.0)
sigma.sq_oum=c(0.45,0.45)
theta0=1.0
theta_oum=c(1.0,2.0)

#make simulated data

#OUM
ktdata_sim_oum_ms<- OUwie.sim(trees_regime_ms[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_oum, sigma.sq=sigma.sq_oum, theta0=theta0, theta=theta_oum)
ktdata_sim_oum_ms<- data.frame(ordered_stomatopod_data$Species,ordered_stomatopod_data$type, ktdata_sim_oum_ms$X)

#running OUwie with simulated data

#oum
oum_BM1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("BM1"),simmap.tree=TRUE))
oum_BMS_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("BMS"),simmap.tree=TRUE))
oum_OU1_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OU1"),simmap.tree=TRUE))
oum_OUM_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OUM"),simmap.tree=TRUE))
oum_OUMV_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OUMV"),simmap.tree=TRUE))
oum_OUMVA_sim_ms <- lapply(1:100, function(x) OUwie(trees_regime_ms[[x]], ktdata_sim_oum_ms,model=c("OUMVA"),simmap.tree=TRUE))

#model selection
#oum
meanaicc_oum_BM1_sim_ms<-t(sapply(oum_BM1_sim_ms,function(x) x$AICc))
meanaicc_oum_BM1_sim_ms
mean(meanaicc_oum_BM1_sim_ms) #64.6109

meanaicc_oum_BMS_sim_ms<-t(sapply(oum_BMS_sim_ms,function(x) x$AICc))
meanaicc_oum_BMS_sim_ms
mean(meanaicc_oum_BMS_sim_ms) #64.6109,66.09303

meanaicc_oum_OU1_sim_ms<-t(sapply(oum_OU1_sim_ms,function(x) x$AICc))
meanaicc_oum_OU1_sim_ms
mean(meanaicc_oum_OU1_sim_ms) #64.6109,66.09303,65.39562

meanaicc_oum_OUM_sim_ms<-t(sapply(oum_OUM_sim_ms,function(x) x$AICc))
meanaicc_oum_OUM_sim_ms
mean(meanaicc_oum_OUM_sim_ms) #64.6109,66.09303,65.39562,57.34175

meanaicc_oum_OUMV_sim_ms<-t(sapply(oum_OUMV_sim_ms,function(x) x$AICc))
meanaicc_oum_OUMV_sim_ms
mean(meanaicc_oum_OUMV_sim_ms) #64.6109,66.09303,65.39562,57.34175,60.10721

meanaicc_oum_OUMVA_sim_ms<-t(sapply(oum_OUMVA_sim_ms,function(x) x$AICc))
meanaicc_oum_OUMVA_sim_ms
mean(meanaicc_oum_OUMVA_sim_ms) #64.6109,66.09303,65.39562,57.34175,60.10721,35.59733 ###wtf???

sim_oum_ms<- c(64.6109,66.09303,65.39562,57.34175,60.10721,35.59733)
akaike.weights(sim_oum_ms)





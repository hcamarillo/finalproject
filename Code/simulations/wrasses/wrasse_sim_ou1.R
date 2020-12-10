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


#simulations for wrasses#####


regime_labrid <- as.matrix(ordered_labrid_data)[, "type"]
regime_labrid ## should list your species, and whether they are shallow or deep

### need to make stochastic character maps of habitat type using the make.simmap function
trees_regime<-make.simmap(pruned_labrid,regime_labrid, nsim=100) ### Make 500 stochastic character maps for habitat depth
trees_regime ### should say "500 phylogenetic trees with mapped discrete characters"

class(trees_regime)


#Simulate an Ornstein-Uhlenbeck model with different state means
#and sigma^2 per selective regime, but same alpha

alpha_ou1=c(0.1,0.1)
sigma.sq_ou1=c(0.9,0.9)
theta0=1.0
theta_ou1=c(1.0,1.0)

#making new data & data frames with simulated data

#OU1
ktdata_sim_ou1<- OUwie.sim(trees_regime[[1]], data=NULL, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha_ou1, sigma.sq=sigma.sq_ou1, theta0=theta0, theta=theta_ou1)

simkt_ou1<- data.frame(ordered_labrid_data$Species,ordered_labrid_data$type, ktdata_sim_ou1$X)

#running OUwie with simulated data
#ou1
ou1_BM1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("BM1"),simmap.tree=TRUE))
ou1_BMS_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("BMS"),simmap.tree=TRUE))
ou1_OU1_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OU1"),simmap.tree=TRUE))
ou1_OUM_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OUM"),simmap.tree=TRUE))
ou1_OUMV_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OUMV"),simmap.tree=TRUE))
ou1_OUMVA_sim <- lapply(1:100, function(x) OUwie(trees_regime[[x]], simkt_ou1,model=c("OUMVA"),simmap.tree=TRUE))

#model selection from simulated data
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

sim_ou1<- c(389.9667,385.6754,367.171,369.2455,365.2558,-1320069)
akaike.weights(sim_ou1)



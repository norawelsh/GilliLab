#WELSH et al 2021
#Clear Workspace
rm(list = ls())

#Install and load all necessary packages
install.packages("Hmisc")

devtools::install_github("kuijjerlab/lionessR")

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("enrichplot")

library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(readxl)
library(column_to_rownames)
library(Hmisc)
library(readxl)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(tidyverse)
library("imputeTS")


#Read in data
#Index are all values (both inflammatory and CNS targets) 
index<- read_excel("C:/Users/norawelsh/Desktop/Dyna July2020/Lioness.Patient.Exp.xlsx")

index<- data.frame(index)
index<- data.frame(index, row.names = 1)

index<-data.matrix(index)

#Info are the patient data
info<- read_excel("C:/Users/norawelsh/Desktop/Dyna July2020/Targets.xlsx")
info<- data.frame(info)

#Disease_nodes are nodes categorized as inflammatory or CNS injury
disease_nodes <- read.csv("C:/Users/norawelsh/Desktop/Dyna July2020/Disease Nodes.csv", header=T,as.is=T)

#Separating data by Diagnosis
#List from Info sheet which data to pull

CIS<- which(info$Diagnosis=="CIS")
CISCon <- which(info$Diagnosis=="CISCON")
RRMS <- which(info$Diagnosis=="RRMS",)
Neg <- which(info$Diagnosis=="NEG")

#Generalized holder, change which group you are analyzing
group_holder<- CISCon


#Analyze Data with similarity matrices

#Pearson Correlations on those data
netgroup_holder <- cor(t(index[,group_holder]))

#Convert adjacency to edgelists, 
melted <- melt(upper.tri(netgroup_holder))
melted <- melted[which(melted$value),]
values <- netgroup_holder[which(upper.tri(netgroup_holder))]
melted <- cbind(melted[,1:2], values)
genes <- row.names(netgroup_holder)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep="_")
tosub <- melted

#Creates a list of row names where the r value is greater than 0.7 or less than -0.7
group_holderedges<-(tosub[which(abs(tosub[,3])>0.7),])
group_holderedges <- row.names(tosub[which(abs(tosub[,3])>0.7),])

#Get P Values from Index 
P<- rcorr(t(index[,group_holder]))
netgroup_holder_Pvalue<- P$P

#Convert adjacency to Edgelist
melted <- melt(upper.tri(netgroup_holder_Pvalue))
melted <- melted[which(melted$value),]
values <- netgroup_holder_Pvalue[which(upper.tri(netgroup_holder_Pvalue))]
melted <- cbind(melted[,1:2], values)
genes <- row.names(netgroup_holder_Pvalue)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep="_")
tosub2 <- melted

#Creates a list of row names where the p value < 0.05
netgroup_holderedges_Pvalues <- (tosub2[which(abs(tosub2[,3])<0.05),])

netgroup_holderedges_Pvalues <- row.names(tosub2[which(abs(tosub2[,3])<0.05),])


#Convert rownames to list
p_list<- as.list(netgroup_holderedges_Pvalues)
r_list<- as.list(group_holderedges)
#Make final list of targets with P Value<0.05 and R value >0.7
Final_list<- Reduce(intersect, list(netgroup_holderedges_Pvalues, group_holderedges))


#Lioness

#model single sample networks based on co expression, so the raw data, using Lioness
#output is Pearson Correlation between every gene for each patient 
cormat <- lioness(index, netFun)
cormat
#Only choose patients from group
COL<- cormat[,c(1,2)]
NUM<- cormat[,group_holder+2]

#Unsorted Correlations of all group_holder Patients
Lioness_net_group_holder<- cbind(COL, NUM)

#paste row.names from adjacency matrix so rows are the same
row.names(Lioness_net_group_holder) <- paste(Lioness_net_group_holder[,1], 
                                             Lioness_net_group_holder[,2], sep="_")

#Choose targets that are the same as list made above and print first two columns 
l_net_group_holder<- Lioness_net_group_holder[which(row.names(Lioness_net_group_holder) %in% Final_list ), 
                            1:ncol(Lioness_net_group_holder)]


#Making Aggregate Network, no row columns so we can analyze
l_net_group_holder_agg<- Lioness_net_group_holder[which(row.names(Lioness_net_group_holder)
                                        %in% Final_list ),3:ncol(Lioness_net_group_holder)]
l_net_group_holder_agg<- as.matrix(l_net_group_holder_agg)

Avg<- rowMeans(l_net_group_holder_agg)
Avg<-as.data.frame(Avg)

Agg_net<-cbind(l_net_group_holder[,1], l_net_group_holder[,2], Avg$Avg)


#Make final aggregate network 
net_group_holder_Lioness <- graph_from_data_frame(d=Agg_net, vertices = disease_nodes, directed=F)


#Make graphs for analysis
#Graphopt
l_graphopt_group_holder_Lioness <- layout_with_graphopt(net_group_holder_Lioness, charge=0.08)


#coloring
#Inflammatory vs CNS Targets
catergory_color <- c("pink","light blue")

#Change graph name to appropriate network
holder <- net_group_holder_Lioness
l_holder <- l_graphopt_group_holder_Lioness
l_holder<- norm_coords(l_holder,ymin = -1, ymax = 1, xmin = -1, xmax = 1)
graph_name <- "CIS Con" #Change this to title of group currently analyzing

#Inflammatory v CNS injury graphopt
V(holder)$color <- catergory_color[V(holder)$category]
E(holder)$color <- ifelse(E(holder)$V3 > 0, 'dark grey','red')
deg<-degree(holder, mode="total")
V(holder)$size <- deg*3
plot(holder,
     edge.arrow.size=0.07,
     vertex.label.font=2,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     main=graph_name, 
     vertex.frame.color = "white",
     rescale=F,
     layout=l_holder*1.2)
legend(x=-1.5,y=-1.5,c("Inflammatory","CNS Injury"),
       pch=21,col="#777777",pt.bg=catergory_color, pt.cex=2,cex=.8,bty="n",ncol=3)


#Community Graph
clp <- cluster_optimal(holder)
class(clp)
V(holder)$color <- type_color[V(holder)$type]
plot(clp, holder,
     edge.arrow.size=0.07,
     edge.color = "grey",
     vertex.label.font=2,
     vertex.size=8,
     vertex.label.cex = 1,
     vertex.label.color = "black",
     main=graph_name, 
     vertex.frame.color = "white",
     rescale=F,
     layout=l_holder*1.2)


#Circle Plot
l_circle_group_holder_Lioness <- layout_in_circle(net_group_holder_Lioness)

holder <- net_group_holder_Lioness
l_holder <- l_circle_group_holder_Lioness
graph_Ame <- "CIS" #Change this
l_holder<- norm_coords(l_holder,ymin = -1, ymax = 1, xmin = -1, xmax = 1)



#Circle Plot black 
plot(holder, edge.arrow.size=0.1,
     vertex.label=NA,
     vertex.size=8,
     vertex.color = "black",
     main=graph_Ame,
     rescale=F,
     layout=l_holder*1.2)

#Circle plot Inflammatory v CNS injury
E(holder)$color <- ifelse(E(holder)$V3 > 0, 'dark grey','red')
V(holder)$color <- catergory_color[V(holder)$category]
plot(holder, edge.arrow.size=0.1,
     vertex.label=NA,
     vertex.size=8,
     main=graph_Ame,
     rescale=F,
     layout=l_holder*1.2)
legend(x=-1.5,y=-1.5,c("Inflammatory","CNS Injury"),
       pch=21,col="#777777",pt.bg=catergory_color, pt.cex=2,cex=.8,bty="n",ncol=3)


#1 degree of freedom from CNS targets
#List of all CNS targets
N_nodes = c("NFL","p231 Tau","Total Tau","AB 1-42","AB 1-40","CHI3L1","GFAP","RAGE","s100b","TREM2","MBP")


N_group_holder <-ego(net_group_holder_Lioness, order=1, nodes = N_nodes, mode="all",mindist=0)
Nnet_group_holder <- induced_subgraph(net_group_holder_Lioness,unlist(N_group_holder))
l_graphopt_Ngroup_holder <- layout_with_graphopt(Nnet_group_holder, charge=.08)

holder <- Nnet_group_holder
l_holder <- l_graphopt_Ngroup_holder
l_holder<- norm_coords(l_holder,ymin = -1, ymax = 1, xmin = -1, xmax = 1)
graph_name <- "group_holder"  #Change this


#Inflammatory v CNS injury, CNS targets
V(holder)$color <- catergory_color[V(holder)$category]
E(holder)$color <- ifelse(E(holder)$V3 > 0, 'dark grey','red')

plot(holder,
     edge.arrow.size=0.07,
     vertex.label.font=2,
     vertex.label.cex =1,
     vertex.label.color = "black",
     main=graph_name, 
     vertex.frame.color = "white",
     rescale=F,
     layout=l_holder*1.2)



#Centrality analysis 
#Closeness Centrality

CC<- as.data.frame(closeness(net_group_holder_Lioness))
CC<- CC[order(CC$`closeness(net_group_holder_Lioness)`), , drop = FALSE]
CC20<- tail (CC, 56)
setDT(CC20, keep.rownames = TRUE)


#Degree of network
deg<- as.data.frame(degree(net_group_holder_Lioness))
deg<- deg[order(deg$`degree(net_group_holder_Lioness)`), , drop = FALSE]
deg20<- tail(deg, 56)
setDT(deg20, keep.rownames = TRUE)

#Intersecting Targets from Closeness and Degree lists
Topgroup_holder<-intersect(deg20$rn, CC20$rn)

#Creating a dataframe that has all targets from list above and value from Closeness measurment
TC_group_holder<- CC20[CC20$rn %in% Topgroup_holder,]
TC_group_holder <-TC_group_holder %>% rename(Target = rn, 
                                            group_holder.Closeness = 'closeness(net_group_holder_Lioness)')

TD_group_holder<- deg20[deg20$rn %in% Topgroup_holder]
TD_group_holder <-TD_group_holder %>% rename(Target = rn, 
                                            group_holder.Degree = 'degree(net_group_holder_Lioness)' )


All_group_holder<- merge(TC_group_holder, TD_group_holder)# Change these so they are named by group
All_group_holder<- All_CIS
All_group_holder<- All_CIS_Con
All_group_holder<- All_Neg
All_group_holder<- All_RRMS_NA

write.csv(All_group_holder,
          file="group_holder.All.csv")  #Change File name so all group names are different
                                        #Repeat this for all groups


#Making plot of all 'Influential' Points in each Network
#Plotting with degree and closeness as variables
#Combine all objects by targets and clean data
test<- list(All_CIS, All_CIS_Con, All_Neg, All_RRMS_NA) %>% 
        reduce(full_join, by = "Target")
test[is.na(test)] = 0
test<- test %>% gather(sample, value, -Target)

#Export to excel
write.csv(test, file= "Paper.Centrality.csv")

#Had to clean data in Excel AND SORT BY TWO VARIABLES

test2<- read_excel("C:/Users/norawelsh/Desktop/Dyna July2020/Paper.Centrality.xlsx")
test2<- data.frame(test2, row.names = 1)

#Example of what data should look like
#       Target Disease.Course    Closeness Degree
#1        CCL1            CIS 0.0004364906      3
#2       CCL11            CIS 0.0004397537      4
#3       CCL13            CIS 0.0003495281      3
#4       CCL19            CIS 0.0004384042      4
#5       CCL26            CIS 0.0004387889      5
  
#Plot Influential targets
ggplot(test2, aes(x=Disease.Course, y=Target, size = Degree, color = Closeness)) +geom_point(alpha = 0.75) +
        labs( x= "Disease Course", y="Targets", color = "Closeness") + 
        scale_color_gradient(low = "red", high = "blue", limits=c(0.0002, .0022)) + theme_bw()



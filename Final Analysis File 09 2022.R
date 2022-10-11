#########################################################
#########################################################
#Final Data analysis for networks###
rm(list = ls())
setwd("C:/Users/norawelsh/Desktop/092022 Networks FINAL/")
library(cutpointr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(rmarkdown)
library(knitr)
library(ggbiplot)
library(gplots)
library(ggpubr)
library(survival)
library(survminer)
library(caret)
library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(readxl)
library(Hmisc)
library(readxl)
library(data.table)
library(reshape2)
library(tidyr)
library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(devtools)
library(pastecs)
library(table1)
library(jpeg)
library(qvalue)
library(ggVennDiagram)
library(gplots)
library(UpSetR)
library(igraph)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condSURV)
#1) Read in Data
#Most recent protein data = ALL.Patients.09.2022.xlsx
index<- read_excel("C:/Users/norawelsh/Desktop/Dyna July2020/ALL.Patients.09.2022.xlsx")
index<- data.frame(index)
index<- data.frame(index, row.names = 1)
index<-data.matrix(index)

#Most recent patient data = ALL.Targets.09.2022.xlsx
info<- read_excel("C:/Users/norawelsh/Desktop/Dyna July2020/ALL.Targets.09.2022.xlsx")
info<- data.frame(info)
#Create new variables for MS v NON MS
info$MS <- ifelse(info$Diagnosis == "NEG", "no MS", "MS")
info$MS_1 <- ifelse(info$Diagnosis == "NEG", 0, 1)
info$BMI<- as.numeric(info$BMI)
info$TTanyactivity<- ifelse(info$TTanyAct =="50", NA, info$TTanyAct)
info$TTMRIActivity<- ifelse(info$TTActivity =="50", NA, info$TTActivity)
info$TTAttack<- ifelse(info$Attack =="50", NA, info$Attack)
info$PatientwithAct<- ifelse(info$TTanyactivity=="NA", "no","yes")
info$PatientwithOCB<- ifelse(info$OCB==0, "no","yes")
info$TTMRIActivity<- as.numeric(info$TTMRIActivity)

#2) Patient exclusion (did not meet inclusion criteria)
#No PPMS
#Brain tumor is confounding
#Follow-up needs to be >12 mo
#CANNOT be on treatment at time of LP
PPMS<- which(info$Diagnosis=="PMS")
tumor<- which(info$Patient=="X0456.2016")
short<- which(info$Follow.Up <11)
question<- which(info$Patient=="X0523.2016")
RX<- which(info$Treatment.first.12.months=="yes-on drug at time of LP")
remove_patients<- c(PPMS, tumor, short, question, RX)
#remove patients from index and info
info<-info[-remove_patients,] 
index<- index[,-remove_patients]


#3) Table1 for patient statistics
#Function to determine P Value
pvalue <- function(x, ...) {
  x <- x[-length(x)]  # Remove "overall" group
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform an ANOVA
    p <- summary(aov(y ~ g))[[1]][["Pr(>F)"]][1]
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

table1(~Sex+Age+Diagnosis+Follow.Up+PatientwithAct+TTanyactivity+TTMRIActivity+TTAttack+
         PatientwithOCB+OCB+Ig.Index|Activity,data=info, 
       render.missing =NULL, extra.col = list('P-Value'=pvalue), extra.col.pos = 4)

#Normal Data Check
ggqqplot(indext$CXCL13) #not normal
h<- data.frame(t(model_index))
ggqqplot(h$CXCL13) #normalized by log

#4) Determine proteins that are differentially expressed from NIND Controls
#Correcting for Age, Sex, and BMI
#Normalize by log(data) since data is not normal
model_index<- data.frame(log(index)) 
design<- model.matrix(~0+MS +Age+Sex+BMI, data = info)
colnames(design) <- c("MS", "NIND", "Age", "Sex", "BMI") #rename columns
contrasts<- makeContrasts(MS-NIND, levels= colnames(design)) #create model matrix

#This relies on assumption of normal distribution, which our data are not
fit<- lmFit(model_index, design) #WLS linear model method
fit2 <- contrasts.fit(fit, contrasts) #based on linear model, estimated coefficients
fit2e <- eBayes(fit2) #given linear fit compute moderated T-statistics
top_ALL<- topTable(fit2e, number = 67) #table table of linear model fit
view(top_ALL)
sig_proteins_pvalue<- top_ALL[which( top_ALL$P.Value<0.05),]
sig_proteins_logFC<- top_ALL[which( top_ALL$logFC>1 |top_ALL$logFC< -1 ),]
Final_sig_list<- Reduce(intersect, list(sig_proteins_pvalue, sig_proteins_logFC))

#Graph of these targets
test<- top_ALL
test$gene_symbol<- row.names(test)
test$diffexpressed <- "NO"
test$diffexpressed[test$logFC > 1 & test$P.Value < 0.05] <- "UP"
test$diffexpressed[test$logFC < -1 & test$P.Value < 0.05] <- "DOWN"

test$delabel <- NA
test$delabel[test$diffexpressed != "NO"] <- test$gene_symbol[test$diffexpressed != "NO"]


ggplot(data=test, aes(x=logFC, y=-log10(P.Value), label = delabel,col=diffexpressed)) + 
  geom_point() +
  geom_text_repel()+
  theme_bw()  +
  theme(text = element_text(size = 17), legend.position = "none")+
  scale_color_manual(values=c("black", "red", "blue")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


#List of significant targets
Final_Targets<- as.list(rownames(Final_sig_list))
view(Final_Targets)
Final_Targets_List<- c("CXCL13","CXCL9","IFNG","CCL22", "NFL", "CCL11","CCL26",
                       "IL4", "IL1B", "IgG1", "CXCL10","IgM","TNFA","IL10","CCL13","IL6",
                       "TREM2")

#5) Graphs of significant targets NIND compared to MS
indext<- as.data.frame(t(model_index))
indext$MS<- info$MS


holder<- "CXCL13"
graph_name<- holder
holder2<- .2


my_comps<- list(c("MS", "no MS"))
mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("MS_NIND_log_", graph_name, ".tiff", sep = ""))

tiff(file=mypath,units="in", width=5, height=4, res=300, compression = 'lzw' )
mytitle = paste("MS", graph_name)
ggboxplot(indext, x = "MS", y = holder,
          color = "MS", palette = "jco")+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width =0.3)+ 
  stat_compare_means(comparisons = my_comps) +
  theme(legend.position="none", plot.title = element_text(size=11),
        axis.title.x = element_blank() )
dev.off()


#6) Determine predictive value of intrathecally produced proteins in activity
#Subset patients and patient info for only MS patients
Neg2 <- which(info$Diagnosis=="NEG")
remove_patients<- c(Neg2)
info_MS<-info[-remove_patients,] 
index_MS<- index[,-remove_patients]

#Ensure only MS patients
table1(~Sex+Age+Diagnosis+Follow.Up+PatientwithAct+TTanyactivity+TTMRIActivity+TTAttack+
         PatientwithOCB+OCB+Ig.Index|Activity,data=info_MS, 
       render.missing =NULL, extra.col = list('P-Value'=pvalue), extra.col.pos = 4)
#Subset data for only significant proteins
Final_Targets_List_1<- c("CXCL13","CXCL9","IFNG","CCL22", "NFL", "CCL11","CCL26", "TREM2")
Final_Targets_List_2<- c("IL4", "IL1B", "IgG1", "CXCL10","IgM","TNFA","IL10","CCL13","IL6")
index_MS<- index_MS[Final_Targets_List_2,]


#Create data frame with pertinent variables (Age, Sex, BMI, and Treatment in first 12 mo)
use_me<- data.frame(t(index_MS))
use_me$Age<- info_MS$Age
use_me$Sex<- factor(info_MS$Sex)
use_me$Sex<- unclass(use_me$Sex)
use_me$Activity<- factor(info_MS$Activity)
use_me$Activity<- unclass(use_me$Activity)
use_me$TX<- factor(info_MS$Treatment.first.12.months)
use_me$TX<- unclass(use_me$TX)
use_me$BMI<- info_MS$BMI

use_me2<- use_me %>% relocate(Activity) # Change which variable to compare

variables <- colnames(use_me2)[2:ncol(use_me2)]
formulas <- list()
for (i in seq_along(variables)) {
  tmp <- combn(variables, i)
  tmp <- apply(tmp, 2, paste, collapse="+TX+Age+Sex+BMI+")
  tmp <- paste0("Activity~", tmp) ####Change this
  formulas[[i]] <- tmp
}

#Loop to make linear models of different combos of Igs to predict age

formulas <- unlist(formulas)
formulas <- sapply(formulas, as.formula)
models <- lapply(formulas, lm, data=use_me2)


#Extract p values from model
#Function to extract P value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

reg_val<- lapply(models, lmp)
reg_val<- as.data.frame(t(reg_val))
reg_val<- as.data.frame(t(reg_val))
sig_rv<- which(reg_val$V1<0.001)

Sig_Act_5<- (reg_val[sig_rv,])
view(Sig_Act_3)
Sig_act_list<- c("CCL11", "IL10", "CXCL9", "TREM2") #Identified two proteins


#7) AUROC and Survival Curve for proteins identified by linear regression

index_MS<- index_MS[Sig_act_list,]
view(index_MS) 
tindex<- data.frame(t(index_MS))
use_me<- data.frame(t(index_MS))
use_me$Age<- info_MS$Age
use_me$Sex<- factor(info_MS$Sex)
use_me$Sex<- unclass(use_me$Sex)
use_me$Activity<- factor(info_MS$Activity)
use_me$Activity<- unclass(use_me$Activity)
use_me$TX<- factor(info_MS$Treatment.first.12.months)
use_me$TX<- unclass(use_me$TX)
use_me$BMI<- info_MS$BMI
use_me$act_status<-ifelse(info_MS$TTanyAct>12,"MSNA","MSA")
use_me$OCB<- info_MS$OCB

use_me2<- use_me %>% relocate(Activity) # Change which variable to compare
use_me2$CCL11<- as.numeric(use_me2$CCL11)
use_me2$IL10<- as.numeric(use_me2$IL10)
use_me2$Age<- as.numeric(use_me2$Age)
use_me2$IgIndex<- as.numeric(info_MS$Ig.Index)




#Change protein of interest here
#"CCL11", "IL10", "CXCL9", "TREM2") 
holder_cp<- use_me2$CXCL10
holder<- "CXCL10"  
graph_name<- holder




#Determine cutoff for each protein 
cp<- cutpointr(use_me2, holder_cp,act_status,
             method = maximize_metric, metric = youden)

#Save statistics for each protein
summary(cp)
sink(paste(holder, ".txt")) ####change####
print(summary(cp))
sink()
plot(cp) #Save if wanted

#Use optimal cutpoint found above to create new binary variable
use_me2$level<- ifelse(holder_cp > cp$optimal_cutpoint ,1, 0) 

use_me2$level<- ifelse(    #use_me2$Sex ==2 &
                          #use_me2$BMI <25 &
                          #use_me2$OCB >3 &
                          #use_me2$CCL11 >39.99 &
                          holder_cp < cp$optimal_cutpoint & 
                          use_me2$Age > 34.5, "MSNA", "MSA"  )
#Compare new variable to Activity to determine predictive content
use_me2$act_status<- factor(use_me2$act_status)
use_me2$level<- factor(use_me2$level)
#Create confusion matrix
example <- confusionMatrix(data=use_me2$level, reference = use_me2$act_status)
example
table <-data.frame(example$table)

plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))



mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Conf.Matr.Age", graph_name, ".tiff", sep = ""))

tiff(file=mypath,units="in", width=5, height=4, res=300, compression = 'lzw' )
ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "green", bad = "red")) +
  theme_bw() +
  xlim(rev(levels(table$Reference)))

dev.off()

#Survival Curve based on cutoffs
#Make new variables
#TTact = time to activity, in first 12 months
#status = anything that will not be used in SC
use_me2$TTAct<- as.numeric(info_MS$TTanyAct)
use_me2$TTAct<- ifelse(use_me2$TTAct>13, 13, use_me2$TTAct)
use_me2$status<-ifelse(info_MS$TTanyAct>12, 0,1)
use_me2$status<- as.numeric(use_me2$status)
fit<- survfit(Surv(TTAct,status) ~ level, data =use_me2)
survdiff(Surv(TTAct,status) ~ level, data =use_me2)

summary(fit)

holder_high<- paste0(holder, "<", round(cp$optimal_cutpoint, digits = 2),"<br/>", 
                     #"+ CXCL10_CXCL9<0.33", "<br/>",
                     #"+IFNG_IL6 >0.746","<br/>",
                     "+ Age <34.5")
holder_low<- paste0(holder,">", round(cp$optimal_cutpoint, digits = 2), "<br/>",
                    #"+ CXCL10_CXCL9>0.33", "<br/>",
                    #"+IFNG_IL6 <0.746","<br/>",
                    "+ Age >34.5")
mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Survival.Curve.Age", graph_name, ".tiff", sep = ""))

tiff(file=mypath,units="in", width=8, height=7, res=300, compression = 'lzw' )
ggsurvplot(fit, data = use_me2, 
           conf.int = TRUE,censor = TRUE,
           pval = TRUE, risk.table = TRUE,
           legend.labs =
             c(holder_high, holder_low),legend = c("none"),
           ggtheme= theme_bw())
dev.off()


#8) Network Analysis based on Activity
#First create edge lists from each group you wish to analyze
view(index)
view(info)
view(Final_Targets_List)
active<- which(info$Activity == "MSA")
list_active<- info$Activity[active]
inactive<- which(info$Activity =="MSNA")
list_inactive<- info$Activity[inactive]
Neg <- which(info$Diagnosis=="NEG")
list_Neg<- info$Patient[Neg]

disease_nodes <- read_excel("C:/Users/norawelsh/Desktop/Dyna July2020/Disease Nodes.xlsx")
disease_nodes<-as.data.frame(disease_nodes[disease_nodes$gene %in% Final_Targets_List,])

index2<- as.data.frame(index[row.names(index) %in% Final_Targets_List,])
index<- index2

#Group analyzing############################CHANGE THESE##############3
group_holder<- inactive
group_list<- as.matrix(list_inactive)
group_index<- index[,group_holder]

#Create edge list for networks using Spearman correlations (data not normal)
netgroup_holder <- cor(t(index[,group_holder]), method = "spearman",
                       use ="na.or.complete")
melted <- melt(upper.tri(netgroup_holder))
melted <- melted[which(melted$value),]
values <- netgroup_holder[which(upper.tri(netgroup_holder))]
melted <- cbind(melted[,1:2], values)
genes <- row.names(netgroup_holder)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep="_")
tosub <- melted

#Creates a list of row names where the r value is greater than 0.5 or less than -0.5
group_holderedges<-(tosub[which(abs(tosub[,3])>0.5),])
group_holderedges <- row.names(tosub[which(abs(tosub[,3])>0.5),])

#P values
P<- rcorr(t(index[,group_holder]), type = "spearman")
netgroup_holder_Pvalue<- P$P
melted <- melt(upper.tri(netgroup_holder_Pvalue))
melted <- melted[which(melted$value),]
values <- netgroup_holder_Pvalue[which(upper.tri(netgroup_holder_Pvalue))]
melted <- cbind(melted[,1:2], values)
genes <- row.names(netgroup_holder_Pvalue)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep="_")
tosub2 <- melted
pvalues<- na.omit(tosub2)
qobj <- qvalue(p=pvalues$values, pi0=1)

#Creates a list of row names where the q value <0.1
q_20<- which(qobj$qvalues<0.2)
netgroup_holderedges_Pvalues <- row.names(tosub2[q_20,])

#Combine Lists
p_list<- as.list(netgroup_holderedges_Pvalues)
r_list<- as.list(group_holderedges)
#Make final list of targets with q value <0.1 and R value >0.5
Final_list<- Reduce(intersect, list(netgroup_holderedges_Pvalues, group_holderedges))
view(Final_list)

#Make list of all signficant correlations of each group you want to compare
NIND_Connections<- Final_list
Non_active_Connection<- Final_list
Active_Connections<- Final_list

#Compare significant edge lists with upset plot/Venn Diagram
x = list('MSA' = Active_Connections, 'MSNA' = Non_active_Connection, 
         'NIND' =NIND_Connections)
mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Venn.Diagram_Activity", ".tiff", sep = ""))
tiff(file=mypath,units="in", width=8, height=7, res=300, compression = 'lzw' )
ggVennDiagram(x, category.names = c("MSA", "MSNA","NIND"),
                  label = "count") +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none", text = element_text(size = 100))
dev.off()

vd<- venn(x)
xy<- attributes(vd)$intersections
sink("activity.vd.txt")
print(xy)
sink()
mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Upset_Activity", ".tiff", sep = ""))
tiff(file=mypath,units="in", width=8, height=7, res=300, compression = 'lzw' )
upset(fromList(x), keep.order = T, point.size = 3.5, line.size = 2, text.scale = 1.75, 
      mainbar.y.label = "Network Correlations", 
      sets.x.label = "Total Correlations")
dev.off()


#9) Create networks from each patient subgroup using Lioness#########
view(NIND_Connections)
view(Active_Connections)
view(Non_active_Connection)

group_holder<- Neg
group_index<- index[,group_holder]
group_final_list<-NIND_Connections
graph_name<- "NIND"
  
  
#Lioness Algorithm
netspear <- function(x, ...) {
  stats::cor(t(x), method="spearman", use ="na.or.complete") 
}
cormat <- lioness(group_index, netspear)
cormat

l_net_group_holder <- as.data.frame(assay(cormat[which(row.names(cormat) %in% group_final_list), ]))
l_net_group_holder<- tibble::rownames_to_column(l_net_group_holder, "Targets")
l_net_group_holder<- l_net_group_holder %>% separate(Targets, c("Tar", "Reg"), sep = "([_])")
row.names(l_net_group_holder) <- paste(l_net_group_holder[,1], l_net_group_holder[,2], sep="_")
Lioness_net_group_holder<- l_net_group_holder
l_net_group_holder_agg<- Lioness_net_group_holder[which(row.names(Lioness_net_group_holder)
                              %in% group_final_list ),3:ncol(Lioness_net_group_holder)]
l_net_group_holder_agg<- as.matrix(l_net_group_holder_agg)
Avg<- rowMeans(l_net_group_holder_agg)
Avg<-as.data.frame(Avg)
Agg_net<-cbind(l_net_group_holder[,1], l_net_group_holder[,2], Avg$Avg)

#Make final aggregate network 
net_group_holder_Lioness <- graph_from_data_frame(d=Agg_net, vertices = disease_nodes, directed=F)


#Make graphs for viewing
#Graphopt
l_graphopt_group_holder_Lioness <- layout_with_graphopt(net_group_holder_Lioness, charge=0.08)

par(mfrow=c(1,1))
catergory_color <- c("pink","light blue")
holder <- net_group_holder_Lioness
l_holder <- l_graphopt_group_holder_Lioness
l_holder<- norm_coords(l_holder,ymin = -1, ymax = 1, xmin = -1, xmax = 1)

#Graphopt
V(holder)$color <- catergory_color[V(holder)$category]
E(holder)$color <- ifelse(E(holder)$V3 > 0, 'dark grey','red')
deg<-degree(holder, mode="total")
V(holder)$size <- deg*2
mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Activity_Network", graph_name, ".tiff", sep = ""))
tiff(file=mypath,units="in", width=8, height=8, res=300, compression = 'lzw' )
plot(holder,
     edge.arrow.size=0.07,
     vertex.label.font=2,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     main=graph_name, adj = 0.75,
     vertex.frame.color = "white",
     rescale=F, repel = T,
     layout=l_holder)
dev.off()


#Community Graph
clp <- cluster_optimal(holder)
class(clp)

mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Activity_Cluster", graph_name, ".tiff", sep = ""))
tiff(file=mypath,units="in", width=8, height=8, res=300, compression = 'lzw' )
plot(clp, holder,
     edge.arrow.size=0.07,
     edge.color = "grey",
     vertex.label.font=2,
     vertex.size=8,
     vertex.label.cex = 0.9,
     vertex.label.color = "black",
     main=graph_name, 
     vertex.frame.color = "white",
     rescale=F,
     layout=l_holder)
dev.off()


#Circle Plot
l_circle_group_holder_Lioness <- layout_in_circle(net_group_holder_Lioness)

holder <- net_group_holder_Lioness
l_holder <- l_circle_group_holder_Lioness
graph_Ame <- graph_name
l_holder<- norm_coords(l_holder,ymin = -1, ymax = 1, xmin = -1, xmax = 1)

#Circle plot Inflammatory v CNS injury
E(holder)$color <- ifelse(E(holder)$V3 > 0, 'dark grey','red')
V(holder)$color <- catergory_color[V(holder)$category]

mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Activity_Circle", graph_name, ".tiff", sep = ""))
tiff(file=mypath,units="in", width=8, height=8, res=300, compression = 'lzw' )
mytitle = paste("Network", graph_name)
plot(holder, edge.arrow.size=0.1,
     vertex.label.font=2,
     vertex.label.cex =1,
     vertex.size=8,
     main=graph_Ame,
     rescale=F,
     layout=l_holder)
dev.off()


#####################################
#10)#Create data frame with pertinent variables from LIONESS to compare based on activity
#(Age, Sex, BMI, and Treatment in first 12 mo)
Neg2 <- which(info$Diagnosis=="NEG")
remove_patients<- c(Neg2)
info_MS<-info[-remove_patients,] 
index_MS<- index[,-remove_patients]

MSA_only_vd_1<- c("CCL26_IFNG","IFNG_IL1B","IFNG_IL4","CCL11_IL6","IgG1_CXCL10", "CCL11_CCL13",
                "CCL26_CCL13", "IFNG_CCL13")
MSA_only_vd_2<- c( "IL1B_CCL13","CCL26_CCL22", "IL1B_CCL22" , "IL10_CCL22", 
                "CCL26_CXCL9", "IL4_CXCL9","CCL13_CXCL9", "IgG1_TNFA")

MSNA_only_vd<- c("IFNG_IL6","IgM_CCL22","CXCL10_CXCL9" ,"CCL22_CXCL9","CXCL9_TNFA","IgG1_NFL")

netspear <- function(x, ...) {
  stats::cor(t(x), method="spearman", use ="na.or.complete") 
}
cormat <- lioness(index_MS, netspear)
cormat
data<- assay(cormat)
data_combo<- data


data2<- data_combo[Act_sig_conns_p,] #Select subset of data

data_cormat<-data.frame(t(data2))
data_cormat$Age<- info_MS$Age
data_cormat$Sex<- factor(info_MS$Sex)
data_cormat$Sex<- unclass(data_cormat$Sex)
data_cormat$Activity<- factor(info_MS$Activity)
data_cormat$Activity<- unclass(data_cormat$Activity)
data_cormat$TX<- factor(info_MS$Treatment.first.12.months)
data_cormat$TX<- unclass(data_cormat$TX)
data_cormat$BMI<- info_MS$BMI

data_cormat2<- data_cormat %>% relocate(Activity) # Change which variable to compare
variables <- colnames(data_cormat2)[2:ncol(data_cormat2)]
formulas <- list()
for (i in seq_along(variables)) {
  tmp <- combn(variables, i)
  tmp <- apply(tmp, 2, paste, collapse="+TX+Age+Sex+BMI+")
  tmp <- paste0("Activity~", tmp) ####Change this
  formulas[[i]] <- tmp
}

#Loop to make linear models of different combos of Igs to predict age

formulas <- unlist(formulas)
formulas <- sapply(formulas, as.formula)
models <- lapply(formulas, lm, data=data_cormat2)


#Extract p values from model
#Function to extract P value
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

reg_val<- lapply(models, lmp)
reg_val<- as.data.frame(t(reg_val))
reg_val<- as.data.frame(t(reg_val))
sig_rv<- which(reg_val$V1<0.0000015)

Combined_sigs<- (reg_val[sig_rv,])

model_test<- lm("Activity"~"IFNG_IL6"+"IgM_CCL22"+"CXCL10_CXCL9" +"TX"+"Age"+"Sex"+"BMI",
              data = data_cormat2)

#MSNA = IGM_CCL22  CCL11_CCL13 IL1B_CCL13
#MSA_1 = IgG1_CXCL10  CCL26_CCL22 CCL13_CXCL9 CCL11_CCL13
#MSA_2 = IFNG_IL6 IGM_CCL22 CXCL10_CXCL9 IL1B_CCL13

#MSNA_combo = CXCL10_CXCL9, CXCL9_TNFA

#COMBINED SIGS = IL4_CCL11 IgM_CCL22 IL1B_CCL13 CXCL10_CXCL9

#####################################
#10)AUROC and Survival Curve for correlations identified by network analysis
Act_sig_conns_p<- c("CXCL10_IgG1","IL4_CCL11","CCL13_IFNG", "CCL11_IL4", "IgM_CCL22",
                    "CCL11_CCL13", "IL1B_CCL13", "CCL11_CCL13", "CXCL10_CXCL9", "CCL26_CCL22",
                    "CCL13_CXCL9", "IgG1_NFL")

MSA_only_vd<- c("CCL26_IFNG","IFNG_IL1B","IFNG_IL4","CCL11_IL6","IgG1_CXCL10", "CCL11_CCL13",
"CCL26_CCL13", "IFNG_CCL13","IL1B_CCL13","CCL26_CCL22", "IL1B_CCL22" , "IL10_CCL22", 
"CCL26_CXCL9", "IL4_CXCL9","CCL13_CXCL9", "IgG1_TNFA")

MSNA_only_vd<- c("IFNG_IL6","IgM_CCL22","CXCL10_CXCL9" ,"CCL22_CXCL9","CXCL9_TNFA","IgG1_NFL")
view(data_combo)

index_Act<- data_combo

view(index_Act) 
use_me<-data.frame(t(index_Act))
use_me$Age<- info_MS$Age
use_me$Sex<- factor(info_MS$Sex)
use_me$Sex<- unclass(use_me$Sex)
use_me$Activity<- factor(info_MS$Activity)
use_me$Activity<- unclass(use_me$Activity)
use_me$TX<- factor(info_MS$Treatment.first.12.months)
use_me$TX<- unclass(use_me$TX)
use_me$BMI<- info_MS$BMI
use_me$OCB<- info_MS$OCB
use_me$act_status<-ifelse(info_MS$TTanyAct>12,"MSNA","MSA")

use_me2<- use_me %>% relocate(Activity) # Change which variable to compare

#COMBINED SIGS = IL4_CCL11(0.8652) IgM_CCL22( 0.1727) IL1B_CCL13( 0.4421) CXCL10_CXCL9

#Change connection of interest here
#IgG1_CXCL10, Age, CXCL10_CXCL9
holder_cp<- use_me2$IgG1_CXCL10
holder<- "IgG1_CXCL10"  
graph_name<- holder




#Determine cutoff for each protein 
cp<- cutpointr(use_me2, holder_cp,act_status,
               method = maximize_metric, metric = youden)

#Save statistics for each protein
summary(cp)
sink(paste(holder, "w.TREM2.txt")) ####change####
print(summary(cp))
sink()
plot(cp) #Save if wanted

#Use optimal cutpoint found above to create new binary variable
use_me2$level<- ifelse(holder_cp < cp$optimal_cutpoint ,"MSNA", "MSA" ) 

use_me2$level<- ifelse(  holder_cp < cp$optimal_cutpoint &  
                        use_me2$Age > 34.5, "MSNA", "MSA"  )
#Change
#Compare new variable to Activity to determine predictive content
use_me2$act_status<- factor(use_me2$act_status)
use_me2$level<- factor(use_me2$level)
#Create confusion matrix
example <- confusionMatrix(data=use_me2$level, reference = use_me2$act_status)
example
table <-data.frame(example$table)

plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))



mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Conf.Matr.Age", 
                          graph_name, ".tiff", sep = ""))

tiff(file=mypath,units="in", width=5, height=4, res=300, compression = 'lzw' )
ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "green", bad = "red")) +
  theme_bw() +
  xlim(rev(levels(table$Reference)))

dev.off()

#Survival Curve based on cutoffs
#Make new variables
#TTact = time to activity, in first 12 months
#status = anything that will not be used in SC
use_me2$TTAct<- as.numeric(info_MS$TTanyAct)
use_me2$TTAct<- ifelse(use_me2$TTAct>13, 13, use_me2$TTAct)
use_me2$status<-ifelse(info_MS$TTanyAct>12, 0,1)
use_me2$status<- as.numeric(use_me2$status)

fit<- survfit(Surv(TTAct,status) ~ level, data =use_me2)
coxph(Surv(TTAct,status) ~ level, data =use_me2)


survdiff(Surv(TTAct,status) ~ level, data =use_me2)


summary(fit)

holder_high<- paste0(holder, ">", round(cp$optimal_cutpoint, digits = 2),"<br/>", 
                     #"+ CXCL10_CXCL9<0.33", "<br/>",
                     #"+IFNG_IL6 >0.746","<br/>",
                     "+ Age <34.5")
holder_low<- paste0(holder,"<", round(cp$optimal_cutpoint, digits = 2), "<br/>",
                    #"+ CXCL10_CXCL9>0.33", "<br/>",
                    #"+IFNG_IL6 <0.746","<br/>",
                    "+ Age >34.5")

mypath <- file.path("C:","/Users/norawelsh/Desktop/", "092022 Networks FINAL",
                    paste("Survival.Curve.Age2",
                          graph_name, ".tiff", sep = ""))

tiff(file=mypath,units="in", width=8, height=7, res=300, compression = 'lzw' )
ggsurvplot(fit, data = use_me2, 
           conf.int = TRUE,censor = TRUE,
           pval = TRUE, risk.table = TRUE,
           legend.labs =
             c(holder_high, holder_low), legend = c("none"),
           ggtheme= theme_bw())
dev.off()








#Title: Ecological interactions without co-occurrence are not evidence of spillover risk: Predicting the impact of landscape conversion on small mammal co-occurrence and the transmission of zoonotic and epizootic diseases in the Atlantic Forest#
#Authors: Costa, APLC; Winck G. R.; Andreazzi, C. S.

# Load require packages
library(Hmsc)
library(MASS)
library(ape)
library(phytools)
library(tidyverse)
library(bipartite)
library(ggplot2)
library(rlang)
library(vioplot)
library(colorspace)
library(vegan)
library(picante)
library(lmPerm)
library(PerformanceAnalytics)

# Load data

phylotree=read.newick("./Filogenias/Smallmammals_maxlike_nf.nwk")
com.pa.f= read.csv("./Data/com.f_pa.csv")
land.feat.na.e.f= read.csv("./Data/land.feat.na.e.f.csv")
TrData= read.csv("./Data/TrData.csv", header=T, na.strings=c("","NA"))
geo= read.csv("./Data/buffer_id.csv", sep = ";")
net.all= read.csv("Data/Net_metrics.csv", sep = ";")
load("./Models/FittedModelElev_5_10_250_final.RData")

# Working on data

land.feat.na.e.f =land.feat.na.e.f %>% remove_rownames %>% column_to_rownames(var="X")
geo =geo %>% remove_rownames %>% column_to_rownames(var="b_id2")
TrData =TrData %>% remove_rownames %>% column_to_rownames(var="X") 
net.all = net.all %>% remove_rownames %>% column_to_rownames(var="X")

str(TrData)
TrData <- TrData %>%
  mutate(across(c("Activity", "TrophicLevel.Kissling", "MainGuild", "Locomotor.Paglia"), as.factor)) %>% 
  mutate(across(c("Comp.cauda", "comp.corpo", "massa.corp"), as.numeric))

# Match Data

corM = match(phylotree$tip.label, colnames(com.pa.f))
com.pa.f= com.pa.f[,corM]
all(phylotree$tip.label  == colnames(com.pa.f))

all(rownames(land.feat.na.e.f) == rownames(com.pa.f))
reorder = match(rownames(com.pa.f), rownames(land.feat.na.e.f))
land.feat.na.e.f= land.feat.na.e.f[reorder,]

all(rownames(geo) == rownames(com.pa.f))
reorder = match(rownames(com.pa.f), rownames(geo))
geo= geo[reorder,]

all(rownames(TrData) == colnames(com.pa.f))
reorder= match(colnames(com.pa.f), rownames(TrData))
TrData= TrData[reorder,]

xy= as.matrix(geo[,4:3])
rownames(xy)= rownames(geo)
plot(xy)

colnames(xy)=c("x-coordinate","y-coordinate")

dist(xy)

############################# Building Model ###################################

rL = HmscRandomLevel(sData=xy)

n=136

studyDesign = data.frame(sample = as.factor(row.names(xy)))

TrFormula.pa = ~ Ca_body + Activity
Xformula.pa= ~ forest_pland+ agriculture_pland+ mosaic_pland+ pasto_pland+forest_np + shdi + Elev_m1

m = Hmsc(Y = com.pa.f, XData = land.feat.na.e.f, XFormula = Xformula.pa,TrData = TrData, TrFormula = TrFormula.pa,  phyloTree = phylotree, studyDesign = studyDesign, ranLevels = list(sample = rL), distr = "probit",YScale=TRUE)

############################## Fitting Model ###################################

nChains = 5
thin = 10
samples = 250
transient = 50*thin
verbose = 0

model_pa = sampleMcmc(m, thin = thin, samples = samples, transient = transient,nChains = nChains, verbose = verbose)

############### Plot Convergence ############
samples_list = c(250)
thin_list = c(10) 
nst = length(thin_list)
nChains = 5

ma = NULL
na = NULL
for (Lst in 1:nst) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  nm = length(list(model_pa))
    mpost = convertToCodaObject(model_pa, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    if(is.null(ma)){
      ma=psrf.beta[,1]
      na = paste0(as.character(thin),",",as.character(samples))
    } else {
      ma = cbind(ma,psrf.beta[,1])
      if(j==1){
        na = c(na,paste0(as.character(thin),",",as.character(samples)))
      } else {
        na = c(na,"")
      }
    }
}

par(mfrow=c(2,1))

vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")

par(mfrow=c(1,1))

###################### Explanatory power ################################
preds.spatial = computePredictedValues(model_pa)
MF.spatial = evaluateModelFit(hM=model_pa, predY=preds.spatial)

mean(MF.spatial$TjurR2)
mean(MF.spatial$AUC)

Perfomance.exp= MF.spatial  %>% ggplot(aes(x=TjurR2, y= AUC)) +
  geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(size= 16,vjust = 0.5, hjust=1,), 
        axis.text.y = element_text(size= 16, vjust = 0.5, hjust=1, face = "bold"), 
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        title = element_text(size=16),
        legend.position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  labs(x= "Tjur’s R2", y="AUC")
Perfomance.exp

############### Predictive power ################
partition = createPartition(model_pa, nfolds = 2, column = "sample")
cvpreds.spatial = computePredictedValues(model_pa, partition=partition,
                                         nParallel = nChains)
cvMF.spatial = evaluateModelFit(hM=model_pa, predY=cvpreds.spatial)

mean(cvMF.spatial$TjurR2)
mean(cvMF.spatial$AUC)

Perfomance.pred= cvMF.spatial  %>% ggplot(aes(x=TjurR2, y= AUC)) +
  geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(size= 16,vjust = 0.5, hjust=1,), 
        axis.text.y = element_text(size= 16, vjust = 0.5, hjust=1, face = "bold"), 
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        title = element_text(size=16),
        legend.position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  labs(x= "Tjur’s R2", y="AUC")
Perfomance.pred

########################## Variation partitioning ##############################

groupnames = c("Forest", "Farming", "Fragmentation", "Diversity", "Elevation")
group = c(1,4,1,3,2,2,2,5) 
VPG = computeVariancePartitioning(model_pa, group = group, groupnames =groupnames)

group.vp= data.frame(VPG$vals)
group.vp = rownames_to_column(group.vp)
group.vp.plot= group.vp %>% group_by(rowname) %>%  gather("Species","Value", 2:53) %>% ungroup()
group.vp.plot$rowname= as.factor(group.vp.plot$rowname)

varplt = group.vp.plot  %>% ggplot(aes(x=Species, y= Value, fill=rowname)) +
  geom_bar(stat = "identity", width=0.5, position = position_stack(reverse = TRUE))+
  theme_classic()+
  scale_fill_manual(values=c("#d099a1", "#c07860", "#af7c16", "#1a5a1f", "#647718", "#cfb08c"), labels=c('Landscape diversity = 8.5', 'Elevation = 16.1','Farming = 29.6', 'Forest = 27.9', 'Fragmentation = 4.6', 'Random:sample = 13.3'))+
  theme(axis.text.x = element_text(size= 16, angle = 90, vjust = 0.5, hjust=1, face = "italic"), 
        axis.text.y = element_text(size= 16, vjust = 0.5, hjust=1, face = "bold"), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        title = element_text(size=16),
        legend.position = "top")+
  guides(fill = guide_legend(nrow = 1))+
  labs(fill='Landscape Variables Means', title= 'Variance Partitioning') 
varplt

# Traits explanation
VP = computeVariancePartitioning(model_pa)

Total.var.traits= VP$R2T$Y #proportion of variation in species occurrence that the traits explain
ExpVar.traits= data.frame(By.covar= VP$R2T$Beta, Total.var=Total.var.traits) #how much of the variation they explain among the responses of the species to their covariates.

# Phylogenetic signal in species niches

mpost = convertToCodaObject(model_pa)
round(summary(mpost$Rho, quantiles = c(0.025, 0.5, 0.975))
      [[2]],2)
round(summary(mpost$Beta, quantiles = c(0.025, 0.5, 0.975))
      [[2]],2)
summary(mpost$Rho)
hist(mpost[["Rho"]][[3]], xlab = "psrf (Omega)")

####################### Cooccurrence of Small mammals ##########################

OmegaCor = computeAssociations(model_pa)
supportLevel = 0.75
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
toPlot = sign(toPlot)

#######################Cooccurrence in function of modules#######################
str(net.all)
# Melt the coocurrence data frame
melted_df <- reshape2::melt(toPlot, id.vars = "row_names")
melted_df2<-melted_df %>% 
  dplyr::rename("Species"=Var1)

melted_df2<- left_join(melted_df2, net.all)
str(melted_df2)

melted_df2<-melted_df2 %>% 
  dplyr::select(Species,Var2, value, Modules)  %>% 
  dplyr::rename("Var1"="Species", "Species"=Var2, Modules_1="Modules")

melted_df2<- left_join(melted_df2, net.all)

melted_df2<-melted_df2 %>% 
  dplyr::select(Species,Var1,value, Modules_1, Modules)  %>% 
  dplyr::rename("Var2"="Species", "Modules_2"=Modules)

melted_group <- melted_df2 %>% 
  dplyr::select(value, Modules_1, Modules_2)  %>% 
  group_by(Modules_1, Modules_2) %>% 
  summarise_all(mean) %>%  
  pivot_wider(names_from = Modules_2, values_from = value) %>% 
  column_to_rownames(var="Modules_1")

melted_group<- as.matrix(melted_group)

# Cooccurrence plot
par(mfrow=c(1,1))
png("./Figures_elevmodel/CorrplotPA_bymod.tiff",  width=20, height=20, units="cm", res=600)

corrplot(melted_group, col.lim = c(-1,1), tl.col= "black", tl.cex=1,cl.cex = 0.5, col=colorRampPalette(c("blue", "white", "red"))(255), type= "lower", diag = TRUE, cl.pos = "r")

dev.off()

############ Phylogenetic distance inside and ouside each module ###############
# Working on data
corM = match(phylotree$tip.label, net.all$Species)
net.all= net.all[corM,]
all(phylotree$tip.label  == net.all$Species)

# Phylogenetic distance 

phydist<-cophenetic(phylotree)
phydist[1:3,1:4] 
net.all= rownames_to_column(net.all)
module_vec= net.all %>% select(Species, Modules)

### Building dataframe
melted_phy <- reshape2::melt(phydist, id.vars = "row_names")

melted_phy<-melted_phy %>% 
  dplyr::rename("Species"=Var1)

melted_phy<- left_join(melted_phy, module_vec)
str(melted_phy)

melted_phy2 <-melted_phy %>% 
  dplyr::select(Species,Var2, value, Modules)  %>% 
  dplyr::rename("Var1"="Species", "Species"=Var2, Module_1="Modules")

melted_phy2<- left_join(melted_phy2, module_vec)

melted_phy2<-melted_phy2 %>% 
  dplyr::select(Species,Var1,value, Module_1, Modules)  %>% 
  dplyr::rename("Var2"="Species", "Module_2"=Modules) %>% filter(Module_1 != "no registry", Module_2 != "no registry")

melted_phy2 = melted_phy2 %>%   mutate(Group = ifelse(Module_1 == Module_2, as.character(Module_1), "outside"))

phy_group <- melted_phy2 %>% unite("Group_module", Module_1, Module_2, sep="_") %>%  unite("Pairwise", Var1, Var2, sep= "_") 

phy_group$sorted_pairwise_species <- apply(phy_group["Pairwise"], 1, function(x) {
  paste(sort(unlist(strsplit(x, "_"))), collapse = "_")
})

# remove duplicates based on the new sorted_pairwise_species column
phy_group <- phy_group %>%
  distinct(sorted_pairwise_species, .keep_all = TRUE) %>%
  select(-sorted_pairwise_species) # Remove the temporary column


# Analyse the difference in pairwise distance inside the module and outiside modules
hist(phy_group$value)
shapiro.test(phy_group$value)

phy.test= lmPerm::lmp(value~ Group, data = phy_group)

summary(phy.test)

# Plot results

Phy_mod_plt = phy_group %>%  ggplot(aes(x = Group, y= value, fill= Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette="Set2")+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(size= 16, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size= 16, vjust = 0.5, hjust=1), 
        axis.title.y = element_text(size= 16, face = "bold"),
        axis.title.x = element_text(size= 16, face = "bold"),
        title = element_text(size=16),
        legend.position = "none")+
  labs(x ='Species pairs', y= 'Pairwise phylogenetic distance') 
Phy_mod_plt

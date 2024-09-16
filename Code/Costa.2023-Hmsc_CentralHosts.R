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
library(devtools)
library(cowplot)
library(Rmisc)
library(colorspace)
library(abind)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(rstatix)
library(FSA)

################ Load data ##############
load("FittedModelElev_5_10_250_final.RData")
net.all= read.csv("./Data/Net_metrics.csv", sep = ";")
central_res = read.csv("./Data/central_res.csv")

c_res <- central_res %>% 
  dplyr::select(rowname, betweenness, residuals, Central_effort) %>% 
  dplyr::rename(Species = rowname) %>% 
  full_join(net.all) %>% 
  mutate(betweenness = case_when(is.na(betweenness) ~ 0,
                                 .default = as.numeric(betweenness)),
         residuals= case_when(is.na(residuals)~0,
                              .default = as.numeric(residuals)),
         Central_effort = case_when(is.na(Central_effort) ~ "unrecorded", .default = as.character(Central_effort)),
         Central_res2 = case_when(betweenness >= 200 ~ "Central", .default = as.character(Central_effort))) %>% 
  dplyr::select(Species, betweenness, residuals, Central_effort, Central_hosts, Central_res2)

############################### Predictions ####################################
head(model_pa$XFormula)

# Make gradients for each landscape metric

Gradient1 = constructGradient(model_pa,focalVariable = "forest_pland")
Gradient2 = constructGradient(model_pa,focalVariable = "agriculture_pland")
Gradient3 = constructGradient(model_pa,focalVariable = "mosaic_pland")
Gradient4 = constructGradient(model_pa,focalVariable = "pasto_pland")
Gradient5 = constructGradient(model_pa,focalVariable = "shdi")
Gradient6 = constructGradient(model_pa,focalVariable = "forest_np")
Gradient7 = constructGradient(model_pa,focalVariable = "Elev_m1")

# Calculating probabilities based on each gradient
predY1 = predict(model_pa, Gradient=Gradient1, expected = TRUE)
predY2 = predict(model_pa, Gradient=Gradient2, expected = TRUE)
predY3 = predict(model_pa, Gradient=Gradient3, expected = TRUE)
predY4 = predict(model_pa, Gradient=Gradient4, expected = TRUE)
predY5 = predict(model_pa, Gradient=Gradient5, expected = TRUE)
predY6 = predict(model_pa, Gradient=Gradient6, expected = TRUE)
predY7 = predict(model_pa, Gradient=Gradient7, expected = TRUE)

############################ Forest_pland data ################################
# Summarizing predictions
NApredY= as.data.frame(Reduce("+",predY1)/length(predY1))
names= colnames(NApredY)

# Support level of prediction
tmp = abind(predY1, along = 3)
xx = Gradient1$XDataNew[, 1]
ngrid = length(xx)

# Selection of species with a positive relationship
p.sup.nat= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] > tmp[1, i, ])
  print(Pr)
  p.sup.nat= rbind(p.sup.nat, Pr)
}
rownames(p.sup.nat)<-names
colnames(p.sup.nat)<- 'Forest_p.pred'
p.sup.nat= rownames_to_column(p.sup.nat)
p.sup.nat=data.frame(p.sup.nat[p.sup.nat$Forest_p.pred>=0.8,])

# Selection of species with a negative relationship
n.sup.nat= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] < tmp[1, i, ])
  print(Pr)
  n.sup.nat= rbind(n.sup.nat, Pr)
}

rownames(n.sup.nat)<-names
colnames(n.sup.nat)<- 'Forest_p.pred'
n.sup.nat= rownames_to_column(n.sup.nat)
n.sup.nat=as.data.frame(n.sup.nat[n.sup.nat$Forest_p.pred>=0.8,])

# Data species by Nat Areas
NApredY.fp= NApredY[,colnames(NApredY) %in% p.sup.nat$rowname, drop=FALSE]
NApredY.fn= NApredY[,colnames(NApredY) %in% n.sup.nat$rowname, drop=FALSE]
NApredY.fp$G.value = Gradient1$XDataNew$forest_pland
NApredY.fn$G.value = Gradient1$XDataNew$forest_pland

PredForestP= NApredY.fp %>% gather("Species", "prob", 1:32) 
str(PredForestP)
PredForestN= NApredY.fn %>% gather("Species", "prob", 1:11) 
str(PredForestN)
PredForestN$type= "Negative"
PredForestP$type= "Positive"

PredForest= full_join(PredForestP,PredForestN)

########################### Agriculture data ###################################

ApredY= as.data.frame(Reduce("+",predY2)/length(predY2))

# Support level of prediction
tmp = abind(predY2, along = 3)
xx = Gradient2$XDataNew[, 1]
ngrid = length(xx)

# Selection of species with a positive relationship

p.sup.ag= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] > tmp[1, i, ])
  print(Pr)
  p.sup.ag= rbind(p.sup.ag, Pr)
}

rownames(p.sup.ag)<-names
colnames(p.sup.ag)<- 'Agr.pred'
p.sup.ag= rownames_to_column(p.sup.ag)
p.sup.ag=as.data.frame(p.sup.ag[p.sup.ag$Agr.pred>=0.8,])

# Selection of species with a negative relationship

n.sup.ag= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] < tmp[1, i, ])
  print(Pr)
  n.sup.ag= rbind(n.sup.ag, Pr)
}

rownames(n.sup.ag)<-names
colnames(n.sup.ag)<- 'Agr.pred'
n.sup.ag= rownames_to_column(n.sup.ag)
n.sup.ag=as.data.frame(n.sup.ag[n.sup.ag$Agr.pred>=0.8,])

# Data species by Agriculture

ApredY.p= ApredY[,colnames(ApredY) %in% p.sup.ag$rowname, drop=FALSE]
ApredY.n= ApredY[,colnames(ApredY) %in% n.sup.ag$rowname, drop=FALSE]

ApredY.p$G.value = Gradient2$XDataNew$agriculture_pland
ApredY.n$G.value = Gradient2$XDataNew$agriculture_pland

PredAgriP= ApredY.p %>% gather("Species", "prob", 1:10) 
PredAgriN= ApredY.n %>% gather("Species", "prob", 1:27) 
PredAgriN$type= "Negative"
PredAgriP$type= "Positive"
PredAgriculture= full_join(PredAgriP,PredAgriN)

############################## Mosaic data #####################################

MpredY= as.data.frame(Reduce("+",predY3)/length(predY3))

# Support level of prediction
tmp = abind(predY3, along = 3)
xx = Gradient3$XDataNew[, 1]
ngrid = length(xx)

# Selection of species with a positive relationship

p.sup.mo= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] > tmp[1, i, ])
  print(Pr)
  p.sup.mo= rbind(p.sup.mo, Pr)
}

rownames(p.sup.mo)<-names
colnames(p.sup.mo)<- 'Mo.pred'
p.sup.mo= rownames_to_column(p.sup.mo)
p.sup.mo=as.data.frame(p.sup.mo[p.sup.mo$Mo.pred>=0.8,])

# Selection of species with a negative relationship

n.sup.mo= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] < tmp[1, i, ])
  print(Pr)
  n.sup.mo= rbind(n.sup.mo, Pr)
}

rownames(n.sup.mo)<-names
colnames(n.sup.mo)<- 'Mo.pred'
n.sup.mo= rownames_to_column(n.sup.mo)
n.sup.mo=as.data.frame(n.sup.mo[n.sup.mo$Mo.pred>=0.8,])

# Data species by Mosaic of uses

MpredY.p= MpredY[,colnames(MpredY) %in% p.sup.ag$rowname, drop=FALSE]
MpredY.n= MpredY[,colnames(MpredY) %in% n.sup.ag$rowname, drop=FALSE]

MpredY.p$G.value = Gradient3$XDataNew$mosaic_pland
MpredY.n$G.value = Gradient3$XDataNew$mosaic_pland

PredMoP= MpredY.p %>% gather("Species", "prob", 1:10) 
str(PredMoP)
PredMoN= MpredY.n %>% gather("Species", "prob", 1:27) 
str(PredMoN)
PredMoN$type= "Negative"
PredMoP$type= "Positive"
PredMosaic= full_join(PredMoP,PredMoN)

############################### Pasture data ###################################

PpredY= as.data.frame(Reduce("+",predY4)/length(predY4))

# Support level of prediction
tmp = abind(predY4, along = 3)
xx = Gradient4$XDataNew[, 1]
ngrid = length(xx)

# Selection of species with a positive relationship

p.sup.pa= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] > tmp[1, i, ])
  print(Pr)
  p.sup.pa= rbind(p.sup.pa, Pr)
}

rownames(p.sup.pa)<-names
colnames(p.sup.pa)<- 'Past.pred'
p.sup.pa= rownames_to_column(p.sup.pa)
p.sup.pa=as.data.frame(p.sup.pa[p.sup.pa$Past.pred>=0.8,])

# Selection of species with a negative relationship

n.sup.pa= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] < tmp[1, i, ])
  print(Pr)
  n.sup.pa= rbind(n.sup.pa, Pr)
}

rownames(n.sup.pa)<-names
colnames(n.sup.pa)<- 'Past.pred'
n.sup.pa= rownames_to_column(n.sup.pa)
n.sup.pa=as.data.frame(n.sup.pa[n.sup.pa$Agr.pred>=0.8,])

# Data species by Pasture

PpredY.p= PpredY[,colnames(PpredY) %in% p.sup.pa$rowname, drop=FALSE]
PpredY.n= PpredY[,colnames(PpredY) %in% n.sup.pa$rowname, drop=FALSE]

PpredY.p$G.value = Gradient4$XDataNew$pasto_pland
PpredY.n$G.value = Gradient4$XDataNew$pasto_pland

PredPP= PpredY.p %>% gather("Species", "prob", 1:15) 
str(PredPP)
PredPP$type= "Positive"

######################### Shannon diversity Data ################################

DpredY= as.data.frame(Reduce("+",predY5)/length(predY5))

# Support level of prediction
tmp = abind(predY5, along = 3)
xx = Gradient5$XDataNew[, 1]
ngrid = length(xx)

# Selection of species with a positive relationship

p.sup.div= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] > tmp[1, i, ])
  print(Pr)
  p.sup.div= rbind(p.sup.div, Pr)
}

rownames(p.sup.div)<-names
colnames(p.sup.div)<- 'Div.pred'
p.sup.div= rownames_to_column(p.sup.div)
p.sup.div=as.data.frame(p.sup.div[p.sup.div$Div.pred>=0.8,])

# Selection of species with a negative relationship

n.sup.div= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] < tmp[1, i, ])
  print(Pr)
  n.sup.div= rbind(n.sup.div, Pr)
}

rownames(n.sup.div)<-names
colnames(n.sup.div)<- 'Div.pred'
n.sup.div= rownames_to_column(n.sup.div)
n.sup.div=as.data.frame(n.sup.div[n.sup.div$Div.pred>=0.8,])

# Data species by Shannon evenness
DpredY.p= DpredY[,colnames(DpredY) %in% p.sup.div$rowname, drop=FALSE]
DpredY.n= DpredY[,colnames(DpredY) %in% n.sup.div$rowname, drop=FALSE]

DpredY.p$G.value = Gradient5$XDataNew$shdi
DpredY.n$G.value = Gradient5$XDataNew$shdi

PredDivP= DpredY.p %>% gather("Species", "prob", 1:13)
PredDivN= DpredY.n %>% gather("Species", "prob", 1:28) 
PredDivN$type= "Negative"
PredDivP$type= "Positive"
PredDiversity= full_join(PredDivP,PredDivN)

####################### Number of forest patches data ##########################

NPpredY= as.data.frame(Reduce("+",predY6)/length(predY6))

# Support level of prediction
tmp = abind(predY6, along = 3)
xx = Gradient6$XDataNew[, 1]
ngrid = length(xx)

# Selection of species with a positive relationship

p.sup.frag= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] > tmp[1, i, ])
  print(Pr)
  p.sup.frag= rbind(p.sup.frag, Pr)
}

rownames(p.sup.frag)<-names
colnames(p.sup.frag)<- 'Frag.pred'
p.sup.frag= rownames_to_column(p.sup.frag)
p.sup.frag=as.data.frame(p.sup.frag[p.sup.frag$Frag.pred>=0.8,])

# Selection of species with a negative relationship

n.sup.frag= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] < tmp[1, i, ])
  print(Pr)
  n.sup.frag= rbind(n.sup.frag, Pr)
}

rownames(n.sup.frag)<-names
colnames(n.sup.frag)<- 'Frag.pred'
n.sup.frag= rownames_to_column(n.sup.frag)
n.sup.frag=as.data.frame(n.sup.frag[n.sup.frag$Frag.pred>=0.8,])

# Data species by number of forest patches
NPpredY.p= NPpredY[,colnames(NPpredY) %in% p.sup.frag$rowname, drop=FALSE]
NPpredY.n= NPpredY[,colnames(NPpredY) %in% n.sup.frag$rowname, drop=FALSE]

NPpredY.p$G.value = Gradient6$XDataNew$forest_np
NPpredY.n$G.value = Gradient6$XDataNew$forest_np

PredNPP= NPpredY.p %>% gather("Species", "prob", 1:14) 
PredNPN= NPpredY.n %>% gather("Species", "prob", 1:23) 
PredNPN$type= "Negative"
PredNPP$type= "Positive"
PredFragmentation= full_join(PredNPP,PredNPN)

############################ Elevation data ####################################

EpredY= as.data.frame(Reduce("+",predY7)/length(predY7))

# Seeing support level of prediction
tmp = abind(predY7, along = 3)
xx = Gradient7$XDataNew[, 1]
ngrid = length(xx)

# Selection of species with a positive relationship

p.sup.elev= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] > tmp[1, i, ])
  print(Pr)
  p.sup.elev= rbind(p.sup.elev, Pr)
}

rownames(p.sup.elev)<-names
colnames(p.sup.elev)<- 'Elev.pred'
p.sup.elev= rownames_to_column(p.sup.elev)
p.sup.elev=as.data.frame(p.sup.elev[p.sup.elev$Elev.pred>=0.8,])

# Selection of species with a negative relationship

n.sup.elev= data.frame(NULL)
for (i in 1:52) {
  Pr = mean(tmp[ngrid, i, ] < tmp[1, i, ])
  print(Pr)
  n.sup.elev= rbind(n.sup.elev, Pr)
}

rownames(n.sup.elev)<-names
colnames(n.sup.elev)<- 'Elev.pred'
n.sup.elev= rownames_to_column(n.sup.elev)
n.sup.elev=as.data.frame(n.sup.elev[n.sup.elev$Elev.pred>=0.8,])

# Data species by Elevation

EpredY.p= EpredY[,colnames(EpredY) %in% p.sup.elev$rowname, drop=FALSE]
EpredY.n= EpredY[,colnames(EpredY) %in% n.sup.elev$rowname, drop=FALSE]

EpredY.p$G.value = Gradient7$XDataNew$Elev_m1
EpredY.n$G.value = Gradient7$XDataNew$Elev_m1

PredEP= EpredY.p %>% gather("Species", "prob", 1:34) 
PredEN= EpredY.n %>% gather("Species", "prob", 1:10)
PredEN$type= "Negative"
PredEP$type= "Positive"
PredElevation= full_join(PredEP,PredEN)

########################## Aggregating data ####################################

PredForest$Gradient= "Forest"
PredAgriculture$Gradient = "Agriculture"
PredMosaic$Gradient = "Mosaic"
PredDiversity$Gradient = "Diversity"
PredFragmentation$Gradient = "Forest NP"
PredPP$Gradient = "Pasture"
PredElevation$Gradient= "Elevation"

Predictions.pa= full_join(PredForest, PredAgriculture)
Predictions.pa= full_join(Predictions.pa, PredMosaic)
Predictions.pa= full_join(Predictions.pa, PredDiversity)
Predictions.pa= full_join(Predictions.pa, PredFragmentation)
Predictions.pa= full_join(Predictions.pa, PredPP)
Predictions.pa= full_join(Predictions.pa, PredElevation)
Predictions.pa= full_join(Predictions.pa, c_res)
Predictions.pa = Predictions.pa %>% dplyr::filter(Species != "caluromys.philander")
unique(Predictions.pa$Gradient)

############################ Gam Model #########################################
central_ef <- central_res %>% 
  select(Species, sample.effort)

Pred_gam = Predictions.pa %>%
  mutate(Gradient = case_when(
    Gradient== "Agriculture" ~ "Farming",
    Gradient == "Pasture"~ "Farming", 
    Gradient == 'Mosaic' ~ "Farming",
    .default = as.character(Gradient))) %>% 
  mutate(Central_hosts = as.factor(Central_hosts)) %>% 
  left_join(central_ef) %>% 
  mutate(sample.effort = case_when(is.na(sample.effort) ~ 0, .default = as.numeric(sample.effort)))

theme_set(theme_classic(base_size = 12, base_family = 'serif') +
            theme(panel.spacing = unit(0, 'lines'),
                  axis.title.y = element_blank(),
                  legend.direction =  "horizontal", 
                  legend.position = "none",
                  legend.title = element_blank()))

pred= unique(Pred_gam$Gradient)

print(pred)

gam_res = NULL

for (i in pred) {
  Pred <- Pred_gam %>% 
    filter(Gradient == paste(i)) %>% 
    filter(prob>= 0.02)
  
  mN= gam(prob ~ s(G.value, by = Central_hosts)+ Central_hosts + s(sample.effort, bs = "re"), data = Pred, family= betar)
  
  sum_n = summary(mN)
  
  print(sum_n)
  
  max_length <- max(length(sum_n$p.coef), length(names(sum_n$p.coeff)), length(sum_n$se))
  
  df_n= data.frame(Coefficient = c(sum_n$p.coef, rep(NA, max_length - length(sum_n$p.coef))),
                   coef= c(names(sum_n$p.coeff), rep(NA, max_length - length(names(sum_n$p.coeff)))),
                   se = sum_n$se, 
                   "p-value" = c(sum_n$p.pv, rep(NA, max_length - length(sum_n$p.pv))), 
                   "Smooth_Gvalue"= c(sum_n$edf, rep(NA, max_length - length(sum_n$edf))), 
                   "p-value_sm" = c(sum_n$s.pv, rep(NA, max_length - length(sum_n$s.pv))), 
                   dev.explain = c(sum_n$dev.expl, rep(NA, max_length - length(sum_n$dev.expl))),
                   Grad= paste0(i))
  
  gam_res=rbind(gam_res, df_n)
  
  var_name <- paste0("gam_", i)
  
  assign(var_name, mN)
  
  a<- plot_predictions(mN, condition = c('G.value', 'Central_hosts'),
                       type = 'link')+
    #coord_cartesian(ylim = c(0, 0.6))+
    scale_color_manual(values=c("#954b3c", "#367500", "#609db1"))
  b<- plot_predictions(mN,
                       condition = c('Central_hosts','Central_hosts'),
                       type = 'response')
  
  var_name <- paste0("gplot_", i)
  
  assign(var_name, a)
  
  var_name <- paste0("Cboxplot_", i)
  
  assign(var_name, b)
}

############################ Plot Results ###################################### 

gplot_Diversity <- gplot_Diversity+
  labs(y = "Occurrence probability",
       x= "Landscape diversity variation")

gplot_Elevation <- gplot_Elevation +
  labs(y = "Occurrence probability",
       x= "Elevation variation")

gplot_Farming <- gplot_Farming +
  labs(y = "Occurrence probability",
       x= "% of Farming variation")

gplot_Forest <- gplot_Forest +
  labs(y = "Occurrence probability",
       x= "% of Forest variation")

gplot_ForestNP <- `gplot_Forest NP`+
  labs(y = "Occurrence probability",
       x= "Number of patches variation")

grid = plot_grid(
  gplot_Forest, gplot_Elevation, gplot_Farming, gplot_ForestNP, gplot_Diversity,
  nrow = 5,
  ncol = 1
)

# Add a common y-axis label
y_label <- ggdraw() + 
  draw_label("Occurrence probability (link scale)", angle = 90, fontfamily  = 'serif')

# Combine the grid and the y-axis label
final_plot1 <- plot_grid(y_label, grid, ncol = 2, rel_widths = c(0.05, 1))

# Display the final plot
print(final_plot1)

png(filename="./Figures_elevmodel/Prob.pred.smooth_gam.tiff", width=10, height=25, units="cm", res=600)
plot(final_plot)
dev.off()

cb_div<-Cboxplot_Diversity+
  coord_cartesian(ylim = c(0, 0.35))+
  labs(y = "",
    x= "")+
  scale_x_discrete(labels = c("Central"= "Central", "Non-Central" = "Non-central", "unrecorded"= "with no registry"))+
  scale_color_manual(values=c("Central" = "#954b3c","Non-Central" = "#367500","unrecorded"=  "#609db1"))+
  theme_bw(base_size = 12, base_family = 'serif') +
  theme(panel.spacing = unit(0, 'lines'),
        legend.direction =  "horizontal", 
        legend.position = "none",
        legend.title = element_blank())

cb_for<-Cboxplot_Forest+
  coord_cartesian(ylim = c(0, 0.35))+
  labs(y = "",
    x= "")+
  scale_x_discrete(labels = c("Central"= "Central", "Non-Central" = "Non-central", "unrecorded"= "with no registry"))+
  scale_color_manual(values=c("Central" = "#954b3c","Non-Central" = "#367500","unrecorded"=  "#609db1"))+
  theme_bw(base_size = 12, base_family = 'serif') +
  theme(panel.spacing = unit(0, 'lines'),
        legend.direction =  "horizontal", 
        legend.position = "none",
        legend.title = element_blank())

cb_el<-Cboxplot_Elevation+
  coord_cartesian(ylim = c(0, 0.35))+
  labs(y = "",
    x= "")+
  scale_x_discrete(labels = c("Central"= "Central", "Non-Central" = "Non-central", "unrecorded"= "with no registry"))+
  scale_color_manual(values=c("Central" = "#954b3c","Non-Central" = "#367500","unrecorded"=  "#609db1"))+
  theme_bw(base_size = 12, base_family = 'serif') +
  theme(panel.spacing = unit(0, 'lines'),
        legend.direction =  "horizontal", 
        legend.position = "none",
        legend.title = element_blank())

cb_far<-Cboxplot_Farming+
  coord_cartesian(ylim = c(0, 0.35))+
  labs(y = "",
    x= "")+
  scale_x_discrete(labels = c("Central"= "Central", "Non-Central" = "Non-central", "unrecorded"= "with no registry"))+
  scale_color_manual(values=c("Central" = "#954b3c","Non-Central" = "#367500","unrecorded"=  "#609db1"))+
  theme_bw(base_size = 12, base_family = 'serif') +
  theme(panel.spacing = unit(0, 'lines'),
        legend.direction =  "horizontal", 
        legend.position = "none",
        legend.title = element_blank())

cb_np<-`Cboxplot_Forest NP`+
  coord_cartesian(ylim = c(0, 0.35))+
  labs(y = "",
    x= "")+
  scale_x_discrete(labels = c("Central"= "Central", "Non-Central" = "Non-central", "unrecorded"= "with no registry"))+
  scale_color_manual(values=c("Central" = "#954b3c","Non-Central" = "#367500","unrecorded"=  "#609db1"))+
  theme_bw(base_size = 12, base_family = 'serif') +
  theme(panel.spacing = unit(0, 'lines'),
        legend.direction =  "horizontal", 
        legend.position = "none",
        legend.title = element_blank())

cb_grid = plot_grid(
  cb_for, cb_el, cb_far, cb_np, cb_div,
  nrow = 5,
  ncol = 1
)

# Add a common y-axis label
y_label2 <- ggdraw() + 
  draw_label("Occurrence probability", angle = 90, fontfamily  = 'serif')

# Combine the grid and the y-axis label
final_plot <- plot_grid(y_label, grid, cb_grid, ncol = 3, rel_widths = c(0.05,1, 1))

# Display the final plot
print(final_plot)

png(filename="./Figures_elevmodel/Prob.pred.cent_gam.tiff", width=18, height=25, units="cm", res=600)
plot(final_plot)
dev.off()

write.csv(gam_res, "./Results_Elevmodel/gam_cent_type.csv")

####################################################################################################################
####################################################################################################################
####################################################################################################################
#   Pipeline
####################################################################################################################
####################################################################################################################
####################################################################################################################

####################################################################################################################
# WORKING DIRECTORY
####################################################################################################################

rm(list=ls())
gc()

setwd("/Users/celiatalavan/Desktop/Taller VICA-2025/2025")

####################################################################################################################
# PACKAGES
####################################################################################################################


# Mandatory packages
required_pkgs <- c("INLA", "spdep", "DCluster", "DT", "openxlsx", "sf","ggplot2")
install_if_missing <- function(pkgs){
  for(pkg in pkgs){
    if(!requireNamespace(pkg, quietly = TRUE)){
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}
# INLA
if(!requireNamespace("INLA", quietly = TRUE)){
  install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable")
}
install_if_missing(setdiff(required_pkgs, "INLA"))

library(INLA)
library(spdep)
library(DCluster)
library(DT)
library(openxlsx)
library(sf)
library(RColorBrewer)
library(ggplot2)

####################################################################################################################
# Step 0: Data
####################################################################################################################

# Shapefile
mapa <- st_read("datos_ejemplo_practico_taller/mapa.shp", quiet = TRUE)

# Data
data_input <- read.xlsx("datos_ejemplo_practico_taller/datos.xlsx",sheet=1)


# Shapefile
# mapa <- st_read("map_example/map_example.shp", quiet = TRUE)

# Datos
# data_input <- read.csv("data_example.csv",header = TRUE, stringsAsFactors = FALSE)


#Check 1: verify the expected columns in the CSV: ID, Obs, population, Exp (adjust if they have different names)
required_cols <- c("ID", "Obs", "population", "Exp")
if(!all(required_cols %in% names(data_input))){
  paste0("Los datos debe contener columnas: ", paste(required_cols, collapse = ", "))
}

# Check 2: Make sure there is an ID field in both the map and the data_input.
# If they are not the same, for example cod_ine names(mapa)[which(names(mapa)%in%"cod_ine")]<-"ID"
shp_id_col <- "ID"
if(!shp_id_col %in% names(mapa)){
  paste0("El mapa no contiene la columna '", shp_id_col, "'. Revisa nombres: ", paste(names(mapa), collapse = ", "))
}


# Check 3: Same IDs column in mapa and data_input
length(intersect(mapa$ID,data_input$ID))
dim(mapa)
dim(data_input)


####################################################
# Creation of the necessary variables and objects
####################################################

# Create variable "region" (sequential integer number) in the data_input
data_input$"region" <- seq_len(nrow(data_input))

# Match by ID: use merge to avoid incorrect reordering
# We will preserve the shapefile order and add columns from data_input

map_df <- as.data.frame(mapa)
union <- merge(map_df, data_input, by = "ID", sort = FALSE, all.x = TRUE)


# Create an index by ID for reordering:
order_idx <- match(map_df$ID, union$ID)
union <- union[order_idx, ]


# Add columns to map
mapa$Observed  <- as.numeric(union$Obs)
mapa$Population <- as.numeric(union$population)
mapa$Expected <- as.numeric(union$Exp)

# Sort data_input to match the order of the map
data_input_ordered <- union[, c(names(data_input))]
# If merge added suffixes, fix them:
names(data_input_ordered)[names(data_input_ordered) == "Obs"] <- "Obs"
data_input <- data_input_ordered

rm(data_input_ordered,union,order_idx,shp_id_col,required_cols,map_df,install_if_missing,required_pkgs)


####################################################################################################################
# Step 1: Spatial neighbours
####################################################################################################################

neigs <- poly2nb(mapa)
x.nb <- neigs

nb2INLA(file="datos_ejemplo_practico_taller/neigs.graph",nb=x.nb)

# Check polygons without neighbours 

plot(st_geometry(mapa), col = "lightblue", border = "black", main = "Map with IDs")
centros <- st_centroid(mapa$geometry[c(1,2)])
text(st_coordinates(centros), labels =mapa$"ID"[c(1,2)], cex = 0.7, col = "red")



##################################################################################################################
# Step 2: Global spatial tests
##################################################################################################################

# Tests of overdispersion: Chi-square test 
chtest <- achisq.test(Obs~offset(log(Exp)),data=data_input,model="multinom",R=999)

chtest

# Tests of overdispersion: Potthoff-Whittinghill's test
pwtest <- pottwhitt.test(Obs~offset(log(Exp)),data=data_input,model="multinom",R=999)

pwtest

# Test global clustering: Moran I Test

longi<-list()
for(i in 1:length(x.nb)){
	
	
	longi[[i]]<-length(x.nb[[i]])
	
	
}
table(unlist(longi))
sum(as.numeric(table(unlist(longi)))) 

col.W <-nb2listw(neighbours=x.nb, zero.policy=TRUE)

resul.moranI.test<-moranI.test(Observed~offset(log(Expected)), 
			as.data.frame(mapa), "negbin", 999,
      listw=col.W, n=length(x.nb), S0=Szero(col.W) )

resul.moranI.test

# Test global clustering: Tango’s test

st_crs(mapa)

spa_units <- st_centroid(mapa)

spa_units_utm <- st_transform(spa_units, 25830) #  Convertir a UTM adecuado (España peninsular → EPSG 25830)

coords <- st_coordinates(spa_units_utm)

dlist <- dnearneigh(coords, 0, Inf)
dlist <- include.self(dlist)

dlist.d <- nbdists(dlist, coords)

col.W.tango <- nb2listw(
  dlist,
  glist = lapply(dlist.d, function(x) exp(-x)),
  style = "C"
) # Create weight (kernel exponencial)

resul.tango.test <- tango.test(
  Obs ~ offset(log(Exp)),
  data_input,
  model = "poisson",
  R = 99,
  listw = col.W.tango,
  zero.policy = TRUE
)

resul.tango.test


##################################################################################################################
# Step 3: Detection of the location of clusters
##################################################################################################################

mapa$"x" <- coords[,1]
mapa$"y" <- coords[,2]


### Prepare data for DCluster
df <- st_drop_geometry(mapa)  # equivalente a as(...,"data.frame")
df <- cbind(df, data_input[,c("Obs","Exp","population","region")])

### Maximum Likelihood Estimation (MLE) of the model
mle <- calculate.mle(df, model="negbin")

### Grid required for opgam (cordenates in meters)
thegrid <- df[,c("x","y")]

### Run opgam
set.seed(10)
knresults <- opgam(
  data=df,
  thegrid=thegrid,
  alpha=0.05,
  iscluster=kn.iscluster,
  fractpop=0.15,
  R=99,
  model="negbin",
  mle=mle
)

### Sort clusters by importance
knresults <- knresults[order(knresults$statistic, decreasing=TRUE),]
knresults$cluster <- seq_len(nrow(knresults))

table_knresults <- knresults[,c("cluster","size","statistic","pvalue")]

table_knresults

### Create final table (results by municipality)
final.result <- data.frame(
  ID=df$ID,
  region=df$region,
  population=df$population,
  Obs=df$Obs,
  Exp=df$Exp,
  SMR=df$Obs/df$Exp,
  x=df$x,
  y=df$y,
  map.index=seq_len(nrow(df)))

### Identify municipalities belonging to each cluster
if(nrow(knresults) > 0){
  clusters <- get.knclusters(df, knresults)

  for(i in seq_along(clusters)){
    final.result[[paste0("cluster_",i)]] <- ""
    final.result[[paste0("cluster_",i)]][clusters[[i]]] <- "cluster"
    final.result[[paste0("cluster_",i)]][clusters[[i]][1]] <- "centre"
  }
}


####################################################################################################################
# Step 4: Weighted Relative Risk estimator
####################################################################################################################

# Relative risk (RR) and Posterior probability (PP)

formula <- Obs ~ f(region,model="bym",graph="datos_ejemplo_practico_taller/neigs.graph")
        
result<-inla(formula, family="poisson",
data=data_input, E=Exp,
control.predictor=list(compute=TRUE, cdf=c(log(1))),
control.inla=list(strategy="simplified.laplace"))   

RR<-result$summary.fitted.values$"0.5quant"
	
PP<- 1-(result$summary.linear.predictor[1:dim(result$summary.linear.predictor)[1],"0cdf"])

final.result$"RR"<-RR
final.result$"lCre"<-result$summary.fitted.values$"0.025quant"
final.result$"uCre"<-result$summary.fitted.values$"0.975quant"
final.result$"PP"<-PP


# Weighted Relative Risk (WRR) 

final.result$"score"<-final.result$"RR"*final.result$"PP"


####################################################################################################################
# Step 5: Absolute risk estimator (DOE)
####################################################################################################################

final.result$"Diff_obs_exp"<-final.result$"Obs"-final.result$"Exp"



####################################################################################################################
# Step 6: Ranking spatial areas
####################################################################################################################

final.result$"score"[final.result$"PP"<0.8 | final.result$"Diff_obs_exp"<0]<-0

final.result$"Diff_obs_exp.new"<-final.result$"Diff_obs_exp"
final.result$"Diff_obs_exp.new"[final.result$"score"==0]<-0

final.result<-final.result[order(final.result$"Diff_obs_exp.new",final.result$"score",decreasing=T),]
final.result$"ranking"<-seq(1,dim(final.result)[1],1)

if(sum(final.result$"score"==0)==length(final.result$"score")){
	
	final.result$"ranking"<-1
	
}
if(sum(final.result$"score"==0)<length(final.result$"score")){

	final.result$"ranking"[final.result$"score"%in%0]<-max(final.result$"ranking"[final.result$"score"!=0])+1

}


final.result[,c("Exp","SMR","RR","lCre","uCre","PP","Diff_obs_exp")]<-round(final.result[,c("Exp","SMR","RR","lCre",
"uCre","PP","Diff_obs_exp")],digits=2)


####################################################################################################################
# Outputs (datos)
####################################################################################################################

complete.results<-final.result[order(final.result$"map.index"),c("ranking","ID","population","Obs","Exp","SMR",
"Diff_obs_exp","RR","lCre","uCre","PP")]

if(dim(knresults)[1]>0){
  
  table.clusters.complete<-as.data.frame(final.result[order(final.result$"map.index"),grep("cluster",names(final.result))])
  names(table.clusters.complete)<-names(final.result)[grep("cluster",names(final.result))]
  table1_clusters<-cbind(complete.results,table.clusters.complete)
  row.names(table1_clusters)<-NULL
  
}

if(dim(knresults)[1]==0){
  
  table1_clusters<-complete.results
  row.names(table1_clusters)<-NULL
  
}

table_results_pipeline<-table1_clusters[order(table1_clusters$"ranking"),]

row.names(table_results_pipeline)<-NULL

table_results_pipeline

####################################################################################################################
# Outputs graficos
####################################################################################################################


# Cluster

clusters <- get.knclusters(mapa, knresults)

sel_clus<-1

mapa_cluster<-mapa

mapa_cluster$cluster <- ""
  mapa_cluster$cluster[clusters[[sel_clus]]] <- "cluster"
  mapa_cluster$cluster[clusters[[sel_clus]][1]] <- "centre"
  mapa_cluster$cluster <- factor(mapa_cluster$cluster)



ggplot(mapa_cluster) +
	  geom_sf(aes(fill = cluster), color="black", size=0.2) +
	  scale_fill_manual(values = c("lightgrey", "orange", "red")) +
	  geom_sf_text(data=mapa_cluster[mapa_cluster$cluster%in%c("centre","cluster"),],aes(label=ID),size=2,fontface="bold") +
	  labs(
	    title = paste0("Cluster_", sel_clus, " — Kulldorff spatial scan (OPGAM)"),
	    fill = "State"
	  ) +
	  theme_void()
	  
	  
	  

		
# Relative Risk
	
plotvar1 <- final.result$RR[order(final.result$map.index)]

plotclr <- c("#298137","#5BAD37","#B4CC2D","#D1E691","#F3F090","#ED9E57",
             "#EF692E","#D73923","#A03227")

breaks_map <- sort(unique(c(0, .67, .77, .91, .95, 1.05, 1.1, 1.3, 1.5, max(plotvar1)+0.1)))

# Cortar en clases
class1 <- cut(plotvar1, breaks = breaks_map, include.lowest = TRUE, right = TRUE)

# Asignar colores a cada nivel
colcode1 <- as.character(class1)
levels_class <- levels(class1)
for(i in seq_along(levels_class)){
  colcode1[class1 == levels_class[i]] <- plotclr[i]
}

# Fusionar con shapefile sf
mapa_sf_plot <- merge(
  mapa, 
  data.frame(ID = final.result$ID, class = class1, color = colcode1),
  by = "ID",
  all.x = TRUE
)

# Crear etiquetas de leyenda manualmente 
make_legend_labels <- function(breaks, counts){
  n <- length(breaks)
  labs <- paste0("(", head(breaks, -1), " - ", tail(breaks, -1), "]")
  labs[1] <- paste0("≤ ", breaks[2])
  labs[n-1] <- paste0("> ", breaks[n-1])
  # añadir conteo N
  labs <- paste(labs, " (N=", counts, ")", sep="")
  return(labs)
}

# Conteo de municipios por color
tabla.N <- table(colcode1)
values.tabla.N <- as.numeric(tabla.N[plotclr])
values.tabla.N[is.na(values.tabla.N)] <- 0

leyenda.text <- make_legend_labels(breaks = breaks_map, counts = values.tabla.N)

# 
ggplot(mapa_sf_plot) +
  geom_sf(aes(fill = class), color = NA) +
  scale_fill_manual(
    values = plotclr,
    drop = FALSE,
    labels = leyenda.text,
    name = "Relative Risk (N)"
  ) +
  #geom_sf_text(aes(label = ID), size = 3) +
  theme_void() +
  ggtitle("Relative Risk Map")

	
# Mapa PP


# Cortes y clasificación 
nivel <- cut(
  final.result$PP[order(final.result$map.index)],
  breaks = c(-1,0.10,0.20,0.40,0.60,0.80,0.90,1),
  include.lowest = TRUE,
  right = TRUE
)

nParts <- as.numeric(nivel)

# Paleta de colores
mypalette <- brewer.pal(7, "RdBu")
#mipaleta <- rev(mypalette)  # invertir paleta
fgs <- mypalette[nParts]

# Fusionar colores y niveles al sf 

mapa_sf_plot <- merge(
  mapa,
  data.frame(ID = final.result$ID, nivel = nivel, color = fgs),
  by = "ID",
  all.x = TRUE
)

# Leyenda Manual
leyenda.text <- c("(0.9-1.0]", "(0.8-0.9]", "(0.6-0.8]", "(0.4-0.6]",
                  "(0.2-0.4]", "(0.1-0.2]", "<=0.1")


ggplot(mapa_sf_plot) +
  geom_sf(aes(fill = nivel), color = NA) +
  scale_fill_manual(
    values = mypalette,
    drop = FALSE,
    labels = leyenda.text,
    name = "Prob(RR>1)"
  ) +
  theme_void() +
  ggtitle("Prob(RR>1) Map")


##################		
# Mapa Ranking
##################		


# Datos a pintar
table_a_pintar <- final.result[final.result$score != 0,c("ranking","ID","population","Obs","Exp","SMR",
 "Diff_obs_exp","RR","lCre","uCre","PP")]
								 
						
 table_a_pintar$level_diff <- cut(
   table_a_pintar$Diff_obs_exp,
   breaks = quantile(table_a_pintar$Diff_obs_exp, probs = seq(0,1,1/6)),
   include.lowest = TRUE
 )

# Fusionar con shapefile 
mapa_sf_plot <- merge(
  mapa, 
  table_a_pintar[, c("ID", "Diff_obs_exp","level_diff")], 
  by = "ID", 
  all.x = TRUE
)


# Paleta
mypalette <- rev(heat.colors(6))


ggplot(mapa_sf_plot) +
  geom_sf(aes(fill = level_diff), color = "black", size = 0.1) +
  scale_fill_manual(
    values = mypalette,
    na.value = "white",     
    name = "Diff Obs - Exp"
  ) +
  geom_sf_text(
    data = mapa_sf_plot[!is.na(mapa_sf_plot$Diff_obs_exp), ],
    aes(label = ID),
    size = 3
  ) +
  theme_void() +
  ggtitle("Diff(Obs - Exp) Map")


####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

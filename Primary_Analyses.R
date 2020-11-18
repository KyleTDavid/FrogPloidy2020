library(spatstat)
library(maptools)
library(raster)
library(ggplot2)
library(sf)
library(rotl)
library(rgdal)
library(viridis)
library(ape)
library(dismo)
library(dplyr)
library(gridExtra)
library(factoextra)
library(phyloclim)
library(geiger)
library(phytools)
library(ggridges)

#setwd("~/FrogPloidy2020")
setwd("~/Desktop/frogs/")

#South American Frog Genera with Diploid and Polyploid Species
Genera <- c("Ceratophrys", "Chiasmocleis", "Odontophrynus", "Phyllomedusa", "Pleurodema")
#South American Frog Polyploid Species
PolySpecies <- c("Ceratophrys aurita", "Ceratophrys joazeirensis", "Ceratophrys ornata",
                 "Chiasmocleis leucosticta", "Odontophrynus americanus",
                 "Phyllomedusa tetraploidea", "Pleurodema cordobae", "Pleurodema bibroni",
                 "Pleurodema kriegi")

#GBIF occurences queried Jan 16 2020
gbif <- read.delim("GBIF.csv", stringsAsFactors = F, row.names = NULL)
gbif <- gbif[gbif$species!="",]

### IUCN ranges downloaded Oct 15 2019###
iucn <- read_sf("ANURA/ANURA.shp", stringsAsFactors = F)
iucn <- iucn[(iucn$genus %in% Genera),]

#get full list of species across vertnet and IUCN
species <- unique(c(gbif$species, iucn$binomial))
resolved_names <- tnrs_match_names(species)

#fix synonyms
gbif[["species"]] <- resolved_names[match(tolower(gbif[['species']]), resolved_names[['search_string']] ) , 'unique_name']
iucn[["species"]] <- resolved_names[match(tolower(iucn[['binomial']]), resolved_names[['search_string']] ) , 'unique_name']

#check discrepancies with AmphibiaWeb
gbif$species <- gsub("Syncope", "Chiasmocleis", gbif$species)
gbif$species <- gsub("Pithecopus", "Phyllomedusa", gbif$species)
gbif$species <- gsub("Physalaemus", "Pleurodema", gbif$species)
gbif$species <- gsub("nordestinus", "nordestina", gbif$species)
gbif$species <- gsub("azureus", "azurea", gbif$species)

#filter taxa that don't appear to be belong to any of the genera
gbif <- gbif[grepl("Ceratophrys|Chiasmocleis|Odontophrynus|Pleurodema|Phyllomedusa", gbif$species, ignore.case = F),]

#filter taxa with both diploid and polyploid members
gbif <- gbif[gbif$species != "Phyllomedusa burmeisteri",]

#assign ploidy
gbif$ploidy <- ifelse(gbif$species %in% PolySpecies, "polyploid", "diploid")
iucn$ploidy <- ifelse(iucn$binomial %in% PolySpecies, "polyploid", "diploid")

#make map of South America
mapSA <- borders("world", regions = c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Ecuador", "French Guiana", "Guyana", "Paraguay", "Peru", "Suriname", "Uruguay", "Venezuela"),
                 fill = "grey90", colour = "grey90")

#plot occurrences and IUCN ranges for diploids and polyploids
ggplot(mapping=aes(fill=ploidy)) + mapSA + geom_sf(data=iucn, aes(alpha=0.5, color=ploidy)) +
  geom_point(data=gbif, mapping=aes(decimalLongitude, decimalLatitude), size=1.5, stroke=1.5, shape=21) + theme_void()

#percent of polyploids in SESA
length(gbif[gbif$decimalLatitude < -25 & gbif$decimalLongitude > -65 & gbif$ploidy=="polyploid",]$species) / 
  length(gbif[gbif$ploidy=="polyploid",]$species)

#percent of diploids in SESA
length(gbif[gbif$decimalLatitude < -25 & gbif$decimalLongitude > -65 & gbif$ploidy=="diploid",]$species) / 
  length(gbif[gbif$ploidy=="diploid",]$species)

#read rasters of worlclim climatic variables
temp <- raster("worldclim/temp.tif")
prec <- raster("worldclim/prec.tif")
alt <- raster("worldclim/alt.tif")
season_temp <- raster("worldclim/season_temp.tif")
season_prec <- raster("worldclim/season_prec.tif")
climate <- stack(temp, prec, alt, season_temp, season_prec)
climate <- mask(climate, cropshape)

#maxent of all polyploid and diploid species collectively 
poly_maxent <- maxent(climate, subset(gbif[gbif$ploidy=="polyploid" & gbif$genus=="Ceratophrys",], select=c("decimalLongitude","decimalLatitude")))
dip_maxent <- maxent(climate, subset(gbif[gbif$ploidy=="diploid" & gbif$genus=="Ceratophrys",], select=c("decimalLongitude","decimalLatitude")))
#predictiabilty scores of polyploids and diploids
poly_predict <- predict(poly_maxent, climate)
dip_predict <- predict(dip_maxent, climate)
#plot diploids - polyploids
dif_predict <- dip_predict-poly_predict
pal <- colorRampPalette(c("#00BFC4","grey90","#FF776D"))
plot(dif_predict, col=pal(100), breaks=seq(-1,1,0.02), main="Ceratophrys")
#read in tree
ultrametric_tree <- read.tree("tree/ultrametric_tree.nh")
ultrametric_tree$tip.label<-gsub("_"," ",ultrametric_tree$tip.label)
#plot tree of just genera level relationships
genus_tree <- keep.tip(ultrametric_tree, ultrametric_tree$tip.label[ultrametric_tree$tip.label %in% c("Ceratophrys aurita", "Odontophrynus americanus", "Chiasmocleis leucosticta", "Phyllomedusa tetraploidea", "Pleurodema bibroni")])
plot(genus_tree, show.tip.label=T)

#function to return grid of predictability scores from a subset of the gbif dataframe
maxent.ag <- function(df) {
  subdf <- subset(df, select=c("decimalLongitude","decimalLatitude"))
  df.maxent <- maxent(climate, subdf)
  df.predict <- predict(df.maxent, climate)
  df.grid <- as(df.predict, 'SpatialGridDataFrame')
  return(df.grid)
}

#run maxent on each polyploid species and return grid of predictability scores
gbif_polyploid <- gbif[gbif$ploidy=="polyploid",]
gbif_polyploid$species <- factor(gbif_polyploid$species)
polyploid_predictlist <- by(gbif_polyploid, gbif_polyploid$species, maxent.ag)
#estimate pairwise range overlap from each maxent result
polyploid_overlap <- niche.overlap(polyploid_predictlist)
#tree with only polyploid taxa
polytree <- keep.tip(ultrametric_tree, ultrametric_tree$tip.label[ultrametric_tree$tip.label %in% dimnames(polyploid_overlap)[[1]]])
#age range correlation
polyploid_ARC <- age.range.correlation(polytree, polyploid_overlap)
polyploid_ARC_df <- as.data.frame(polyploid_ARC$age.range.correlation)

#run maxent on each diploid species and return grid of predictability scores
gbif_diploid <- gbif[gbif$ploidy=="diploid",]
gbif_diploid$species <- factor(gbif_diploid$species)
diploid_predictlist <- by(gbif_diploid, gbif_diploid$species, maxent.ag)
#estimate pairwise range overlap from each maxent result
diploid_overlap <- niche.overlap(diploid_predictlist)
#tree with only diploid taxa
ditree <- keep.tip(ultrametric_tree, ultrametric_tree$tip.label[ultrametric_tree$tip.label %in% dimnames(diploid_overlap)[[1]]])
#age range correlation
diploid_ARC <- age.range.correlation(ditree, diploid_overlap)
diploid_ARC_df <- as.data.frame(diploid_ARC$age.range.correlation)

#wwf biomes and ecoregions
wwf <- read_sf("wwf/wwf_terr_ecos.shp")
wwf <- st_crop(wwf, xmin = -100, xmax = -30, ymin = -60, ymax = 15)
wwf <-wwf[c("geometry", "BIOME", "ECO_NAME")]
wwf_data <- st_intersection(st_as_sf(gbif, coords = c("decimalLongitude", "decimalLatitude"), crs=4326), wwf)
wwf_data <- wwf_data[c("species", "ploidy", "BIOME", "ECO_NAME")]

#kg climate data
kg_map <- raster("Map_KG-Global/KG_1986-2010.grd")
kg_map <- crop(kg_map, extent(-100, -30, -60, 15))
points <- SpatialPoints(data.frame(decimallongitude=gbif$decimalLongitude, decimallatitude=gbif$decimalLatitude), proj4string = kg_map@crs)
kg_climate <- data.frame(kg=extract(kg_map, points), ploidy=gbif$ploidy)

#count occurences in each KG climate
kg <- data.frame(ecoregion=character(), dip_count=integer(), poly_count=integer(),
                 dip_percent=integer(), poly_percent=integer())
for (k in unique(kg_climate$kg)) {
  dip_count <- length(kg_climate[kg_climate$kg==k & kg_climate$ploidy=="diploid",]$kg)
  poly_count <- length(kg_climate[kg_climate$kg==k & kg_climate$ploidy=="polyploid",]$kg)
  dip_percent <- dip_count/length(kg_climate[kg_climate$ploidy=="diploid",]$kg)
  poly_percent <- poly_count/length(kg_climate[kg_climate$ploidy=="polyploid",]$kg)
  kg <- rbind(kg, data.frame(kg=k, dip_count=dip_count, poly_count=poly_count,
                             dip_percent=dip_percent, poly_percent=poly_percent))
}

#count occurrences in each ecoregion
ecoregions<-data.frame(ecoregion=character(), dip_count=integer(), poly_count=integer(),
                       dip_percent=integer(), poly_percent=integer())
for (eco in unique(wwf_data$ECO_NAME)) {
  dip_count <- length(wwf_data[wwf_data$ECO_NAME==eco & wwf_data$ploidy=="diploid",]$ECO_NAME)
  poly_count <- length(wwf_data[wwf_data$ECO_NAME==eco & wwf_data$ploidy=="polyploid",]$ECO_NAME)
  dip_percent <- dip_count/sum(length(wwf_data[wwf_data$ploidy=="diploid",]$ECO_NAME))
  poly_percent <- poly_count/sum(length(wwf_data[wwf_data$ploidy=="polyploid",]$ECO_NAME))
  ecoregions <- rbind(ecoregions, data.frame(ecoregion=eco, dip_count=dip_count, poly_count=poly_count,
                                             dip_percent=dip_percent, poly_percent=poly_percent))
}

#count occurrences in each biome
biomes<-data.frame(biome=character(), dip_count=integer(), poly_count=integer(),
                   dip_percent=integer(), poly_percent=integer())
for (eco in unique(wwf_data$BIOME)) {
  dip_count <- length(wwf_data[wwf_data$BIOME==eco & wwf_data$ploidy=="diploid",]$BIOME)
  poly_count <- length(wwf_data[wwf_data$BIOME==eco & wwf_data$ploidy=="polyploid",]$BIOME)
  dip_percent <- dip_count/sum(length(wwf_data[wwf_data$ploidy=="diploid",]$BIOME))
  poly_percent <- poly_count/sum(length(wwf_data[wwf_data$ploidy=="polyploid",]$BIOME))
  biomes <- rbind(biomes, data.frame(BIOME=eco, dip_count=dip_count, poly_count=poly_count,
                                     dip_percent=dip_percent, poly_percent=poly_percent))
}

#test for similarity of climates, biomes, and ecoregions between polyploids and diploids
chisq.test(kg %>% select(dip_count,poly_count))
chisq.test(biomes %>% select(dip_count,poly_count))
chisq.test(ecoregions %>% select(dip_count,poly_count))

#read rasters of SEDAC anthropogenic variables
Pfert <- raster("SEDAC/pfertilizer_samerica.tif")
Nfert <- raster("SEDAC/nfertilizer_samerica.tif")
fert <- Pfert + Nfert
Pmanur <- raster("SEDAC/pmanure_samerica.tif")
Nmanur <- raster("SEDAC/nmanure_samerica.tif")
manur <- Pmanur + Nmanur
cropland <- raster("SEDAC/sa_cropland.tif")
pasture <- raster("SEDAC/sa_pasture.tif")
pest <- raster("SEDAC/pesticide_samerica.tif")

#append occurence dataframe with the values from each quantitative variable at the same lat/long
points <- SpatialPoints(data.frame(decimallongitude=gbif$decimalLongitude, decimallatitude=gbif$decimalLatitude))
gbif$temp <- extract(temp, points)
gbif$prec <- extract(prec, points)
gbif$alt <- extract(alt, points)
gbif$season_temp <- extract(season_temp, points)
gbif$season_prec <- extract(season_prec, points)
gbif$fert <- extract(fert, points)
gbif$manur <- extract(manur, points)
gbif$cropland <- extract(cropland, points)
gbif$pasture <- extract(pasture, points)
gbif$pest <- extract(pest, points)

gbif <- as.data.frame(gbif)

#function that accepts a variable name as a string from gbif dataframe, will plot violin plots
#of variable distribution between polyploids and diploids within and across genera, will also
#report wilcox pval
violin_plotter <- function(var) {
  all <- ggplot(gbif, mapping=aes(x=ploidy, y=gbif[[var]], fill=ploidy)) + geom_violin() + geom_boxplot(width=0.1, outlier.size=-1) +
    theme_void() + scale_y_continuous() + 
    theme_void() + theme(legend.position = "none") + theme(axis.text.y = element_text(), panel.grid.major = element_line()) + 
    coord_cartesian(ylim=c(min(na.omit(gbif[[var]])), max(na.omit(gbif[[var]])))) + 
    annotate("text", x=1.5, y=mean(na.omit(gbif[[var]])), size=5, label=as.character(signif(wilcox.test(gbif[gbif$ploidy=="polyploid",][[var]], gbif[gbif$ploidy=="diploid",][[var]])$p.value, 3)))
  Ce <- ggplot(gbif[gbif$genus=="Ceratophrys",], mapping=aes(x=ploidy, y=gbif[gbif$genus=="Ceratophrys",][[var]], fill=ploidy)) +
    geom_violin() + geom_boxplot(width=0.1, outlier.size=-1) + theme_void() + theme(legend.position = "none") +
    coord_cartesian(ylim=c(min(na.omit(gbif[[var]])), max(na.omit(gbif[[var]])))) +
    annotate("text", x=1.5, y=mean(na.omit(gbif[gbif$genus=="Ceratophrys",][[var]])), size=5, label=as.character(signif(wilcox.test(gbif[gbif$ploidy=="polyploid" & gbif$genus =="Ceratophrys" ,][[var]], gbif[gbif$ploidy=="diploid" & gbif$genus=="Ceratophrys",][[var]])$p.value, 3)))
  Ch <- ggplot(gbif[gbif$genus=="Chiasmocleis",], mapping=aes(x=ploidy, y=gbif[gbif$genus=="Chiasmocleis",][[var]], fill=ploidy)) +
    geom_violin() + geom_boxplot(width=0.1, outlier.size=-1) + theme_void() + theme(legend.position = "none") +
    coord_cartesian(ylim=c(min(na.omit(gbif[[var]])), max(na.omit(gbif[[var]])))) +
    annotate("text", x=1.5, y=mean(na.omit(gbif[gbif$genus=="Chiasmocleis",][[var]])), size=5, label=as.character(signif(wilcox.test(gbif[gbif$ploidy=="polyploid" & gbif$genus =="Chiasmocleis" ,][[var]], gbif[gbif$ploidy=="diploid" & gbif$genus=="Chiasmocleis",][[var]])$p.value, 3)))
  O <- ggplot(gbif[gbif$genus=="Odontophrynus",], mapping=aes(x=ploidy, y=gbif[gbif$genus=="Odontophrynus",][[var]], fill=ploidy)) +
    geom_violin() + geom_boxplot(width=0.1, outlier.size=-1) + theme_void() + theme(legend.position = "none") +
    coord_cartesian(ylim=c(min(na.omit(gbif[[var]])), max(na.omit(gbif[[var]])))) +
    annotate("text", x=1.5, y=mean(na.omit(gbif[gbif$genus=="Odontophrynus",][[var]])), size=5, label=as.character(signif(wilcox.test(gbif[gbif$ploidy=="polyploid" & gbif$genus =="Odontophrynus" ,][[var]], gbif[gbif$ploidy=="diploid" & gbif$genus=="Odontophrynus",][[var]])$p.value, 3)))
  Ph <- ggplot(gbif[gbif$genus=="Phyllomedusa",], mapping=aes(x=ploidy, y=gbif[gbif$genus=="Phyllomedusa",][[var]], fill=ploidy)) +
    geom_violin() + geom_boxplot(width=0.1, outlier.size=-1) +
    theme_void() + theme(legend.position = "none") +
    coord_cartesian(ylim=c(min(na.omit(gbif[[var]])), max(na.omit(gbif[[var]])))) +
    annotate("text", x=1.5, y=mean(na.omit(gbif[gbif$genus=="Phyllomedusa",][[var]])), size=5, label=as.character(signif(wilcox.test(gbif[gbif$ploidy=="polyploid" & gbif$genus =="Phyllomedusa" ,][[var]], gbif[gbif$ploidy=="diploid" & gbif$genus=="Phyllomedusa",][[var]])$p.value, 3)))
  Pl <- ggplot(gbif[gbif$genus=="Pleurodema",], mapping=aes(x=ploidy, y=gbif[gbif$genus=="Pleurodema",][[var]], fill=ploidy)) +
    geom_violin() + geom_boxplot(width=0.1, outlier.size=-1) + theme_void() + theme(legend.position = "none") +
    coord_cartesian(ylim=c(min(na.omit(gbif[[var]])), max(na.omit(gbif[[var]])))) +
    annotate("text", x=1.5, y=mean(na.omit(gbif[gbif$genus=="Pleurodema",][[var]])), size=5, label=as.character(signif(wilcox.test(gbif[gbif$ploidy=="polyploid" & gbif$genus =="Pleurodema" ,][[var]], gbif[gbif$ploidy=="diploid" & gbif$genus=="Pleurodema",][[var]])$p.value, 3)))
  grid.arrange(all, Ce, Ch, O, Ph, Pl, nrow=1)
}

#PCA
frog.pca <- prcomp(na.omit(gbif[,c(52:61)]), center = T, scale. = T)
frogdf <- na.omit(gbif[c(9 , 51:61)])
#PCA plot
grid.arrange(
  fviz_pca_biplot(frog.pca, col.ind = frogdf$ploidy, label = "var", addEllipses = T, ellipse.level=0.68, col.var = "grey25", pointshape=19) + theme(legend.position = "none"),
  fviz_pca_biplot(frog.pca, col.ind = frogdf$genus, label = "var", addEllipses = T, ellipse.level=0.68, pointshape=19, palette = c("#B873E0", "#F49A40", "#7676D3", "#76DD73", "#8BEFEF"),col.var = "grey25") + theme(legend.position = "none"),
  nrow=1)

#phylogenetic anova between diploid and polyploid species for a environmental variable (in this case temperature seasonality)
gbif_ag <- aggregate(list(
  temp=gbif$temp, prec=gbif$prec, alt=gbif$alt, season_temp=gbif$season_temp, season_prec=gbif$season_prec, fert=gbif$fert, manur=gbif$manur, cropland=gbif$cropland, pasture=gbif$pasture, pest=gbif$pest
), by=list(species=gbif$species, genus=gbif$genus, ploidy=gbif$ploidy), mean, na.rm=T)
gbif_ag$species <- gsub(" ", "_", gbif_ag$species)
x <- gbif_ag$season_temp
names(x) <- unlist(gbif_ag$species)
grop <- factor(gbif_ag$ploidy)
names(grop) <- unlist(gbif_ag$species)
aov.phylo(x~grop, ultrametric_tree, nsim = 1000) 

#phylogenetic PCA and plots
gbif_ag <- gbif_ag[gbif_ag$species %in% ultrametric_tree$tip.label,]
rownames(gbif_ag) <- gbif_ag$species
sub_gbif = subset(gbif_ag,select=-c(species,genus,ploidy))
res <- phyl.pca(ultrametric_tree, na.omit(sub_gbif))
tmp <- data.frame(res$S)
plot(tmp$PC1, tmp$PC2)
frog <- merge(tmp, gbif_ag, by.x=0, by.y='species')
l <- data.frame(res$L)
pcaploidy <- ggplot(frog) + geom_point(aes(PC1, PC2, color=ploidy)) + stat_ellipse(level = 0.68, aes(PC1, PC2, color=ploidy)) +
  geom_segment(data=l, aes(0,0,xend=PC1*2000,yend=PC2*2000), arrow=arrow()) +
  annotate("text", x = l$PC1*2000, y = l$PC2*2000, label = rownames(l)) + theme_classic() + theme(legend.position = "none")
pcagenus <- ggplot(frog) + geom_point(aes(PC1, PC2, color=genus)) + stat_ellipse(level = 0.68, aes(PC1, PC2, color=genus)) + scale_color_manual(values=c("#B873E0", "#F49A40", "#7676D3", "#76DD73", "#8BEFEF")) +
  geom_segment(data=l, aes(0,0,xend=PC1*2000,yend=PC2*2000), arrow=arrow()) +
  annotate("text", x = l$PC1*2000, y = l$PC2*2000, label = rownames(l)) + theme_classic() + theme(legend.position = "none") 
grid.arrange(
  pcaploidy,
  pcagenus,
  nrow=1)

#temperature seasonality ridge plots
gbif$species <- gsub(" ", "_", gbif$species)
small_samples <- c("Ceratophrys_joazeirensis", "Phyllomedusa_ayeaye", "Phyllomedusa_centralis", "Phyllomedusa_araguari", "Chiasmocleis_lacrimae")
gbif <- gbif[!(gbif$species %in% small_samples),]
tree <- drop.tip(ultrametric_tree, c("Sooglossus_thomasseti", small_samples))
gbif$factor <- factor(gbif$species, levels = ultrametric_tree$tip.label)
ggplot(gbif, aes(x = season_temp, y = factor, fill = stat(x), color=ploidy)) +
  geom_density_ridges_gradient(size=0.25) +
  scale_fill_viridis_c(option='plasma') +
  theme_classic() +
  theme(legend.position = 'none')






#used to get pesticide raster by taking the mean between high and low estimates for all 2015 variables
files <- list.files(path="SEDAC/ferman-v1-pest-chemgrids-geotiff/ferman-v1-pest-chemgrids_geotiff/ApplicationRate/GEOTIFF", pattern="*_2015_*H.tif", full.names=TRUE, recursive=FALSE)
for (f in files){
  r <- raster(f)
  croper <- crop(r, extent(-100, -30, -60, 15))
  croper[croper<0] <- 0
  high_estimate <- high_estimate + croper }
high_estimate[high_estimate==0] <- NA
pest_stack <- stack(low_estimate, high_estimate)
pest_mean <- calc(pest_stack, fun = mean)
writeRaster(pest_mean, "SEDAC/pesticide_samerica.tif", format="GTiff")







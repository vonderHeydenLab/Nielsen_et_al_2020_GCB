################################################################################
# Scripts to perfrom species distribution models as described in
# Nielsen et al. (2021) 'Distinct inter- and intra-specific vulnerability of 
# coastal species to global change'
#
# Adapted and/or written by ES Nielsen, Stellenbosch University, South Africa
#
# Code is provided as is, without support 
################################################################################

#CP=cyclograpsus punctatus
#PA=parechinus angulosus
#SG=scutellastra granularis

################################################################################
####### Code to run Gradient Forest models (on each species individually) ######
################################################################################

library(RStoolbox)
library(gradientForest)
library(rasterVis)
library(raster)
library(tidyverse)
library(sdmpredictors)
library(vegan)
library(latticeExtra)


########## Download Environmental Layers #######################

######################## DOWNLOAD #############################
######################## CURRENT ##############################

a.contp <- load_layers( layercodes = c("WC_bio1", "WC_bio5", "WC_bio7", "WC_bio12", "WC_bio15") , equalarea=FALSE, rasterstack=TRUE)
environment.contp <- load_layers( layercodes = c("BO2_salinityrange_bdmin", "BO2_salinitymean_bdmin", "BO2_temprange_bdmin", "BO2_tempmean_bdmin"))
chl <- load_layers("BO2_chlomean_ss")
a<-extent(-180, 180, -90, 90)
o<-extent(-180, 180, -90, 90)
extent_list<-list(a, o)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)
reso<-res(a.contp)
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=a.contp@crs)
a.c.r.2 <-resample(a.contp, s2, method="ngb")
o.c.r.2 <-resample(environment.contp, s2, method="ngb")
chl.2 <-resample(chl, s2, method="ngb")
contemp.r.2=stack(a.c.r.2, o.c.r.2, chl.2)
writeRaster(contemp.r.2, "10env.curr.tif", format="GTiff", overwrite=TRUE)

SA.cp.ext <- extent(16, 34.08333, -34.91667, -25)
SA.pa.ext <- extent(14.7, 32.3, -37, -25)

setwd("~/Desktop/PhD_stuffies/CH2/NEW_SDMS/Cpunctatus")
MyBinCP_curr <- raster::stack("Cpunctatus/proj_current/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

cont.env.c <- crop(contemp.r.2, MyBinCP_curr)
cont.env.CP <- mask(cont.env.c, MyBinCP_curr)
cont.env.CP <- crop(cont.env.CP, SA.cp.ext)

cont.env.c <- crop(contemp.r.2, MyBinCP_curr)
cont.env.CP <- mask(cont.env.c, MyBinCP_curr)
cont.env.PA <- crop(cont.env.CP, SA.pa.ext)

cont.env.c <- crop(contemp.r.2, MyBinCP_curr)
cont.env.CP <- mask(cont.env.c, MyBinCP_curr)
cont.env.SG <- crop(cont.env.CP, SA.pa.ext)

setwd("~/Desktop/PhD_stuffies/CH3/new_GF/CP")
env_trns <-as.data.frame(cont.env.CP, row.names=NULL, optional=FALSE, xy=FALSE, na.rm=TRUE)
env_trns.xy <-as.data.frame(cont.env.CP, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=TRUE)


######################## DOWNLOAD #############################
######################## FUTURE ###############################

############################ 2050 RCP layers ##########################

##RCP45 2050
mr.50.45.t <- load_layers( layercodes = c("WC_bio1_mr45_2050", "WC_bio5_mr45_2050", "WC_bio7_mr45_2050", "WC_bio12_mr45_2050", "WC_bio15_mr45_2050") , equalarea=FALSE, rasterstack=TRUE)
cc.50.45.t <- load_layers( layercodes = c("WC_bio1_cc45_2050", "WC_bio5_cc45_2050", "WC_bio7_cc45_2050", "WC_bio12_cc45_2050", "WC_bio15_cc45_2050") , equalarea=FALSE, rasterstack=TRUE)
he.50.45.t <- load_layers( layercodes = c("WC_bio1_he45_2050", "WC_bio5_he45_2050", "WC_bio7_he45_2050", "WC_bio12_he45_2050", "WC_bio15_he45_2050") , equalarea=FALSE, rasterstack=TRUE)

t.50.45 <- overlay(mr.50.45.t, cc.50.45.t, he.50.45.t, fun=mean)

o.50.45 <- load_layers( layercodes = c("BO2_RCP45_2050_salinityrange_bdmin", "BO2_RCP45_2050_salinitymean_bdmin", "BO2_RCP45_2050_temprange_bdmin", "BO2_RCP45_2050_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
chl.50.45 <- load_layers("BO2_RCP45_2050_chlomean_bdmin")
a.mid.2 <-resample(t.50.45, s2, method="ngb")
o.mid.2 <-resample(o.50.45, s2, method="ngb")
chl.mid.2 <-resample(chl.50.45, s2, method="ngb")
env.50.45 <- stack(a.mid.2, o.mid.2, chl.mid.2)


##RCP85 CCSM4
mr.50.85.t <- load_layers( layercodes = c("WC_bio1_mr85_2050", "WC_bio5_mr85_2050", "WC_bio7_mr85_2050", "WC_bio12_mr85_2050", "WC_bio12_mr85_2050") , equalarea=FALSE, rasterstack=TRUE)
cc.50.85.t <- load_layers( layercodes = c("WC_bio1_cc85_2050", "WC_bio5_cc85_2050", "WC_bio7_cc85_2050", "WC_bio12_cc85_2050", "WC_bio12_cc85_2050") , equalarea=FALSE, rasterstack=TRUE)
he.50.85.t <- load_layers( layercodes = c("WC_bio1_he85_2050", "WC_bio5_he85_2050", "WC_bio7_he85_2050", "WC_bio12_he85_2050", "WC_bio12_he85_2050") , equalarea=FALSE, rasterstack=TRUE)

t.50.85 <- overlay(mr.50.85.t, cc.50.85.t, he.50.85.t, fun=mean)

o.50.85 <- load_layers( layercodes = c("BO2_RCP85_2050_salinityrange_bdmin", "BO2_RCP85_2050_salinitymean_bdmin", "BO2_RCP85_2050_temprange_bdmin", "BO2_RCP85_2050_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
chl.50.85 <- load_layers("BO2_RCP85_2050_chlomean_bdmin")
a.mid.2 <-resample(t.50.85, s2, method="ngb")
o.mid.2 <-resample(o.50.85, s2, method="ngb")
chl.mid.2 <-resample(chl.50.85, s2, method="ngb")
env.50.85 <- stack(a.mid.2, o.mid.2, chl.mid.2)

env.50 <- merge(env.50.45, env.50.85)

############################ 2070-2100 RCP layers ##########################

## RCP 45
mr.90.45.t <- load_layers( layercodes = c("WC_bio1_mr45_2070", "WC_bio5_mr45_2070", "WC_bio7_mr45_2070", "WC_bio12_mr45_2070", "WC_bio15_mr45_2070") , equalarea=FALSE, rasterstack=TRUE)
cc.90.45.t <- load_layers( layercodes = c("WC_bio1_cc45_2070", "WC_bio5_cc45_2070", "WC_bio7_cc45_2070", "WC_bio12_cc45_2070", "WC_bio15_cc45_2070") , equalarea=FALSE, rasterstack=TRUE)
he.90.45.t <- load_layers( layercodes = c("WC_bio1_he45_2070", "WC_bio5_he45_2070", "WC_bio7_he45_2070", "WC_bio12_he45_2070", "WC_bio15_he45_2070") , equalarea=FALSE, rasterstack=TRUE)

t.90.45 <- overlay(mr.90.45.t, cc.90.45.t, he.90.45.t, fun=mean)

o.90.45 <- load_layers( layercodes = c("BO2_RCP45_2100_salinityrange_bdmin", "BO2_RCP45_2100_salinitymean_bdmin", "BO2_RCP45_2100_temprange_bdmin", "BO2_RCP45_2100_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
chl.90.45 <- load_layers("BO2_RCP45_2100_chlomean_bdmin")
a.mid.2 <-resample(t.90.45, s2, method="ngb")
o.mid.2 <-resample(o.90.45, s2, method="ngb")
chl.mid.2 <-resample(chl.90.45, s2, method="ngb")
env.90.45 <- stack(a.mid.2, o.mid.2, chl.mid.2)


## RCP 85
mr.90.85.t <- load_layers( layercodes = c("WC_bio1_mr85_2070", "WC_bio5_mr85_2070", "WC_bio7_mr85_2070", "WC_bio12_mr85_2070", "WC_bio15_mr85_2070") , equalarea=FALSE, rasterstack=TRUE)
cc.90.85.t <- load_layers( layercodes = c("WC_bio1_cc85_2070", "WC_bio5_cc85_2070", "WC_bio7_cc85_2070", "WC_bio12_cc85_2070", "WC_bio15_cc85_2070") , equalarea=FALSE, rasterstack=TRUE)
he.90.85.t <- load_layers( layercodes = c("WC_bio1_he85_2070", "WC_bio5_he85_2070", "WC_bio7_he85_2070", "WC_bio12_he85_2070", "WC_bio15_he85_2070") , equalarea=FALSE, rasterstack=TRUE)

t.90.85 <- overlay(mr.90.85.t, cc.90.85.t, he.90.85.t, fun=mean)

o.90.85 <- load_layers( layercodes = c("BO2_RCP85_2100_salinityrange_bdmin", "BO2_RCP85_2100_salinitymean_bdmin", "BO2_RCP85_2100_temprange_bdmin", "BO2_RCP85_2100_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
chl.90.85 <- load_layers("BO2_RCP85_2100_chlomean_bdmin")
a.mid.2 <-resample(t.90.85, s2, method="ngb")
o.mid.2 <-resample(o.90.85, s2, method="ngb")
chl.mid.2 <-resample(chl.90.85, s2, method="ngb")
env.90.85 <- stack(a.mid.2, o.mid.2, chl.mid.2)

env.90 <- merge(env.90.45, env.90.85)

writeRaster(env.50.45, "env.50.45.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.50.85, "env.50.85.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.90.45, "env.90.45.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.90.85, "env.90.85.tif", format="GTiff", overwrite=TRUE)

writeRaster(env.50, "env.50.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.90, "env.90.tif", format="GTiff", overwrite=TRUE)

# Set species extents
SA.ext <- extent(5, 45, -40, -10)
SA.cp.ext <- extent(10, 40, -37, -22)
SA.pa.ext <- extent(10, 40, -37, -19)
SA.sg.ext <- extent(10, 40, -37, -12)

##################### Import saved raster layers from above
##################

setwd("~/Desktop/PhD_stuffies/CH3/GF")
env.50 <- 'env.50.tif'
env.50=stack(env.50)

env.90 <- 'env.90.tif'
env.90=stack(env.90)

#crop by CP
env.50.CP <- crop(env.50, MyBinCP_curr)
env.50.CP <- mask(env.50.CP, MyBinCP_curr)
env.50.CP <- crop(env.50.CP, SA.cp.ext)

env.90.CP <- crop(env.90, MyBinCP_curr)
env.90.CP <- mask(env.90.CP, MyBinCP_curr)
env.90.CP <- crop(env.90.CP, SA.cp.ext)

stack.50.CP <- stack(cont.env.CP, env.50.CP)
env_trns_50 <- raster::extract(stack.50.CP, 1:ncell(stack.50.CP), df = TRUE)
env_trns_50 <- na.omit(env_trns_50)
env_cont <- as.data.frame(env_trns_50[,1:11])
env_f.50 <- as.data.frame(env_trns_50[c(1,12:21)])


stack.90.CP <- stack(cont.env.CP, env.90.CP)
env_trns_90 <- raster::extract(stack.90.CP, 1:ncell(stack.90.CP), df = TRUE)
env_trns_90 <- na.omit(env_trns_90)
env_f.90 <- as.data.frame(env_trns_90[c(1,12:21)])


env_fut.xy <-as.data.frame(stack.50.CP, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=TRUE)
env_fut.xy <- env_fut.xy[,1:2]

#crop by PA/SG

env.50.CP <- crop(env.50, MyBinCP_curr)
env.50.CP <- mask(env.50.CP, MyBinCP_curr)
env.50.PA <- crop(env.50.CP, SA.pa.ext)

env.90.CP <- crop(env.90, MyBinCP_curr)
env.90.CP <- mask(env.90.CP, MyBinCP_curr)
env.90.PA <- crop(env.90.CP, SA.pa.ext)

stack.50.PA <- stack(cont.env.PA, env.50.PA)
env_trns_50 <- raster::extract(stack.50.PA, 1:ncell(stack.50.PA), df = TRUE)
env_trns_50 <- na.omit(env_trns_50)
env_cont <- as.data.frame(env_trns_50[,1:11])
env_f.50 <- as.data.frame(env_trns_50[c(1,12:21)])


stack.90.PA <- stack(cont.env.PA, env.90.PA)
env_trns_90 <- raster::extract(stack.90.PA, 1:ncell(stack.90.PA), df = TRUE)
env_trns_90 <- na.omit(env_trns_90)
env_f.90 <- as.data.frame(env_trns_90[c(1,12:21)])


env_fut.xy <-as.data.frame(stack.50.PA, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=TRUE)
env_fut.xy <- env_fut.xy[,1:2]


######################## DOWNLOAD #############################
######################## INPUTS ###############################

# read in data file with minor allele freqs & env/space variables
## NB: make sure csv is separated by "," and not ";" - can find replace in text editor
gfData <- read.csv("cp.env.GF.csv")
envGF <- gfData[,4:13] # get climate & MEM variables

## Calculate PCNMs to account for geographic space

coord <- gfData[,c("x","y")]
pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
keep <- round(length(which(pcnm$value > 0))/2)
pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
pcnm.keep <- as.data.frame(pcnm.keep)

envGF <- cbind(envGF, pcnm.keep)

# build individual SNP datasets
SNPs_neut <- gfData[,grep("merged",colnames(gfData))] # neutral SNPs
SNPs_out <- gfData[,grep("outlier",colnames(gfData))] # outlier SNPs

# GRADIENT FOREST MODELING -----------------------------------------------------

maxLevel <- log2(0.368*nrow(envGF)/2) #account for correlations, see ?gradientForest 

gfNeut <- gradientForest(cbind(envGF, SNPs_neut), predictor.vars=colnames(envGF),
                         response.vars=colnames(SNPs_neut), ntree=500, mtry=NULL,transform=NULL,
                         maxLevel=maxLevel, trace=T)

# Fit gf models for GI5 SNPs
gfOuts <- gradientForest(cbind(envGF, SNPs_out), predictor.vars=colnames(envGF),
                         response.vars=colnames(SNPs_out), ntree=500, mtry=NULL,transform=NULL,
                         maxLevel=maxLevel, trace=T)

# plot output, see ?plot.gradientForest
type = "O"
pdf("Neut.GF.plot.pdf")
plot(gfNeut, plot.type=type)
dev.off()

pdf("Outs.GF.plot.pdf")
plot(gfOuts, plot.type=type)
dev.off()

#Get R2 variable importance (per env var)
R2.imp <- as.data.frame(gfNeut$imp.rsq)
write.table(R2.imp, file="CP.R2.imp.txt", sep="\t", quote = FALSE)

#Get R2 per allele
R2perspp <- as.data.frame(gfNeut$result)
write.table(R2perspp, file="CP.R2perspp.txt", sep="\t", quote = FALSE)

## Here need to extract importance numbers to plot separately (as heatmap):
Neut.imp <- as.data.frame(gfNeut$overall.imp)
Neut.r2<- as.data.frame(gfNeut$overall.imp2)

### Mapping heatmap of R2 weighted importance
library(ggplot2)
r2_heatmap <- read_excel("~/Desktop/PhD_stuffies/CH3/new_GF/r2.heatmap.xlsx")
heatmap <- ggplot(data = r2_heatmap, mapping = aes(x = Loci, y = Variable, fill = Imp)) + geom_tile() + scale_fill_gradient(name = "R2 weighted importance", low = "#FFFFFF", high = "#012345")
heatmap

# Plotting combined cumulative importance --------------------------------------------
#Can combine Neut and Outs GFs and plot combined cumulative importance
#These show cumulative
#change in abundance of individual species, where changes occur on the gradient, and the species
#changing most on each gradient

f12 <- combinedGradientForest(Neutral=gfNeut,Outlier=gfOuts)
prednames <- c("Trange", "SSSrange", "SSSmean", "SSTrange", "SSTmean")
pdf("GF.comb.imp.plot.pdf")
plot(f12,plot.type="Cumulative.Importance", sort=FALSE, prednames=prednames, show.weights=TRUE)
dev.off()

# Plotting biological space in 2 dimensions --------------------------------------------
### Plotting bi-plot
#Different coordinate positions in the biplot represent differing compositions, as associated
#with the predictors.

#CP
env_trns <- raster::extract(cont.env.CP, 1:ncell(cont.env.CP), df = TRUE)
#SG PA
env_trns <- raster::extract(cont.env.PA, 1:ncell(cont.env.CP), df = TRUE)

env_trns <- na.omit(env_trns)
colnames(env_trns)[2] <- "Trange"
colnames(env_trns)[3] <- "SSSrange"
colnames(env_trns)[4] <- "SSSmean"
colnames(env_trns)[5] <- "SSTrange"
colnames(env_trns)[6] <- "SSTmean"


colnames(env_trns)[2] <- "Tmean"
colnames(env_trns)[3] <- "Tmax"
colnames(env_trns)[4] <- "Trange"
colnames(env_trns)[5] <- "Pmean"
colnames(env_trns)[6] <- "Prange"
colnames(env_trns)[7] <- "SSSrange"
colnames(env_trns)[8] <- "SSSmean"
colnames(env_trns)[9] <- "SSTrange"
colnames(env_trns)[10] <- "SSTmean"
colnames(env_trns)[11] <- "CHLmean"

predNeut <- predict(gfNeut, env_trns[,-1])  #note the removal of the cell ID column with [,-1])
predOuts <- predict(gfOuts, env_trns[,-1])


pca <- prcomp(predNeut, center=TRUE, scale.=FALSE)

a1 <- pca$x[, 1]
a2 <- pca$x[, 2]
a3 <- pca$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255

nvs <- dim(pca$rotation)[1]
vec <- c("Trange", "SSSrange", "SSSmean", "SSTrange", "SSTmean")
lv <- length(vec)
vind <- rownames(pca$rotation) %in% vec
scal <- 40

xrng <- range(pca$x[, 1], pca$rotation[, 1]/scal) *
  +     + 1.1

yrng <- range(pca$x[, 2], pca$rotation[, 2]/scal) *
  +     + 1.1

pdf("PA.neut.biplot.labs.pdf")
pdf("PA.neut.biplot.no.labs.pdf")

pdf("PA.outs.biplot.labs.pdf")
pdf("PA.outs.biplot.no.labs.pdf")

plot((pca$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 8, col = rgb(r, g, b, max = 255), asp = 1)
points(pca$rotation[!vind, 1:2]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), pca$rotation[vec, + 1]/scal, pca$rotation[vec, 2]/scal, length = 0.1)

jit <- 0.0015
text(pca$rotation[vec, 1]/scal + jit * sign(pca$rotation[vec,
                                                         + 1]), pca$rotation[vec, 2]/scal + jit * sign(pca$rotation[vec,
                                                                                                                    + 2]), labels = vec)

dev.off()

# Mapping spatial genetic variation --------------------------------------------
###### Scripts adapted from Keller & Fitzpatrick et al. (2015)- Ecology Letters


###### functions to support mapping #####
# builds RGB(red,green,blue) raster from transformed environment
# snpPreds = dataframe of transformed variables from gf or gdm model
# rast = a raster mask to which RGB values are to be mapped
# cellNums = cell IDs to which RGB values should be assigned
pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}

# Function to map difference between spatial genetic predictions
# predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
# predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
# rast = a raster mask to which Procrustes residuals are to be mapped
# mapCells = cell IDs to which Procrustes residuals values should be assigned
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}



# OK, on to mapping. Script assumes:
# (1) a dataframe named "env_trns" containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask, named "mask" of the study region to which the RGB data will be written

#create mask of curr layer
mask<-cont.env.CP$X5env.curr.8

# make all values in mask =1
mask[]<-as.numeric(mask[]>0)
#could also do this: mask <- all.s[[1]]/all.s[[1]] - 1

# transform env using gf models, see ?predict.gradientForest
predNeut <- predict(gfNeut, env_trns[,-1]) # remove cell column before transforming (which we just added so 6th col)
predOuts <- predict(gfOuts, env_trns[,-1])

# map continuous variation - reference SNPs
refRGBmap <- pcaToRaster(predNeut, mask, env_trns$ID)
plotRGB(refRGBmap)
writeRaster(refRGBmap, "CP.neut.rgb.tif", format="GTiff", overwrite=TRUE)

# map continuous variation - outlier SNPs
outRGBmap <- pcaToRaster(predOuts, mask, env_trns$ID)
plotRGB(outRGBmap)
writeRaster(outRGBmap, "CP.outs.rgb.tif", format="GTiff", overwrite=TRUE)

# Difference between maps (GI5 and reference) 
diffGI5 <- RGBdiffMap(predNeut, predOuts, rast=mask, mapCells=env_trns$ID)
plot(diffGI5[[2]])
writeRaster(diffGI5[[2]], "CP.n.o.diff.tif", format="GTiff", overwrite=TRUE)

# Calculate and map "genetic offset" under climate change ----------------------
# Script assumes:
# (1) a dataframe of transformed env. variables for CURRENT climate 
# (e.g., predGI5 from above).
#
# (2) a dataframe named env_trns_future containing extracted raster data of 
# env. variables for FUTURE a climate scenario, same structure as env_trns


colnames(env_cont)[2] <- "Trange"
colnames(env_cont)[3] <- "SSSrange"
colnames(env_cont)[4] <- "SSSmean"
colnames(env_cont)[5] <- "SSTrange"
colnames(env_cont)[6] <- "SSTmean"

colnames(env_f.50)[2] <- "Trange"
colnames(env_f.50)[3] <- "SSSrange"
colnames(env_f.50)[4] <- "SSSmean"
colnames(env_f.50)[5] <- "SSTrange"
colnames(env_f.50)[6] <- "SSTmean"

colnames(env_f.90)[2] <- "Trange"
colnames(env_f.90)[3] <- "SSSrange"
colnames(env_f.90)[4] <- "SSSmean"
colnames(env_f.90)[5] <- "SSTrange"
colnames(env_f.90)[6] <- "SSTmean"


# first transform FUTURE env. variables
predNeut <- predict(gfNeut, env_cont[,-1]) 
predOuts <- predict(gfOuts, env_cont[,-1])

#do sequentially for 50 then for 90
projNeut <- predict(gfNeut, env_f.50[,-1])
projNeut <- predict(gfNeut, env_f.90[,-1])

projOuts <- predict(gfOuts, env_f.50[,-1])
projOuts <- predict(gfOuts, env_f.90[,-1])

# calculate euclidean distance between current and future genetic spaces 
## NB change to number of predictor variables used (in example there's 7)
genOffsetNeut.50 <- sqrt(((projNeut[,1]-predNeut[,1])^2)+((projNeut[,2]-predNeut[,2])^2)
                         +((projNeut[,3]-predNeut[,3])^2)+((projNeut[,4]-predNeut[,4])^2)
                         +((projNeut[,5]-predNeut[,5])^2)+((projNeut[,6]-predNeut[,6])^2))

genOffsetOuts.50 <- sqrt((projOuts[,1]-predOuts[,1])^2+(projOuts[,2]-predOuts[,2])^2
                         +(projOuts[,3]-predOuts[,3])^2+(projOuts[,4]-predOuts[,4])^2
                         +(projOuts[,5]-predOuts[,5])^2+((projNeut[,6]-predNeut[,6])^2))



genOffsetNeut.90 <- sqrt(((projNeut[,1]-predNeut[,1])^2)+((projNeut[,2]-predNeut[,2])^2)
                         +((projNeut[,3]-predNeut[,3])^2)+((projNeut[,4]-predNeut[,4])^2)
                         +((projNeut[,5]-predNeut[,5])^2)+((projNeut[,6]-predNeut[,6])^2))

genOffsetOuts.90 <- sqrt((projOuts[,1]-predOuts[,1])^2+(projOuts[,2]-predOuts[,2])^2
                         +(projOuts[,3]-predOuts[,3])^2+(projOuts[,4]-predOuts[,4])^2
                         +(projOuts[,5]-predOuts[,5])^2+((projNeut[,6]-predNeut[,6])^2))


# We then need to save the above as dataframe, bind with env_ DFs, and plot
#Nuetral SNPS
genOffsetNeut.50 <- as.data.frame(norm.genOffsetNeut.50)
gen.off.n50 <- cbind(env_fut.xy,genOffsetNeut.50)
gen.off.n50.ras <- rasterFromXYZ(gen.off.n50)
projection(gen.off.n50.ras) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(gen.off.n50.ras)
writeRaster(gen.off.n50.ras, "CP.gen.off.n50.tif", format="GTiff", overwrite=TRUE)

genOffsetNeut.90 <- as.data.frame(norm.genOffsetNeut.90)
gen.off.n90 <- cbind(env_fut.xy,genOffsetNeut.90)
gen.off.n90.ras <- rasterFromXYZ(gen.off.n90)
projection(gen.off.n90.ras) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(gen.off.n90.ras)
writeRaster(gen.off.n90.ras, "CP.gen.off.n90.tif", format="GTiff", overwrite=TRUE)

#Outlier SNPS
genOffsetOuts.50 <- as.data.frame(genOffsetOuts.50)
gen.off.o50 <- cbind(env_fut.xy,genOffsetOuts.50)
gen.off.o50.ras <- rasterFromXYZ(gen.off.o50)
projection(gen.off.o50.ras) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(gen.off.o50.ras)
writeRaster(gen.off.o50.ras, "CP.gen.off.o50.tif", format="GTiff", overwrite=TRUE)

genOffsetOuts.90 <- as.data.frame(genOffsetOuts.90)
gen.off.o90 <- cbind(env_fut.xy,genOffsetOuts.90)
gen.off.o90.ras <- rasterFromXYZ(gen.off.o90)
projection(gen.off.o90.ras) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(gen.off.o90.ras)
writeRaster(gen.off.o90.ras, "CP.gen.off.o90.tif", format="GTiff", overwrite=TRUE)



################################################################################
#### Code to run Species Distribution models (on each species individually) ####
################################################################################


##### FORECASTING WITH BIOMOD DATA (at ~10km resolution) #####

LIB <- c("rgbif", "biomod2", "ggplot2", "gridExtra", "knitr", "raster", 
         "ade4", "rworldmap", "cleangeo", "maptools", "rasterVis", "rgdal", "sdmpredictors", "usdm")
for(i in LIB) { install.packages(i, repos="http://ftp.sun.ac.za/ftp/pub/mirrors/cran.za.r/") ; library(i, character.only=T) }
for(i in LIB) { library(i, character.only=T) }

############################ Present Day layers ##########################

a.curr <- load_layers( layercodes = c("WC_bio5", "WC_bio6") , equalarea=FALSE, rasterstack=TRUE)
o.curr <- load_layers( layercodes = c("BO2_salinitymean_ss", "BO2_tempmean_ss") , equalarea=FALSE, rasterstack=TRUE)
a<-extent(-180, 180, -90, 90)
o<-extent(-180, 180, -90, 90)
extent_list<-list(a, o)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)
reso<-res(a.curr)
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=a.curr@crs)

a.c.2 <-resample(a.curr, s2, method="ngb")
o.c.2 <-resample(o.curr, s2, method="ngb")
curr <- stack(a.c.2, o.c.2)
writeRaster(curr, "4env.curr.tif", format="GTiff", overwrite=TRUE)


############################ 2070-2100 RCP layers ##########################


##RCP45 MIROC
a.mid <- load_layers( layercodes = c("WC_bio5_mr45_2050", "WC_bio6_mr45_2050") , equalarea=FALSE, rasterstack=TRUE)
o.mid <- load_layers( layercodes = c("BO2_RCP45_2050_salinitymean_ss", "BO2_RCP45_2050_tempmean_ss") , equalarea=FALSE, rasterstack=TRUE)
a.mid.2 <-resample(a.mid, s2, method="ngb")
o.mid.2 <-resample(o.mid, s2, method="ngb")
mr.50.45 <- stack(a.mid.2, o.mid.2)

writeRaster(mr.50.45, "mr.45.50.tif", format="GTiff", overwrite=TRUE)

##RCP45 CCSM4
a.mid <- load_layers( layercodes = c("WC_bio5_cc45_2050", "WC_bio6_cc45_2050") , equalarea=FALSE, rasterstack=TRUE)
o.mid <- load_layers( layercodes = c("BO2_RCP45_2050_salinitymean_bdmin", "BO2_RCP45_2050_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
a.mid.2 <-resample(a.mid, s2, method="ngb")
o.mid.2 <-resample(o.mid, s2, method="ngb")
cc.50.45 <- stack(a.mid.2, o.mid.2)

writeRaster(cc.50.45, "cc.45.50.tif", format="GTiff", overwrite=TRUE)

##RCP85 MIROC
miroc.a.85 <- load_layers( layercodes = c("WC_bio5_mr85_2050", "WC_bio6_mr85_2050") , equalarea=FALSE, rasterstack=TRUE)
o.85 <- load_layers( layercodes = c("BO2_RCP85_2050_salinitymean_bdmin", "BO2_RCP85_2050_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
a.85.2 <-resample(miroc.a.85, s2, method="ngb")
o.85.2 <-resample(o.85, s2, method="ngb")
mr.50.85 <- stack(a.85.2, o.85.2)

writeRaster(mr.50.85, "mr.85.50.tif", format="GTiff", overwrite=TRUE)

##RCP85 CCSM4
miroc.a.85 <- load_layers( layercodes = c("WC_bio5_cc85_2050", "WC_bio6_cc85_2050") , equalarea=FALSE, rasterstack=TRUE)
o.85 <- load_layers( layercodes = c("BO2_RCP85_2050_salinitymean_bdmin", "BO2_RCP85_2050_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
a.85.2 <-resample(miroc.a.85, s2, method="ngb")
o.85.2 <-resample(o.85, s2, method="ngb")
cc.50.85 <- stack(a.85.2, o.85.2)

writeRaster(cc.50.85, "cc.50.85.tif", format="GTiff", overwrite=TRUE)

############################ 2070-2100 RCP layers ##########################

## RCP 45 MIROC
miroc.R2080.t <- load_layers( layercodes = c("WC_bio5_mr45_2070", "WC_bio6_mr45_2070") , equalarea=FALSE, rasterstack=TRUE)
R2080.sst <- load_layers( layercodes = c("BO2_RCP45_2100_salinitymean_bdmin", "BO2_RCP45_2100_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
R2080.sst.2 <-resample(R2080.sst, s2, method="ngb")
R2080.t.2 <-resample(miroc.R2080.t, s2, method="ngb")
mr.90.45 <- stack(R2080.t.2, R2080.sst.2)

writeRaster(mr.90.45, "mr.45.90.tif", format="GTiff", overwrite=TRUE)

## RCP 45 CCSM
miroc.R2080.t <- load_layers( layercodes = c("WC_bio5_cc45_2070", "WC_bio6_cc45_2070"), equalarea=FALSE, rasterstack=TRUE)
R2080.sst <- load_layers( layercodes = c( "BO2_RCP45_2100_salinitymean_bdmin", "BO2_RCP45_2100_tempmean_bdmin"), equalarea=FALSE, rasterstack=TRUE)
R2080.sst.2 <-resample(R2080.sst, s2, method="ngb")
R2080.t.2 <-resample(miroc.R2080.t, s2, method="ngb")
cc.90.45 <- stack(R2080.t.2, R2080.sst.2)

writeRaster(cc.90.45, "cc.45.90.tif", format="GTiff", overwrite=TRUE)

## RCP 8.5 MIROC
miroc.R2080.t <- load_layers( layercodes = c("WC_bio5_mr85_2070", "WC_bio6_mr85_2070") , equalarea=FALSE, rasterstack=TRUE)
R2080.sst <- load_layers( layercodes = c("BO2_RCP85_2100_salinitymean_bdmin", "BO2_RCP85_2100_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
R2080.sst.2 <-resample(R2080.sst, s2, method="ngb")
R2080.t.2 <-resample(miroc.R2080.t, s2, method="ngb")
mr.90.85 <- stack(R2080.t.2, R2080.sst.2)

writeRaster(mr.90.85, "mr.85.90.tif", format="GTiff", overwrite=TRUE)

## RCP 8.5 CCSM4
miroc.R2080.t <- load_layers( layercodes = c("WC_bio5_cc85_2070", "WC_bio6_cc85_2070") , equalarea=FALSE, rasterstack=TRUE)
R2080.sst <- load_layers( layercodes = c("BO2_RCP85_2100_salinitymean_bdmin", "BO2_RCP85_2100_tempmean_bdmin") , equalarea=FALSE, rasterstack=TRUE)
R2080.sst.2 <-resample(R2080.sst, s2, method="ngb")
R2080.t.2 <-resample(miroc.R2080.t, s2, method="ngb")
cc.90.85 <- stack(R2080.t.2, R2080.sst.2)

writeRaster(cc.90.85, "cc.85.90.tif", format="GTiff", overwrite=TRUE)

t.50.45 <- overlay(mr.50.45, cc.50.45, fun=mean)
t.50.85 <- overlay(mr.50.85, cc.50.85, fun=mean)
t.90.45 <- overlay(mr.90.45, cc.90.45, fun=mean)
t.90.85 <- overlay(mr.90.85, cc.90.85, fun=mean)

writeRaster(t.50.45, "fut.50.45.tif", format="GTiff", overwrite=TRUE)
writeRaster(t.50.85, "fut.50.85.tif", format="GTiff", overwrite=TRUE)
writeRaster(t.90.45, "fut.90.45.tif", format="GTiff", overwrite=TRUE)
writeRaster(t.90.85, "fut.90.85.tif", format="GTiff", overwrite=TRUE)

#crop to species
SA.cp.ext <- extent(10, 40, -37, -22)
SA.pa.ext <- extent(10, 40, -37, -19)
SA.sg.ext <- extent(10, 40, -37, -12)

##crop rasters
curr.4 <- crop(curr, SA.cp.ext)
curr.4=stack(curr.4)
plot(curr.4)


#2050
t.50.45 <- 'fut.50.45.tif'
t.50.45=stack(t.50.45)

t.50.45 <- crop(t.50.45, SA.cp.ext)
t.50.45=stack(t.50.45)  
names(t.50.45) <- c("X4env.curr.1", "X4env.curr.2 ", "X4env.curr.3", "X4env.curr.4")
plot(t.50.45)

t.50.85 <- 'fut.50.85.tif'
t.50.85=stack(t.50.85)

t.50.85 <- crop(t.50.85, SA.cp.ext)
t.50.85=stack(t.50.85)  
names(t.50.85) <- c("X4env.curr.1", "X4env.curr.2 ", "X4env.curr.3", "X4env.curr.4")
plot(t.50.85)

#2090
t.90.45 <- 'fut.90.45.tif'
t.90.45=stack(t.90.45)

t.90.45 <- crop(t.90.45, SA.cp.ext)
t.90.45=stack(t.90.45)  
names(t.90.45) <- c("X4env.curr.1", "X4env.curr.2 ", "X4env.curr.3", "X4env.curr.4")
plot(t.90.45)

t.90.85 <- 'fut.90.85.tif'
t.90.85=stack(t.90.85)

t.90.85 <- crop(t.90.85, SA.cp.ext)
t.90.85=stack(t.90.85)  
names(t.90.85) <- c("X4env.curr.1", "X4env.curr.2 ", "X4env.curr.3", "X4env.curr.4")
plot(t.90.85)

t.50.45 <- stack(t.50.45)
t.50.85 <- stack(t.50.85)
t.90.45 <- stack(t.90.45)
t.90.45 <- stack(t.90.45)


########################################################################
#####################      Biomod formatting     ######################
########################################################################

#convert the column named "PRESENCE" to a character class
myResp<-as.numeric(Crabs_xy$Cpunctatus)

myRespName <- 'Cpunctatus'

#create presence and pseudo absences
SPC_PresAbs <- BIOMOD_FormatingData(resp.var = myResp,
                                    expl.var = 4env.curr,
                                    resp.xy = Crabs_xy[,c('x', 'y')],
                                    resp.name = myRespName,
                                    PA.nb.rep = 3,
                                    PA.nb.absences = 15000,
                                    PA.strategy = 'random') 

SPC_PresAbs
plot(SPC_PresAbs)


# Set modeling options
MySpc_options <- BIOMOD_ModelingOptions(
  GLM = list( type = 'quadratic', interaction.level = 1 ),
  GBM = list( n.trees = 1000 ),
  GAM = list( algo = 'GAM_mgcv' ) )

# Set up models
MySpc_models <- BIOMOD_Modeling( data = SPC_PresAbs,
                                 models = c("GLM","GAM", "GBM", "RF","MARS", "FDA"),
                                 models.options = MySpc_options,
                                 NbRunEval = 10,
                                 DataSplit = 70,
                                 VarImport = 3,
                                 models.eval.meth=c('TSS','ROC'),
                                 do.full.models = F )

########################################################################
##################        Evaluate models          ####################
########################################################################

# Get models evaluation scores
MyModels_scores <- get_evaluations(MySpc_models)
dim(MyModels_scores)
dimnames(MyModels_scores)

#Graphically see model scores
models_scores_graph(MySpc_models, by = "models" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(MySpc_models, by = "cv_run" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(MySpc_models, by = "data_set" , metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))

# The predictive accuracy of the models are good when the AUC (area under the curve- here ROC) 
#  0.8 and TSS  0.65.`RF`/ 'GLM' models seems to be the most accurate ones on average,
# followed by `GBM`then `GAM` 
# You can also view scores numerically...

MyModels_scores["ROC","Testing.data",,,]
MyModels_scores["TSS","Testing.data",,,]


# Variable importance
# The higher a score is, the more important is the variable. We will visualize this as a barplot

MyModels_var_import <- get_variables_importance(MySpc_models)
MyModels_var_import
dimnames(MyModels_var_import)

# Average variable importance by algorithm
mVarImp <- apply(MyModels_var_import, c(1,2), median) 
mVarImp <- apply(mVarImp, 2, function(x) x*(1/sum(x))) # standardize the data
mVarImp 
#Visualize this as a bar plot
VarImpBarPlot <- barplot(mVarImp, legend.text=row.names(mVarImp), xlim=c(0,7))


# Response curves
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# To analyze how each environmental variable influences modelled probability of species presence,
# we will use an evaluation procedure proposed by Elith et al.(2005). 

#A plot of these predictions allows visualisation of the modeled response(y-axis) 
# to the given variable (x-axis),conditional to the other variables being held constant.

# We have first to name and load the produced models.
MySpc_glm <- BIOMOD_LoadModels(MySpc_models, models='GLM')
MySpc_gam <- BIOMOD_LoadModels(MySpc_models, models='GAM')
MySpc_gbm <- BIOMOD_LoadModels(MySpc_models, models='GBM')
MySpc_rf  <- BIOMOD_LoadModels(MySpc_models, models='RF')

glm_eval_strip <- biomod2::response.plot2(
  models  = MySpc_glm, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

gam_eval_strip <- biomod2::response.plot2(
  models  = MySpc_gam, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

gbm_eval_strip <- biomod2::response.plot2(
  models  = MySpc_gbm, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

rf_eval_strip <- biomod2::response.plot2(
  models  = MySpc_rf, Data = get_formal_data(MySpc_models,'expl.var'), 
  show.variables= get_formal_data(MySpc_models,'expl.var.names'),
  do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE,
  display_title = FALSE, data_species = get_formal_data(MySpc_models,'resp.var'))

# Map model prediction on the current South African climate
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_proj_current <- BIOMOD_Projection( modeling.output = MySpc_models,
                                                new.env = 4env.curr,
                                                proj.name = "current",
                                                binary.meth = "ROC",
                                                output.format = ".img",
                                                do.stack = FALSE )

#A list of all the models that were just executed
MySpc_models_proj_current

# Plot and compare the maps for the potential current distribution projected by the different 
# models.
plot(MySpc_models_proj_current,  str.grep="PA1_RUN1")


########################################################################
#####################   Ensemble modelling       ######################
########################################################################

#Run ensemble models
MySpc_ensemble_models <- BIOMOD_EnsembleModeling( modeling.output = MySpc_models,
                                                  chosen.models ='all',
                                                  em.by = 'all',  #combine all models
                                                  eval.metric = 'all',
                                                  eval.metric.quality.threshold = c(0.55,0.8),
                                                  models.eval.meth = c('TSS','ROC'),
                                                  prob.mean = FALSE,
                                                  prob.cv = TRUE, #coefficient of variation across predictions
                                                  committee.averaging = TRUE,
                                                  prob.mean.weight = TRUE,
                                                  VarImport = 0 )

#check scores
MySpc_ensemble_models_scores <- get_evaluations(MySpc_ensemble_models)
MySpc_ensemble_models_scores

# Ensemble model forecasts
# ...............................................
MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_current,
  binary.meth = "ROC",  #make binary predictions (pres/abs) based on ROC score
  output.format = ".img",
  do.stack = FALSE )

# The projections for current conditions are stored in the 'proj_current' directory. 
list.files(paste(SDM1, "/proj_current/individual_projections", sep=""))
get_projected_models(MySpc_ensemble_models_proj_current)

#Plot them all--- this may take a while---
plot(MySpc_ensemble_models_proj_current)
plot(MySpc_ensemble_models_proj_current, str.grep="EMcaByTSS")
plot(MySpc_ensemble_models_proj_current, str.grep="EMwmeanByTSS")


########################################################################
##############  Model projection with past conditions  ################
########################################################################


# Run the projections
# ...............................................
MySpc_models_proj_50.45 <- BIOMOD_Projection( modeling.output = MySpc_models,
                                          new.env = t.50.45,
                                          proj.name = "t.50.45",
                                          binary.meth = c("ROC"),
                                          output.format = ".img",
                                          do.stack = FALSE)

MySpc_ensemble_models_proj_50.45 <- BIOMOD_EnsembleForecasting( 
  EM.output = MySpc_ensemble_models,
  projection.output = MySpc_models_proj_50.45,
  binary.meth = "ROC",
  output.format = ".img",
  do.stack = FALSE,
  build.clamping.mask=F)

# Repeat above for other future timepoints, and for each species. 

############################## Plotting future SDMs ############################

library(rasterVis)
library(gridExtra)
library(rgdal)
library(tmap)
library(viridisLite)

#Upload outputs

MyBinCP_curr <- raster::stack("Cpunctatus/Cpunctatus/proj_current/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_50.45 <- raster::stack("Cpunctatus/Cpunctatus/proj_CP.50.85/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_90.45 <- raster::stack("Cpunctatus/Cpunctatus/proj_CP.90.85/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_50.85 <- raster::stack("Cpunctatus/Cpunctatus/proj_CP.50.85/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinCP_90.85 <- raster::stack("Cpunctatus/Cpunctatus/proj_CP.90.85/individual_projections/Cpunctatus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

MyBinPA_curr <- raster::stack("Pangulosus2/proj_current/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_50.45 <- raster::stack("Pangulosus2/proj_PA.50.85/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_90.45 <- raster::stack("Pangulosus2/proj_PA.90.85/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_50.85 <- raster::stack("Pangulosus2/proj_PA.50.85/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinPA_90.85 <- raster::stack("Pangulosus2/proj_PA.90.85/individual_projections/Pangulosus_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

MyBinSG_curr <- raster::stack("Sgranularis/Sgranularis/proj_current/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_50.45 <- raster::stack("Sgranularis/Sgranularis/proj_SG.50.85/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_90.45 <- raster::stack("Sgranularis/Sgranularis/proj_SG.90.85/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_50.85 <- raster::stack("Sgranularis/Sgranularis/proj_SG.50.85/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
MyBinSG_90.85 <- raster::stack("Sgranularis/Sgranularis/proj_SG.90.85/individual_projections/Sgranularis_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

#Merge RCPs so we have one map per time point, and aggregate for easy viewing
CP.50 <- merge(MyBinCP_50.45, MyBinCP_50.85)
CP.90 <- merge(MyBinCP_90.45, MyBinCP_90.85)
CP.ens <- stack(MyBinCP_curr, CP.50, CP.90)
CP.ens.agg <- aggregate(CP.ens, fact=3)
names(CP.ens.agg) <- c('Current', '2050','2070/2100')


PA.50 <- merge(MyBinPA_50.45, MyBinPA_50.85)
PA.90 <- merge(MyBinPA_90.45, MyBinPA_90.85)
PA.ens <- stack(MyBinPA_curr, PA.50, PA.90)
PA.ens.agg <- aggregate(PA.ens, fact=3)
names(PA.ens.agg) <- c('Current', '2050','2070/2100')


SG.50 <- merge(MyBinSG_50.45, MyBinSG_50.85)
SG.90 <- merge(MyBinSG_90.45, MyBinSG_90.85)
SG.ens <- stack(MyBinSG_curr, SG.50, SG.90)
SG.ens.agg <- aggregate(SG.ens, fact=3)
names(SG.ens.agg) <- c('Current', '2050','2070/2100')

#Extend PA and CP to SG extent for stacking
CP.ens.agg <- extend(CP.ens.agg, SG.ens.agg)
PA.ens.agg <- extend(PA.ens.agg, SG.ens.agg)

# Do for RCP 45 and 85 separate
PA.45 <- stack(MyBinPA_curr, MyBinPA_50.45, MyBinPA_90.45)
PA.85 <- stack(MyBinPA_curr, MyBinPA_50.85, MyBinPA_90.85)
PA.45.agg <- aggregate(PA.45, fact=3)
PA.85.agg <- aggregate(PA.85, fact=3)
names(PA.45.agg) <- c('Current', '2050','2070/2100')
names(PA.85.agg) <- c('Current', '2050','2070/2100')

SG.45 <- stack(MyBinSG_curr, MyBinSG_50.45, MyBinSG_90.45)
SG.85 <- stack(MyBinSG_curr, MyBinSG_50.85, MyBinSG_90.85)
SG.45.agg <- aggregate(SG.45, fact=3)
SG.85.agg <- aggregate(SG.85, fact=3)
names(SG.45.agg) <- c('Current', '2050','2070/2100')
names(SG.85.agg) <- c('Current', '2050','2070/2100')

#Extend PA and CP to SG extent for stacking
CP.45.agg <- extend(CP.45.agg, SG.45.agg)
PA.45.agg <- extend(PA.45.agg, SG.45.agg)

CP.85.agg <- extend(CP.85.agg, SG.85.agg)
PA.85.agg <- extend(PA.85.agg, SG.85.agg)

#Stack species
spp.raster <- stack(CP.ens.agg, PA.ens.agg, SG.ens.agg)

spp.raster.45 <- stack(CP.45.agg, PA.45.agg, SG.45.agg)
spp.raster.85 <- stack(CP.85.agg, PA.85.agg, SG.85.agg)

#Upload background map
setwd("~/Desktop/PhD_stuffies/SDMS/SDM.shps")
map1 <- readOGR("Africa.shp")
SA.sg.ext <- extent(10, 40, -37, -12)
map4 <- crop(map1, SA.sg.ext)
map5 <- tm_shape(map4) + tm_borders(col = "lightgrey")+ tm_fill(col = "white")


#Set colors
vir <- viridis(16, direction = -1)

#Plotting with LatticeExtra
setwd("~/Desktop/PhD_stuffies/CH3/SDMs_june")
pdf("3spp.fut.SDMs.pdf")
levelplot(spp.raster, col.regions=vir, layout=c(3, 3)) + layer(sp.polygons(map4, fill='white', alpha=0.2))
dev.off()

#Plotting with tmap
setwd("~/Desktop/PhD_stuffies/CH3/SDMs_june")
pdf("3spp.fut.SDMs.pdf")
map5 + tm_shape(spp.raster) + tm_raster(style = "cont", palette = vir)+tm_legend(outside=TRUE) 
dev.off()

################################################################################
######## Code to create convex hulls (code adapted from Mark Miller) ##########
################################################################################

library(FD)
library(ade4)
library(ggplot2)
library(raster)
library(cowplot)
library("ggpubr")


# upload current layers (MARSPEC)
curr <- 'contemp.4vars.tif'
curr <- stack(curr)

# upload past layers
cc6 <- 'ccsm.6.crop.tif'
cc6 <- stack(cc6)

mr6 <- 'miroc.6.crop.tif'
mr6 <- stack(mr6)

p.6 <- merge(cc6, mr6)

cc21 <- 'ccsm.21.tif'
cc21 <- stack(cc21)

mr21 <- 'miroc.nc.21.tif'
mr21 <- stack(mr21)

p.21 <- merge(cc21, mr21)

# upload future layers
cc.50.45 <- 'cc.50.45.tif'
cc.50.45=stack(cc.50.45)

cc.50.85 <- 'cc.50.85.tif'
cc.50.85=stack(cc.50.85)

cc.90.45 <- 'cc.90.45.tif'
cc.90.45=stack(cc.90.45)

t.90.85 <- 'fut.90.85.tif'
t.90.85=stack(t.90.85)

SA.sg.ext <- extent(10, 40, -37, -12)
t.50.45.c <- crop(t.50.45, SA.sg.ext)
t.50.85.c <- crop(t.50.85, SA.sg.ext)
t.90.45.c <- crop(t.90.45, SA.sg.ext)
t.90.85.c <- crop(t.90.85, SA.sg.ext)

t.50 <- merge(t.50.45.c, t.50.85.c)
t.90 <- merge(t.90.45.c, t.90.85.c)

p6.enf <- extract(x=p.6, y=cbind(curr_hull_xy$X, curr_hull_xy$Y))
p6.enf.sea <- extract(x=p.6, y=cbind(curr_sea$X, curr_sea$Y))
p6.enf.land <- extract(x=p.6, y=cbind(curr_land$X, curr_land$Y))

write.table(p6.enf, file="MH.enf.1.txt", sep="\t", quote = FALSE)
write.table(p6.enf.land, file="MH.enf.L.txt", sep="\t", quote = FALSE)
write.table(p6.enf.sea, file="MH.enf.S.txt", sep="\t", quote = FALSE)


# create Euclidean matrix
env.mat <- dist(hulls_input[,1:4], method = "euclidean", diag = FALSE, upper=FALSE)

# principal coord analyses of Gower Dist trait data, 4 axes
pc2<-dudi.pco(d = env.mat, scannf = FALSE, nf = 4) 

pc2.dfs <- data.frame(pc2$li, hulls_input) # combine PCoA axes with trait/site data

#make global hull of first 2 PCs
glob_hull<-pc2.dfs[chull(pc2.dfs$A1, pc2.dfs$A2),]

# code to setup ggplot enviroment
ppp <- ggplot() + coord_fixed() +
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") +
  geom_vline(xintercept=0, col="darkgrey")

# plot global hull
ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+theme_bw()

# make regional hulls
pc2_subset <- pc2.dfs[ which(pc2.dfs$Time=='CURR'),]
reg_hull<-pc2_subset[chull(pc2_subset$A1, pc2_subset$A2),]
Curr <- ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+
  geom_polygon(data=reg_hull,aes(x=A1,y=A2),alpha=0.08, fill='green',colour="black")+
  geom_point(data=pc2_subset, aes(x=A1, y=A2), colour='green')+theme_bw()

pc2_subset <- pc2.dfs[ which(pc2.dfs$Time=='MH'),]
reg_hull<-pc2_subset[chull(pc2_subset$A1, pc2_subset$A2),]
MH <- ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+
  geom_polygon(data=reg_hull,aes(x=A1,y=A2),alpha=0.08, fill='blue',colour="black")+
  geom_point(data=pc2_subset, aes(x=A1, y=A2), colour='blue')+theme_bw()

pc2_subset <- pc2.dfs[ which(pc2.dfs$Time=='LGM'),]
reg_hull<-pc2_subset[chull(pc2_subset$A1, pc2_subset$A2),]
LGM <- ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+
  geom_polygon(data=reg_hull,aes(x=A1,y=A2),alpha=0.08, fill='purple',colour="black")+
  geom_point(data=pc2_subset, aes(x=A1, y=A2), colour='purple')+theme_bw()

pc2_subset <- pc2.dfs[ which(pc2.dfs$Time=='2050'),]
reg_hull<-pc2_subset[chull(pc2_subset$A1, pc2_subset$A2),]
f.50 <- ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+
  geom_polygon(data=reg_hull,aes(x=A1,y=A2),alpha=0.08, fill='orange',colour="black")+
  geom_point(data=pc2_subset, aes(x=A1, y=A2), colour='orange')+theme_bw()

pc2_subset <- pc2.dfs[ which(pc2.dfs$Time=='2090'),]
reg_hull<-pc2_subset[chull(pc2_subset$A1, pc2_subset$A2),]
f.90 <- ppp+
  geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
  geom_point(data=pc2.dfs, aes(x=A1, y=A2), colour='grey70')+
  geom_polygon(data=reg_hull,aes(x=A1,y=A2),alpha=0.08, fill='red',colour="black")+
  geom_point(data=pc2_subset, aes(x=A1, y=A2), colour='red')+theme_bw()

pdf("hull.plots.pdf")
plot_grid(LGM, MH, Curr, f.50, f.90, label_size = 12, ncol = 2)
dev.off()


### Plot PCA with loadings
#with dudi.pca
Y.pca <- dudi.pca(hulls_input[,1:4], scannf=F, nf=4)
pdf("dudi.pdf")
scatter(Y.pca,
        posieig = "none", # Hide the scree plot
        clab.row = 0      # Hide row labels
)
dev.off()

groups <- as.factor(hulls_input$Time)
pdf("dudi.groups.pdf")
fviz_pca_ind(Y.pca,
             col.ind = groups, # color by groups
             palette = c("red", "orange", "green", "blue", "purple"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = F
)
dev.off()


pca1.dfs <- data.frame(Y.pca$li, hulls_input)
ppp + geom_point(data=pca1.dfs, aes(x=Axis1, y=Axis2, col=Time))

#with ggfortify
library(ggfortify)
df <- hulls_input[c(1, 2, 3, 4)]

pdf("dudi.gg.lablels.pdf")
autoplot(prcomp(df), data = hulls_input, colour = 'Time',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE,  loadings.label.size = 3)
dev.off()

pdf("dudi.gg.no.lablels.pdf")
autoplot(prcomp(df), data = hulls_input, colour = 'Time',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = F)
dev.off()


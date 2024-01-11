## Import libraries and packages

library(sp)
library(rgdal)
library(raster)
library(ggplot2)
library(viridis)
library(rasterVis)
library(ggplot2)
library(tidyr)
library(ClusterR)
library(cluster)
library(e1071)
library(tidyverse)
library(broom)
library(pROC) 

## Function to convert a rasterbrick to a NDVI
## Input: brick_img = a rasterbrick image, k = band corresponding to NIR, i = band corresponding to Red 
## Output =  A dataframe of the NDVI index for each pixel in the image
VI <- function(brick_img, k, i) {
  
    bk <- brick_img[[k]] ## Extract NIR rasterlayer from rasterbrick
    bi <- brick_img[[i]] ## Extract Red rasterlayer from rasterbrick
    vi <- (bk - bi) / (bk + bi) ## Calculate NDVI

    return(vi) ## Return NDVI for each pixel
  
}

## Function to convert a UTM coordinate system to a Longitude/Latitude system 
## Input: file1 = a tif file representing sentinal 1 data, file2 = a tif file representing sentinal 2 data
## Output =  A new rasterbrick object that re-projects file to to a Longitude/Latitude coordinate system and resamples the pixels to ensure both images hav ethe same number of cells
UTM_to_LongLat <- function(file1, file2) {
  
  setwd("~/Documents/Edinburgh/Dissertation/Space/SDS_Project_2022_FillingTheGap/") ## Set the working directory
  
  s1_raster <- brick(paste("s1/", file1, sep="")) ## Convert file1 into a rasterbrick object
  s2_raster <- brick(paste("s2/", file2, sep="")) ## Convert file2 into a rasterbrick object
  
  reprojectedData2 <- projectRaster(s2_raster, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ") ## Reproject file2 
  
  new_raster <- resample(reprojectedData2, s1_raster, method='bilinear') ## Resample file2
  
  return(new_raster) ## Return new rasterbrick object
  
}


## Function to convert two rasterbrick objects into a dataframe with the VV, VH and NDVI values for each cell
## Input: brick1 = A raster brick object represnting sentinal 1 data, brick2 = A raster brick object represnting sentinal 2 data
## Output =  A new dataframe containing VV, VH, NDVI and vegetation classification for the image 
bricks_to_df <- function(brick1, brick2) {
  
  b1 <- brick1[[1]] ## Get raster layer representing VV index
  b1_df <- as.data.frame(b1, xy = T) ## Convert VV index into a dataframe
  b2 <- brick1[[2]] ## Get raster layer representing VH index
  b2_df <- as.data.frame(b2, xy = T) ## Convert VH index into a dataframe
  
  NDVI <- VI(brick2, 4, 3) ## Generate NDVI index from sentinal 2 rasterbrick
  rdf_ndvi <- as.data.frame(NDVI, xy = T) #Convert NDVI index rasterlayer into a dataframe
  
  df <- data.frame("NDVI" = rdf_ndvi$layer,"VV" = b1_df$VV, "VH" = b2_df$VH) ## Merge VV, VH and NDVI index dataframes
  df$Veg <- rep(1, nrow(df)) ## Storage vector to store veg classification 
  df$Veg[df$NDVI < 0.4] <- 0 ## Generate classification for vegetation based on NDVI index 
  df$Veg[is.na(df$NDVI)] <- NA ## Convert all NAs for NDVI index into NAs for veg classification 
  
  return(df) ## Return dataframe 
  
}


## Function to impute NAs from a dataframe based on surrounding pixels
## Input: df = A dataframe that contains information about the VV, VH and Index index. Usually an output from the above function
## Output =  A new dataframe with imputed VV, VH and NDVI values
impute <- function(df) {
  
  missing <- which(is.na(df$VV)) ## Vector of rows for which the VV index is missing
  
  for (i in missing) { ## Loop through the elements of the vector to calculate the value to impute 
    
    impute <- c(df$VV[i-5], df$VV[i-4], df$VV[i-3], df$VV[i-2], df$VV[i-1], df$VV[i+1], df$VV[i+2], df$VV[i+3], df$VV[i+4], df$VV[i+5]) ## Generate a vector of the surrounding 10 pixels 
    
    df$VV[i] <- mean(impute, na.rm = TRUE) ## Calculate the mean of the surrounding 10 pixels, ignoring any NA values 
    
  }
  
  missing2 <- which(is.na(df$VH)) ## Vector of rows for which the VH index is missing
  
  for (i in missing2) { ## Loop through the elements of the vector to calculate the value to impute 
    
    impute <- c(df$VH[i-5], df$VH[i-4], df$VH[i-3], df$VH[i-2], df$VH[i-1], df$VH[i+1], df$VH[i+2], df$VH[i+3], df$VH[i+4], df$VH[i+5]) ## Generate a vector of the surrounding 10 pixels 
    
    df$VH[i] <- mean(impute, na.rm = TRUE) ## Calculate the mean of the surrounding 10 pixels, ignoring any NA values
    
  }
  
  missing3 <- which(is.na(df$NDVI)) ## Vector of rows for which the VH index is missing
  
  for (i in missing3) { ## Loop through the elements of the vector to calculate the value to impute 
    
    impute <- c(df$NDVI[i-5], df$NDVI[i-4], df$NDVI[i-3], df$NDVI[i-2], df$NDVI[i-1], df$NDVI[i+1], df$NDVI[i+2], df$NDVI[i+3], df$NDVI[i+4], df$NDVI[i+5]) ## Generate a vector of the surrounding 10 pixels 
    
    df$NDVI[i] <- mean(impute, na.rm = TRUE) ## Calculate the mean of the surrounding 10 pixels, ignoring any NA values
    
  }
  
  
  return(df) ## Dataframe with imputed values 
  
}

## Function to generate the red and blue correlation coefficient between two sentinal datasets  
## Input: brick1 = A rasterbrick representing a sentinal 1 dataset, brick2 = A rasterbrick representing a Sentinal 2 dataset
## Output =  A vector of length 2 where the first item represents the red correlation between the two sentinal images and the second item represents the blue correlation between the two sentinal images 
correl <- function(brick1, brick2){
  
  s1_df <- as.data.frame(brick1[[1]], xy = T) ## Convert the sentinal 1 VV band image to a dataframe
  s2_df <- as.data.frame(brick2[[3]], xy = T) ## Convert the sentinal 2 red band image to a dataframe
  red_corr <- cor(s1_df$VV, s2_df$B4, use="complete.obs") ## Generate a correlation coefficient between VV polarization index and red spectral band
  
  s1_df <- as.data.frame(brick1[[2]], xy = T) ## Convert the sentinal 1 VH band image to a dataframe
  s2_df <- as.data.frame(brick2[[2]], xy = T) ## Convert the sentinal 1 blue band image to a dataframe
  blue_corr <- cor(s1_df$VH, s2_df$B3, use="complete.obs") ## Generate a correlation coefficient between VH polarization index and blue spectral band
  
  corr <- c(red_corr, blue_corr) ## Return vector of correlation coefficients
  
}

## Function to synthesise clouds by replacing cells with NA values 
## Input: brick = A rasterbrick object  representing a Sentinal 2 image, extent: an object of class extent corresponding to the area of a cloud to be synthesied
## Output =  A new rasterbrick with the area represented by the extent replaced with NAs
Simulate_Clouds <- function(brick, extent) {
  
  b1 <- brick[[1]] ## Convert the rasterbrick into a rasterlayer for the first band
  b2 <- brick[[2]] ## Convert the rasterbrick into a rasterlayer for the second band
  b3 <- brick[[3]] ## Convert the rasterbrick into a rasterlayer for the third band
  b4 <- brick[[4]] ## Convert the rasterbrick into a rasterlayer for the fourth band
  
  b1_df <- as.data.frame(b1, xy=TRUE) ## Convert the rasterlayer into a df for the first band
  b2_df <- as.data.frame(b2, xy=TRUE) ## Convert the rasterlayer into a df for the second band
  b3_df <- as.data.frame(b3, xy=TRUE) ## Convert the rasterlayer into a df for the thrid band
  b4_df <- as.data.frame(b4, xy=TRUE) ## Convert the rasterlayer into a df for the fourth band
  
  b1_df$B2[b1_df$x > extent[1] & b1_df$x < extent[2] & b1_df$y < extent[4] & b1_df$y > extent[3]] <- NA ## Replace all picels within the extent boudaries with NA for first band
  b2_df$B3[b1_df$x > extent[1] & b1_df$x < extent[2] & b1_df$y < extent[4] & b1_df$y > extent[3]] <- NA ## Replace all picels within the extent boudaries with NA for second band
  b3_df$B4[b1_df$x > extent[1] & b1_df$x < extent[2] & b1_df$y < extent[4] & b1_df$y > extent[3]] <- NA ## Replace all picels within the extent boudaries with NA for thrid band
  b4_df$B8[b1_df$x > extent[1] & b1_df$x < extent[2] & b1_df$y < extent[4] & b1_df$y > extent[3]] <- NA ## Replace all picels within the extent boudaries with NA for fourth band
  
  raster1 <- rasterFromXYZ(b1_df, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ") ## Convert new dataframes with synthesised NA values/clouds to raster layer for band 1
  raster2 <- rasterFromXYZ(b2_df, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ") ## Convert new dataframes with synthesised NA values/clouds to raster layer for band 2
  raster3 <- rasterFromXYZ(b3_df, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ") ## Convert new dataframes with synthesised NA values/clouds to raster layer for band 3
  raster4 <- rasterFromXYZ(b4_df, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ") ## Convert new dataframes with synthesised NA values/clouds to raster layer for band 4
  
  brick <- brick(raster1, raster2, raster3, raster4) ## Combine rasterlayers to form new rasterbrick
  
  return(brick) ## Return rasterbrick
  
}

## CONVERT .TIF IMAGES TO RASTERBRICKS

setwd("~/Documents/Edinburgh/Dissertation/Space/SDS_Project_2022_FillingTheGap/") ## Set workingdirector to where all rasterbricks are located

s1_raster_1 <- brick('s1/s1_aoi_1.tif') ## Convert first location image from tif to rasterbrick
s2_raster_1 <- UTM_to_LongLat('s1_aoi_1.tif', 's2_aoi_1.tif') ## Convert first location image from UTM coordinates to Longitude Latutude coordinates
s1_raster_1 <- aggregate(s1_raster_1, 3,fun=mean) ## Aggregate first image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_1 <- aggregate(s2_raster_1, 3,fun=mean) ## Aggregate first image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_1, s1_raster_1) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_2 <- brick('s1/s1_aoi_2.tif') ## Convert second location image from tif to rasterbrick
s2_raster_2 <- UTM_to_LongLat('s1_aoi_2.tif', 's2_aoi_2.tif') ## Convert second location image from UTM coordinates to Longitude Latutude coordinates
s1_raster_2 <- aggregate(s1_raster_2, 3,fun=mean) ## Aggregate second image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_2 <- aggregate(s2_raster_2, 3,fun=mean) ## Aggregate secondf image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_2, s1_raster_2) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_3 <- brick('s1/s1_aoi_3.tif') ## Convert thrid location image from tif to rasterbrick
s2_raster_3 <- UTM_to_LongLat('s1_aoi_3.tif', 's2_aoi_3.tif') ## Convert thrid location image from UTM coordinates to Longitude Latutude coordinates
s1_raster_3 <- aggregate(s1_raster_3, 3,fun=mean) ## Aggregate third image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_3 <- aggregate(s2_raster_3, 3,fun=mean) ## Aggregate thrid image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_3, s1_raster_3) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_4 <- brick('s1/s1_aoi_4.tif') ## Convert fourth location image from tif to rasterbrick
s2_raster_4 <- UTM_to_LongLat('s1_aoi_4.tif', 's2_aoi_4.tif') ## Convert fourth location image from UTM coordinates to Longitude Latutude coordinates
s1_raster_4 <- aggregate(s1_raster_4, 3,fun=mean) ## Aggregate fourth image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_4 <- aggregate(s2_raster_4, 3,fun=mean) ## Aggregate fourth image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_4, s1_raster_4) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_5 <- brick('s1/s1_aoi_5.tif') ## Convert fifth location image from tif to rasterbrick
s2_raster_5 <- UTM_to_LongLat('s1_aoi_5.tif', 's2_aoi_5.tif') ## Convert fifth location image from UTM coordinates to Longitude Latutude coordinates
s1_raster_5 <- aggregate(s1_raster_5, 3,fun=mean) ## Aggregate fifth image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_5 <- aggregate(s2_raster_5, 3,fun=mean) ## Aggregate fifth image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_5, s1_raster_5) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_6 <- brick('s1/s1_aoi_6.tif') ## Convert sixth location image from tif to rasterbrick
s2_raster_6 <- UTM_to_LongLat('s1_aoi_6.tif', 's2_aoi_6.tif') ## Convert sixth location image from UTM coordinates to Longitude Latutude coordinates
s1_raster_6 <- aggregate(s1_raster_6, 3,fun=mean) ## Aggregate sixth image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_6 <- aggregate(s2_raster_6, 3,fun=mean) ## Aggregate sixth image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_5, s1_raster_5) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_7 <- brick('s1/s1_aoi_7.tif') ## Convert seventh location image from tif to rasterbrick
s2_raster_7 <- brick('s2/s2_aoi_7.tif') ## Convert seventh location image from tif to rasterbrick
s1_raster_7 <- aggregate(s1_raster_7, 3,fun=mean) ## Aggregate seventh image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_7 <- aggregate(s2_raster_7, 3,fun=mean) ## Aggregate seventh image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_7, s1_raster_7) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_8 <- brick('s1/s1_aoi_8.tif') ## Convert eighth location image from tif to rasterbrick
s2_raster_8 <- brick('s2/s2_aoi_8.tif') ## Convert eighth location image from tif to rasterbrick
s1_raster_8 <- aggregate(s1_raster_8, 3,fun=mean) ## Aggregate eighth image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_8 <- aggregate(s2_raster_8, 3,fun=mean) ## Aggregate eighth image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_8, s1_raster_8) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_9 <- brick('s1/s1_aoi_9.tif') ## Convert ninth location image from tif to rasterbrick
s2_raster_9 <- brick('s2/s2_aoi_9.tif') ## Convert ninth location image from tif to rasterbrick
s1_raster_9 <- aggregate(s1_raster_9, 3,fun=mean) ## Aggregate ninth image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_9 <- aggregate(s2_raster_9, 3,fun=mean) ## Aggregate ninth image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_9, s1_raster_9) ## Check Sentinal 1 and 2 rasterbricks have same features

s1_raster_10 <- brick('s1/s1_aoi_10.tif') ## Convert tenth location image from tif to rasterbrick
s2_raster_10 <- brick('s2/s2_aoi_10.tif') ## Convert tenth location image from tif to rasterbrick
s1_raster_10 <- aggregate(s1_raster_10, 3,fun=mean) ## Aggregate tenth image Sentinal 1 rasterbrick to lower resolution and reduce cells
s2_raster_10 <- aggregate(s2_raster_10, 3,fun=mean) ## Aggregate tenth image Sentinal 2 rasterbrick to lower resolution and reduce cells
compareRaster(s2_raster_10, s1_raster_10) ## Check Sentinal 1 and 2 rasterbricks have same features

## CORRELATION ANALYSIS

img1_corr <- correl(s1_raster_1, s2_raster_1) ## Get correlations between red and blue bands for image 1
img2_corr <- correl(s1_raster_2, s2_raster_2) ## Get correlations between red and blue bands for image 2
img3_corr <- correl(s1_raster_3, s2_raster_3) ## Get correlations between red and blue bands for image 3
img4_corr <- correl(s1_raster_4, s2_raster_4) ## Get correlations between red and blue bands for image 4
img5_corr <- correl(s1_raster_5, s2_raster_5) ## Get correlations between red and blue bands for image 5
img6_corr <- correl(s1_raster_6, s2_raster_6) ## Get correlations between red and blue bands for image 6
img7_corr <- correl(s1_raster_7, s2_raster_7) ## Get correlations between red and blue bands for image 7
img8_corr <- correl(s1_raster_8, s2_raster_8) ## Get correlations between red and blue bands for image 8
img9_corr <- correl(s1_raster_9, s2_raster_9) ## Get correlations between red and blue bands for image 9
img10_corr <- correl(s1_raster_10, s2_raster_10) ## Get correlations between red and blue bands for image 10

dummy <- as.data.frame(rbind(img1_corr, img2_corr, img3_corr, img4_corr, img5_corr, img6_corr, img7_corr, img8_corr, img9_corr, img10_corr)) ## create dataframe of all correlations

corr_df <- data.frame(dummy, "Image" = c("Spain/Portugal", "Spain/Portugal", "Spain/Portugal", "Spain/Portugal", "Ireland", "Ireland", "California", "Argentina", "Austrailia","Madagascar")) ## Add location of image to correlatiom dataframe

## Plot scatter plot of corellations for images 
ggplot(corr_df, aes(x=V1, y=V2, color = Image)) +
  geom_point() + 
  geom_text(label=corr_df$Image, size = 3, hjust = 0, vjust = 1.5) + stat_ellipse() + xlab("Red correlation") + ylab("Blue Correlation") + ggtitle("Correlation Coefficients Between Sentinal 1 & 2 bands")

## GENERATE DATAFRAMES TO BUILD MODEL ON

df_1 <- bricks_to_df(s1_raster_1, s2_raster_1) ## Convert first image into a dataframe with VV, VH, NDVI and Veg classification
df_2 <- bricks_to_df(s1_raster_2, s2_raster_2) ## Convert second image into a dataframe with VV, VH, NDVI and Veg classification
df_3 <- bricks_to_df(s1_raster_3, s2_raster_3) ## Convert thrid image into a dataframe with VV, VH, NDVI and Veg classification
df_4 <- bricks_to_df(s1_raster_4, s2_raster_4) ## Convert fourth image into a dataframe with VV, VH, NDVI and Veg classification
df_5 <- bricks_to_df(s1_raster_5, s2_raster_5) ## Convert fifth image into a dataframe with VV, VH, NDVI and Veg classification
df_6 <- bricks_to_df(s1_raster_6, s2_raster_6) ## Convert sixth image into a dataframe with VV, VH, NDVI and Veg classification
df_7 <- bricks_to_df(s1_raster_7, s2_raster_7) ## Convert seventh image into a dataframe with VV, VH, NDVI and Veg classification
df_8 <- bricks_to_df(s1_raster_8, s2_raster_8) ## Convert eighth image into a dataframe with VV, VH, NDVI and Veg classification
df_9 <- bricks_to_df(s1_raster_9, s2_raster_9) ## Convert ninth image into a dataframe with VV, VH, NDVI and Veg classification
df_10 <- bricks_to_df(s1_raster_10, s2_raster_10) ## Convert tenth image into a dataframe with VV, VH, NDVI and Veg classification

## MISSINGNESS ANALYSED

image1 <- t(data.frame("image_1" = sapply(df_1, function(x) sum(is.na(x))))) ## Generate dataframe for location 1 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image2 <- t(data.frame("image_2" = sapply(df_2, function(x) sum(is.na(x))))) ## Generate dataframe for location 2 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image3 <- t(data.frame("image_3" = sapply(df_3, function(x) sum(is.na(x))))) ## Generate dataframe for location 3 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image4 <- t(data.frame("image_4" = sapply(df_4, function(x) sum(is.na(x))))) ## Generate dataframe for location 4 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image5 <- t(data.frame("image_5" = sapply(df_5, function(x) sum(is.na(x))))) ## Generate dataframe for location 5 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image6 <- t(data.frame("image_6" = sapply(df_6, function(x) sum(is.na(x))))) ## Generate dataframe for location 6 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image7 <- t(data.frame("image_7" = sapply(df_7, function(x) sum(is.na(x))))) ## Generate dataframe for location 7 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image8 <- t(data.frame("image_8" = sapply(df_8, function(x) sum(is.na(x))))) ## Generate dataframe for location 8 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image9 <- t(data.frame("image_9" = sapply(df_9, function(x) sum(is.na(x))))) ## Generate dataframe for location 9 of number of missing pixels for each column, NDVI, VV, VH and Veg classification
image10 <- t(data.frame("image_10" = sapply(df_10, function(x) sum(is.na(x))))) ## Generate dataframe for location 10 of number of missing pixels for each column, NDVI, VV, VH and Veg classification

missing_analysis <- data.frame(rbind(image1, image2, image3, image4, image5, image6, image7, image8, image9, image10), "image" = c(1, 2,3, 4, 5, 6, 7, 8, 9, 10)) ## Create a full dataframe of the missingness in the images 

missing_analysis_2 <- gather(missing_analysis, Band, Total, NDVI:VH) #Create long format of the missing dataframe

## Side by side barplot of the number of missing datapoints in each image location 
plot <- ggplot(missing_analysis_2, aes(as.factor(image), Total, fill=Band))
plot <- plot + geom_bar(stat = "identity", position = 'dodge') + ggtitle("Missing pixel values by band category") + xlab("Sentinal Image") + ylab("Missing Data Points")
plot

## IMPUTE TRAINING DATA SET

df_1 <- impute(df_1) ## Impute missing pixels for dataframe corresponding to the first image
df_6 <- impute(df_6) ## Impute missing pixels for dataframe corresponding to the sixth image
df_9 <- impute(df_9) ## Impute missing pixels for dataframe corresponding to the ninth image
df_10 <- impute(df_10) ## Impute missing pixels for dataframe corresponding to the tenth image

df_1$Region <- "Forest/Grassland" ## Generate new variable for terrain for location 1
df_6$Region <- "Forest/Grassland" ## Generate new variable for terrain for location 6
df_9$Region <- "Tundra/Desert" ## Generate new variable for terrain for location 9
df_10$Region <- "Tundra/Desert" ## Generate new variable for terrain for location 110

## BASELINE MODEL 

df_10_basline <- mean(df_10$NDVI) ## Generate mean NDVI for the basline model

## LOOK INTO THE NUMBER OF CELLS FOR EACH RASTER BRICK

bricks <- list(s2_raster_1, s2_raster_2, s2_raster_3, s2_raster_4, s2_raster_5, s2_raster_6, s2_raster_7, s2_raster_8, s2_raster_9, s2_raster_10) ## Generate list of raster bricks for sentinal 2 to loop through
cells <- rep(0, 10) ## Create storage vector to store number of cells for each image
j <- 1 ## Create counter variable

## Loop through list of raster images and look at number of cells in each raster brick
for (i in bricks) {
  cells[j] <- ncell(i) ## Get number of cells
  j <- j + 1 ## Increment counter 
}

cell_df <- data.frame("Image" = c("Spain/Portugal 1", "Spain/Portugal 2", "Spain/Portugal 3", "Spain/Portugal 4", "Ireland 1", "Ireland 2", "California", "Argentina", "Austrailia","Madagascar"), "Cells" = cells) ## Create full dataframe of number of cells in each image

## Plot barchart of number of cells in each image
p <- ggplot(cell_df, aes(x = Image, y = Cells)) + geom_col() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + ggtitle("Image Cell Count") + ylab("Number of Cells")
p

## NDVI INDEX DISTRIBUTION GRAPHS

NDVI_dfs <- list(df_1$NDVI, df_2$NDVI, df_3$NDVI, df_4$NDVI, df_5$NDVI, df_6$NDVI, df_7$NDVI,df_8$NDVI, df_9$NDVI, df_10$NDVI) ## Create dataframe of just NDVI vectors for each image
country <- list("Spain/Portugal 1", "Spain/Portugal 2", "Spain/Portugal 3", "Spain/Portugal 4", "Ireland 1", "Ireland 2", "California","Argentina", "Austrailia","Madagascar") ## Create list of location of each image to help with plot titles
j <- 1 ## Counter variable 

par(mfrow = c(2, 5)) ## Plot dimensions
## Look through the list of NDVI vectors to create a distribution of NDVIs for each image
for (i in NDVI_dfs) {
  hist(i,
       main = country[j],
       xlab = "NDVI",
       ylab= "Frequency",
       col = "aquamarine3",
       xlim = c(-0.5, 1),
       breaks = 30,
       xaxt = 'n')
  axis(side = 1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))
  j <- j + 1 ## Increment counter variable
}


## SIMULATING CLOUDS FOR TEST DATASETS

par(mfrow = c(2, 2)) ## Define plot areas

## Location 2
b1 <- s2_raster_2[[1]] ## Convert first band of raster brick to raster layer
plot(b1) ## Draw image
e <- drawExtent() ## User to draw first cloud   
e1 <- drawExtent() ## User to draw second cloud  
e2 <- drawExtent() ## User to draw third cloud  

## Location 5
b1 <- s2_raster_5[[1]] ## Convert first band of raster brick to raster layer
plot(b1) ## Draw image
e3 <- drawExtent() ## User to draw first cloud   
e4 <- drawExtent() ## User to draw second cloud  
e5 <- drawExtent() ## User to draw third cloud  

## Location 7
b1 <- s2_raster_7[[1]] ## Convert first band of raster brick to raster layer
plot(b1) ## Draw image
e6 <- drawExtent() ## User to draw first cloud   
e7 <- drawExtent() ## User to draw second cloud  
e8 <- drawExtent() ## User to draw third cloud  

## Location 8
b1 <- s2_raster_8[[1]] ## Convert first band of raster brick to raster layer
plot(b1) ## Draw image
e9 <- drawExtent() ## User to draw first cloud   
e10 <- drawExtent() ## User to draw second cloud  
e11 <- drawExtent() ## User to draw third cloud  

s2_raster_2_new <- Simulate_Clouds(s2_raster_2, e) ## Convert first the extent for image 2 to NAs
s2_raster_2_new <- Simulate_Clouds(s2_raster_2_new, e1) ## Convert second the extent for image 2 to NAs
s2_raster_2_new <- Simulate_Clouds(s2_raster_2_new, e2) ## Convert third the extent for image 2 to NAs

RGB <- stack(list(s2_raster_2_new[[3]], s2_raster_2_new[[2]], s2_raster_2_new[[1]])) ## Create raster stack of RBG image with clouds to ensure clouds are correctly defined
plotRGB(RGB, axes = TRUE, stretch = "lin", main = "Spain RGB colour composite with clouds") ## Draw image


s2_raster_5_new <- Simulate_Clouds(s2_raster_5, e3) ## Convert first the extent for image 5 to NAs
s2_raster_5_new <- Simulate_Clouds(s2_raster_5_new, e4) ## Convert second the extent for image 5 to NAs
s2_raster_5_new <- Simulate_Clouds(s2_raster_5_new, e5) ## Convert third the extent for image 5 to NAs

RGB <- stack(list(s2_raster_5_new[[3]], s2_raster_5_new[[2]], s2_raster_5_new[[1]])) ## Create raster stack of RBG image with clouds to ensure clouds are correctly defined
plotRGB(RGB, axes = TRUE, stretch = "lin", main = "Ireland RGB colour composite with clouds") ## Draw image

s2_raster_7_new <- Simulate_Clouds(s2_raster_7, e6) ## Convert first the extent for image 7 to NAs
s2_raster_7_new <- Simulate_Clouds(s2_raster_7_new, e7) ## Convert second the extent for image 7 to NAs
s2_raster_7_new <- Simulate_Clouds(s2_raster_7_new, e8) ## Convert third the extent for image 7 to NAs
  
RGB <- stack(list(s2_raster_7_new[[3]], s2_raster_7_new[[2]], s2_raster_7_new[[1]])) ## Create raster stack of RBG image with clouds to ensure clouds are correctly defined   
plotRGB(RGB, axes = TRUE, stretch = "lin", main = "California RGB colour composite with clouds") ## Draw image


s2_raster_8_new <- Simulate_Clouds(s2_raster_8, e9) ## Convert first the extent for image 8 to NAs
s2_raster_8_new <- Simulate_Clouds(s2_raster_8_new, e10) ## Convert second the extent for image 8 to NAs
s2_raster_8_new <- Simulate_Clouds(s2_raster_8_new, e11) ## Convert third the extent for image 8 to NAs

RGB <- stack(list(s2_raster_8_new[[3]], s2_raster_8_new[[2]], s2_raster_8_new[[1]])) ## Create raster stack of RBG image with clouds to ensure clouds are correctly defined   
plotRGB(RGB, axes = TRUE, stretch = "lin", main = "Argentina RGB colour composite with clouds") ## Draw image

## CREATE TEST IMAGE DATASETS

df_2 <- bricks_to_df(s1_raster_2, s2_raster_2) ## Convert second image into a dataframe with VV, VH, NDVI and Veg classification
df_5 <- bricks_to_df(s1_raster_5, s2_raster_5) ## Convert third image into a dataframe with VV, VH, NDVI and Veg classification
df_7 <- bricks_to_df(s1_raster_7, s2_raster_7) ## Convert seventh image into a dataframe with VV, VH, NDVI and Veg classification
df_8 <- bricks_to_df(s1_raster_8, s2_raster_8) ## Convert eighth image into a dataframe with VV, VH, NDVI and Veg classification

df_2 <- impute(df_2) ## Impute missing pixels for dataframe corresponding to the second location 
df_5 <- impute(df_5) ## Impute missing pixels for dataframe corresponding to the fifth location 
df_7 <- impute(df_7) ## Impute missing pixels for dataframe corresponding to the seventh location 
df_8 <- impute(df_8) ## Impute missing pixels for dataframe corresponding to the eighth location 


df_2_new <- bricks_to_df(s1_raster_2, s2_raster_2_new) ## Convert second image with synthesized clouds into a dataframe with VV, VH, NDVI and Veg classification
df_5_new <- bricks_to_df(s1_raster_5, s2_raster_5_new) ## Convert fifth image with synthesized clouds into a dataframe with VV, VH, NDVI and Veg classification
df_7_new <- bricks_to_df(s1_raster_7, s2_raster_7_new) ## Convert seventh image with synthesized clouds into a dataframe with VV, VH, NDVI and Veg classification
df_8_new <- bricks_to_df(s1_raster_8, s2_raster_8_new) ## Convert eighth image with synthesized clouds into a dataframe with VV, VH, NDVI and Veg classification


df_2$Region <- "Forest/Grassland" ## Create new variable to describe terrain 
df_5$Region <- "Forest/Grassland" ## Create new variable to describe terrain 
df_7$Region <- "Tundra/Desert" ## Create new variable to describe terrain 
df_8$Region <- "Tundra/Desert" ## Create new variable to describe terrain 
df_2_new$Region <- "Forest/Grassland" ## Create new variable to describe terrain 
df_5_new$Region <- "Forest/Grassland" ## Create new variable to describe terrain 
df_7_new$Region <- "Tundra/Desert" ## Create new variable to describe terrain 
df_8_new$Region <- "Tundra/Desert" ## Create new variable to describe terrain 
df_2_new$AcVeg <- df_2$Veg ## Add actual vegetation classification to the synthesised cloud dataframe for location 2
df_5_new$AcVeg <- df_5$Veg ## Add actual vegetation classification to the synthesised cloud dataframe for location 5
df_7_new$AcVeg <- df_7$Veg ## Add actual vegetation classification to the synthesised cloud dataframe for location 7
df_8_new$AcVeg <- df_8$Veg ## Add actual vegetation classification to the synthesised cloud dataframe for location 8
df_2_new$Image <- 2 ## Create indicator variable designed to easily identify location when combining dataframes
df_5_new$Image <- 5 ## Create indicator variable designed to easily identify location when combining dataframes
df_7_new$Image <- 7 ## Create indicator variable designed to easily identify location when combining dataframes
df_8_new$Image <- 8 ## Create indicator variable designed to easily identify location when combining dataframes
df_2$Image <- 2 ## Create indicator variable designed to easily identify location when combining dataframes
df_5$Image <- 5 ## Create indicator variable designed to easily identify location when combining dataframes
df_7$Image <- 7 ## Create indicator variable designed to easily identify location when combining dataframes
df_8$Image <- 8 ## Create indicator variable designed to easily identify location when combining dataframes

## KMEANS CLUSTERING MODEL

par(mfrow = c(1, 1)) ## define plot area 

cluster_df_2 <- data.frame(df_2$VV, df_2$VV) ## Store VV and VH polarizations in a datafrome for location 2
set.seed(99) ## Set the seed so results are not entirely random each time
kmncluster_2 <- kmeans(cluster_df_2, centers = 2, iter.max = 1000,
                     nstart = 5, algorithm = "Lloyd") ## Run clustering algorithm for location 2

NDVI_2 <- VI(s2_raster_2, 4, 3) # Create raster of NDVI values
knr_2 <- NDVI_2 ## copy the raster
values(knr_2) <- kmncluster_2$cluster ## Replace raster values with classification from k means
kmncluster_2$cluster[kmncluster_2$cluster == 2] <- 0 ## Change cluster object from 0 to 1 

check_df_2 <- data.frame(df_2_new$Veg, df_2$Veg, kmncluster_2$cluster) ## Create dataframe to manipulate and store all synthesized clouds
check_df_2 <- check_df_2[is.na(check_df_2$df_2_new.Veg),] ## Remove all pixels apart from synthesised clouds

Loss <- rep(0,4) ## Create storage vector for 0-1 loss function 
Loss[1] <- mean(abs(check_df_2$df_2.Veg - check_df_2$kmncluster_2.cluster)) ## Store losses for second location  

cluster_df_5 <- data.frame(df_5$VV, df_5$VV) ## Store VV and VH polarizations in a datafrome for location 5
set.seed(99) ## Set the seed so results are not entirely random each time
kmncluster_5 <- kmeans(cluster_df_5, centers = 2, iter.max = 1000,
                       nstart = 5, algorithm = "Lloyd") ## Run clustering algorithm for location 5

NDVI_5 <- VI(s2_raster_5, 4, 3) # Create raster of NDVI values
knr_5 <- NDVI_5 ## copy the raster
values(knr_5) <- kmncluster_5$cluster ## Replace raster values with classification from k means


kmncluster_5$cluster[kmncluster_5$cluster == 2] <- 0 ## Change cluster object from 0 to 1 

check_df_5 <- data.frame(df_5_new$Veg, df_5$Veg, kmncluster_5$cluster)  ## Create dataframe to manipulate and store all synthesized clouds
check_df_5 <- check_df_5[is.na(check_df_5$df_5_new.Veg),] ## Remove all pixels apart from synthesised clouds

Loss[2] <- mean(abs(check_df_5$df_5.Veg - check_df_5$kmncluster_5.cluster)) ## Store losses for fifth location  

 
cluster_df_7 <- data.frame(df_7$VV, df_7$VV) ## Store VV and VH polarizations in a datafrome for location 7
set.seed(99) ## Set the seed so results are not entirely random each time
kmncluster_7 <- kmeans(cluster_df_7, centers = 2, iter.max = 1000,
                       nstart = 5, algorithm = "Lloyd") ## Run clustering algorithm for location 7

NDVI_7 <- VI(s2_raster_7, 4, 3) # Create raster of NDVI values
knr_7 <- NDVI_7 ## copy the raster
values(knr_7) <- kmncluster_7$cluster ## Replace raster values with classification from k means


check_df_7 <- data.frame(df_7_new$Veg, df_7$Veg, kmncluster_7$cluster-1) ## Create dataframe to manipulate and store all synthesized clouds
check_df_7 <- check_df_7[is.na(check_df_7$df_7_new.Veg),] ## Remove all pixels apart from synthesised clouds
 
Loss[3] <- mean(abs(check_df_7$df_7.Veg - check_df_7$kmncluster_7.cluster)) ## Store losses for seventh location  

cluster_df_8 <- data.frame(df_8$VV, df_8$VV) ## Store VV and VH polarizations in a datafrome for location 8
set.seed(99) ## Set the seed so results are not entirely random each time
kmncluster_8 <- kmeans(cluster_df_8, centers = 2, iter.max = 1000,
                       nstart = 5, algorithm = "Lloyd") ## Run clustering algorithm for location 8

NDVI_8 <- VI(s2_raster_8, 4, 3) # Create raster of NDVI values
knr_8 <- NDVI_8 ## copy the raster
values(knr_8) <- kmncluster_8$cluster ## Replace raster values with classification from k means


check_df_8 <- data.frame(df_8_new$Veg, df_8$Veg, kmncluster_8$cluster-1) ## Create dataframe to manipulate and store all synthesized clouds
check_df_8 <- check_df_8[is.na(check_df_8$df_8_new.Veg),] ## Remove all pixels apart from synthesised clouds

Loss[4] <- mean(abs(check_df_8$df_8.Veg - (check_df_8$kmncluster_8.cluster-1))) ## Store losses for eighth location  


## PLOT K MEANS TEST IMAGES 

par(mfrow = c(1, 2)) ## Set plot dimensions 

## Location 2
veg_2 <- NDVI_2 ## Create dummy raster to store actual vegetation classifications 
values(veg_2) <- df_2$Veg ## Replace values to actual vegetation classification
plot(knr_2, main = "Spain Cluster Classification", col = viridis_pal(option = "G")(2)) ## Draw cluster classification image plot
plot(veg_2, main = "Spain Actual Classification", col = viridis_pal(option = "G", direction  = -1)(2)) ## Draw actual classification image plot

## Location 5
veg_5 <- NDVI_5 ## Create dummy raster to store actual vegetation classifications 
values(veg_5) <- df_5$Veg # Replace values to actual vegetation classification
plot(knr_5, main = "Ireland Cluster Classification", col = viridis_pal(option = "G")(2)) ## Draw cluster classification image plot
plot(veg_5, main = "Ireland Actual Classification", col = viridis_pal(option = "G", direction = -1)(2)) ## Draw actual classification image plot

## Location 7
veg_7 <- NDVI_7 ## Create dummy raster to store actual vegetation classifications 
values(veg_7) <- df_7$Veg # Replace values to actual vegetation classification
plot(knr_7, main = "California Cluster Classification", col = viridis_pal(option = "G")(2)) ## Draw cluster classification image plot
plot(veg_7, main = "California Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot

## Location 8
veg_8 <- NDVI_8 ## Create dummy raster to store actual vegetation classifications 
values(veg_8) <- df_8$Veg # Replace values to actual vegetation classification
plot(knr_8, main = "Argentina Cluster Classification", col = viridis_pal(option = "G")(2)) ## Draw cluster classification image plot
plot(veg_8, main = "Argentina Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot

## LOGISTIC REGRESSION

logistic_df <- rbind(df_10, df_1, df_6, df_9) ## Create dataframe of all training datapoints

model2 <- glm(Veg ~ VV + VH + as.factor(Region), data = logistic_df, family= "binomial") ## Build logistic regression model
summary(model2) ## Sumarise the model coefficients and scoring metrics

probabilities <- predict(model2, type = "response") ## Predict model probabilities for training dataset
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg") ## Specify cutoff of 0.5 for probabilities and create a vector with predicted binary classification 

mydata <- logistic_df[,c("VV", "VH")] ## Create new dataframe with VV and VH polarization inputs
predictors <- colnames(mydata) ## Create a storage vector for names of predictor variables 
mydata <- mydata %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)
mydata <- mydata[sample(nrow(mydata), 1000, replace=FALSE), ] ## Transform data to get the logit of probabilities and group by predictor variables

## Plot assumption of linearity graph 
ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() + 
  facet_wrap(~predictors, scales = "free_y")


plot(model2, which = 4, id.n = 3) ## Plot outliers graph showing cooks distance

## Plot graph showing residuals vs fitted values 
model.data <- augment(model2) %>% 
  mutate(index = 1:n()) 
model.data <- as.data.frame(model.data)
model.data <- model.data[sample(nrow(model.data), 500000, replace=FALSE), ]
ggplot(model.data, aes(index, .std.resid)) + 
  geom_point(aes(color = Veg), alpha = .5) +
  theme_bw()

ROC_df <- data.frame("veg" = logistic_df$Veg, "predicted" = probabilities) ## Create a dataframe with actual and predicted classifications 

ROC <- roc(veg ~ predicted, data = ROC_df) ## Create ROC curve dataset
plot(ROC, main= "ROC curve") ## Plot first model ROC curve

final_predict <- rbind(df_2_new, df_5_new, df_7_new, df_8_new) ## Create dataframe to store synthesised clouds, their actual vegetation index and predicted classification 
final_predict <- final_predict[is.na(final_predict$NDVI),] ## Remove all apart from NAs

new_predict <- predict(model2, final_predict, type = "response") ## Predict test dataset probabilities
final_predict$prob <- new_predict
final_predict$predicted <- new_predict ## Add newly predicted probabilities to the test dataset for comparison
final_predict <- final_predict[,c("AcVeg", 'predicted', "Image")] ## Remove unwanted columns
final_predict$predicted[final_predict$predicted > 0.9] <- 1 ## Classify predicted vegetation based on cutoff probability of 0.9
final_predict$predicted[final_predict$predicted <= 0.9] <- 0 ## Everything less than 0.9 is classified as non-vegetated 
final_predict$Loss <- abs(final_predict$AcVeg - final_predict$predicted) ## Store loss function vector 
final_predict$brierloss <- (final_predict$AcVeg - final_predict$prob)^2


Loss_logistic <- rep(0, 4)  ## Create storage vector for logistic regression 0-1 loss function 
Loss_logistic[1] <- mean(final_predict$Loss[final_predict$Image == 2]) ## Loss for test image 1
Loss_logistic[2] <- mean(final_predict$Loss[final_predict$Image == 5], na.rm = TRUE) ## Loss for test image 2
Loss_logistic[3] <- mean(final_predict$Loss[final_predict$Image == 7]) ## Loss for test image 3
Loss_logistic[4] <- mean(final_predict$Loss[final_predict$Image == 8]) ## Loss for test image 4

Brier_logistic <- rep(0, 4) ## Create storage vector for brier score
Brier_logistic[1] <- mean(final_predict$brierloss[final_predict$Image == 8]) ## Brier score for test image 1
Brier_logistic[2] <- mean(final_predict$brierloss[final_predict$Image == 5]) ## Brier score for test image 2
Brier_logistic[3] <- mean(final_predict$brierloss[final_predict$Image == 7]) ## Brier score for test image 3
Brier_logistic[4] <- mean(final_predict$brierloss[final_predict$Image == 2]) ## Brier score for test image 4

## LOGISTIC REGRESSION PLOTS

final_predict_2 <- rbind(df_2, df_5, df_7, df_8) ## Create new dataframe to use to predict all test images to recreate an image
new_predict_2 <- predict(model2, final_predict_2, type = "response") ## Predict probabilities based on input test dataset
new_predict_2[new_predict_2 > 0.5] <- 1 ## Classify test dataset probabilities based on a cutoff of 0.5
new_predict_2[new_predict_2 <= 0.5] <- 0 ## Anything less than 0.5 we classify as non-vegetated 
images_df <- data.frame("Predicted" = new_predict_2, "Image" = final_predict_2$Image) ## Create new datafram to also store the image number 

veg_2_log <- NDVI_2 ## Create dummy raster to store logistic regression vegetation classifications for image 2
values(veg_2_log) <- images_df[images_df$Image == 2, c("Predicted")] ## Replace values with those of the predicted classification from logistic regression for image 2
veg_5_log <- NDVI_5 ## Create dummy raster to store logistic regression vegetation classifications for image 5
values(veg_5_log) <- images_df[images_df$Image == 5, c("Predicted")] ## Replace values with those of the predicted classification from logistic regression for image 5
veg_7_log <- NDVI_7 ## Create dummy raster to store logistic regression vegetation classificationsfor image 7
values(veg_7_log) <- images_df[images_df$Image == 7, c("Predicted")] ## Replace values with those of the predicted classification from logistic regression for image 7
veg_8_log <- NDVI_8 ## Create dummy raster to store logistic regression vegetation classifications for image 8
values(veg_8_log) <- images_df[images_df$Image == 8, c("Predicted")] ## Replace values with those of the predicted classification from logistic regression for image 8

par(mfrow = c(1, 2)) ## Define plot areas 

plot(veg_2_log, main = "Spain Logistic Regression Classification", col = viridis_pal(option = "G")(2)) ## Draw logistic regression classification image plot for image 2
plot(veg_2, main = "Spain Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 2

plot(veg_5_log, main = "Ireland Logistic Regression Classification", col = viridis_pal(option = "G")(2)) ## Draw logistic regression classification image plot for image 5
plot(veg_5, main = "Ireland Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 5

plot(veg_7_log, main = "California Logistic Regression Classification", col = viridis_pal(option = "G")(2)) ## Draw logistic regression classification image plot for image 7
plot(veg_7, main = "California Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 7

plot(veg_8_log, main = "Argentina Logistic Regression Classification", col = viridis_pal(option = "G")(2)) ## Draw logistic regression classification image plot for image 8
plot(veg_8, main = "Argentina Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 8


## NEURAL NETWORK TEST DATASET

setwd("~/Documents/Edinburgh/Dissertation/Space/") ## Set working directory 

neural_predicted_vals <- read.csv(file = 'pred_vals.csv') ## Read csv file for neural network predicted classification 
neural_predicted_probs <- read.csv(file = 'pred_prob_2.csv') ## Read csv file for neural network predicted probs


final_predict_3 <- rbind(df_2_new, df_5_new, df_7_new, df_8_new) ## Create new test dataset dataframe with VV, VH, NDVI and actual vegeation classifications
final_predict_3$predicted <- neural_predicted_vals$X0 ## Add the neural network classification to the test dataset
final_predict_3$prob <- neural_predicted_probs$X0 ## Add the neural network probs to the test dataset
final_predict_3$brierloss <- (final_predict_3$AcVeg - final_predict_3$prob)^2

veg_2_neural <- NDVI_2 ## Create dummy raster to store neural network vegetation classifications for image 2
values(veg_2_neural) <- final_predict_3[final_predict_3$Image == 2, c("predicted")] ## Replace values with those of the predicted classification from neural network for image 2
veg_5_neural <- NDVI_5 ## Create dummy raster to store neural network vegetation classifications for image 5
values(veg_5_neural) <- final_predict_3[final_predict_3$Image == 5, c("predicted")] ## Replace values with those of the predicted classification from neural network for image 5
veg_7_neural <- NDVI_7 ## Create dummy raster to store neural network vegetation classifications for image 7
values(veg_7_neural) <- final_predict_3[final_predict_3$Image == 7, c("predicted")] ## Replace values with those of the predicted classification from neural network for image 7
veg_8_neural <- NDVI_8 ## Create dummy raster to store neural network vegetation classifications for image 8
values(veg_8_neural) <- final_predict_3[final_predict_3$Image == 8, c("predicted")] ## Replace values with those of the predicted classification from neural network for image 8

par(mfrow = c(1, 2)) ## Define plot areas 

plot(veg_2_neural, main = "Spain Neural Network Classification", col = viridis_pal(option = "G")(2)) ## Draw neural network regression classification image plot for image 2
plot(veg_2, main = "Spain Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 2

plot(veg_5_neural, main = "Ireland Neural Network Classification", col = viridis_pal(option = "G")(2)) ## Draw neural network regression classification image plot for image 5
plot(veg_5, main = "Ireland Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 5

plot(veg_7_neural, main = "California Logistic Regression Classification", col = viridis_pal(option = "G")(2)) ## Draw neural network regression classification image plot for image 7
plot(veg_7, main = "California Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 7

plot(veg_8_neural, main = "Argentina Logistic Regression Classification", col = viridis_pal(option = "G")(2)) ## Draw neural network regression classification image plot for image 8
plot(veg_8, main = "Argentina Actual Classification", col = viridis_pal(option = "G")(2)) ## Draw actual classification image plot for image 8

final_predict_3$Loss <- abs(final_predict_3$AcVeg - final_predict_3$predicted) ## Generate the loss function for predicted images 

Loss_neural <- rep(0, 4) ## Create storage vector for losses of neural network
Loss_neural[1] <- mean(final_predict_3$Loss[final_predict_3$Image == 2]) ## Loss for test image 1
Loss_neural[2] <- mean(final_predict_3$Loss[final_predict_3$Image == 5]) ## Loss for test image 2
Loss_neural[3] <- mean(final_predict_3$Loss[final_predict_3$Image == 7]) ## Loss for test image 3
Loss_neural[4] <- mean(final_predict_3$Loss[final_predict_3$Image == 8]) ## Loss for test image 4

Brier_neural <- rep(0, 4) ## Create storage vector for brier score
Brier_neural[1] <- mean(final_predict_3$brierloss[final_predict_3$Image == 2]) ## Brier score for test image 1
Brier_neural[2] <- mean(final_predict_3$brierloss[final_predict_3$Image == 5]) ## Brier score for test image 2
Brier_neural[3] <- mean(final_predict_3$brierloss[final_predict_3$Image == 7]) ## Brier score for test image 3
Brier_neural[4] <- mean(final_predict_3$brierloss[final_predict_3$Image == 8]) ## Brier score for test image 4


veg_amt_total <- length(final_predict_3[final_predict_3$AcVeg == 1, c("AcVeg")])/length(final_predict_3[, c("AcVeg")]) ## See how will the baseline model would have predicted images by finding the percentage of vegetated areas

error_final <- data.frame("Neural Network" = Loss_neural, "Logistic_Regression"=Loss_logistic, "Clustering"=Loss, "Image" = c("Spain","Ireland","Austrailia","Madagascar"))
Error_final_long <- gather(error_final, Band, Loss, Neural.Network:Clustering) #Create long format of the missing dataframe

## Create plot comparing loss function across the 4 locations for each of the threee  regression technqiues
## Add lines to the model to highlight random and baseline model
ggplot(data=Error_final_long, aes(x=Image, y=Loss, fill=Band)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Loss function for synthesised clouds") + 
  geom_abline(slope=0, intercept=0.5,  col = "lightblue",lty=2) +
  annotate("text", label = "Ramdom Model", x = 1, y = 0.52) + 
  geom_abline(slope=0, intercept=veg_amt_total,  col = "lightblue",lty=2) +
  annotate("text", label = "Mean Baseline", x = 1, y = veg_amt_total + 0.02) 

brier_final <- data.frame("Neural Network" = Brier_neural, "Logistic Regression"= Brier_logistic, "Image" = c("Spain","Ireland","Austrailia","Madagascar"))
brier_final_long <- gather(brier_final, Band, Loss, Neural.Network:Logistic.Regression) #Create long format of the missing dataframe

ggplot(data=brier_final_long, aes(x=Image, y=Loss, fill=Band)) +
  geom_bar(stat="identity", position=position_dodge()) + ggtitle("Brier Loss function for synthesised clouds") 

## CREATE CSV FILES TO EXPORT FOR USE IN PYTHON

dfs <- list(df_1, df_2, df_3, df_4, df_5, df_6, df_7, df_8, df_9, df_10) ## Create a list of dataframes to loop through 
 
setwd("~/Documents/Edinburgh/Dissertation/Space/CSV") ## Set working directory 

j = 1 ## Create counter for file naming 

## Loop through dataframes containig VV, VH, Region, NDVI index and classification and export each to CSV file
for (i in dfs) {
  write.csv(i,paste("df", j, ".csv", sep=""), row.names = FALSE) ## Create csv file
  print(paste("done", j))
  j = j + 1
}

dfs_new <- list(df_2_new, df_5_new, df_7_new, df_8_new) ## Create a list of dataframes with synthesised clouds to loop through 

setwd("~/Documents/Edinburgh/Dissertation/Space/CSV") ## Set working directory

j = 1 ## Create counter for file naming 

## Loop through dataframes containig VV, VH, Region, NDVI index and classification for synthesised cloud images and export each to CSV file
for (i in dfs_new) {
  write.csv(i,paste("df_new", j, ".csv", sep=""), row.names = FALSE)
  print(paste("done", j))
  j = j + 1
}

dfs_logistic <- list(logistic_df, final_predict_2, final_predict) ## Create a list of dataframes to loop through, these are bespoke and maybe used in neural network training 

setwd("~/Documents/Edinburgh/Dissertation/Space/CSV") ## Set working directory 

j = 1 ## Create counter for file naming 

## Loop through dataframes for neural network training and prediction to write to CSV file
for (i in dfs_logistic) {
  write.csv(i,paste("df_logistic", j, ".csv", sep=""), row.names = FALSE) ## Write a csv file
  print(paste("done", j))
  j = j + 1
}

## Plot of VH Polarization Vs VV Polarization Split by Vegetation Classification
p <- ggplot(logistic_df, aes(VV, VH, color = Veg))+
  geom_point() + ggtitle("VH Polarization Vs VV Polarization Split by Vegetation Classification")


## SVM

SVM_10 <- as.data.frame(df_10, xy = T) ## Create a datafram of just 1 image to test SVM algorithm on
SVM_10 <- SVM_10[,c("VV", "VH", "Veg")] ## Reduce size of dataframe to just three variables 

fit <- svm(factor(Veg) ~ ., data = SVM_10, scale = FALSE, kernel = "radial") ## Fit a SVM

## Attempt at converting the SVM results, however, processing time is too long

#xgrid = expand.grid(x = SVM_10_samp$x, y = SVM_10_samp$y)
#ygrid = predict(fit, xgrid)

#(xgrid, col = as.numeric(ygrid), pch = 20, cex = .5)
#points(SVM_10_samp[,1:2], col = SVM_10_samp$Veg, pch = 19)

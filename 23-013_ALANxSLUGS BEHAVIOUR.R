
# REMOVE ALL OBJECT FROM R-WORKING SPACE
rm(list=ls()) 

# EMPTY THE GARBAGE
gc()

# MAKE A NEW LOCAL DIRECTORY ON YOUR MACHINE IN WHICH R WILL LOOK FOR THIS ANALYSIS AND WILL EXPORT GRAPHS & CO.
DIRECTORY = "C:/Users/vgrog/OneDrive/00_AGROSCOPE - ALAN PROJECT/06_ANALYSIS/2023/"
HOME = paste0(DIRECTORY, "01_ALAN x HERBIVORY/OUTPUT")
rm(DIRECTORY)

# CHECK IF THAT DIRECTORY EXIST
file.exists(HOME)

# MAKE HOME OUR WORK DIRECTORY, THE DEFAULT FOLDER IN WHICH R WILL LOOK
setwd(HOME)

##############################################################################################################

# I NEVER REMEMBER WHICH PACKAGES I USE. I ALWAYS LOAD THEM ALL....

library("stringr")
# library("stringi")
# library("tidyverse")
library("ggplot2")
library("ggpubr")
# library("data.table")
library("dplyr")
# library("plyr")           # LOAD MULTIPLE CSV TABLES
# library("readr")
library("lubridate")
library("wesanderson")  # Color palettes 1
library("RColorBrewer") # Color palettes 2
# library("maptools")     # Import .shp and work with it
# library("spatstat")     # Import .shp and work with it
# library("rgdal")        # Deal with .shp
# library("rgeos")        # Deal with .shp
# library("installr")
# library("terra")
# library("raster")
library("gridExtra")
library("grid")
library("gridExtra")    # MERGE DIFFERENT GGPLOT
library("ggcorrplot")   # PEARSON CORRELATION PLOT
library("reshape2")
# library("rmarkdown")
# library("ggfortify")    # GLM OUTPUT WITH GGPLOT
library("pscl")         # ZERO INFLATED MODELS
# library("lmtest")       # LOGLIKELIHOOD TEST
# library("MASS")         # GLM.NB
library("lme4")         # glmer
# library("nlme")         # lme
# library("glmmADMB")
library("GLMMadaptive")
# library("FactoMineR")   # Compute MCA
# library("factoextra")   # Visualize MCA
# library("piecewiseSEM")
# library("lavaan")
# library("semEff")
library("ggpattern")
library("rlist")
library("Hmisc")
library("piecewiseSEM")
library("multcompView")
library("knitr")
library("arm")
library("stargazer")
library("png")
library("patchwork")
library("gtsummary")
library("broom.mixed")
library("car")
# library("ggstats")
library("gt")
library("webshot2")
library("sjPlot")

library("imager")
library("imagefx")
library("jpeg")
library("rasterImage")

library("EBImage")
library("exiftoolr")

##############################################################################################################

# Things to change when loading a new table
#------------------------------------------

# 1 - Image directory (01.1)
# 2 - Adjust Grid over first picture (01.1)
# 3 - Change directory for output picture (01.3)
# 4 - Change file name of the output picture (01.3)
# 5 - Change Data in the SLUG BEHAV for loop (01.3)
        # Lit vs Dark
        # N slugs
        # Week number
# 6 - Change output name for SLUG BEHAV.csv for loop file (01.4)

# 01.1 ADJUST THE GRID
######################

# IMPORT JPEG
#------------

# IDENTIFY WHERE THE IMAGES ARE LOCATED
img.dir <- 'C:/Users/vgrog/OneDrive/00_AGROSCOPE - ALAN PROJECT/05_DATA/2023/2023_H_SLUGS/2024_BEHAVIOUR/ILLUMINATED/RAW/WEEK04/'

# LIST ALL THE FILES WITHIN DIRECTORY
img.files <- list.files(img.dir)

# CHOSE THE FIRST FILE FROM THE DIRECTORY
cur.file = img.files[1]

# Read in the JPEG image
cur.img <- readJPEG(paste0(img.dir, cur.file))

# OVERLAY THE THREE DIMENSIONS
cur.img.sum <- cur.img[,,1] + cur.img[,,2] + cur.img[,,3]
cur.img.bw <- cur.img.sum / max(cur.img.sum)

# DOWNSIZE RESOLUTION
cur.img.bw = resize(cur.img.bw, w = 200, filter = "none")

#-------------------------------------------------------------------------------------------------------------

# GRID THE JPEG
#--------------

# Rotate the image so the origin is located in the bottom right
img.rot <- t(apply(cur.img.bw, 2 ,rev))
rm(cur.img.sum, cur.img.bw)

# Define the domain of this image
xdim <- nrow(img.rot)
ydim <- ncol(img.rot)

# Choose the number of grids in the x and y direction
num.xgrid = 50
num.ygrid = 40

# Set limits
lims.xs.min <- 16
lims.xs.max <- 36
lims.ys.min <- 10
lims.ys.max <- 40

grid.limits <- c(lims.xs.min, lims.xs.max, lims.ys.min, lims.ys.max)

# Find the grid xs and ys
grid.xs <- seq(1, xdim, length.out = num.xgrid + 1)[lims.xs.min:lims.xs.max]
grid.ys <- seq(1, ydim, length.out = num.ygrid + 1)[lims.ys.min:lims.ys.max]

## Plot the grid locations

image(1:xdim,
      1:ydim,
      img.rot,
      col = gray.colors(100),
      asp = 1)

abline(v = grid.xs, col = 'red', lwd = 0.5)
abline(h = grid.ys, col = 'red', lwd = 0.5)

rm(img.rot, cur.img, grid.xs, grid.ys, num.xgrid, num.ygrid, xdim, ydim, cur.file)
rm(lims.xs.max, lims.xs.min, lims.ys.max, lims.ys.min)

##############################################################################################################

# 01.2 TRACK MOTIONS IN GRID
############################

# Initilize a list to hold all the shifts in each grid over every image sequence pair
shift.list <- list()

# Create an empty dataframe to collect data
SLUG.BEHAV <- data.frame(
    
    WEEK = NA,
    LIGHT.TREA = NA,
    DATE = NA,
    TIME = NA,
    N_CELLS = NA,
    N_SLUGS = NA,
    BEHAV.INDEX = NA,
    x = NA,
    y = NA)

# Loop over all the images

ii <- 1 # first pair of images

while(ii < length(img.files)) { 
    
    # Empty data frame
    i.SLUG.BEHAV <- SLUG.BEHAV[1,]
  
    print(ii)

    # Read in the current 2 images
    img1.org <- readJPEG(paste0(img.dir, img.files[ii]))
    img2.org <- readJPEG(paste0(img.dir, img.files[ii + 1]))
    
    # install_exiftool()
    dat <- exif_read(paste0(img.dir, img.files[ii + 1]))
    DateTime  <- dat$DateTimeOriginal

    # Overlay RGB layers
    img1.org <- img1.org[,,1] + img1.org[,,2] + img1.org[,,3]
    img1.org <- img1.org / max(img1.org)
    
    img2.org <- img2.org[,,1] + img2.org[,,2] + img2.org[,,3]
    img2.org <- img2.org / max(img2.org)

    # Rotate them 90 degrees clockwise
    img1 <- t(apply(img1.org, 2, rev))
    img2 <- t(apply(img2.org, 2, rev))

    # Set image dimensions
    xdim <- nrow(img1)
    ydim <- ncol(img1)

    # Discritize the raw grid
    num.xgrid = 50
    num.ygrid = 40

    # Find the grid xs and ys
    grid.xs <- seq(1, xdim, length.out = num.xgrid + 1)[grid.limits[1]:grid.limits[2]]
    grid.ys <- seq(1, ydim, length.out = num.ygrid + 1)[grid.limits[3]:grid.limits[4]]

    # Create a matrix holding the x and y shifts and the correlation value for reach grid cell
    shift.mat = matrix(NA, nrow = (length(grid.xs) - 1) * (length(grid.ys) - 1), ncol = 5)

    # For loop over grids to find the movement between frames in each grid (one loop = one cell)
    k <- 1
    i <- 1
    
    while(i < length(grid.ys)) {

        j <- 1
        
        while( j < length(grid.xs)) {

            # Define the x and y indices for the current grid region
            x.inds <- grid.xs[j]:grid.xs[j + 1]
            y.inds <- grid.ys[i]:grid.ys[i + 1]

            # Find the cur grid for each image [WHY MINUS MEAN?]
            cur.grid1 <- img1[x.inds, y.inds] - mean(img1[x.inds, y.inds])
            cur.grid2 <- img2[x.inds, y.inds] - mean(img2[x.inds, y.inds])
            
            # Sometimes it is helpful to window the gridded regions
            gaus.xdim = nrow(cur.grid1)
            gaus.ydim = ncol(cur.grid1)
            gaus.sig <- 10

            wind.gaus = build.gaus(gaus.xdim, gaus.ydim, gaus.sig)
            cur.grid1 = cur.grid1 * wind.gaus
            cur.grid2 = cur.grid2 * wind.gaus
            
            # Track the movement between grids using...
              
              # [1] The CROSS correlation function
              cur.corr <- xcorr3d(cur.grid1, cur.grid2)

            # OR track the movement using...
            
              # [2] The PHASE correlation function 
              # cur.corr <- pcorr3d(cur.grid1,cur.grid2)

            # Save the current current correlation components
            cur.xys = c(mean(x.inds), mean(y.inds) )
            cur.sc = c(cur.corr$max.shifts, cur.corr$max.corr)
            shift.mat[k,] <- c(cur.xys, cur.sc)

            k = k + 1
            j = j + 1
        }
        
        i = i + 1
    }
    
    # Save the shift matrix in the shift list
    shift.list[[ii]] <- shift.mat

    ii = ii + 1
    
# }

rm(cur.corr, cur.grid1, cur.grid2, shift.mat)
rm(cur.sc, cur.xys, dat, wind.gaus)
rm(img1.org, img2.org, x.inds, y.inds)
rm(gaus.sig, gaus.xdim, gaus.ydim)

##############################################################################################################

# 01.3 PLOT GRID MOTION
#######################

num.xgrid <- length(grid.xs) - 1
num.ygrid <- length(grid.ys) - 1

# Convert each matix in the shift list into a data frame
shift.list <- lapply(shift.list, as.data.frame)

# Create a matrix of all the shift values in each list element
shift.xs.unlist = unlist(lapply(shift.list, "[", 3))
all.shift.xs <- matrix(shift.xs.unlist, nrow = num.xgrid * num.ygrid)

shift.ys.unlist = unlist(lapply(shift.list, "[", 4))
all.shift.ys = matrix(shift.ys.unlist, nrow = num.xgrid * num.ygrid)

rm(shift.xs.unlist, shift.ys.unlist)

# Take the row averages for all the shifts in the x and y direction
avg.shift.xs <- all.shift.xs[,(ii-1)]
avg.shift.ys <- all.shift.ys[,(ii-1)]

rm(all.shift.xs, all.shift.ys)

# Take the x and y grid locations from an example shift element
avg.xs = shift.list[[(ii-1)]][,1]
avg.ys = shift.list[[(ii-1)]][,2]

# Find the shifts that are greater than some minimum tolerance (this is to avoid warnings when plotting arrows later on)
t <- 10
# pos.shift.inds <- which(abs(avg.shift.xs) > t | abs(avg.shift.ys) > t)
pos.shift.inds <- which(sqrt(abs(avg.shift.xs)^2 + abs(avg.shift.ys)^2) > t)

# Limit the average shift xs, ys to the positive indices
pos.shift.xs <- avg.shift.xs[pos.shift.inds]
pos.shift.ys <- avg.shift.ys[pos.shift.inds]
pos.xs <- avg.xs[pos.shift.inds]
pos.ys <- avg.ys[pos.shift.inds]

# rm(avg.shift.xs, avg.xs, avg.shift.ys, avg.ys)

  # Plot

  OUTPUT.PIC <- "C:/Users/vgrog/OneDrive/00_AGROSCOPE - ALAN PROJECT/05_DATA/2023/2023_H_SLUGS/2024_BEHAVIOUR/ILLUMINATED/EDITED/WEEK04/"

  jpeg(file = paste0(OUTPUT.PIC, "LIT04.", ii, ".jpeg"), res = 100)

    mar.org = par()$mar
    
    # Set the margins for this particular plot
    par(mar = c(0, 0, 0, 0))
    
    # Plot an example image
    image(
      1:xdim,
      1:ydim,
      img2,
      col = gray.colors(20),
      useRaster = TRUE,
      asp = 1,
      axes = FALSE,
      xlab = '',
      ylab = '')
    
    # Return the margins to the users original values
    par(mar = mar.org)
    
    # Define how to scale the arrows (for easier visualization)
    arrow.scale = 5
    
    # Define the arrow locations
    x0 = pos.xs
    y0 = pos.ys
    x1 = x0 + (pos.shift.xs * arrow.scale)
    y1 = y0 + (pos.shift.ys * arrow.scale)
    
    # Plot the motion indicating the average shift between frames
    arrows(x0 = x0, y0 = y0, x1 = x1, y1 = y1, length = 0.05, lwd = 1, angle = 25, col = 'red')
    
  dev.off()
  
  # Calculte shift
  i.shift.list <- shift.list[[ii - 1]]
  i.shift.list.pos <- i.shift.list[(sqrt(abs(i.shift.list$V3)^2 + abs(i.shift.list$V4)^2)) > t,]
  i.shift <- sqrt(abs(i.shift.list.pos$V3)^2 + abs(i.shift.list.pos$V4)^2)
  
  ifelse(
    
    # ARGUMENT
    length(i.shift) > 0,
    
    # IF YES - FILL SLUB BEHAV
    {
    
    # Repete n times number of lines 
    i.SLUG.BEHAV <- i.SLUG.BEHAV[rep(seq_len(nrow(i.SLUG.BEHAV)), length(i.shift)), ]    
      
    i.SLUG.BEHAV$WEEK <- "04"
    i.SLUG.BEHAV$LIGHT.TREA <- "LIT"
    i.SLUG.BEHAV$DATE <- substr(DateTime, 1, 10)
    i.SLUG.BEHAV$TIME <- substr(DateTime, 12, 16)
    i.SLUG.BEHAV$N_CELLS <- num.xgrid * num.ygrid
    i.SLUG.BEHAV$N_SLUGS <- 11
    i.SLUG.BEHAV$BEHAV.INDEX <- i.shift
    i.SLUG.BEHAV$x <- pos.xs
    i.SLUG.BEHAV$y <- pos.ys
    SLUG.BEHAV <- rbind(SLUG.BEHAV, i.SLUG.BEHAV)},
  
    # IF NO
    {
    
    i.SLUG.BEHAV$WEEK <- "04"
    i.SLUG.BEHAV$LIGHT.TREA <- "LIT"
    i.SLUG.BEHAV$DATE <- substr(DateTime, 1, 10)
    i.SLUG.BEHAV$TIME <- substr(DateTime, 12, 16)
    i.SLUG.BEHAV$N_CELLS <- num.xgrid * num.ygrid
    i.SLUG.BEHAV$N_SLUGS <- 11
    i.SLUG.BEHAV$BEHAV.INDEX <- NA
    i.SLUG.BEHAV$x <- NA
    i.SLUG.BEHAV$y <- NA
    SLUG.BEHAV <- rbind(SLUG.BEHAV, i.SLUG.BEHAV)})
  
rm(img1, img2)

}

# REMOVE FIRST LINE FULL OF NA's
SLUG.BEHAV <- SLUG.BEHAV[-1,]

rm(arrow.scale, grid.limits, grid.xs, grid.ys, i, ii, k, j, mar.org, num.xgrid, num.ygrid)
rm(OUTPUT.PIC, pos.shift.inds, pos.shift.xs, pos.shift.ys, pos.xs, pos.ys, x0, x1, y0, y1, xdim, ydim)
rm(i.shift.list, i.shift.list.pos, i.SLUG.BEHAV, shift.list)
rm(avg.shift.xs, avg.shift.ys, avg.xs, avg.ys, i.shift, t, DateTime)

##############################################################################################################

# 01.4 SAVE SLUG BEHAV
######################

# SAVE
DATADIR <- "C:/Users/vgrog/OneDrive/00_AGROSCOPE - ALAN PROJECT/05_DATA/2023/2023_H_SLUGS/2024_BEHAVIOUR/SHIFT TABLES/"
write.csv(SLUG.BEHAV, 
          paste0(
            DATADIR, 
            "LIT.",
            "WEEK04.",
            "csv"))

rm(DATADIR, SLUG.BEHAV, img.dir, img.files)

# END
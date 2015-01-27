# code to generate figures for paper

paperDir <- getwd()

setwd('~/research/jmac/composition')
source("config")

# TMP
TMP <- TRUE
if(TMP){
productVersion <- 0.2
outputDir="/var/tmp/paciorek/paleon/comp_0.2/output"
dataDir = '/var/tmp/paciorek/paleon/comp_0.2/data'
}
source(file.path(codeDir, "plot.R"))
source(file.path(codeDir, "set_domain.R"))

require(ncdf4)
require(ggplot2)
require(maptools)
require(rgdal)
require(raster)

usShp <- readShapeLines(file.path(dataDir, 'us_alb.shp'), proj4string=CRS('+init=epsg:3175'))
usShp@data$id <- rownames(usShp@data)
usFortified <- fortify(usShp, region='id')

region = t(as.matrix(raster(file.path(dataDir, 'paleonDomain.tif'))))
water = t(as.matrix(raster(file.path(dataDir, 'water.tif'))))
# t()  manipulates matrix so plots correctly W-E and N-S in R

# region[region %in% c(2,3,5,6,11,12)] <- NA
water[water == 100] <- NA
mask = is.na(region)
maskWater = is.na(water)


# raw stuff would need more work to show raw data
# western data/results

load(file.path(dataDir, paste0('data_western_', runID, '.Rda')))
coordsWest <- coord

raw_west <- matrix(0, nrow = nCells, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$cell[data$taxon == p])
  raw_west[as.numeric(names(tbl)) , p] <- tbl
}
total <- rowSums(raw_west)
raw_west[total == 0] = NA
total[total == 0] <- 1
raw_west <- raw_west / total
dimnames(raw_west)[[2]] <- gsub("/", "ZZZ", taxa$taxonName)  # otherwise both / and " " become "." so can't distinguish when I substitute back in for "."

# for distinguishing no trees from no data
west_presence <- read.csv(file.path(dataDir, 'western.csv'))
west_presence <- west_presence[west_presence$x <= easternLimitOfWesternDomain, ]
west_presence$data <- as.numeric(rowSums(west_presence[ , c(3:ncol(west_presence))]) > 0)
# west_presence$data = as.numeric(rowSums(west_presence[ , c(3:ncol(west_presence))]) == west_presence$no.tree)
# get in same order as data for fitting (from far N, W to E in rows)
tmp <- matrix(west_presence$data, length(westernDomainY), length(westernDomainX))
west_presence <- c(t(tmp[rev(seq_along(westernDomainY)), ]))

finalNcdfName <- paste0('PLScomposition_western_', productVersion, '-release.nc')

ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))
test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nSamples <- dim(test)[3]

if(nCells != prod(dim(test)[1:2]))
  stop("nCells does not match first dimension of netCDF file.")

preds <- array(0, c(nCells, nTaxa, nSamples))
#dimnames(preds)[[2]] <- taxa$taxonName
for(p in 1:nTaxa) 
  try(preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1)))


attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

pmWest <- apply(preds, c(1, 2), 'mean')
psdWest <- apply(preds, c(1, 2), 'sd')



# eastern data/results

load(file.path(dataDir, paste0('data_eastern_', runID, '.Rda')))
load(file.path(dataDir, paste0('intersection_eastern_', productVersion, '.Rda')))

nTowns <- dim(townCellOverlap)[1]
easternDataDir <- 'eastern'
ohioDataDir <- 'ohio'
eastern_townships <- readOGR(file.path(dataDir, easternDataDir), paste0(easternVersionID, 'polygons_v', easternVersion))
ohio_townships <- readOGR(file.path(dataDir, ohioDataDir), paste0('OH', ohioVersionID, 'polygons_v', ohioVersion))
ohio_townships <- spTransform(ohio_townships, CRSobj=CRS('+init=epsg:3175'))  # transform to Albers

raw_east <- matrix(0, nrow = nTowns, ncol = nTaxa)
for(p in 1:nTaxa) {
  tbl <- table(data$town[data$taxon == p])
  raw_east[as.numeric(names(tbl)) , p] <- tbl
}
raw_east <- raw_east / rowSums(raw_east)
attributes(raw_east)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

#for(p in seq_len(nTaxa))
#  raw_east[ , p] <- cut(raw_east[ , p], propBreaks, include.lowest = TRUE, labels = FALSE)
  

finalNcdfName <- paste0('PLScomposition_eastern_', productVersion, '-release.nc')

ncdfPtr <- nc_open(file.path(outputDir, finalNcdfName))
test <- ncvar_get(ncdfPtr, "Oak", c(1, 1, 1), c(-1, -1, -1))

nSamples <- dim(test)[3]

nCells <- m1*m2

if(nCells != prod(dim(test)[1:2]))
  stop("nCells does not match first dimension of netCDF file.")


preds <- array(0, c(nCells, nTaxa, nSamples))
for(p in 1:nTaxa) 
  try(preds[ , p, ] <- ncvar_get(ncdfPtr, taxa$taxonName[p], c(1, 1, 1), c(-1, -1, -1)))


attributes(preds)$dimnames[[2]] <- gsub("/", "ZZZ", taxa$taxonName) 

pmEast <- apply(preds, c(1, 2), 'mean')
psdEast <- apply(preds, c(1, 2), 'sd')

## combine across grids


pmFull <- matrix(NA, c(xRes*yRes), ncol(pmEast))
dimnames(pmFull)[[2]] <- dimnames(pmEast)[[2]]
fullTmp <- psdFull <- rawFull <- pmFull

ids <- matrix(1:(xRes*yRes), nrow = xRes, ncol = yRes)
eastSubset <- c(ids[easternDomainX, easternDomainY])
westSubset <- c(ids[westernDomainX, westernDomainY])

eastOnly <- dimnames(pmEast)[[2]]
eastOnly <- eastOnly[!(eastOnly %in% dimnames(pmWest)[[2]])]

pmFull[eastSubset, ] <- pmEast
fullTmp[westSubset, dimnames(pmWest)[[2]]] <- pmWest
fullTmp[!(region %in% c(6,12,11, 5, 2, 3)), ] <- NA
pmFull[!is.na(fullTmp)] <- fullTmp[!is.na(fullTmp)]
pmFull[region %in% c(5, 12), eastOnly] <- NA

psdFull[eastSubset, ] <- psdEast
fullTmp[westSubset, dimnames(psdWest)[[2]]] <- psdWest
fullTmp[!(region %in% c(6,12,11, 5, 2, 3)), ] <- NA
psdFull[!is.na(fullTmp)] <- fullTmp[!is.na(fullTmp)]
psdFull[region %in% c(5, 12), eastOnly] <- NA

pmFull[mask, ] <- NA
psdFull[mask, ] <- NA
rawFull[mask, ] <- NA

coordFull <- expand.grid(X = xGrid, Y = rev(yGrid))

## get raw data into polygons

raw_west_poly <- rasterToPolygons(raster(ncols = length(westernDomainX), nrows = length(westernDomainY),
                                       xmn = xRange[1], xmx = easternLimitOfWesternDomain + gridSpacing/2,
                         ymn = yRange[1], ymx = yRange[2]),
                      fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)
west_fort <- fortify(raw_west_poly)
names(west_fort)[1:2] <- c('X', 'Y')
west_presence_fort <- west_fort
if(TMP) {
  tmp <- dimnames(raw_west)[[2]]
  NAs <- rep(as.numeric(NA), nrow(raw_west))
  raw_west <- cbind(raw_west[ , 1], NAs, raw_west[ , c(2:7)], NAs, raw_west[ , c(8:21)])
  dimnames(raw_west)[[2]] <- c(tmp[1], 'Atlantic white cedar', tmp[2:7], 'Chestnut', tmp[8:21])
}
west_fort$id <- as.numeric(west_fort$id)
west_fort <- cbind(west_fort, raw_west[west_fort$id, ])

west_presence_fort$id <- as.numeric(west_presence_fort$id)
west_presence_fort <- cbind(west_presence_fort, west_presence[west_presence_fort$id])
names(west_presence_fort)[ncol(west_presence_fort)] <- "data"

northeast_fort <- fortify(eastern_townships)
northeast_fort$id <- as.numeric(northeast_fort$id)
ohio_fort <- fortify(ohio_townships)
ohio_fort$id <- as.numeric(ohio_fort$id)

northeast_idMap <- data.frame(id = sort(eastern_townships@data$ID), town = seq_along(eastern_townships))
ohio_idMap <- data.frame(id = sort(ohio_townships@data$ID), town = length(eastern_townships) + seq_along(ohio_townships))

northeast_fort <- merge(northeast_fort, northeast_idMap, by.x = 'id', by.y = 'id', all.x = TRUE, all.y = FALSE)
northeast_fort <- northeast_fort[order(as.numeric(northeast_fort$id), northeast_fort$order), ]

ohio_fort <- merge(ohio_fort, ohio_idMap, by.x = 'id', by.y = 'id', all.x = TRUE, all.y = FALSE)
ohio_fort <- ohio_fort[order(as.numeric(ohio_fort$id), ohio_fort$order), ]

east_fort <- rbind(northeast_fort, ohio_fort)
east_fort <- cbind(east_fort, raw_east[east_fort$town,])
east_fort$id <- east_fort$town + max(west_fort$id)
names(east_fort)[2:3] <- c('X', 'Y')
east_fort$town <- NULL
east_fort <- east_fort[ , names(west_fort)]

east_presence_fort <- east_fort
east_presence_fort$data <- as.numeric(rowSums(east_presence_fort[ , 8:ncol(east_presence_fort)], na.rm = TRUE) > 0)

fort_presence <- rbind(west_presence_fort, east_presence_fort[ , c(1:7, ncol(east_presence_fort))])

fort <- rbind(west_fort, east_fort)
#ggplot(fort, aes(long, lat, group = id)) + geom_polygon(aes(fill = Oak))

names(fort) <- gsub("ZZZ", "/", names(fort))
taxon_dat_long <- melt(fort, c('X', 'Y', 'order', 'hole', 'piece', 'group', 'id' ))

  
#########################################################################################
## plot the data, plots in large grid
#########################################################################################

focalTaxa <- c('Beech', 'Cherry', 'Chestnut', 'Elm', 'Hemlock', 'Oak', 'Pine')
focalTaxa <- focalTaxa[-2]

propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

figHgt = 16
figWth = 22
figHgtIndiv = 16*.5
figWthIndiv = 22*.5

psdFull[psdFull > .25] = 0.25
psdBreaks = c(0, 0.01, 0.03, 0.05, 0.075, 0.10, 0.15, 0.2, 0.25)

make_areal_map <- function(data, variables = NULL, breaks, legendName = 'Raw proportions', map_data = usShp, facet = TRUE, col = terrain.colors, zero_color = terrain.colors(40)[39], reverse_colors = TRUE, print = TRUE, ncol = 5, legend = TRUE) {
  make_map <- function(p) {
    
    nm <- data$variable[1]
    nm <- gsub("ZZZ", "/", nm)
    nm <- gsub("\\.", " ", nm)

    col <- col(length(breaks)-1)
    if(reverse_colors) col <- rev(col)
    if(legend) guide <- "legend" else guide <- "none"
    col[1] <- zero_color
    d <- ggplot(data, aes(X, Y, group = id)) + geom_polygon(aes(fill = as.factor(value))) + 
      scale_fill_manual(values = col, labels = breaklabels, name = legendName, guide = guide) +
        theme(strip.text.x = element_text(size = 16), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size = 16), legend.title = element_text(size = 16))
    d <- add_map_albers(plot_obj = d, map_data = map_data, dat = data)
    if(facet) {
      d <- d + facet_wrap(~variable, ncol = ncol)
    } else {
      d <- d + ggtitle(nm)
    }
    d <- theme_clean(d)
    return(d)
  }

  if(!is.null(variables))
    data <- data[data$variable %in% variables, ]
  data$value <- cut(data$value, breaks, include.lowest = TRUE, labels = FALSE)
  
  breaklabels <- apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                       function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })

  if(facet) {
    d <- make_map(p)
    if(print) print(d)
  } else {
    for ( p in variables ) {
      d <- make_map(p)
      if(print) print(d)
    }
  }
  invisible(d)
}

presence_long <- melt(fort_presence, c('X', 'Y', 'order', 'hole', 'piece', 'group', 'id'))
presence_long <- presence_long[presence_long$value == 1, ]

col = gray(c(.2,.5,.8))
d <- ggplot(presence_long, aes(X, Y, group = id)) + geom_polygon(aes(fill = as.factor(value))) + 
  scale_fill_manual(values = col, guide = FALSE)
d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = presence_long)
d <- theme_clean(d)

pdf('fig1.pdf', width=4, height=3)
print(d)
dev.off()


figs <- list()
length(figs) <- length(focalTaxa) * 3
cnt <- 1
for( taxon in focalTaxa ) {
  figs[[cnt]] <- make_areal_map(data = taxon_dat_long, variables = taxon, breaks = propBreaks, legendName = 'raw proportions', map_data = usFortified, facet = FALSE, ncol = 1, legend = FALSE) + theme(plot.margin = unit(rep(0,4), 'lines'))
  figs[[cnt + 1]] <- make_veg_map(data = pmFull[ , taxon, drop = FALSE], breaks = propBreaks, coords = coordFull, legendName = 'fitted proportions', map_data = usFortified, facet = FALSE, ncol = 1, legend = FALSE) + theme(plot.margin = unit(rep(0,4), 'lines'))
  figs[[cnt + 2]] <- make_veg_map(data = psdFull[ , taxon, drop = FALSE], breaks = psdBreaks, coords = coordFull, legendName = 'standard error', map_data = usFortified, facet = FALSE, ncol = 1, legend = FALSE) + theme(plot.margin = unit(rep(0,4), 'lines'))
  cnt <- cnt + 3
}

setwd(paperDir)

pdf('fig2.pdf', width=7, height=10)
do.call(grid.arrange, c(figs, nrow = length(focalTaxa), ncol = 3))
dev.off()

scaling = .5
ix=1
iy=seq(0,1,len=300)
iz=matrix(iy,nc=300)
pdf(file='legendRaw.pdf',width=5/scaling,height=.7/scaling,paper="special")
par(fig=c(0,1,0,1),mai=c(0.6,.2,0,.2),mgp=c(1.8,.7,0),cex.axis=1.5)
cols <- rev(terrain.colors(length(propBreaks)-1))
cols[1] <- terrain.colors(40)[39]
image(iy,ix,t(iz), yaxt = "n", xlab = "", xaxt='n', 
            ylab = "", col = cols, bty= 'n')
axis(1, at = seq(0,1, len = length(propBreaks)), labels = as.character(propBreaks),
     cex.axis = 2, tick = FALSE)
box()
dev.off()

ix=1
iy=seq(0,1,len=300)
iz=matrix(iy,nc=300)
pdf(file='legendProp.pdf',width=5/scaling,height=.7/scaling,paper="special")
par(fig=c(0,1,0,1),mai=c(0.6,.2,0,.2),mgp=c(1.8,.7,0),cex.axis=1.5)
image(iy,ix,t(iz), yaxt = "n", xlab = "", xaxt='n', 
            ylab = "", col = rev(terrain.colors(length(propBreaks)-1)))
axis(1, at = seq(0,1, len = length(propBreaks)), labels = as.character(propBreaks),
     cex.axis = 2, tick = FALSE)
box()
dev.off()

ix=1
iy=seq(0,1,len=300)
iz=matrix(iy,nc=300)
pdf(file='legendSD.pdf',width=5/scaling,height=.7/scaling,paper="special")
par(fig=c(0,1,0,1),mai=c(0.6,.2,0,.2),mgp=c(1.8,.7,0),cex.axis=1.5)
image(iy,ix,t(iz), yaxt = "n", xlab = "", xaxt='n', 
            ylab = "", col = rev(terrain.colors(length(psdBreaks)-1)))
axis(1, at = seq(0,1, len = length(psdBreaks)), labels = as.character(psdBreaks),
     cex.axis = 2, tick = FALSE)
box()
dev.off()


  
# Fig 3

  # take plot.full code with appropriate # cols passed to make_veg_map



  

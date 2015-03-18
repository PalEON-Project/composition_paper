# code to generate figures for paper

paperDir <- getwd()

setwd('~/research/jmac/composition')
source("config")

source(file.path(codeDir, "plot.R"))
source(file.path(codeDir, "prep_plot_data.R"))
source(file.path(codeDir, "set_domain.R"))

require(ncdf4)
require(ggplot2)
require(maptools)
require(rgdal)
require(raster)
require(gridExtra)

# for distinguishing no trees from no data
west_presence <- read.csv(file.path(dataDir, 'western.csv'))
west_presence <- west_presence[west_presence$x <= easternLimitOfWesternDomain, ]
west_presence$data <- as.numeric(rowSums(west_presence[ , c(3:ncol(west_presence))]) > 0)
# get in same order as data for fitting (from far N, W to E in rows)
tmp <- matrix(west_presence$data, length(westernDomainY), length(westernDomainX))
west_presence <- c(t(tmp[rev(seq_along(westernDomainY)), ]))

west_presence_fort <- west_fort
west_presence_fort$id <- as.numeric(west_presence_fort$id)
west_presence_fort <- cbind(west_presence_fort, west_presence[west_presence_fort$id])
names(west_presence_fort)[ncol(west_presence_fort)] <- "data"

east_presence_fort <- east_fort
east_presence_fort$data <- as.numeric(rowSums(east_presence_fort[ , 8:ncol(east_presence_fort)], na.rm = TRUE) > 0)

fort_presence <- rbind(west_presence_fort[ , c(1:7, ncol(west_presence_fort))], east_presence_fort[ , c(1:7, ncol(east_presence_fort))])

  
#########################################################################################
## Fig. 1
#########################################################################################

presence_long <- melt(fort_presence, c('X', 'Y', 'order', 'hole', 'piece', 'group', 'id'))
presence_long <- presence_long[presence_long$value == 1, ]

col = gray(c(.2,.5,.8))
d <- ggplot(presence_long, aes(X, Y, group = id)) + geom_polygon(aes(fill = as.factor(value))) + 
  scale_fill_manual(values = col, guide = FALSE) 
d <- add_map_albers(plot_obj = d, map_data = usFortified, dat = presence_long)
d <- d + scale_y_continuous(limits = c(70000, 1490000)) 
d <- theme_clean(d)

setwd(paperDir)

pdf('fig1.pdf', width=4.85, height=3)
print(d)
dev.off()

###########################################################################

focalTaxa <- c('Beech', 'Cherry', 'Chestnut', 'Elm', 'Hemlock', 'Oak', 'Pine')
focalTaxa <- focalTaxa[-2]

propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

figHgt = 16
figWth = 22
figHgtIndiv = 16*.5
figWthIndiv = 22*.5

psd[psd > .25] = 0.25
psdBreaks = c(0, 0.01, 0.03, 0.05, 0.075, 0.10, 0.15, 0.2, 0.25)




figs <- list()
length(figs) <- length(focalTaxa) * 3
cnt <- 1
for( taxon in focalTaxa ) {
  figs[[cnt]] <- make_areal_map(data = taxon_dat_long, variables = taxon, breaks = propBreaks, legendName = 'raw proportions', map_data = usFortified, facet = FALSE, ncol = 1, legend = FALSE, title = FALSE) + theme(plot.margin = unit(rep(0,4), 'lines'))
  figs[[cnt + 1]] <- make_veg_map(data = pm[ , taxon, drop = FALSE], breaks = propBreaks, coords = coord, legendName = 'fitted proportions', map_data = usFortified, facet = FALSE, ncol = 1, legend = FALSE) + theme(plot.margin = unit(rep(0,4), 'lines'))
  figs[[cnt + 2]] <- make_veg_map(data = psd[ , taxon, drop = FALSE], breaks = psdBreaks, coords = coord, legendName = 'standard error', map_data = usFortified, facet = FALSE, ncol = 1, legend = FALSE, col = heat.colors, title = FALSE) + theme(plot.margin = unit(rep(0,4), 'lines'))
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
par(fig=c(0,1,0,1),mai=c(0.6,.2,0,.3),mgp=c(1.8,.7,0),cex.axis=1.5)
cols <- rev(terrain.colors(length(propBreaks)-1))
cols[1] <- terrain.colors(40)[39]
image(iy,ix,t(iz), yaxt = "n", xlab = "", xaxt='n', 
            ylab = "", col = cols, bty= 'n')
axis(1, at = seq(0,1, len = length(propBreaks)), labels = as.character(propBreaks),
     cex.axis = 1.8, tick = FALSE)
box()
dev.off()

ix=1
iy=seq(0,1,len=300)
iz=matrix(iy,nc=300)
pdf(file='legendProp.pdf',width=5/scaling,height=.7/scaling,paper="special")
par(fig=c(0,1,0,1),mai=c(0.6,.2,0,.3),mgp=c(1.8,.7,0),cex.axis=1.5)
image(iy,ix,t(iz), yaxt = "n", xlab = "", xaxt='n', 
            ylab = "", col = rev(terrain.colors(length(propBreaks)-1)))
axis(1, at = seq(0,1, len = length(propBreaks)), labels = as.character(propBreaks),
     cex.axis = 1.8, tick = FALSE)
box()
dev.off()

ix=1
iy=seq(0,1,len=300)
iz=matrix(iy,nc=300)
pdf(file='legendSD.pdf',width=5/scaling,height=.7/scaling,paper="special")
par(fig=c(0,1,0,1),mai=c(0.6,.2,0,.3),mgp=c(1.8,.7,0),cex.axis=1.5)
image(iy,ix,t(iz), yaxt = "n", xlab = "", xaxt='n', 
            ylab = "", col = rev(heat.colors(length(psdBreaks)-1)))
axis(1, at = seq(0,1, len = length(psdBreaks)), labels = as.character(psdBreaks),
     cex.axis = 1.8, tick = FALSE)
box()
dev.off()


  
# Fig 3

# take plot.full code with appropriate # cols passed to make_veg_map

propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

figHgt = 11
figWth = 8.5

pdf(file = 'fig3.pdf', height = figHgt, width = figWth)
make_veg_map(data = pm, breaks = propBreaks, coords = coord, legendName = 'fitted proportions', map_data = usFortified, facet = TRUE, ncol = 3)
dev.off()


  

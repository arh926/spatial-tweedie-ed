
# usa_shape <- readRDS("GADM_2.8_USA_adm2.rds") # you shuld get this from raster::getData function
# ct_shape <- subset(usa_shape,NAME_1=="Connecticut")

sp_plot <- function(col.seq.length=NULL,
                    col.text = NULL,
                    data_frame=NULL,
                    shape=NULL, # shape file like shown above for CT
                    zlim=NULL,
                    lab=NULL){
  col.br <- colorRampPalette(RColorBrewer::brewer.pal(col.seq.length,col.text))
  surf <- MBA::mba.surf(data_frame,no.X=300,no.Y=300,extend=T)$xyz.est # , sp=T
  # sp::proj4string(surf) <- sp::proj4string(shape)
  # surf <- surf[!is.na(over(surf,shape))[,1],]
  # surf <- as.image.SpatialGridDataFrame(surf)
  if(is.null(zlim)) zlim <- range(surf[["z"]],na.rm=T)
  fields::image.plot(surf, xaxs="i", yaxs="i",col=rev(col.br(100)), axes=T,
                     zlim=zlim,
                     #xlim=shape@bbox["x",],
                     #ylim=shape@bbox["y",],
                     xlab=latex2exp::TeX("Longitude$\\degree$"),
                     ylab=latex2exp::TeX("Latitude$\\degree$"),
                     legend.lab = lab)
  # plot(shape, add=T)
  # points(data_frame[,(1:2)], pch="+", cex=0.1)
  contour(surf,add=T) #lwd = 0.1,labcex = 0.1
}

# Test
# require(sp)
# require(spdep)
# require(rgeos)
# require(rgdal)
# require(sf)
# require(MBA)
# generate spatial effect
# blockwise
# pdf("/Users/aritrahalder/Desktop/speff.pdf")


# sp_plot(11,"Spectral",cbind(ct_zip[,c("longitude","latitude")],speff),ct_shape)
# plot(ct_shape)
# points(ct_zip[,c("longitude","latitude")], pch="+", col=col,cex=.5)


# sp_plot(11,"Spectral",cbind(ct_zip[,c("longitude","latitude")],speff),ct_shape)



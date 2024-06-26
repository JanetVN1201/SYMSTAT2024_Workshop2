library(leaflet)
library(leafpop)
library(leaflet.extras)
library(htmltools)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(lattice)
library(htmlwidgets)
library(widgetframe)
library(shiny)
library(sp)
library(wallace)
library(mapdata)
library(raster)
library(terra)
require(devtools)
library(fields)
library(dplyr)
#install_github("eblondel/cleangeo")
library(cleangeo)
library(inlabru)
library(fmesher)
library(INLAspacetime)
par(mfrow= c(1,1))

Data1 <- readRDS("~/Documents/GitHub/RSG_Projects/Data/Survey data of RSP/Helicopter and incidental sightings/incidental_sightings_dugongs.rds")

Data2 <- Data1

coords <- cbind(as.matrix(Data1$Long), as.matrix(Data1$Lat) )

d_icon = "https://cdn-icons-png.flaticon.com/512/8652/8652530.png"

plot1 <- leaflet() %>% 
  setView(lng = median(Data2$Long), lat = median(Data2$Lat), zoom = 8) %>%
  addProviderTiles(providers$OpenStreetMap, group = "Open Street Map") %>% 
  # Add additional basemap layers
  addProviderTiles(providers$Esri.WorldImagery, group = "ESRI World Imagery") %>% 
  addProviderTiles(providers$Esri.OceanBasemap, group = "ESRI Ocean Basemap") %>% 
  # Add a User-Interface (UI) control to switch layers
  addLayersControl(
    baseGroups = c("Open Street Map","ESRI World Imagery","ESRI Ocean Basemap"),
    options = layersControlOptions(collapsed = FALSE)) 

frameWidget(plot1, width = "100%", height = "500")

plot2 <- leaflet() %>% 
  setView(lng = mean(Data2$Long), lat = mean(Data2$Lat), zoom = 8) %>%
  addProviderTiles(providers$OpenStreetMap, group = "Open Street Map") %>% 
  # Add additional basemap layers
  addProviderTiles(providers$Esri.WorldImagery, group = "ESRI World Imagery") %>% 
  addProviderTiles(providers$Esri.OceanBasemap, group = "ESRI Ocean Basemap") %>% 
  # Add a User-Interface (UI) control to switch layers
  addLayersControl(
    baseGroups = c("Open Street Map","ESRI World Imagery","ESRI Ocean Basemap"),
    options = layersControlOptions(collapsed = FALSE)) %>%
  addMarkers(lng = Data2$Long, lat = Data2$Lat, group = "Dugong",
             icon = 
               list(iconUrl = d_icon,
                    iconSize = c(30,30)))

frameWidget(plot2, width = "100%", height = "500")

#Use Red sea bathymetry to get the barriers
ras1= raster("~/Documents/GitHub/RSG_Projects/Data/Survey data of RSP/Helicopter and incidental sightings/RSPbathym100m.tif")

rs_bath <- projectRaster(ras1, crs = CRS("+proj=longlat +datum=WGS84"))

#summary(values(rs_bath))
rs_bath1 <- rs_bath
rs_bath1[rs_bath1>-2] <- NA #Sand parts - common for all

rs_bath1A <- rast(rs_bath1)
plot(rs_bath1A)

poly_sea <- as.polygons(rs_bath1A > -Inf)
plot(poly_sea, border = "blue", add = TRUE, lwd = 1)
poly_sea1 <- as(poly_sea, "Spatial")

#Clean - maybe not needed
if (T){
rr <- clgeo_CollectionReport(poly_sea1)
summary(rr)
issues <- rr[rr$type == NA,]
nv <- clgeo_SuspiciousFeatures(rr)
my_poly_sea <- poly_sea1
poly_sea_c <- clgeo_Clean(my_poly_sea)
}
poly_sea2 <- poly_sea_c

plot4 <- plot2 %>%
  addPolygons(data = poly_sea2,
              group = "Red Sea",
              color = "blue",
              weight = 0.5)

frameWidget(plot4, width = "100%", height = "500")

kmproj <- CRS("+proj=utm +zone=39 ellps=WGS84 +units=km")
poly.RS_l = spTransform(poly_sea2, kmproj)

xx <- spsample(SpatialPolygons(poly_sea2@polygons), type = "random", n = 25)

coords <- rbind(coords, xx@coords) #Adding some extra coordinatess

n = nrow(coords)
Data2 <- data.frame(rep(NA, n))
Data2$Long<- coords[,1]
Data2$Lat <- coords[,2]
Data2 <- Data2[,-1]
Data2$y <- c(rep(1,25), rep(0,25))

coords_p <- spTransform(SpatialPoints(coords,proj4string = CRS("+proj=longlat +datum=WGS84")), kmproj)

mesh2 <- inla.mesh.2d(boundary = poly.RS_l, max.edge = c(2, 10), offset = c(10, 10), cutoff = c(1))
mesh2$crs <- kmproj

nv = mesh2$n
library(fmesher)
water.tri = fm_contains(x = poly.RS_l, y = mesh2,
                              type = "centroid", ignore.CRS = TRUE)
num.tri = length(mesh2$graph$tv[, 1])
barrier.tri = setdiff(1:num.tri, water.tri)

plot(mesh2)
points(coords_p, pch = 8, cex = 0.8, col = c(rep(2, 25), rep(5, 25)))


bmodel <- barrierModel.define(
  mesh = mesh2, 
  barrier.triangles = barrier.tri,
  prior.range = c(1, 0.01), ## P(range < 1) = 0.01
  prior.sigma = c(20, 0.01), ## P(sigma > 1) = 0.01
  range.fraction = 0.1)

model <- y ~ Intercept(1) + 
  field(coords_p, model = bmodel)

result <- bru(
  model, data = Data2, family = "binomial")

pred_map <- rasterize(x = mesh2$loc[,1:2], y = rs_bath1p1, field = result$summary.random$field$mean, fun = mean)
pred_map <- projectRaster(pred_map, crs = CRS("+proj=longlat +datum=WGS84"))
pal <- colorNumeric("viridis", domain = c(-0.5,2), na.color = "transparent")

raster::plot(rs_bath1p1)
plot8 <- plot2 %>%
  leaflet::addRasterImage(pred_map, colors = pal, opacity = 1) %>%
  addScaleBar(position = c("bottomleft"))%>%
  leaflet::addLegend("bottomright",
                     pal = pal,
                     values = seq(-0.5,2,by = 0.05))

frameWidget(plot8, width = "100%", height = "500")

pred_sd = rasterize(x = mesh2$loc[,1:2], y = rs_bath1p1, field = result$summary.random$field$sd, fun = mean)
pred_sd <- projectRaster(pred_sd, crs = CRS("+proj=longlat +datum=WGS84"))
pal <- colorNumeric("viridis", domain = c(0,2), na.color = "transparent")

plot9 <- plot2 %>%
  leaflet:::addRasterImage(pred_map, colors = pal, opacity = 0.8) %>%
  addScaleBar(position = c("bottomleft"))%>%
  leaflet::addLegend("bottomright",
                     pal = pal,
                     values = seq(0,2, by = 0.01))

frameWidget(plot9, width = "100%", height = "500")


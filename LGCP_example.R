library(sf)
library(splancs)
library(fmesher)
library(INLA)
library(inlabru)
library(ggplot2)
library(gridExtra)

## load Burkitt data from the splancs package
data('burkitt')

ls()

head(burkitt)

par(mar = c(0,0,0,0))
plot(burbdy, bty = 'n', type = "l", asp = 1)
with(burkitt, points(x, y, cex = age/10))

## data as sf
bkt <- st_as_sf(burkitt, coords = c("x", "y"))
bnd <- st_sf(st_sfc(st_polygon(list(burbdy))))

## mesh
mesh <- fm_mesh_2d(
    boundary = bnd, 
    max.edge = c(5, 15),
    cutoff = 2)

## visualize
ggplot() + theme_minimal() + 
    geom_sf(data = bnd) +
    gg(mesh) +
    geom_sf(data = bkt) 

## The spatial model
spde <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(5, 0.01),
    prior.sigma = c(0.3, 0.01)
)

## model components 
components <- ~  
    Intercept(1) +        ## "inlabru" way of doing it
    spatial(geometry,     ## index on the geometry
            model = spde) ## the actual model definition

## likelihood object
lhood <- like(
    geometry ~ .,         ## locations are the outcome
    family = "cp",
    data = bkt,
    domain = list(geometry = mesh),
    samplers = bnd)

## fit the model
fit <- bru(
    components,
    lhood)

fit$cpu.used
fit$bru_timings

## summaries
fit$summary.hyperpar

## prepare for visualization
post.range <- spde.posterior(fit, name = "spatial", what = "range")
post.sigma <- spde.posterior(fit, name = "spatial", what = "variance")
post.corr <- spde.posterior(fit, name = "spatial", what = "matern.correlation")

## visualize
grid.arrange(
    plot(post.range),
    plot(post.sigma),
    plot(post.corr),
    nrow = 1)

fit$summary.fixed

nrow(burkitt) / st_area(bnd)
exp(-6.61)

## the bounding box
st_bbox(bnd)
bb <- matrix(st_bbox(bnd), 2)
bb

apply(bb, 1, diff)


## setup a grid
grid <- fm_pixels(
    mesh = mesh,
    dims = c(95, 182),
    mask = bnd,
    format = "sf")

str(grid)

## inlabru::predict()
##  drawn (Monte Carlo: independent) samples
## from the model parameters posterior fitted by INLA
## and compute functions from these
pred <- predict(
    fit, 
    grid, 
    ~ data.frame(
    lambda = exp(spatial),
    loglambda = spatial
    )
)

str(pred)

## visualize
ggplot() + theme_minimal() +
    geom_sf(data = pred$loglambda,            
            aes(color = mean)) +
    scale_color_distiller(
        palette = "Spectral"
    ) +
    geom_sf(data = bkt, cex = 0.1) +
  theme(legend.title = element_blank(),
        legend.key.height = unit(1, "null")) +
  guides(colour = guide_colourbar(position = "right"))

## see a complete example at
## https://inlabru-org.github.io/inlabru/articles/2d_lgcp.html

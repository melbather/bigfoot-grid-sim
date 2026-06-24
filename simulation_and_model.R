# Code to simulate covariate information across a study area

# Libaries
library(dplyr)
library(spatstat.geom)
library(spatstat.model)
library(mvtnorm)

# Pull in BFRO dataset and clean
bfro_data <- read.csv("bfro_reports_geocoded.csv")

bfro_data_clean <- bfro_data |> 
  filter(
    !(county == "Cowlitz County" & date == "2006-07-05"),
    !(county == "Flathead County" & date == "2012-09-16"),
    !(county == "Idaho County" & date == "2013-12-04"),
    !is.na(latitude),
    !is.na(longitude)
    ) |> 
  mutate(
    rain = grepl("rain", precip_type),
    snow = grepl("snow", precip_type)
  )

# Create study window
W <- owin(
    xrange = range(bfro_data_clean$longitude),
    yrange = range(bfro_data_clean$latitude)
  )

# ppp object
bigfoot_pp <- ppp(
  x = bfro_data_clean$longitude,
  y = bfro_data_clean$latitude,
  window = W
)

# Create a grid
nx <- 40
ny <- 40

xgrid <- seq(W$xrange[1], W$xrange[2], length.out = nx)
ygrid <- seq(W$yrange[1], W$yrange[2], length.out = ny)

grid <- expand.grid(x = xgrid, y = ygrid)
n_grid <- nrow(grid)

# Create spatial covariance matrix
D <- as.matrix(dist(grid))
Sigma_spatial <- exp(-D / 5) + diag(1e-6, n_grid)

weather <- matrix(
  c(
    1.0, -0.4,  0.5,
   -0.4,  1.0, -0.6,
    0.5, -0.6,  1.0
  ),
  nrow = 3,
  byrow = TRUE
)

Sigma_full <- kronecker(weather, Sigma_spatial)

# Latent Gaussian field
latent <- rmvnorm(
  n = 1,
  mean = rep(0, 3 * n_grid),
  sigma = Sigma_full
)

latent <- as.numeric(latent)

z_temp <- latent[1:n_grid]
z_precip <- latent[(n_grid + 1):(2 * n_grid)]
z_visibility <- latent[(2 * n_grid + 1):(3 * n_grid)]

# Helper function for rescaling
rescale_to_range <- function(z, new_min, new_max) {
  z_scaled <- (z - min(z)) / (max(z) - min(z))
  new_min + z_scaled * (new_max - new_min)
}

# Create covariates
temp_vals <- rescale_to_range(z_temp, 25, 85)
precip_vals <- rescale_to_range(z_precip, 0, 0.25)
vis_vals <- rescale_to_range(z_visibility, 2, 10)

# Rain more likely when precipitation is high
rain_vals <- as.numeric(precip_vals > 0.10)

# Snow more likely when precipitation is high and temperature is low
snow_vals <- as.numeric(precip_vals > 0.10 & temp_vals < 35)

# Convert into spatstat images
temperature_mid <- im(
  matrix(temp_vals, nrow = nx, ncol = ny),
  xcol = xgrid,
  yrow = ygrid
)

precip_intensity <- im(
  matrix(precip_vals, nrow = nx, ncol = ny),
  xcol = xgrid,
  yrow = ygrid
)

visibility <- im(
  matrix(vis_vals, nrow = nx, ncol = ny),
  xcol = xgrid,
  yrow = ygrid
)

rain <- im(
  matrix(rain_vals, nrow = nx, ncol = ny),
  xcol = xgrid,
  yrow = ygrid
)

snow <- im(
  matrix(snow_vals, nrow = nx, ncol = ny),
  xcol = xgrid,
  yrow = ygrid
)

# IPP model
bigfoot_ipp <- ppm(
  bigfoot_pp ~ temperature_mid +
       precip_intensity +
       rain +
       snow +
       visibility,
  covariates = list(
    temperature_mid = temperature_mid,
    precip_intensity = precip_intensity,
    rain = rain,
    snow = snow,
    visibility = visibility
  )
)

summary(bigfoot_ipp)

# Plot predicted intensity
pred <- predict(bigfoot_ipp)
plot(pred)
plot(bigfoot_pp, add = TRUE, pch = 16, cex = 0.3)
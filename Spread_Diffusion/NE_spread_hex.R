### Title
### Description 
### Date

### ### Libraries
library(sf)
library(ggplot2)
library(patchwork)
library(spData)
library(terra)
library(raster)
library(tidyverse)


### Landscape
# Data source - https://www.weather.gov/gis/Counties
US_county_outlines <- st_read(
  "Spread_Diffusion/Data/USCountyOutlines/c_05mr24.shp")

# subset counties (vectors)
NEcounties <- US_county_outlines[US_county_outlines$STATE == "MA" |
                                   US_county_outlines$STATE == "NY" |
                                   US_county_outlines$STATE == "PA" |
                                   US_county_outlines$STATE == "VT" |
                                   US_county_outlines$STATE == "NJ" |
                                   US_county_outlines$STATE == "CT",]

# try and deal with same county name issue
NEcounties$Count_state <- paste0(NEcounties$STATE, NEcounties$COUNTYNAME)

## As one big region (should make below faster..)
# Combine all polygons into one
NE_single_polygon <- st_union(NEcounties)

# Convert the result into an sf object
NE_single_polygon <- st_sf(geometry = NE_single_polygon)

## Set up hex grid
plot(st_make_grid(NEcounties, cellsize = .2, square = FALSE))
plot(NEcounties, add = TRUE)

hex <- st_make_grid(NEcounties, cellsize = .2, square = FALSE)
# Will need to change size later...

# Convert hex grid to an sf object
hex_sf <- st_sf(geometry = hex)

# cut to state boarder.. 

# Assign CRS to hex_sf (assume NEcounties' CRS is correct and use it for hex_sf)
st_crs(hex_sf) <- st_crs(NEcounties)


# Check and fix invalid geometries in hex_sf
hex_sf <- st_make_valid(hex_sf)

# Find intersections between hexagonal grid and NEcounties
NE_hex_grid <- st_intersection(NE_single_polygon, hex_sf) ## This takes forever
# This looks a litte weird because not full hex on edges, but will still work the same


# Plot results
plot(st_geometry(hex_sf), col = 'lightblue', main = 'Intersection of Hex Grid and NEcounties')
plot(st_geometry(NEcounties), add=TRUE, border = 'red')
plot(st_geometry(NE_hex_grid), add = TRUE, col = 'green', border = 'black', lwd = 2)


### ID each cell with number and empty prev 
NE_hex_grid$value <- seq_along(NE_hex_grid$geometry)
NE_hex_grid_df <- as.data.frame(NE_hex_grid$value)
NE_hex_grid_df$prev <- c(0)
NE_hex_grid$prev <- NE_hex_grid_df$prev


### Identify the boundaries 

mymat_adj <- st_intersects(NE_hex_grid, NE_hex_grid)

bordervec <- c()
for(i in 1:nrow(mymat_adj)){
  adj <- (unlist(mymat_adj[i])) 
  if(length(adj) != 7){
    bordervec <- c(bordervec, i)
  }
  
}



NE_bordercells <- NE_hex_grid[bordervec, ]
## Will have to deal with this later - looks kinda wonky but I think it's right
# Probably best going through GIS to ID water body hard boundaries... 






### Identify which PA counties have CWD 
# choosing a few to start with, can correct later 

PA_counties_data <- readxl::read_excel("Spread_Diffusion/Data/PACountyPrev_9_8_24.xlsx", col_names = TRUE) # assembled from PA heat map - 7/8/24

PAcounties <- US_county_outlines[US_county_outlines$STATE == "PA" ,]
PAsubset_Bedford <- PAcounties[PAcounties$COUNTYNAME == "Bedford",]
PAsubset_Somerset <- PAcounties[PAcounties$COUNTYNAME == "Somerset",]

# Ensure valid geometry for the subset county
PAsubset_Bedford <- st_make_valid(PAsubset_Bedford)
PAsubset_Somerset <- st_make_valid(PAsubset_Somerset)

PAgrid <- st_intersection(PAcounties, NE_hex_grid)


### For each county with non zero prevalence.. 
for(i in unique(PA_counties_data$County)){
  ## Pull out the county as a geo subset
  PAsubset_Bedford <- PAcounties[PAcounties$COUNTYNAME == i,]
  #Make it a spatial object and intersect with hex grid
  PAsubset_Bedford <- st_make_valid(PAsubset_Bedford)
  Bed <- st_intersection(PAsubset_Bedford, NE_hex_grid)
  
  # Add prevalence 
  IDcountyprev <- PA_counties_data[PA_counties_data$County == i, 4]
  Bed$prev <- rep(unlist(IDcountyprev), length(Bed$value))
  Bed <- st_make_valid(Bed)
  # Get the area
  Bed$area <- st_area(Bed)
  
  # Differentiate between first object and all others
  if(i == "Bedford"){
    combined <- Bed
  } else{
    # Combine all counties into one dataframe
    combined <- bind_rows(combined, Bed)
  }
  
}





# Ensure combined contains geometry and area
combined <- combined %>%
  mutate(area = as.numeric(st_area(geometry)))  # Ensure area is numeric


# Compute the area of each intersected polygon
combined$area <- st_area(combined)
NE_hex_grid$area <- st_area(NE_hex_grid)

# Compute total area for each hexagon
total_area_df <- NE_hex_grid %>%
  st_set_geometry(NULL) %>%
  group_by(value) %>%
  summarize(total_area = sum(area, na.rm = TRUE), .groups = 'drop') %>% 
  mutate(total_area = as.vector(total_area))

# Calculate proportional area and weighted prevalence
aggregated <- combined %>%
  st_set_geometry(NULL) %>%
  left_join(total_area_df, by = "value") %>%
  mutate(proparea = area / total_area) %>%
  group_by(value) %>%
  summarize(weightedprev = sum(proparea * prev, na.rm = TRUE), .groups = 'drop') %>% 
  mutate(weightedprev = as.vector(weightedprev)) #drop units

# Add the aggregated prevalence to the original hexagonal grid
NE_hex_grid2 <- NE_hex_grid %>%
  left_join(aggregated, by = "value") %>%
  mutate(prev = replace_na(weightedprev, 0)) %>%
  select(-weightedprev)

# Optionally, plot the results to check alignment
ggplot() +
  geom_sf(data = NE_hex_grid2, aes(fill = prev)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Hexagonal Grid with Aggregated Prevalence Values")










Bed <- st_intersection(PAsubset_Bedford, NE_hex_grid)
somer <- st_intersection(PAsubset_Somerset, NE_hex_grid)


## Add prevalence levels
IDcountyprev <- PA_counties_data[PA_counties_data$County == "Somerset", 4]
somer$prev <- rep(unlist(IDcountyprev), length(somer$value))

# doing this mannually, later automate
IDcountyprev <- PA_counties_data[PA_counties_data$County == "Bedford", 4]
Bed$prev <- rep(unlist(IDcountyprev), length(Bed$value))


plot(Bed$prev)
plot(somer$prev, add = TRUE)
plot(Bed$geometry); plot(Bed[, "prev"], add = TRUE)
plot(somer[, "prev"], add = TRUE)

## Get area
Bed$area <- st_area(Bed)
somer$area <- st_area(somer)

ave_height <- st_join(Bed, somer)
#, by = c("STATE", "CWA","COUNTYNAME","FIPS","TIME_ZONE","FE_AREA","LON","LAT","geometry") )

st_touches(Bed, somer)





combined_polygon <- st_union(Bed, somer)

### Going to need to make repro on this one.. not super straigt forward... 
# could do in a loop and run through each value? 



c("STATE", "CWA","COUNTYNAME","FIPS","TIME_ZONE","FE_AREA","LON","LAT","geometry")

a = st_sf(a = 1:3,
          geom = st_sfc(st_point(c(1,1)), st_point(c(2,2)), st_point(c(3,3))))
b = st_sf(a = 11:14,
          geom = st_sfc(st_point(c(10,10)), st_point(c(2,2)), st_point(c(2,2)), st_point(c(3,3))))

me <- st_join(a,b) %>% aggregate(list(.$a.x), mean)









names(PA_counties) <- "PA_Count"
PA_counties$newname <- paste0("PA", PA_counties$PA_Count)

## Get subset of PA counties as hex cells and set prevalence?

startingrast <- NE_count_rast == "Warren"
for(i in (PA_counties$newname)){
  pullcount <- (NE_count_rast == i)
  startingrast <- startingrast + pullcount
}






### ### repro example to get avg prev values from hex split across counties 

library(sf)
library(ggplot2)

# Define a simple bounding box as an sfc object
sfc <- st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))

# Create a hexagonal grid
hex <- st_make_grid(sfc, cellsize = 0.1, square = FALSE)

# Convert the grid to an sf object
hex_sf <- st_sf(geometry = hex)

# Add a 'value' column with example values
hex_sf$value <- seq_along(hex_sf$geometry)  # Ensure this column is numeric and non-NA

# Create polygons that touch but do not overlap
# Coordinates adjusted to ensure only touching along one side

# Polygon 1
x1 <- c(0.1, 0.1, 0.2, 0.2, 0.1)
y1 <- c(0.1, 0.4, 0.4, 0.1, 0.1)

# Polygon 2, shifted to touch the side of Polygon 1 without overlapping
x2 <- c(0.2, 0.2, 0.3, 0.3, 0.2)
y2 <- c(0.1, 0.4, 0.4, 0.1, 0.1)

# Polygon 3, shifted to touch the side of Polygon 2 without overlapping
x3 <- c(0.3, 0.3, 0.4, 0.4, 0.3)
y2 <- c(0.1, 0.4, 0.4, 0.1, 0.1)

# Assign the vertices to polygons
poly1 <- st_polygon(list(cbind(x1, y1)))
poly2 <- st_polygon(list(cbind(x2, y2)))
poly3 <- st_polygon(list(cbind(x3, y2)))

# Create sf objects for the polygons
poly1_sf <- st_sf(geometry = st_sfc(poly1))
poly2_sf <- st_sf(geometry = st_sfc(poly2))
poly3_sf <- st_sf(geometry = st_sfc(poly3))

# Assign a common CRS to polygons and grid
crs <- st_crs(hex_sf)
st_crs(poly1_sf) <- crs
st_crs(poly2_sf) <- crs
st_crs(poly3_sf) <- crs

# Optionally, plot the results to check alignment
ggplot() +
  geom_sf(data = hex_sf, fill = "lightgrey") +
  geom_sf(data = poly1_sf, fill = "blue", alpha = 0.5) +
  geom_sf(data = poly2_sf, fill = "red", alpha = 0.5) +
  geom_sf(data = poly3_sf, fill = "green", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Polygons Touching Without Overlap")



### Okay so now set up the same


Bed <- st_intersection(poly1_sf, hex_sf)
somer <- st_intersection(poly2_sf, hex_sf)
new <- st_intersection(poly3_sf, hex_sf)

# Set prev value 
Bed$prev <- .1
somer$prev <- .4
new$prev <- .3

# get area
Bed$area <- st_area(Bed)
somer$area <- st_area(somer)
new$area <- st_area(new)

plot(Bed$geometry); plot(somer$geometry, add = TRUE)


### Start from here, weighted average is finicky... !!!!!!!!!!


combined_polygon <- st_union(Bed, somer)
# Combine the results
combined <- rbind(Bed, somer)


# Combine Bed and somer into one data frame
# combined <- bind_rows(Bed, somer)
# 
# # Aggregate prevalence values for each hexagon using unique IDs
# aggregated <- combined %>%
#   st_set_geometry(NULL) %>%
#   group_by(id) %>%
#   summarize(prev = sum(prev, na.rm = TRUE), .groups = 'drop')
# 
# # Add the aggregated prevalence to the original hexagonal grid
# hex_sf <- hex_sf %>%
#   left_join(aggregated, by = "id") %>%
#   mutate(prev = replace_na(prev, 0))
# 
# # Optionally, plot the results to check alignment
# ggplot() +
#   geom_sf(data = hex_sf, aes(fill = prev)) +
#   scale_fill_viridis_c() +
#   theme_minimal() +
#   labs(title = "Hexagonal Grid with Aggregated Prevalence Values")




# # Combine Bed and somer into one data frame
# combined <- bind_rows(Bed, somer)
# 
# # Compute the area of each intersected polygon
# combined$area <- st_area(combined)
# hex_sf$area <- st_area(combined)

# # Aggregate prevalence values for each hexagon using weighted average based on area
# aggregated <- combined %>%
#   st_set_geometry(NULL) %>%
#   group_by(id) %>%
#   summarize(
#     total_area = sum(area, na.rm = TRUE),
#     proparea = ((area / total_area) * 100),
#     weightedprev = sum(proparea * prev),
#     .groups = 'drop'
#   )
# 
# # Add the aggregated prevalence to the original hexagonal grid
# hex_sf <- hex_sf %>%
#   left_join(aggregated, by = "id") %>%
#   mutate(prev = replace_na(weighted_prev, 0)) %>%
#   select(-weighted_prev)
# 
# # Optionally, plot the results to check alignment
# ggplot() +
#   geom_sf(data = hex_sf, aes(fill = prev)) +
#   scale_fill_viridis_c() +
#   theme_minimal() +
#   labs(title = "Hexagonal Grid with Aggregated Prevalence Values")
# 





# 
# # Initialize prevalence column in hex_sf
# hex_sf$prev <- 0
# 
# # Calculate prev for each hexagon
# for (i in 1:nrow(hex_sf)) {
#   hex_geom <- hex_sf$geometry[i]
#   
#   
# }



library(sf)
library(dplyr)

# Assuming 'combined' is already defined as per previous context
# and contains the 'area' and 'prev' columns

# Combine Bed and somer into one data frame
combined <- bind_rows(Bed, somer, new)


# Ensure combined contains geometry and area
combined <- combined %>%
  mutate(area = as.numeric(st_area(geometry)))  # Ensure area is numeric



# Convert the grid to an sf object
hex_sf <- st_sf(geometry = hex)

# Add a 'value' column with example values
hex_sf$id <- seq_along(hex_sf$geometry)

# Compute the area of each intersected polygon
combined$area <- st_area(combined)
hex_sf$area <- st_area(hex_sf)

# Compute total area for each hexagon
total_area_df <- hex_sf %>%
  st_set_geometry(NULL) %>%
  group_by(id) %>%
  summarize(total_area = sum(area, na.rm = TRUE), .groups = 'drop')

# Calculate proportional area and weighted prevalence
aggregated <- combined %>%
  st_set_geometry(NULL) %>%
  left_join(total_area_df, by = "id") %>%
  mutate(proparea = area / total_area) %>%
  group_by(id) %>%
  summarize(weightedprev = sum(proparea * prev, na.rm = TRUE), .groups = 'drop')

# aggregated <- hex_sf %>%
#   st_set_geometry(NULL) %>%
#   left_join(total_area_df, by = "id") %>%
#   left_join(combined, by = "id") %>%
#   mutate(proparea = area / total_area) %>%
#   group_by(id) %>%
#   summarize(weightedprev = sum(proparea * prev, na.rm = TRUE), .groups = 'drop')


# Add the aggregated prevalence to the original hexagonal grid
hex_sf <- hex_sf %>%
  left_join(aggregated, by = "id") %>%
  mutate(prev = replace_na(weightedprev, 0)) %>%
  select(-weightedprev)

# Optionally, plot the results to check alignment
ggplot() +
  geom_sf(data = hex_sf, aes(fill = prev)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Hexagonal Grid with Aggregated Prevalence Values")



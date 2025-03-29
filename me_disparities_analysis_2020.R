# Read in data
dat = read.csv("~/Documents/ali-data/me_disparities_analysis_2020.csv")

# Create binned exposures
## Bin based on straight-line access
dat$STRAIGHT_PROXIMITY_BINNED = cut(x = dat$STRAIGHT_PROXIMITY, 
                                    breaks = c(-0.1, seq(1, 5, by = 1), 16))

## Bin based on map-based access
dat$MAP_PROXIMITY_BINNED = cut(x = dat$MAP_PROXIMITY, 
                               breaks = c(-0.1, seq(1, 5, by = 1), 16))

## Bin based on estimated % food insecure
dat$EST_FOOD_INSECURITY_BINNED = cut(x = dat$EST_FOOD_INSECURITY, 
                                     breaks = c(0, seq(10, 20, by = 2.5), 32))

## Bin based on lower-bound % food insecure
dat$LB_FOOD_INSECURITY_BINNED = cut(x = dat$LB_FOOD_INSECURITY, 
                                    breaks = c(0, seq(10, 20, by = 2.5), 32))

## Bin based on upper-bound % food insecure
dat$UB_FOOD_INSECURITY_BINNED = cut(x = dat$UB_FOOD_INSECURITY, 
                                    breaks = c(0, seq(10, 20, by = 2.5), 32))
# Replace factor labels 
dat = dat |> 
  dplyr::mutate(
    STRAIGHT_PROXIMITY_BINNED = factor(x = STRAIGHT_PROXIMITY_BINNED, 
                                       levels = c("(-0.1,1]", "(1,2]", "(2,3]", 
                                                  "(3,4]", "(4,5]", "(5,16]"),
                                       labels = c("Up to 1 Mile", "Between 1-2 Miles", "Between 2-3 Miles", 
                                                  "Between 3-4 Miles", "Between 4-5 Miles", "More than 5 Miles")), 
    MAP_PROXIMITY_BINNED = factor(x = MAP_PROXIMITY_BINNED, 
                                  levels = c("(-0.1,1]", "(1,2]", "(2,3]", 
                                             "(3,4]", "(4,5]", "(5,16]"),
                                  labels = c("Up to 1 Mile", "Between 1-2 Miles", "Between 2-3 Miles", 
                                             "Between 3-4 Miles", "Between 4-5 Miles", "More than 5 Miles")), 
    EST_FOOD_INSECURITY_BINNED = factor(x = EST_FOOD_INSECURITY_BINNED, 
                                        levels = c("(0,10]", "(10,12.5]", "(12.5,15]", 
                                                   "(15,17.5]", "(17.5,20]", "(20,32]"),
                                        labels = c("Up to 10%", "Between 10-12.5%", "Between 12.5-15%", 
                                                   "Between 15-17.5%", "Between 17.5-20%", "More than 20%")), 
    LB_FOOD_INSECURITY_BINNED = factor(x = LB_FOOD_INSECURITY_BINNED, 
                                       levels = c("(0,10]", "(10,12.5]", "(12.5,15]", 
                                                  "(15,17.5]", "(17.5,20]", "(20,32]"),
                                       labels = c("Up to 10%", "Between 10-12.5%", "Between 12.5-15%", 
                                                  "Between 15-17.5%", "Between 17.5-20%", "More than 20%")), 
    UB_FOOD_INSECURITY_BINNED = factor(x = UB_FOOD_INSECURITY_BINNED, 
                                       levels = c("(0,10]", "(10,12.5]", "(12.5,15]", 
                                                  "(15,17.5]", "(17.5,20]", "(20,32]"),
                                       labels = c("Up to 10%", "Between 10-12.5%", "Between 12.5-15%", 
                                                  "Between 15-17.5%", "Between 17.5-20%", "More than 20%"))
  )

# Replace patient IDs with 1:N and remove identifiable columns 
dat = dat |> 
  dplyr::mutate(ID = 1:dplyr::n()) |> 
  dplyr::select(ID, DIABETES, 
                MAP_PROXIMITY, STRAIGHT_PROXIMITY, 
                MAP_PROXIMITY_BINNED, STRAIGHT_PROXIMITY_BINNED, 
                EST_FOOD_INSECURITY, LB_FOOD_INSECURITY, UB_FOOD_INSECURITY, 
                EST_FOOD_INSECURITY_BINNED, LB_FOOD_INSECURITY_BINNED, UB_FOOD_INSECURITY_BINNED)


dat |> 
  write.csv("~/Documents/ME_Disparities_HSORM/me_disparities_analysis_2020.csv", row.names = FALSE)

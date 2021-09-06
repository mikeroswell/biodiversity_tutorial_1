

# We need to write a function that calculates a cumulative of species sampled.
# Currently, `mobsim::calc_PIE` calculates each quadrat individually. We need to
# write a function that computes ENS_PIE for each _collection_ of quadrats.

ens_PIE <- function(community_data){

  # First we're going to see how many plots are in the community data.frame. We will use this to determine how long we want our ENS_PIE vector to be, and how many times we need to repeat our calculation of the ENS_PIE metric.
  rows <- nrow(community_data)
  ens_PIE_value <- vector(length = rows)
  species_data <- community_data[1,]
  ens_PIE_value[1] <- calc_PIE(species_data, ENS = TRUE)

  for(i in 2:rows){

    # Here is where we diverge from the function as it was written. We are now
    # looping through the community dataset, making a cumulative count of the
    # species in the dataset that increases with respect to area/quadrat #. This
    # approximates but is not exactly the same as the "nested" sampling in the
    # paper. We then calculate the ENS_PIE metric on this cumulative dataset and
    # return the list that contains these values.
    species_data <- species_data + community_data[i,]
    ens_PIE_value[i] <- calc_PIE(species_data, ENS = TRUE)

  }

  return(ens_PIE_value)

}


# Plotting function

contr <- control_45
Tx <- reduce_density_45
x <- x
metric <- "rich"
meth = "accum"

plotter <- function(contr
                    , Tx
                    , x
                    , metric = c("rich", "ENS_PIE")
                    , meth = c("accum", "effect")
                    , title){
  # avoid "exp" and other named functions as variable names
  if(meth == "accum"){
    plot(x = log(x), y = if(metric == "rich"){contr$sac.richness} else{contr$ENS_PIE}
         , ylim = c(0, 50), col = "blue"
         , pch = 17
         , xlab = "Log(Number of Quadrats)"
         , ylab = ifelse(metric == "rich", "Species Richness", "ENS_PIE")
         , main = title)
    # points(log(x), contr, col = "blue", pch = 17)
    points(log(x)
           , y = if(metric == "rich"){Tx$sac.richness}else{Tx$ENS_PIE}
           , col = "red", pch = 15)
    # points(x, exp, col = "red", pch = 15)
    legend("topright", inset = .01
           , legend=c("Control", "Treatment"),
           col=c("blue", "red")
           , pch = c(17, 15),
           box.lty  = 0)

  }
  if(meth == "effect"){
    effect_size = if(metric == "rich"){log(Tx$sac.richness) - log(contr$sac.richness)}
      else{log(Tx$ENS_PIE) - log(contr$ENS_PIE)}
    plot(log(x), effect_size, ylim = c(-0.5, 3), col = "black", pch = 17, xlab = "Log(Number of Quadrats)", ylab = "Log Ratio Effect Size", main = title)
  }
}

# In this tutorial, we will use the package "mobsim" to generate biodiversity at
# scale and replicate a similar methodology to that used in the Knight and Chase
# (2013) paper! We will start things off with a toy model that gives us a sense
# of what we will be working through, then we will write four functions that
# will first give us the tools to create several simulated communities under
# different treatments. The second function will modify our community data so we
# can calculate the cumulative ENS_PIE metric as our survey area increases. Our
# third function will allow us to run these simulations many times to maximize
# the robustness of our dataset.  Finally, our last method will allow us to
# easily create a nice set of plots to display our data in a visually appealing
# way. There is a lot here, so it's possible we won't finish in class. Now with
# all that out of the way, let's get started!

# Install packages: "mobsim", "mobr", "vegan"
#install.package("mobsim")
#install.package("mobr")
#install.package("vegan") # MR: note it's good practice to leave the space after
# the comment hash. Improves readability and lets you uncomment with the text
# editor (Rstudio)

# Load the required packages:
library(mobsim)
library(mobr)
library(vegan)


# First we will simulate a SAD (Species Abundance Distribution) with an upper
# bound of 10 species in the species pool, 320 individuals, a log-normal
# distribution model, and a relatively even species distribution.
sad_lnorm <- sim_sad(s_pool = 10 # number of species
                     , n_sim = 320 # number of individuals
                     , sad_type = "lnorm" # theoretical distribution of relative abundance
                     , sad_coef = list(cv_abund = .5)) # parameters for theoretical distribution

# simulate a species abundance distribution to create a simulated community. The
# SAD follows a log-normal distribution, a continuous curve modelling species'
# relative abundances. To fit a community dataset to this curve, we convert a
# continuous distribution to a discrete dataset (integer number of species with
# integer abundances; you can't have fractions of an individual!)

# This is just a test run to get us used to simulating community data, so the
# numbers we are using are small (a simulated species pool of 10 species, drawn
# from a community of 320 individuals). We will repeatedly simulate a larger
# community later on.

# The coefficient of variation for species abundance (cv_abund) determines the
# size of the "spread" of species from high abundance to low abundance. A low
# number corresponds to a more even community, while a higher number as we will
# use later corresponds to a more uneven community.

# MR: is the paragraph above true (evenness monotonically decreasing with
# cv_abund?)

# Our `cv_abund` here of '.5' represents a moderately even community.


# `sim_sad()` simply returns species abundances. To add biological realism, we
# will also generate locations for each individual


# Let's add some coordinates to our SAD to simulate a spatial distribution
# pattern. We will demonstrate both a random and an aggregated spatial
# distribution. This will allow us to sample our community similar to how a
# "real" community would have individuals spread out over a physical landscape,
# our simulated community will have individuals spread out over a coordinate
# plane, in either a randomized or an aggregated distribution.

# Here we will add spatially random coordinates to our SAD, in a grid space of 2 x 2
practice_sad_coords <- sim_poisson_coords(sad_lnorm, xrange = c(0, 2), yrange = c(0,2))
# Alternately, here we will introduce interspecies aggregation for a different
# resultant community landscape. We will similarly place our (now aggregated)
# SAD over a grid space of 2 x 2
practice_sad_coords_agg <- sim_thomas_coords(sad_lnorm, xrange = c(0, 2), yrange = c(0, 2))


# Let's take a look at what our simulated communities looks like! Note the differences!

plot(practice_sad_coords) # The non-aggregated community
plot(practice_sad_coords_agg) # The aggregated community

# MR: there's really strong intraspecific aggregation here but it looks like no
# interspecific. Is that right? Intended? Might add something about how this is
# ecologically realistic





# Now we want to divide our community into a grid and sample quadrats

# mobsim_sample_quadratst returns a data.frame. We will sample the data in 0.25
# x 0.25 quadrats, then display a plot of the sampling, and sample in a grid
# format. NOTE: This sampling method does not perform the nested sampling that
# the paper mentions.
# What might be possible consequences of this?

# MR: not sure people will know what this question is asking.
practice_sample_data <- sample_quadrats(practice_sad_coords,
                               n_quadrats = 64,
                               quadrat_area = .0625,
                               plot = TRUE,
                               method = "grid",
                               x0 = 0,
                               y0 = 0,
                               delta_x = .25,
                               delta_y = .25)

# MR: style thing: I like to hit return before the comma, this makes things more
# modular/easier to comment out

# Let's do the same again for our aggregated community.
practice_sample_data_agg <- sample_quadrats(practice_sad_coords_agg,
                                        n_quadrats = 64,
                                        quadrat_area = .0625,
                                        plot = TRUE,
                                        method = "grid",
                                        x0 = 0,
                                        y0 = 0,
                                        delta_x = .25,
                                        delta_y = .25)


# MR: pretty sure people will be lost around here. Can you help them see what
# sample_quadrats does?
# the plot you made is a great start.




# Now let's make a Species accumulation curve from the sampled quadrats we just simulated.
practice_sac <- specaccum(practice_sample_data$spec_dat, method = "collector")
practice_sac_agg <- specaccum(practice_sample_data_agg$spec_dat, method = "collector")

# And here we will again plot the results.
# x <- log(1:64) #this seemsed weird
x <- 1:64
plot(x, practice_sac$richness, col = "blue", pch = 17)
points(x, practice_sac_agg$richness, col = "red", pch = 15)

# Why are there jumps in these curves?
# What differences are there between our two SACs? Do these changes match what
# we would expect from the paper's predictions?

# MR: can you point to a specific prediction for people to compare with?


###############################################
## Checkpoint 1 ###############################
###############################################


# Now let's turn what we just did into something a little more representative of
# the paper's methods:

# We will write a function that takes parameters to easily repeat the simulation
# MANY times under different initial conditions (treatments).

# MR: question above seemed pedantic, I would assume that people can appreciate
# sampling over many instances... maybe you can ask more specifically what kind
# of variation you're *not* interested in.

# We will set our initial conditions to default to ones that are similar to
# those used in the paper's control, and much larger than our initial
# simulation. We can (and will) change these in the future, but why might we
# want to set default values for our function?

# MR: I think there's so much going on re: simulation that asking questions
# about coding style is probably going to be distracting here.

treatment <- function(s_pool_input = 45
                      , n_sim_input = 28880
                      , cv_abund_input = 0.5
                      , agg = FALSE){

  # Here we are inputting our parameters into our simulated SAD
  sad_lnorm <- sim_sad(s_pool = s_pool_input
                       , n_sim = n_sim_input
                       , sad_type = "lnorm"
                       , sad_coef = list(cv_abund = cv_abund_input))

  # If we changed the parameter for aggregation to be TRUE instead of FALSE, we
  # want to run our simulation to reflect this. Again, why might we want to take
  # this approach, and make the following conditional statement? Can you think
  # of a better approach? If so, please share!

  if(agg == TRUE){
    sad_coords_agg <- sim_thomas_coords(abund_vec = sad_lnorm
                                        , xrange = c(0, 19), yrange = c(0, 19))
    # A 19 x 19 coordinate space will give us 361 quadrats, very close to the
    # 360 mentioned in the paper

    # Now to sample our community...
    sample_data_agg <- sample_quadrats(sad_coords_agg,
                                   n_quadrats = 361,
                                   quadrat_area = 1,
                                   plot = FALSE, # Why did we set plot = FALSE this time?
                                   method = "grid",
                                   x0 = 0,
                                   y0 = 0,
                                   delta_x = 1,
                                   delta_y = 1)
    return(sample_data_agg)

    # We will now repeat the same steps for the 'else' portion of our
    # conditional, but with a non-aggregated community. Can you think of a way
    # to re-write this using fewer lines of code?

    # MR: I think again, there's so much going on w.r.t. the simulation that
    # this is going to be distracting. But, you might (now that you observed
    # some inefficiency) do the cleanup yourself.
  }else{
    sad_coords <- sim_poisson_coords(abund_vec = sad_lnorm,
                                     xrange = c(0, 19),
                                     yrange = c(0, 19))

    sample_data <- sample_quadrats(sad_coords,
                                   n_quadrats = 361,
                                   quadrat_area = 1,
                                   plot = FALSE,
                                   method = "grid",
                                   x0 = 0,
                                   y0 = 0,
                                   delta_x = 1,
                                   delta_y = 1)
    return(sample_data)
  }


}

###############################################
## Checkpoint 2 ###############################
###############################################

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


###############################################
## Checkpoint 3 ###############################
###############################################

# We've now written two functions! One that makes it really easy to quickly
# simulate communities with different initial conditions, and a second that
# calculates the ENS_PIE metric that respects a cumulative counting of both
# individuals sampled and area sampled. Our next step will be to write a
# function that allows us to easily set our initial conditions and run that
# simulation many times.

# We are re-stating the 'default' conditions that we established in our earlier
# function 'treatment' with the exception of one new one, 'loop.' What would be
# the consequences of not restating these default conditions? What is the logic
# for setting loop to default to '1'?

# Hint 1: Think about how 'line 168' might run if we did not restate the
# 'default' conditions. Hint 2: Think about what might happen if we set the
# default of 'loop' to 10000 repetitions and forgot to change it at a later
# point.

# MR Sorry to mess up line reference.
run_loop <- function(s_pool_input_treatment = 45,
                     n_sim_input_treatment = 28880, # this is a weird number, where does it come from?
                     cv_abund_input_treatment = 0.9,
                     agg_treatment = FALSE,
                     loop = 1){

  sample_sim_plots <- treatment(s_pool_input = s_pool_input_treatment,
                                n_sim_input = n_sim_input_treatment,
                                cv_abund_input = cv_abund_input_treatment,
                                agg = agg_treatment)

  sac <- specaccum(sample_sim_plots$spec_dat, method = "collector")
  ENS_PIE <- ens_PIE(sample_sim_plots$spec_dat)

  for(i in 2:loop){

    sample_sim_plots <- treatment(s_pool_input = s_pool_input_treatment,
                                  n_sim_input = n_sim_input_treatment,
                                  cv_abund_input = cv_abund_input_treatment,
                                  agg = agg_treatment)

    sac$richness <- sac$richness + specaccum(sample_sim_plots$spec_dat
                                             , method = "collector")$richness
    ENS_PIE <- ENS_PIE + ens_PIE(sample_sim_plots$spec_dat)


  }

  sac$richness <- sac$richness / (loop)
  ENS_PIE <- ENS_PIE / (loop)

  results <- data.frame(sac$richness, ENS_PIE)

  return(results)
}

# Let's now run some simulations! Once you are sure you can run the simulations
# without needing to further debug your code, feel free to change the default
# loop to something more interesting... (let's try 10 reps per treatment).
# Notice the difference in time this takes to complete. If you are feeling
# really confident you might try for 100 reps!

# MR: Here, i don't think people will know what to be confident about! Maybe
# again, focus on the ecological stuff rather than the programming stuff.

# Baseline - Pre-treatment
control <- run_loop(loop = 10)
experimental <- run_loop(s_pool_input_treatment = 30, loop = 10)

# Run Simulations:

# Treatment a - Reduce individuals by 50%
control_a <- run_loop(n_sim_input_treatment = 14440, loop = 10)
experimental_a <- run_loop(s_pool_input_treatment = 30, n_sim_input_treatment = 14440, loop = 10)

# Treatment b - Reduce SAD, and add rare species to the community
control_b <- run_loop(cv_abund_input_treatment = 10, loop = 10)
experimental_b <- run_loop(s_pool_input_treatment = 30, cv_abund_input_treatment = 10, loop = 10)

# Treatment c - Dramaticaly reduce the SAD, and make the community highly uneven
control_c <- run_loop(cv_abund_input_treatment = 100, loop = 10)
experimental_c <- run_loop(s_pool_input_treatment = 30, cv_abund_input_treatment = 100, loop = 10)

# Treatment d - Introduce interspecies aggregation into the community
control_d <- run_loop(agg_treatment = TRUE, loop = 10)
experimental_d <- run_loop(s_pool_input_treatment = 30, agg_treatment = TRUE, loop = 10)

###############################################
## Checkpoint 4 ###############################
###############################################
# Now to top it all off, let's turn our data into some pretty graphs! Since we
# have a lot of data to visually represent, let's make a function that
# simplifies making our plots to a single call.

# Let's start with making some 'global' variables that we will use within our
# plotting function. Why would we want these variables to be 'global' in scope
# rather than 'local'? What would change if we made them 'local' instead? How
# might this be desirable in some circumstances but not in others? Say if we
# wanted to make our code more 'general'? It's OKAY if you don't know what
# 'global' or 'local' mean, but try to figure it out from context if you can.

x <- log(1:361)
y1 <- control
y2 <- experimental


# The rest of this function is fairly straightforward. We will be using a
# conditional to determine if we want to plot a SAC or the "log ratio effect
# size" of the ENS_PIE metric with respect to increasing area. We will represent
# our control and experimental treatments differently to make our results more
# visually 'striking'.

plotter <- function(contr, exp, sac_plot = TRUE, title){
  # avoid "exp" and other named functions as variable names
  if(sac_plot == TRUE){

    plot(x, y1$sac.richness, ylim = c(0, 50), col = "blue"
         , pch = 17, xlab = "Log(Number of Quadrats)", ylab = "Species Richness"
         , main = title)
    points(x, contr, col = "blue", pch = 17)
    points(x, y2$sac.richness,  col = "red", pch = 15)
    points(x, exp, col = "red", pch = 15)


  }else{

    y_PIE <- log(y1$ENS_PIE) - log(contr)
    treatment_PIE <- log(y2$ENS_PIE) - log(exp)
    plot(x, y_PIE, ylim = c(-0.5, 3), col = "blue", pch = 17, xlab = "Log(Number of Quadrats)", ylab = "Log Ratio Effect Size", main = title)
    points(x, treatment_PIE, col = "red", pch = 15)

  } ########### I'm still having trouble with the log ratio effect ENS displaying weird results, hence the wide x-axis span. This is really the last thing I'd like to clean up/troubleshoot. Otherwise I'm feeling pretty good about where things sit.
  # I think you're not actually taking the log of x!

}

# NOTE: Run the two plot groups separately if you plan to copy-paste the code.
# Plot group 1: Species Accumulation Curves
SAC_plot.par <- par(mfrow=c(2, 2))

ya1 <- control_a$sac.richness # To my eye this is clunky.
# You already have a plotter function so why keep renaming things
# instead of just feeding them to your function?
ya2 <- experimental_a$sac.richness
plotter(ya1, ya2, title = "Treatment A")

yb1 <- control_b$sac.richness
yb2 <- experimental_b$sac.richness
plotter(yb1, yb2, title = "Treatment B")

yc1 <- control_c$sac.richness
yc2 <- experimental_c$sac.richness
plotter(yc1, yc2, title = "Treatment C")

yd1 <- control_d$sac.richness
yd2 <- experimental_d$sac.richness
plotter(yd1, yd2, title = "Treatment D")

legend("bottomright", inset = .01, legend=c("Control", "Experimental"),
       col=c("blue", "red"), pch = c(17, 15),
       box.lty=0)
par(SAC_plot.par)

# Plot group 2: log ratio effect size for ENS_PIE and sampling scale

PIE_plot.par <- par(mfrow=c(2, 2)) #nice!
ya1 <- control_a$ENS_PIE
ya2 <- experimental_a$ENS_PIE
plotter(ya1, ya2, FALSE, "Treatment A")

yb1 <- control_b$ENS_PIE
yb2 <- experimental_b$ENS_PIE
plotter(yb1, yb2, FALSE, "Treatment B") ####### I'm unclear why there is a substantial gap between the "control" group and the "experimental" group. Additionally, why these numbers are so different from those reached in the paper.

yc1 <- control_c$ENS_PIE
yc2 <- experimental_c$ENS_PIE
plotter(yc1, yc2, FALSE, "Treatment C") ####### Similar to above, but gap between control and experimental seems to be smaller. The overall pattern of ENS_PIE being inversely related to community evenness decreases holds, which I suppose is the most important thing.

yd1 <- control_d$ENS_PIE
yd2 <- experimental_d$ENS_PIE
plotter(yd1, yd2, FALSE, "Treatment D") ####### This last one is especially weird. It seems to follow an exponential decay curve instead of how the paper's curve drops off... Thoughts?

legend("topright", inset = .01, legend=c("Control", "Treatment"),
       col=c("blue", "red"), pch = c(17, 15),
       box.lty=0)
par(PIE_plot.par)

##----------------------------------------------------------

####### Any thoughts on finishing and wrapping it up neatly? Is there anything I should cut you think? Anything worth mentioning that I missed?

##----------------------------------------------------------
###############################################
######### THE END!!! ##########################
###############################################


# CONGRATULATIONS! You have now written several functions to run a community
# simulation numerous times that models a variety of communities under different
# treatments, analyzed the data generated from those simulations, and displayed
# the data neatly (hopefully with minimal headaches)! By using functions to
# simplify repeatedly calling the same (or very nearly so) blocks of code, you
# are witnessing firsthand one of the most powerful tools that programming gives
# us. The power of automation! This allows us to assign the drudgery of carrying
# out the same task many times to the computer, both freeing us up to work on
# other problems and (at least in theory) reducing the potential for mistakes.
# We covered a lot of concepts and code today, so no sweat if a lot of this was
# new for you or you struggled with some of it. I hope you found this tutorial
# both engaging and valuable, and feel free to email me any questions or
# comments you have about the code at chacowhitaker@gmail.com. Happy coding!

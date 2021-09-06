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

# some bulky functions
source("ChaseKnight_helper.R")


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
    sad_coords <- sim_thomas_coords(abund_vec = sad_lnorm
                                        , xrange = c(0, 19), yrange = c(0, 19))
    # A 19 x 19 coordinate space will give us 361 quadrats, very close to the
    # 360 mentioned in the paper


  }
  else{

    # We will now repeat the same steps for the 'else' portion of our
    # conditional, but with a non-aggregated community.
    sad_coords <- sim_poisson_coords(abund_vec = sad_lnorm,
                                     xrange = c(0, 19),
                                     yrange = c(0, 19))


  }

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


###############################################
## Checkpoint 2 ###############################
###############################################



# We've now written a function that makes it easy to quickly
# simulate communities with different initial conditions, Our next step will be to write a
# function that allows us to easily set our initial conditions and run that
# simulation many times.

run_loop <- function(s_pool_input_treatment = 45,
                     n_sim_input_treatment = 28880, # this is a weird number, where does it come from?
                     cv_abund_input_treatment = 0.5,
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

  sac$richness <- sac$richness / (loop) # take average, richness already summed across iterations
  ENS_PIE <- ENS_PIE / (loop) # take average

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
control_45 <- run_loop(loop = 10)
control_30 <- run_loop(s_pool_input_treatment = 30, loop = 10)

# Run Simulations:

# Treatment a - Reduce individuals by 50%
reduce_density_45 <- run_loop(n_sim_input_treatment = 14440, loop = 10) # In the MS, isn't control always identical?
reduce_density_30 <- run_loop(s_pool_input_treatment = 30, n_sim_input_treatment = 14440, loop = 10)

# Treatment b - Reduce evenness; now there are very rare species in the community
uneven_45 <- run_loop(cv_abund_input_treatment = 1, loop = 10) # MR: I modified this from 10 to 1
uneven_30 <- run_loop(s_pool_input_treatment = 30, cv_abund_input_treatment = 1, loop = 10)

# Treatment c - Dramatically reduce evenness
dramatic_45 <- run_loop(cv_abund_input_treatment = 10, loop = 10)
dramatic_30 <- run_loop(s_pool_input_treatment = 30, cv_abund_input_treatment = 10, loop = 10)

# Treatment d - Introduce interspecies aggregation into the community
agg_45 <- run_loop(agg_treatment = TRUE, loop = 10)
agg_30<- run_loop(s_pool_input_treatment = 30, agg_treatment = TRUE, loop = 10)

###############################################
## Checkpoint 3 ###############################
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

x <- 1:361

# The rest of this function is fairly straightforward. We will be using a
# conditional to determine if we want to plot a SAC or the "log ratio effect
# size" of the ENS_PIE metric with respect to increasing area. We will represent
# our control and experimental treatments differently to make our results more
# visually 'striking'.


# NOTE: Run the two plot groups separately if you plan to copy-paste the code.
# Plot group 1: Species Accumulation Curves
SAC_plot.par <- par(mfrow=c(2, 2))

# effect of reducing density
plotter(control_45
        , reduce_density_45
        , x = x
        , meth = "accum", metric = "rich", title = "density reduction, 45 spp")

plotter(control_45
        , reduce_density_45
        , x = x
        , meth = "effect", metric = "rich", title = "density reduction, 45 spp")

plotter(control_45
        , reduce_density_45
        , x = x
        , meth = "accum"
        , metric = "ENS_PIE", title = "density reduction, 45 spp")

plotter(control_45, reduce_density_45,
        x = x,
        meth = "effect"
        , metric = "ENS_PIE"
        , title ="density reduction, 45 spp")


# effect of intraspecific aggregation
plotter(control_45
        , agg_45
        , x = x
        , meth = "accum", metric = "rich", title = "density reduction, 45 spp")

# y-axis truncated

plotter(control_45
        , agg_45
        , x = x
        , meth = "effect", metric = "rich", title = "density reduction, 45 spp")

plotter(control_45
        , agg_45
        , x = x
        , meth = "accum"
        , metric = "ENS_PIE", title = "density reduction, 45 spp")

# y-axis truncated
plotter(control_45
        , agg_45,
        x = x,
        meth = "effect"
        , metric = "ENS_PIE"
        , title ="density reduction, 45 spp")


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

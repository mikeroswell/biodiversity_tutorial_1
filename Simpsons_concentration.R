# This is an R script/ tutorial to explore Simpson's concentration. There are 4 "Checkpoints" in it and you will work through it in pairs, regrouping when most groups have hit each checkpoint. So make sure you let Michael/Dan know each time you hit a checkpoint.

# you do not need to submit this tutorial, or any going forward. You are welcome to for feedback at any point!

# There are lots of comments, and many contain questions or directions. I suggest you do your best to add something in response to these questions and directions. This is supposed to be interactive.

###########################################################

# Always good to list packages at the top
# I'll load the required packages with `library()` which returns an error if you haven't installed the package already. Most R packages are on CRAN and can be downloaded and installed with `install.packages()` (and an internet connection!) #See below for an example

library(ggplot2) #if you don't have ggplot2 installed, I'd recommend just uncommenting and running
# install.packages("tidyverse")
# you could also uncomment and run
# install.packages("ggplot2")
# the first install.packages would mean you do not separately need to install these, the second would
library(dplyr)
library(purrr)

# In many of the papers we'll read in the first few weeks, there is a lot of attention to variations on Simpson's index. Simpson (1949) described a measure of "true" diversity (i.e. the diversity of a entire assemblage) that could be estimated in an unbiased manner from a sample of any size. He described it in terms of the "concentration" of abundance in a few species (low diversity) or a smaller concentration (and higher diversity). The measure is "the probability that two individuals chosen at random and independently from the population will be found to belong to the same group." For our purposes we'll usually imagine group to mean species (although this metric is key in genetic diversity (where group represents "allele", and in many other applications)).

# We'll review a tiny bit of probability theory.

# If you have two independent events, occurring with probabilities a and b, what is the probability of a and b both occurring? You can write an answer in a comment below:


# For many problems, it's easier to think about the probability of a and b *not* both happening. What is the formula for this?

# Finally, what is the probability of one of two disjoint events (that is, if one happens, the other doesn't), a, and b?


# Ok. Now we can think about the probability of selecting two individuals of the same species (where we assume that the only thing affecting which species we see is that species' relative abundance)

# we'll begin with a really simple system, with just two species.

a <- 0.5 #relative abundance of species a
b <- 1-a # relative abundance of species b

# probability that if we sample two individuals, they will both be "a"

a*a

# what is the probability that both will be "b'?

# what is the probability that both individuals will be of the the same "group" (i.e., either both a or both b?)



###############################################
## Checkpoint 1 ###############################
#############################################

# Simpson's concentration is the generalization of this probability (i.e. for an arbitrary number of species / groups)

# here is a function to compute Simpson's concentration from relative abundances

simpson = function(x){
  if(sum(x) !=1 ){stop("relative abundances must sum to 1")} #this returns an error if relative abundances don't add up
  sum(x^2)
}

# let's examine how our function works with some test data

rich <- 10

test.data <- data.frame(
  species = 1:rich # give each species a name, using an integer value for now
  , abundance = rep(100, rich)
)

# add relative abundances

test.data <- test.data %>% mutate(rel_abundance = abundance/sum(abundance))

test.data

simpson(test.data$rel_abundance) #pretty boring when all species have the same abundance!


# give each species a random abundance instead
test.data$abundance<- sample(1:1000, rich, replace =T)

# recompute relative abundance
test.data <- test.data %>% mutate(rel_abundance = abundance/sum(abundance))

test.data
simpson(test.data$rel_abundance)

# it went up! Is this community more or less diverse than the one where all the species had the same abundance?


######################################
## checkpoint 2 ######################
####################################

# what happens here?
simpson(test.data$abundance)

# Side note: What does the function set.seed() do in R? What might it do for us here?

# Simpson was concerned about estimating the divesity of a whole community based on a sample from it. What if we pretend that test.data represents our community. Let's sample from it

# here's a function from an R package I'm working on that uses R's sample command to do this. Of course in the R package I've had to say what all the things in here do. If you're curious, I invite you to comment this further! And if you're a programming/R/C++ guru and see an easy way to make this way more efficient, let me know and maybe we can find some fun collaborations.

sample_finite <- function(ab, size = sum(ab)){
  inds <- unlist(lapply(1:length(ab), function(x){
    rep(x, ab[x])
  }))
  sam <- sample(inds, size = size, replace = FALSE)
  ss <- unlist(lapply(1:length(ab), function(y){
    length(which(sam == y))
  }))
  return(ss)
}

# we'll use it to add some sample data to our test.data
sample_size <- 100

test.data <- test.data %>% mutate(sample_abundane = sample_finite(abundance, size = sample_size))

test.data

# let's look at what happens if we compute simpson's index on the sample abundances

simpson(test.data$sample_abundance/sum(test.data$sample_abundance))

# whoops! Why did it return an error?

# I made a typo! there is a missing "c" when i added sample_abundance to test.data. Please fix it and rerun (and the sample_abundane column will still be there).

##########################################################
### Checkpoint 3 #########################################
##########################################################

# was the the Simspons's index higher for the sample abundances, or the true abundances?

# in my research I'm often interested by how indices like this behave. Let's examine whether there is a more general answer to the question above.

# I learned to program initially using "for" loops. Now I just use "map" and especially "map_dfr" and I'm much happier now. My ability to navigate and troubleshoot multidimensional arrays is just not very good. I can often figure out what's working and not in a data.frame, though :-)

sample_sizes <- floor(10^seq(1,3, 0.5)) # when simulating/testing, you can save a lot of time/computation by exploring parameter space in a more log/exponential way (the difference between 10 and 20 might really matter, I doubt the difference between 567 and 767 is going to matter a whole lot, though)

trials <- 100 # for each sample size, how many times do we want to test

simpson_bias <- map_dfr(sample_sizes, function(ss){
  map_dfr(1:trials, function(trial){
    sample_abundance = sample_finite(test.data$abundance, ss)
    data.frame(sample_size = ss
               , naive_simpsons = simpson(sample_abundance/sum(sample_abundance)))
  })
})

# the function "head" is great, there are lots of other ways to see a small snapshot of a bigger dataset (glimpse, head, tail, table, summary, str, names, dim, class...). When I'm writing code, though, I often just look at the whole object (or whatever R shows me) in the console. Not that I am endorsing this.

# did our simulation work?
simpson_bias

# let's make a quick graph of what's going on here.

simpson_bias %>% ggplot(
  aes(sample_size, naive_simpsons)) +
  geom_jitter(alpha = 0.3) +
  geom_hline(yintercept = simpson(test.data$rel_abundance)) +
  theme_classic()

########################################
## Checkpoint 4 #######################
######################################


# Simpson gave an unbiased estimator for his index. Here's a function that should compute it for _abundance_ data (not relative abundances!)

simpson_corrected <- function(x){
  if(!is.integer(x)){stop("Simpson's corrected estimator requires integer abundances")}
  sum(x*(x-1))/(sum(x)*(sum(x)-1))
}

# let's just make sure it gives us reasonable values. This is often how I test my code on the fly (I simply ask myself, does it give a close-looking answer)
simpson_corrected(test.data$sample_abundance)

# Go head and edit the code above that created the object "simpson_bias" to compute both the naive, uncorrected Simpson's estimator (we did this before) as well as the new, corrected estimator, for the samples we created before.

# Plot both on the same axes, but give the corrected and naive estimates different colors.

# couple quick things before you go!

# As communities get more diverse, does Simpsons' index go up or down?
# Scroll down once you've decided.













# The fact that Simpson's concentration is smaller for more diverse communities is at the very least annoying, given that people use it as a measure of diversity. Two variations of Simpson's index avoid this problem.

# the Gini-Simpson index is the probability that two individuals randomly selected are _not_ from the same species, or 1-Simpson. This is, properly speaking, a probability of interspecific encounter. Plugging in Simpson's estimator does give an unbiased estimate for samples of any size.

# the reciprocal of Simpson's index is commonly used. I like to call this "Hill-Simpson diversity" (We'll dig into Hill diversities). Hill diversities are also known as "effective numbers of species."

# A bunch of people working on biodiversity research, including Jonathan Chase, one of the authors of the paper we're reading for next week, call the reciprocal of Simpson's index (Specificially, I think they usually mean this in terms of the naive estimator) "ENS-PIE", for "the Effective Number of Species based on the Probability of Interspecific Encounter" Although this is really verbose, confusing terminology, I think the rest of the paper is really clear!

# As far as I know, no truly unbiased estimator for Hill-Simpson diversity is known. However, plugging in Simpson's estimator (then taking the reciprocal) works pretty well, and Chao and Jost 2015 proposed some general Hill diversity estimators that have very low bias for Hill-Simpson diversity.

###########################
### THE END ###############
###########################



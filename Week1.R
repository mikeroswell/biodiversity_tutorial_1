# Always good to list packages at the top
# I'll load the required packages with `library()` which returns an error if you haven't installed the package already. Most packages are on CRAN and can be installed with `install.packages()` See below for an example

library(ggplot2) #if you don't have ggplot2 installed, I'd recommend just uncommenting and running
# install.packages("tidyverse")
# you could also uncomment and run
# install.packages("ggplot2")
library(dplyr)
library(purrr)


# # in R, if a package is loaded you don't need to reference it explicitly. However, it can be a good defensive programming technique to do that all the time. The syntax packageA::functionA tells R to use a function called `functionA` from a package called `packageA`
# # there are plenty of ways to install packages from github repositories, but devtools comes in handy eventually
# devtools::install_github("mikeroswell/MeanRarity")
# library(MeanRarity)

# In many of the papers we'll read in the first few weeks, there is a lot of attention to variations on Simpson's index. Simpson (1949) described a measure of "true" diversity (i.e. the diversity of a entire assemblage) that could be estimated in an unbiased manner from a sample of any size. He described it in terms of the "concentration" of abundance in a few species (low diversity) or a smaller concentration (and higher diversity). The measure is "the probability that two individuals chosen at random and independnetly from the population will be found to belong ot the same group." For our purposes we'll usually imagine group to mean species (although this metric is key in genetic diversity (where group represents "allele", and in many other applications)).

# We'll review a tiny bit of probability theory.

# If you have two independent events, occurring with probabilities a and b, what is the probability of a and b both occurring? You can write an answer in a comment below:


# For many problems, it's easier to think about the probability of a and b *not* both happening. What is the formula for this?

# Finally, what is the probability of one of two disjoint events (that is, if one happens, the other doesn't), a, and b?


# Ok. Now we can think about the probability of selecting two individuals of the same species (where we assume that the only thing affecting which species we see is that species' relative abundance)

# we'll begin with a really simple system, with just two species.

a <- 0.5 #relative abundance of species a
b <- 1-a # relative abundance of species b

# probability that if we sample two individusls, they will both be "a"

a*a

# what is the probability that both will be "b'?

# what is the probability that both individuals will be the same (i.e., either both a or both b?)

# Simpson's concentration is the generalization of this probability (i.e. for an arbitrary number of species)

# here is a function to compute Simpson's concentration from relativea abundances

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


# what happens here?
simpson(test.data$abundance)

# Side note? What does the function set.seed() do in R? What might it do for us here?

# Simpson was concerned about estimating the diveristy of a whole community based on a sample from it. What if we pretend that test.data represents our community. let's sample from it

# here's a function from an R package I'm working on that uses R's sample command to do this

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

# I made a typo! there is a missing "c" when i added sample_abundance to test.data. you can fix it and rerun (and the sample_abundane column will still be there) run again.

# was the the simspons's index higher for the sample abundances, or the true abundances?

# in my research I'm often interested by how indices like this behave. Let's examine whether there is a more general answer to the question above.

# I learned to program initially using "for" loops. Now I just use "map" and especially "map_dfr" and I'm much happier now. My ability to navigate and troubleshoot multidimensional arrays is just not very good. I can often figure out what's working and not in a data.frame, though :-)

sample_sizes <- floor(10^seq(1,3, 0.5)) # when simulating/testing, you can save a lot of time/computation by exploring parameter space in a more log/exponential way (the difference between 10 and 20 might really matter, I doubt the difference between 567 and 767 is going to matter a whole lot, though)

trials <- 100

simpson_bias <- map_dfr(sample_sizes, function(ss){
  map_dfr(1:trials, function(trial){
    sample_abundance = sample_finite(test.data$abundance, ss)
    data.frame(sample_size = ss
               , naive_simpsons = simpson(sample_abundance/sum(sample_abundance)))
  })
})

simpson_bias

# let's make a quick graph of what's going on here.

simpson_bias %>% ggplot(
  aes(sample_size, naive_simpsons)) +
  geom_jitter(alpha = 0.3) +
  geom_hline(yintercept = simpson(test.data$rel_abundance)) +
  theme_classic()


# Simpson gave an unbiased estimator for his index. Here's a function that should compute it for _abundance_ data (not relative abundances!)

simpson_corrected <- function(x){
  if(!is.integer(x)){stop("Simpson's corrected estimator requires integer abundances")}
  sum(x*(x-1))/(sum(x)*(sum(x)-1))
}

# let's just make sure it gives us reasonable values. This is often how I test my code on the fly (does it give a close-looking answer?)
simpson_corrected(test.data$sample_abundance)

# Go head and

# edit the code above to compute both the naive, uncorrected Simpson's estimator (we did this before) as well as the new, corrected estimator, for the samples we created before. Plot both on the same axes, but give the corrected and naive estimates different colors.




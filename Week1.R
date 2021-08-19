# Always good to list packages at the top
# I'll load the required packages with `library()` which returns an error if you haven't installed the package already. Most packages are on CRAN and can be installed with `install.packages()` See below for an example

library(ggplot2) #if you don't have ggplot2 installed, I'd recommend just uncommenting and running
# install.packages("tidyverse")
# you could also uncomment and run
# install.packages("ggplot2")
library(dplyr)
library(purrr)
# in R, if a package is loaded you don't need to reference it explicitly. However, it can be a good defensive programming technique to do that all the time. The syntax packageA::functionA tells R to use a function called `functionA` from a package called `packageA`
# there are plenty of ways to install packages from github repositories, but devtools comes in handy eventually
devtools::install_github("mikeroswell/MeanRarity")
library(MeanRarity)

# In many of the papers we'll read in the first few weeks, there is a lot of attention to variations on Simpson's index. Simpson (1949) described a measure of "true" diversity (i.e. the diversity of a entire assemblage) that could be estimated in an unbiased manner from a sample of any size. He described it in terms of the "concentration" of abundance in a few species (low diversity) or a smaller concentration (and higher diversity). The measure is "the probability that two individuals chosen at random and independnetly from the population will be found to belong ot the same group." For our purposes we'll usually imagine group to mean species (although this metric is key in genetic diversity (where group represents "allele", and in many other applications)).

# We'll review a tiny bit of probability theory.

# If you have two independent events, occurring with probabilities a and b, what is the probability of a and b both occurring? You can write an answer in a comment below:


# For many problems, it's easier to think about the probability of a and b *not* both happening. What is the formula for this?

# Ok, let's imagine a really simple case. Let's imagine we have a community with two species with equal abundance. What is the probability that if you grabbed two random individuals from this community, they would both be from the same species?
# compute this with the base R function `dbinom()`


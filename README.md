
# MIRT4FC

<!-- badges: start -->
<!-- badges: end -->

The objective of MIRT4FC is to efficiently implement various forced-choice models using the istem algorithm.Currently, it includes Thurstone's Item Response Theory (TIRT, Brown et al., 2011) Model, Multi-Unidimensional Pairwise Preference Two Parameter Logistic Model (MUPP-2PLM, Morillo et al., 2016), Multi-Unidimensional Pairwise Preference Generalized Graded Unfolding Model (MUPP-GGUM ,Stark et al., 2005) and Generalized Graded Unfolding-RANK Model (GGUM-RANK, Lee et al., 2018), and we plan to continue updating and adding new models in the future. In addition to item parameter estimation capabilities, our R package also offers the ability to estimate ability parameters using MAP, EAP, and MLE methods. It can generate simulated response matrices, calculate standard errors (SE) for both ability and item parameters, and include a set of empirical data.

## Reference manual


Reference manual can be found in MIRT4FC.pdf

## Installation

The development version of MIRT4FC can be installed in the following way:

``` r
install.packages("devtools")
library("devtools")
install_github("xpy-ECNU/MIRT4FC")
1
library(MIRT4FC)
```
## Example

This is a simple example that demonstrates the entire process of simulating the latent trait parameters of 1000 participants in 6 dimensions, each with 10 statements and 20 triplet blocks, generating a response matrix, and finally estimating project parameters:

``` r
#################Example1: A simulation example based on the MUPP-2PL model###############
########Set simulation information
library (MIRT4FC)
D <- 6                           	# Dimension
nitem.per.dim <- 10                    # Iems number per dimension
nblock <- D * nitem.per.dim / 3      		# Blocks number
set.seed(123456)                     	# Set random seed
# Simulate block-item-demension correspondence table
BID <- data.frame (Block = rep (1:nblock,each=3),
                  Item=rep (1:3, nblock),
                  Dim=c(combn(D,3) [, sample(choose(D,3), nblock,replace = TRUE)]))
# Simulate item parameter truth value
item.par <- data.frame (a = seq_len (D * nitem.per.dim))
item.par <- within (item.par, {
  a <- runif (D*nitem.per.dim,0.7,3)
  b <- rnorm (D*nitem.per.dim)
  d <- a*b
})
item.par$d <- c (t (aggregate (item.par$d, by=list (BID$Block), function(x)x-mean(x)) [, -1]))
N <- 1000                          	# Sample number
v <- matrix (0.5, D, D)                   	# Intertrait correlation
diag (v) <- 1
# Simulate latent trait parameter truth value
theta <- mvnfast::rmvn (N, seq(-1, 1, length.out = D),sigma = v)

######## Generate a simulated dataset
Y <- data.sim (item.par, theta, BID, blocksize = 3, res = 'rank')

######## Item parameter estimation
fit <- StEM (Y, BID, maxitr = 100, blocksize = 3, res = 'rank', fix.sigma = TRUE, cores = 1)

##################Example2:an empirical example for the 2PL-RANK model####################
##for the paper ***A 2PLM-RANK Multidimensional Forced-choice Model*** and its Fast Estimation##
##########################################################################################
######## Read dataset
library (MIRT4FC)
Y <- data("MAP_data")
######## Item parameter estimation
fit <- StEM (Y, BID, maxitr = 150, blocksize = 3, res = 'rank', fix.sigma = TRUE, cores = 1)

```


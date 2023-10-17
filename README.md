
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
install_github("SAL-Lab-ECNU/MIRT4FC")
1
library(MIRT4FC)
```
## Example1

This is a simple example that demonstrates the entire process of simulating the latent trait parameters of 1000 participants in 6 dimensions, each with 10 statements and 20 triplet blocks, generating a response matrix, and finally estimating item parameters:

``` r
#################A simulation example based on the MUPP-2PL model###############
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
```

## Example2

This is a empirical example for the paper "A 2PLM-RANK Multidimensional Forced-choice Model and its Fast Estimation":
``` r
#########################an empirical example for the 2PL-RANK model######################
######## Read dataset
library (MIRT4FC)
Y <- data("MAP_data")
BID <- data.frame(
  Block = rep(1:88,each=3),
  Item = rep(1:3,times=88),
  Dim = c(11,	1,	14,	15,	13,	18,	6,	2,	7,	1,	21,	24,	23,	9, 22,	5,	3,	8,	6,	18,	20,	9,	7,	5,	19,	1,	9,	14,	 22,	10,	23,	21,	8,	3,	20,	10,	22,	3,	17,	4,	23,	24,
           8,	4,	24,	21,	12,	10,	9,	11,	8,	17,	12,	2,	15,	11,	14,  24,	22,	15,	13,	16,	14,	22,	4,	14,	10,	7,	6,	14,	19,  6,	4,	13,	3,	14,	15,	2,	8,	3,	11,	18,	23,	20,	24,	15,
           22,	3,	20,	1,	14,	21,	16,	18,	4,	13,	16,	18,	5,	1, 8,	23,	2,	24,	11,	19,	23,	15,	12,	11,	10,	20,	9,	21,	10,	 4,	16,	7,	2,	3,	12,	16,	10,	6,	13,	16,	21,	16,	17,	20,
           10,	19,	13,	7,	5,	15,	11,	23,	24,	8,	6,	11,	2,	19,	15,  17,	20,	18,	9,	7,	12,	5,	9,	7,	22,	17,	24,	16,	6,	17, 13,	23,	7,	15,	17,	5,	8,	19,	14,	18,	3,	12,	22,	4,	5,	
           21,	13,	1,	23,	9,	13,	11,	22,	9,	19,	21,	8,	21,	12,	6,	16,	1,	23,	9,	1,	19,	2,	6,	11,	18,	24,	10,	7,	18,	17,	 5,	7,	3,	22,	4,	2,	8,	20,	17,	15,	8,	3,	24,	12,	10,	22,
           2,	11,	23,	19,	19,	13,	17,	6,	20,	24,	9,	17,	5,	20,	12,	6,  19,	18,	16,	15,	21,	7,	5,	1,	18,	2,	4,	14,	1,	13,	12,	 16,	12,	2,	20,	4,	5,	10,	4,	1,	21,	14,	3))
######## Item parameter estimation
fit <- StEM (Y, BID, maxitr = 150, blocksize = 3, res = 'rank', fix.sigma = TRUE, cores = 1)

```


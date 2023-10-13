
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
library(MIRT4FC)
```
## Example

This is a simple illustration that demonstrates how to resolve a typical issue:

``` r
library(MIRT4FC)
#dimension
D <- 6
#number of blocks for each dimension
nitem.per.dim <- 10
#block number
nblock <- D * nitem.per.dim / 3
set.seed(123456)
#Analog true value
item.par <- data.frame(a=seq_len(D*nitem.per.dim))
item.par <- within(item.par,{
 a <- runif(D*nitem.per.dim,0.7,3)
 b <- rnorm(D*nitem.per.dim)
 d <- a*b
})
#Columns are "Block", "Item" and "Dimensions"
BID <- data.frame(Block=rep(1:nblock,each=3),
                 Item=rep(1:3,nblock),
                 Dim=c(combn(D,3)[,sample(choose(D,3),nblock,replace = TRUE)]))
#Put restrictions in place so that the total of the intercepts inside a block equals 0
item.par$d <- c(t(aggregate(item.par$d,by=list(BID$Block),function(x)x-mean(x))[,-1]))
#number of sample
N <- 1000
#True values of correlations between trait dimensions
v <- matrix(0.5,D,D)
diag(v) <- 1
#True values of traits
theta <- mvnfast::rmvn(N,seq(-1,1,length.out = D),sigma = v)
#Simulate a response matrix
Y <- data.sim(item.par,theta,BID,blocksize=3,res='rank')
#Applying the iStEM algorithm to estimate parameters 
x <- StEM(Y,BID,maxitr = 100,blocksize=3,res='rank',fix.sigma = TRUE,cores=1)
```

